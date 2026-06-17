//! Columnar **Parquet** transcript cache (replaces fastVEP's bincode cache).
//!
//! One row per transcript: queryable scalar columns (chrom/start/gene/biotype/…)
//! plus a `model` JSON column for faithful `Transcript` reconstruction. Portable,
//! inspectable, `read_parquet`-able from DuckDB — and DuckDB-free here (uses the
//! `arrow`/`parquet` crates, not libduckdb). See docs/DESIGN.md §5.

use arrow::array::{
    Array, BinaryArray, BooleanArray, Int64Array, Int8Array, RecordBatch, StringArray,
};
use arrow::datatypes::{DataType, Field, Schema};
use fastvep_genome::Transcript;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;
use std::error::Error;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub(crate) fn cache_path(gff3: &str) -> PathBuf {
    PathBuf::from(format!("{gff3}.transcripts.parquet"))
}

/// Cache is usable if it exists and is at least as new as the GFF3.
pub(crate) fn is_fresh(cache: &Path, gff3: &Path) -> bool {
    match (
        cache.metadata().and_then(|m| m.modified()),
        gff3.metadata().and_then(|m| m.modified()),
    ) {
        (Ok(c), Ok(g)) => c >= g,
        _ => false,
    }
}

/// Transcripts serialized per Arrow row group. Bounds the writer's *transient*
/// memory to one batch's worth of serialized + Arrow copies instead of
/// materializing a copy of the **entire** ~280k-transcript model at once.
const BATCH_ROWS: usize = 8192;

/// Default zstd level for the cache's `model` column — a sane space/build-time
/// tradeoff. Callers pass an explicit level via `vep_load_cache(gff, fasta,
/// distance, zstd_level)`; 1 = fastest build / largest file, up to 22 = smallest
/// file / slowest build.
pub(crate) const DEFAULT_ZSTD: i32 = 3;

fn build_batch(schema: &Arc<Schema>, chunk: &[Transcript]) -> Result<RecordBatch, Box<dyn Error>> {
    let mut tid = Vec::with_capacity(chunk.len());
    let (mut chrom, mut gid, mut sym, mut bt, mut model) = (
        Vec::with_capacity(chunk.len()),
        Vec::with_capacity(chunk.len()),
        Vec::with_capacity(chunk.len()),
        Vec::with_capacity(chunk.len()),
        Vec::<Vec<u8>>::with_capacity(chunk.len()),
    );
    let (mut start, mut end) = (Vec::with_capacity(chunk.len()), Vec::with_capacity(chunk.len()));
    let mut strand = Vec::with_capacity(chunk.len());
    let mut canon = Vec::with_capacity(chunk.len());
    for t in chunk {
        tid.push(t.stable_id.to_string());
        chrom.push(t.chromosome.to_string());
        start.push(t.start as i64);
        end.push(t.end as i64);
        strand.push(t.strand.as_int());
        gid.push(t.gene.stable_id.to_string());
        sym.push(t.gene.symbol.as_deref().unwrap_or("").to_string());
        bt.push(t.biotype.to_string());
        canon.push(t.canonical);
        // Compact binary model (matches fastVEP's bincode cache): ~2× smaller than
        // the old JSON text and far cheaper to (de)serialize. The scalar columns
        // above stay plain so the cache is still inspectable via `read_parquet`.
        model.push(bincode::serialize(t)?);
    }
    let model_refs: Vec<&[u8]> = model.iter().map(|v| v.as_slice()).collect();
    Ok(RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(tid)),
            Arc::new(StringArray::from(chrom)),
            Arc::new(Int64Array::from(start)),
            Arc::new(Int64Array::from(end)),
            Arc::new(Int8Array::from(strand)),
            Arc::new(StringArray::from(gid)),
            Arc::new(StringArray::from(sym)),
            Arc::new(StringArray::from(bt)),
            Arc::new(BooleanArray::from(canon)),
            Arc::new(BinaryArray::from(model_refs)),
        ],
    )?)
}

pub(crate) fn save(
    transcripts: &[Transcript],
    path: &Path,
    zstd_level: i32,
) -> Result<(), Box<dyn Error>> {
    let level = zstd_level.clamp(1, 22);
    let schema = Arc::new(Schema::new(vec![
        Field::new("transcript_id", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end_pos", DataType::Int64, false),
        Field::new("strand", DataType::Int8, false),
        Field::new("gene_id", DataType::Utf8, false),
        Field::new("gene_symbol", DataType::Utf8, false),
        Field::new("biotype", DataType::Utf8, false),
        Field::new("canonical", DataType::Boolean, false),
        Field::new("model", DataType::Binary, false),
    ]));
    let props = WriterProperties::builder()
        .set_max_row_group_row_count(Some(BATCH_ROWS))
        .set_compression(Compression::ZSTD(ZstdLevel::try_new(level)?))
        .build();
    let mut w = ArrowWriter::try_new(File::create(path)?, schema.clone(), Some(props))?;
    // Stream one row-group batch at a time: each batch's serialized + Arrow copies
    // are dropped before the next is built, so peak memory is O(BATCH_ROWS).
    for chunk in transcripts.chunks(BATCH_ROWS) {
        w.write(&build_batch(&schema, chunk)?)?;
    }
    w.close()?;
    Ok(())
}

pub(crate) fn load(path: &Path) -> Result<Vec<Transcript>, Box<dyn Error>> {
    let reader = ParquetRecordBatchReaderBuilder::try_new(File::open(path)?)?.build()?;
    let mut out = Vec::new();
    for batch in reader {
        let batch = batch?;
        let col = batch
            .column_by_name("model")
            .ok_or("transcript cache: missing 'model' column")?;
        // New caches store the model as a bincode BLOB; older caches stored JSON
        // text. Dispatch on the Arrow type so a stale JSON cache still loads
        // (back-compat) without forcing a rebuild.
        match col.data_type() {
            DataType::Binary => {
                let m = col
                    .as_any()
                    .downcast_ref::<BinaryArray>()
                    .ok_or("transcript cache: 'model' is not BLOB")?;
                for i in 0..m.len() {
                    out.push(bincode::deserialize(m.value(i))?);
                }
            }
            DataType::Utf8 => {
                let m = col
                    .as_any()
                    .downcast_ref::<StringArray>()
                    .ok_or("transcript cache: 'model' is not VARCHAR")?;
                for i in 0..m.len() {
                    out.push(serde_json::from_str(m.value(i))?);
                }
            }
            other => return Err(format!("transcript cache: unsupported 'model' type {other:?}").into()),
        }
    }
    Ok(out)
}
