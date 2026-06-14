//! Columnar **Parquet** transcript cache (replaces fastVEP's bincode cache).
//!
//! One row per transcript: queryable scalar columns (chrom/start/gene/biotype/…)
//! plus a `model` JSON column for faithful `Transcript` reconstruction. Portable,
//! inspectable, `read_parquet`-able from DuckDB — and DuckDB-free here (uses the
//! `arrow`/`parquet` crates, not libduckdb). See docs/DESIGN.md §5.

use arrow::array::{Array, BooleanArray, Int64Array, Int8Array, RecordBatch, StringArray};
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

pub(crate) fn save(transcripts: &[Transcript], path: &Path) -> Result<(), Box<dyn Error>> {
    let mut tid = Vec::with_capacity(transcripts.len());
    let (mut chrom, mut gid, mut sym, mut bt, mut model) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut start, mut end) = (Vec::new(), Vec::new());
    let mut strand = Vec::new();
    let mut canon = Vec::new();
    for t in transcripts {
        tid.push(t.stable_id.to_string());
        chrom.push(t.chromosome.to_string());
        start.push(t.start as i64);
        end.push(t.end as i64);
        strand.push(t.strand.as_int());
        gid.push(t.gene.stable_id.to_string());
        sym.push(t.gene.symbol.as_deref().unwrap_or("").to_string());
        bt.push(t.biotype.to_string());
        canon.push(t.canonical);
        model.push(serde_json::to_string(t)?);
    }
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
        Field::new("model", DataType::Utf8, false),
    ]));
    let batch = RecordBatch::try_new(
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
            Arc::new(StringArray::from(model)),
        ],
    )?;
    // zstd: the verbose JSON `model` column compresses ~10×.
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::try_new(3)?))
        .build();
    let mut w = ArrowWriter::try_new(File::create(path)?, schema, Some(props))?;
    w.write(&batch)?;
    w.close()?;
    Ok(())
}

pub(crate) fn load(path: &Path) -> Result<Vec<Transcript>, Box<dyn Error>> {
    let reader = ParquetRecordBatchReaderBuilder::try_new(File::open(path)?)?.build()?;
    let mut out = Vec::new();
    for batch in reader {
        let batch = batch?;
        let model = batch
            .column_by_name("model")
            .ok_or("transcript cache: missing 'model' column")?
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or("transcript cache: 'model' is not VARCHAR")?;
        for i in 0..model.len() {
            out.push(serde_json::from_str(model.value(i))?);
        }
    }
    Ok(out)
}
