//! Scan-driven scalar UDFs — the idiomatic, performant path (DESIGN.md §3.2).
//!
//! `vep_load_cache(gff3, fasta)` builds the engine once and stores it in a global
//! (set once, read by all worker threads). `vep_consequence(chrom,pos,ref,alt)`
//! is then a pure scalar over rows DuckDB scans (`read_vcf`, `read_parquet`, any
//! source) — DuckDB owns the scan; the kernel only computes per variant. Returns
//! a JSON array of per-transcript consequences; `UNNEST(from_json(...))` it.

use crate::vep::engine::{build_context, AnnotatedRow, EngineContext, DEFAULT_DISTANCE};
use duckdb::arrow::array::{Array, ArrayRef, AsArray, RecordBatch, StringArray};
use duckdb::arrow::datatypes::{DataType, Int64Type};
use duckdb::vscalar::{ArrowFunctionSignature, VArrowScalar};
use std::error::Error;
use std::sync::{Arc, OnceLock, RwLock};

/// Load-once engine, shared across DuckDB worker threads.
static CACHE: OnceLock<RwLock<Option<Arc<EngineContext>>>> = OnceLock::new();
fn cache() -> &'static RwLock<Option<Arc<EngineContext>>> {
    CACHE.get_or_init(|| RwLock::new(None))
}

fn rows_to_json(rows: &[AnnotatedRow]) -> String {
    let arr: Vec<serde_json::Value> = rows
        .iter()
        .map(|r| {
            serde_json::json!({
                "transcript_id": r.transcript_id,
                "gene_id": r.gene_id,
                "gene_symbol": r.gene_symbol,
                "biotype": r.biotype,
                "canonical": r.canonical,
                "consequence": r.consequence,
                "impact": r.impact,
                "amino_acids": r.amino_acids,
                "codons": r.codons,
                "protein_pos": r.protein_pos,
                "alt": r.alt,
            })
        })
        .collect();
    serde_json::Value::Array(arr).to_string()
}

/// `vep_load_cache(gff3, fasta) -> VARCHAR` — build + store the engine.
pub struct VepLoadCache;

impl VArrowScalar for VepLoadCache {
    type State = ();

    fn invoke(_: &Self::State, input: RecordBatch) -> Result<ArrayRef, Box<dyn Error>> {
        let gff3 = input.column(0).as_string::<i32>();
        let fasta = input.column(1).as_string::<i32>();
        let gff3v = gff3.value(0);
        let fastav = if fasta.is_null(0) || fasta.value(0).is_empty() {
            None
        } else {
            Some(fasta.value(0))
        };
        let ctx = build_context(gff3v, fastav, DEFAULT_DISTANCE)?;
        *cache().write().map_err(|_| "cache lock poisoned")? = Some(Arc::new(ctx));
        let n = input.num_rows();
        Ok(Arc::new(StringArray::from(vec!["loaded"; n])))
    }

    fn signatures() -> Vec<ArrowFunctionSignature> {
        vec![ArrowFunctionSignature::exact(
            vec![DataType::Utf8, DataType::Utf8],
            DataType::Utf8,
        )]
    }
}

/// `vep_consequence(chrom, pos, ref, alt) -> VARCHAR(json)` — per-variant kernel.
pub struct VepConsequence;

impl VArrowScalar for VepConsequence {
    type State = ();

    fn invoke(_: &Self::State, input: RecordBatch) -> Result<ArrayRef, Box<dyn Error>> {
        let guard = cache().read().map_err(|_| "cache lock poisoned")?;
        let ctx = guard
            .as_ref()
            .ok_or("vep_consequence: call vep_load_cache(gff3, fasta) first")?;

        let chrom = input.column(0).as_string::<i32>();
        let pos = input.column(1).as_primitive::<Int64Type>();
        let refa = input.column(2).as_string::<i32>();
        let alta = input.column(3).as_string::<i32>();

        let out: Vec<String> = (0..input.num_rows())
            .map(|i| {
                let rows = ctx.annotate_variant(
                    chrom.value(i),
                    pos.value(i) as u64,
                    refa.value(i),
                    alta.value(i),
                );
                rows_to_json(&rows)
            })
            .collect();
        Ok(Arc::new(StringArray::from(out)))
    }

    fn signatures() -> Vec<ArrowFunctionSignature> {
        vec![ArrowFunctionSignature::exact(
            vec![
                DataType::Utf8,
                DataType::Int64,
                DataType::Utf8,
                DataType::Utf8,
            ],
            DataType::Utf8,
        )]
    }
}
