//! `vep_transcripts()` — the loaded gene model as a first-class SQL relation.
//!
//! Reads the resident `EngineContext` (set by `vep_load_cache`), so the
//! range-join annotation path can do `JOIN vep_transcripts() t ON …` and read
//! transcripts straight from the in-memory model — one source of truth, no
//! re-reading the parquet cache and no double-resident gene model. It is also the
//! building block for relational gene-model queries (panels, biotype filters,
//! the supplementary-annotation joins). See docs/kernel-algorithm.md §6 / §8.

use crate::vep::consequence::loaded_engine;
use duckdb::core::{DataChunkHandle, Inserter, LogicalTypeHandle, LogicalTypeId};
use duckdb::vtab::{BindInfo, InitInfo, TableFunctionInfo, VTab};
use std::error::Error;
use std::sync::atomic::{AtomicUsize, Ordering};

const VECTOR_SIZE: usize = 2048;

pub struct VepTranscripts;

struct TxRow {
    transcript_id: String,
    chrom: String,
    start: i64,
    end: i64,
    strand: i64,
    gene_id: String,
    gene_symbol: String,
    biotype: String,
    canonical: bool,
    coding: bool,
}

pub struct TxBind {
    rows: Vec<TxRow>,
}

pub struct TxInit {
    cursor: AtomicUsize,
}

impl VTab for VepTranscripts {
    type BindData = TxBind;
    type InitData = TxInit;

    fn bind(bind: &BindInfo) -> Result<Self::BindData, Box<dyn Error>> {
        let varchar = || LogicalTypeHandle::from(LogicalTypeId::Varchar);
        let bigint = || LogicalTypeHandle::from(LogicalTypeId::Bigint);
        let boolean = || LogicalTypeHandle::from(LogicalTypeId::Boolean);
        bind.add_result_column("transcript_id", varchar());
        bind.add_result_column("chrom", varchar());
        bind.add_result_column("start", bigint());
        bind.add_result_column("end_pos", bigint());
        bind.add_result_column("strand", bigint());
        bind.add_result_column("gene_id", varchar());
        bind.add_result_column("gene_symbol", varchar());
        bind.add_result_column("biotype", varchar());
        bind.add_result_column("canonical", boolean());
        bind.add_result_column("coding", boolean());

        let ctx = loaded_engine().ok_or("vep_transcripts: call vep_load_cache(gff3, fasta) first")?;
        let rows = ctx
            .transcript_iter()
            .map(|t| TxRow {
                transcript_id: t.stable_id.to_string(),
                chrom: t.chromosome.to_string(),
                start: t.start as i64,
                end: t.end as i64,
                strand: t.strand.as_int() as i64,
                gene_id: t.gene.stable_id.to_string(),
                gene_symbol: t.gene.symbol.as_deref().unwrap_or("").to_string(),
                biotype: t.biotype.to_string(),
                canonical: t.canonical,
                coding: t.is_coding(),
            })
            .collect();
        Ok(TxBind { rows })
    }

    fn init(_: &InitInfo) -> Result<Self::InitData, Box<dyn Error>> {
        Ok(TxInit {
            cursor: AtomicUsize::new(0),
        })
    }

    fn func(
        func: &TableFunctionInfo<Self>,
        output: &mut DataChunkHandle,
    ) -> Result<(), Box<dyn Error>> {
        let bind = func.get_bind_data();
        let init = func.get_init_data();
        let start = init.cursor.load(Ordering::Relaxed);
        let total = bind.rows.len();
        let n = (total - start).min(VECTOR_SIZE);
        let rows = &bind.rows[start..start + n];

        macro_rules! vc {
            ($i:expr, $f:ident) => {{
                let v = output.flat_vector($i);
                for (i, r) in rows.iter().enumerate() {
                    v.insert(i, r.$f.as_str());
                }
            }};
        }
        macro_rules! ic {
            ($i:expr, $f:ident) => {{
                let mut v = output.flat_vector($i);
                let s = unsafe { v.as_mut_slice::<i64>() };
                for (i, r) in rows.iter().enumerate() {
                    s[i] = r.$f;
                }
            }};
        }
        macro_rules! bc {
            ($i:expr, $f:ident) => {{
                let mut v = output.flat_vector($i);
                let s = unsafe { v.as_mut_slice::<bool>() };
                for (i, r) in rows.iter().enumerate() {
                    s[i] = r.$f;
                }
            }};
        }
        vc!(0, transcript_id);
        vc!(1, chrom);
        ic!(2, start);
        ic!(3, end);
        ic!(4, strand);
        vc!(5, gene_id);
        vc!(6, gene_symbol);
        vc!(7, biotype);
        bc!(8, canonical);
        bc!(9, coding);

        init.cursor.store(start + n, Ordering::Relaxed);
        output.set_len(n);
        Ok(())
    }

    fn parameters() -> Option<Vec<LogicalTypeHandle>> {
        Some(vec![])
    }
}
