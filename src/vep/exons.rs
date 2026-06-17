//! `vep_exons()` — the loaded gene model's EXONS as a first-class SQL relation
//! (one row per exon, keyed by the transcript ordinal `tx_idx`). This is the SoA
//! interval data the structural-consequence sweep classifies against (exon vs
//! intron vs splice vs UTR vs CDS), and the start of killing the opaque `model`
//! BLOB by exposing the model's nested structure as typed columns. The `cds_start`
//! /`cds_end` are the transcript's genomic coding bounds (0 for non-coding), so a
//! single per-exon row carries everything the region kernel needs.

use crate::vep::consequence::loaded_engine;
use duckdb::core::{DataChunkHandle, Inserter, LogicalTypeHandle, LogicalTypeId};
use duckdb::vtab::{BindInfo, InitInfo, TableFunctionInfo, VTab};
use std::error::Error;
use std::sync::atomic::{AtomicUsize, Ordering};

const VECTOR_SIZE: usize = 2048;

pub struct VepExons;

struct ExonRow {
    tx_idx: i64,
    chrom: String,
    exon_start: i64,
    exon_end: i64,
    strand: i64,
    cds_start: i64,
    cds_end: i64,
}

pub struct ExBind {
    rows: Vec<ExonRow>,
}

pub struct ExInit {
    cursor: AtomicUsize,
}

impl VTab for VepExons {
    type BindData = ExBind;
    type InitData = ExInit;

    fn bind(bind: &BindInfo) -> Result<Self::BindData, Box<dyn Error>> {
        let varchar = || LogicalTypeHandle::from(LogicalTypeId::Varchar);
        let bigint = || LogicalTypeHandle::from(LogicalTypeId::Bigint);
        bind.add_result_column("tx_idx", bigint());
        bind.add_result_column("chrom", varchar());
        bind.add_result_column("exon_start", bigint());
        bind.add_result_column("exon_end", bigint());
        bind.add_result_column("strand", bigint());
        bind.add_result_column("cds_start", bigint());
        bind.add_result_column("cds_end", bigint());

        let ctx = loaded_engine().ok_or("vep_exons: call vep_load_cache(gff3, fasta) first")?;
        let mut rows = Vec::new();
        // transcript_iter yields transcripts in tx_idx (ordinal) order.
        for (tx_idx, t) in ctx.transcript_iter().enumerate() {
            let cds_start = t.coding_region_start.unwrap_or(0) as i64;
            let cds_end = t.coding_region_end.unwrap_or(0) as i64;
            let strand = t.strand.as_int() as i64;
            for e in &t.exons {
                rows.push(ExonRow {
                    tx_idx: tx_idx as i64,
                    chrom: t.chromosome.to_string(),
                    exon_start: e.start as i64,
                    exon_end: e.end as i64,
                    strand,
                    cds_start,
                    cds_end,
                });
            }
        }
        Ok(ExBind { rows })
    }

    fn init(_: &InitInfo) -> Result<Self::InitData, Box<dyn Error>> {
        Ok(ExInit {
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

        macro_rules! ic {
            ($i:expr, $f:ident) => {{
                let mut v = output.flat_vector($i);
                let s = unsafe { v.as_mut_slice::<i64>() };
                for (i, r) in rows.iter().enumerate() {
                    s[i] = r.$f;
                }
            }};
        }
        ic!(0, tx_idx);
        {
            let v = output.flat_vector(1);
            for (i, r) in rows.iter().enumerate() {
                v.insert(i, r.chrom.as_str());
            }
        }
        ic!(2, exon_start);
        ic!(3, exon_end);
        ic!(4, strand);
        ic!(5, cds_start);
        ic!(6, cds_end);

        init.cursor.store(start + n, Ordering::Relaxed);
        output.set_len(n);
        Ok(())
    }

    fn parameters() -> Option<Vec<LogicalTypeHandle>> {
        Some(vec![])
    }
}
