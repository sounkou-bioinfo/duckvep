//! `vep_cds_seq()` — the translateable (CDS) sequence per coding transcript, as a
//! SQL relation `(tx_idx, cds_seq)`. Only populated when `vep_load_cache` was
//! given a FASTA (the engine builds `translateable_seq` from it). This is the
//! packed-sequence half of the SoA cache: the codon kernel reads it (nt2/nt4
//! packed) to classify synonymous/missense/stop on the CDS bucket only.

use crate::vep::consequence::loaded_engine;
use duckdb::core::{DataChunkHandle, Inserter, LogicalTypeHandle, LogicalTypeId};
use duckdb::vtab::{BindInfo, InitInfo, TableFunctionInfo, VTab};
use std::error::Error;
use std::sync::atomic::{AtomicUsize, Ordering};

const VECTOR_SIZE: usize = 2048;

pub struct VepCdsSeq;

struct SeqRow {
    tx_idx: i64,
    cds_seq: String,
}

pub struct SeqBind {
    rows: Vec<SeqRow>,
}
pub struct SeqInit {
    cursor: AtomicUsize,
}

impl VTab for VepCdsSeq {
    type BindData = SeqBind;
    type InitData = SeqInit;

    fn bind(bind: &BindInfo) -> Result<Self::BindData, Box<dyn Error>> {
        bind.add_result_column("tx_idx", LogicalTypeHandle::from(LogicalTypeId::Bigint));
        bind.add_result_column("cds_seq", LogicalTypeHandle::from(LogicalTypeId::Varchar));
        let ctx = loaded_engine().ok_or("vep_cds_seq: call vep_load_cache(gff3, fasta) first")?;
        let mut rows = Vec::new();
        for (tx_idx, t) in ctx.transcript_iter().enumerate() {
            if let Some(seq) = &t.translateable_seq {
                rows.push(SeqRow {
                    tx_idx: tx_idx as i64,
                    cds_seq: seq.clone(),
                });
            }
        }
        Ok(SeqBind { rows })
    }

    fn init(_: &InitInfo) -> Result<Self::InitData, Box<dyn Error>> {
        Ok(SeqInit {
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
        let n = (bind.rows.len() - start).min(VECTOR_SIZE);
        let rows = &bind.rows[start..start + n];
        {
            let mut v = output.flat_vector(0);
            let s = unsafe { v.as_mut_slice::<i64>() };
            for (i, r) in rows.iter().enumerate() {
                s[i] = r.tx_idx;
            }
        }
        {
            let v = output.flat_vector(1);
            for (i, r) in rows.iter().enumerate() {
                v.insert(i, r.cds_seq.as_str());
            }
        }
        init.cursor.store(start + n, Ordering::Relaxed);
        output.set_len(n);
        Ok(())
    }

    fn parameters() -> Option<Vec<LogicalTypeHandle>> {
        Some(vec![])
    }
}
