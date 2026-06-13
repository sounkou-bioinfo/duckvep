//! `vep_annotate(vcf, gff3 := ..., fasta := ..., distance := 5000)` — end-to-end
//! variant effect prediction as a DuckDB table function.
//!
//! This is the DuckDB frontend: it parses parameters and marshals the engine's
//! [`AnnotatedRow`]s into output vectors. All compute lives in [`super::engine`]
//! (DuckDB-free), so the CLI/C-API frontends share it without linking libduckdb.
//! One output row per (variant × transcript × allele).

use crate::vec_util::fill_string_list;
use crate::vep::engine::{annotate, AnnotatedRow, DEFAULT_DISTANCE};
use duckdb::core::{DataChunkHandle, Inserter, LogicalTypeHandle, LogicalTypeId};
use duckdb::vtab::{BindInfo, InitInfo, TableFunctionInfo, VTab};
use std::error::Error;
use std::sync::atomic::{AtomicUsize, Ordering};

const VECTOR_SIZE: usize = 2048;

pub struct VepAnnotate;

pub struct VepBind {
    rows: Vec<AnnotatedRow>,
}

pub struct VepInit {
    cursor: AtomicUsize,
}

impl VTab for VepAnnotate {
    type BindData = VepBind;
    type InitData = VepInit;

    fn bind(bind: &BindInfo) -> Result<Self::BindData, Box<dyn Error>> {
        let varchar = || LogicalTypeHandle::from(LogicalTypeId::Varchar);
        let bigint = || LogicalTypeHandle::from(LogicalTypeId::Bigint);
        bind.add_result_column("chrom", varchar());
        bind.add_result_column("pos", bigint());
        bind.add_result_column("ref", varchar());
        bind.add_result_column("alt", varchar());
        bind.add_result_column("gene_id", varchar());
        bind.add_result_column("gene_symbol", varchar());
        bind.add_result_column("transcript_id", varchar());
        bind.add_result_column("biotype", varchar());
        bind.add_result_column("canonical", LogicalTypeHandle::from(LogicalTypeId::Boolean));
        bind.add_result_column("consequence", LogicalTypeHandle::list(&varchar()));
        bind.add_result_column("impact", varchar());
        bind.add_result_column("amino_acids", varchar());
        bind.add_result_column("codons", varchar());
        bind.add_result_column("protein_pos", bigint());

        let vcf_path = bind.get_parameter(0).to_string();
        let gff3 = bind
            .get_named_parameter("gff3")
            .map(|v| v.to_string())
            .filter(|s| !s.is_empty())
            .ok_or("vep_annotate requires gff3 := '<path>'")?;
        let fasta = bind
            .get_named_parameter("fasta")
            .map(|v| v.to_string())
            .filter(|s| !s.is_empty());
        let distance = bind
            .get_named_parameter("distance")
            .and_then(|v| v.to_string().parse::<u64>().ok())
            .unwrap_or(DEFAULT_DISTANCE);

        let rows = annotate(&vcf_path, &gff3, fasta.as_deref(), distance)?;
        Ok(VepBind { rows })
    }

    fn init(_: &InitInfo) -> Result<Self::InitData, Box<dyn Error>> {
        Ok(VepInit {
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

        macro_rules! varchar_col {
            ($idx:expr, $field:ident) => {{
                let v = output.flat_vector($idx);
                for (i, r) in rows.iter().enumerate() {
                    v.insert(i, r.$field.as_str());
                }
            }};
        }
        varchar_col!(0, chrom);
        {
            let mut v = output.flat_vector(1);
            let s = unsafe { v.as_mut_slice::<i64>() };
            for (i, r) in rows.iter().enumerate() {
                s[i] = r.pos;
            }
        }
        varchar_col!(2, reference);
        varchar_col!(3, alt);
        varchar_col!(4, gene_id);
        varchar_col!(5, gene_symbol);
        varchar_col!(6, transcript_id);
        varchar_col!(7, biotype);
        {
            let mut v = output.flat_vector(8);
            let s = unsafe { v.as_mut_slice::<bool>() };
            for (i, r) in rows.iter().enumerate() {
                s[i] = r.canonical;
            }
        }
        fill_string_list(output, 9, rows, |r| &r.consequence);
        varchar_col!(10, impact);
        varchar_col!(11, amino_acids);
        varchar_col!(12, codons);
        {
            let mut v = output.flat_vector(13);
            {
                let s = unsafe { v.as_mut_slice::<i64>() };
                for (i, r) in rows.iter().enumerate() {
                    s[i] = r.protein_pos.unwrap_or(0);
                }
            }
            for (i, r) in rows.iter().enumerate() {
                if r.protein_pos.is_none() {
                    v.set_null(i);
                }
            }
        }

        init.cursor.store(start + n, Ordering::Relaxed);
        output.set_len(n);
        Ok(())
    }

    fn parameters() -> Option<Vec<LogicalTypeHandle>> {
        Some(vec![LogicalTypeHandle::from(LogicalTypeId::Varchar)])
    }

    fn named_parameters() -> Option<Vec<(String, LogicalTypeHandle)>> {
        let varchar = || LogicalTypeHandle::from(LogicalTypeId::Varchar);
        Some(vec![
            ("gff3".to_string(), varchar()),
            ("fasta".to_string(), varchar()),
            (
                "distance".to_string(),
                LogicalTypeHandle::from(LogicalTypeId::Bigint),
            ),
        ])
    }
}
