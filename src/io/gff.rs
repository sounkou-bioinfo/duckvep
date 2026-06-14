//! `read_gff_transcripts(gff3)` — the GFF importer: parse a GFF3 gene model and
//! emit one flat row per transcript, matching the cache `transcripts` schema
//! (the same shape the Ensembl MySQL sync produces). Build the runtime cache as
//! a native DuckDB table:
//!   `CREATE TABLE transcripts AS SELECT * FROM read_gff_transcripts('x.gff3');`
//! The cache is a `.duckdb` file read with DuckDB core only at runtime — no
//! Parquet, no MySQL, no extensions (docs/DESIGN.md §1c, §5).
//!
//! Relational, not nested: exons come from a companion `read_gff_exons` so the
//! cache mirrors Ensembl's schema and the importer stays flat.

use crate::vec_util::fill_string_list;
use duckdb::core::{DataChunkHandle, Inserter, LogicalTypeHandle, LogicalTypeId};
use duckdb::vtab::{BindInfo, InitInfo, TableFunctionInfo, VTab};
use fastvep_cache::gff::parse_gff3;
use fastvep_genome::Transcript;
use std::error::Error;
use std::fs::File;
use std::sync::atomic::{AtomicUsize, Ordering};

const VECTOR_SIZE: usize = 2048;

struct TxRow {
    transcript_id: String,
    chrom: String,
    start: i64,
    end: i64,
    strand: i64,
    biotype: String,
    gene_id: String,
    gene_symbol: String,
    canonical: bool,
    coding: bool,
    tsl: Option<i64>,
    appris: String,
    flags: Vec<String>,
}

pub struct ReadGffTranscripts;

pub struct GffBind {
    rows: Vec<TxRow>,
}

pub struct GffInit {
    cursor: AtomicUsize,
}

fn load_transcripts(path: &str) -> Result<Vec<Transcript>, Box<dyn Error>> {
    let f = File::open(path)?;
    let trs = if path.ends_with(".gz") || path.ends_with(".bgz") {
        parse_gff3(flate2::read::MultiGzDecoder::new(f))?
    } else {
        parse_gff3(f)?
    };
    Ok(trs)
}

impl VTab for ReadGffTranscripts {
    type BindData = GffBind;
    type InitData = GffInit;

    fn bind(bind: &BindInfo) -> Result<Self::BindData, Box<dyn Error>> {
        let varchar = || LogicalTypeHandle::from(LogicalTypeId::Varchar);
        let bigint = || LogicalTypeHandle::from(LogicalTypeId::Bigint);
        bind.add_result_column("transcript_id", varchar());
        bind.add_result_column("chrom", varchar());
        bind.add_result_column("start", bigint());
        bind.add_result_column("end_pos", bigint());
        bind.add_result_column("strand", bigint());
        bind.add_result_column("biotype", varchar());
        bind.add_result_column("gene_id", varchar());
        bind.add_result_column("gene_symbol", varchar());
        bind.add_result_column("canonical", LogicalTypeHandle::from(LogicalTypeId::Boolean));
        bind.add_result_column("coding", LogicalTypeHandle::from(LogicalTypeId::Boolean));
        bind.add_result_column("tsl", bigint());
        bind.add_result_column("appris", varchar());
        bind.add_result_column("flags", LogicalTypeHandle::list(&varchar()));

        let path = bind.get_parameter(0).to_string();
        let rows = load_transcripts(&path)?
            .iter()
            .map(|t| TxRow {
                transcript_id: t.stable_id.to_string(),
                chrom: t.chromosome.to_string(),
                start: t.start as i64,
                end: t.end as i64,
                strand: t.strand.as_int() as i64,
                biotype: t.biotype.to_string(),
                gene_id: t.gene.stable_id.to_string(),
                gene_symbol: t.gene.symbol.as_deref().unwrap_or("").to_string(),
                canonical: t.canonical,
                coding: t.is_coding(),
                tsl: t.tsl.map(|v| v as i64),
                appris: t.appris.clone().unwrap_or_default(),
                flags: t.flags.clone(),
            })
            .collect();
        Ok(GffBind { rows })
    }

    fn init(_: &InitInfo) -> Result<Self::InitData, Box<dyn Error>> {
        Ok(GffInit {
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
        vc!(5, biotype);
        vc!(6, gene_id);
        vc!(7, gene_symbol);
        bc!(8, canonical);
        bc!(9, coding);
        {
            let mut v = output.flat_vector(10);
            {
                let s = unsafe { v.as_mut_slice::<i64>() };
                for (i, r) in rows.iter().enumerate() {
                    s[i] = r.tsl.unwrap_or(0);
                }
            }
            for (i, r) in rows.iter().enumerate() {
                if r.tsl.is_none() {
                    v.set_null(i);
                }
            }
        }
        vc!(11, appris);
        fill_string_list(output, 12, rows, |r| &r.flags);

        init.cursor.store(start + n, Ordering::Relaxed);
        output.set_len(n);
        Ok(())
    }

    fn parameters() -> Option<Vec<LogicalTypeHandle>> {
        Some(vec![LogicalTypeHandle::from(LogicalTypeId::Varchar)])
    }
}
