//! gen-hard-variants — generate a synthetic corpus of PATHOLOGICAL hard variants
//! for the consequence classes we target, placed at *exact* offsets relative to
//! transcript features (which random sampling almost never hits). Pure Rust:
//! parquet for coordinates, `fastvep_cache::FastaReader` (noodles-backed) for ref
//! bases — no bash/R/vcfppR/Python.
//!
//! Usage: gen-hard-variants <exons.parquet> <transcripts.parquet> <fasta> <chrom> <out.vcf>

use anyhow::{Context, Result};
use arrow::array::{Array, BooleanArray, Int64Array, RecordBatch, StringArray};
use fastvep_cache::fasta::FastaReader;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::Write;

type Rec = (i64, String, String, String); // (pos, id/tag, ref, alt)

fn cstr<'a>(b: &'a RecordBatch, n: &str) -> &'a StringArray {
    b.column(b.schema().index_of(n).unwrap())
        .as_any()
        .downcast_ref()
        .unwrap()
}
fn ci64<'a>(b: &'a RecordBatch, n: &str) -> &'a Int64Array {
    b.column(b.schema().index_of(n).unwrap())
        .as_any()
        .downcast_ref()
        .unwrap()
}
fn cbool<'a>(b: &'a RecordBatch, n: &str) -> &'a BooleanArray {
    b.column(b.schema().index_of(n).unwrap())
        .as_any()
        .downcast_ref()
        .unwrap()
}

fn seq(fa: &FastaReader, chrom: &str, s: i64, e: i64) -> Option<String> {
    if s < 1 {
        return None;
    }
    let v = fa.fetch(chrom, s as u64, e as u64).ok()?;
    let s = String::from_utf8_lossy(&v).into_owned();
    if s.is_empty() || s.contains('N') {
        None
    } else {
        Some(s)
    }
}

/// VCF anchor-base deletion of genomic [s, e] inclusive (ref = anchor..e, alt = anchor).
fn del(fa: &FastaReader, chrom: &str, r: &mut BTreeSet<Rec>, s: i64, e: i64, tag: &str) {
    if s - 1 < 2 {
        return;
    }
    if let (Some(rf), Some(an)) = (seq(fa, chrom, s - 1, e), seq(fa, chrom, s - 1, s - 1)) {
        r.insert((s - 1, tag.into(), rf, an));
    }
}
/// VCF insertion of `b` after position p (ref = base@p, alt = base@p + b).
fn ins(fa: &FastaReader, chrom: &str, r: &mut BTreeSet<Rec>, p: i64, b: &str, tag: &str) {
    if let Some(an) = seq(fa, chrom, p, p) {
        r.insert((p, tag.into(), an.clone(), format!("{an}{b}")));
    }
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() < 5 {
        anyhow::bail!("usage: gen-hard-variants <exons.parquet> <transcripts.parquet> <fasta> <chrom> <out.vcf>");
    }
    let (exons_p, tx_p, fasta_p, chrom, out_p) = (&a[0], &a[1], &a[2], a[3].clone(), &a[4]);

    // Coding, MANE/canonical transcripts on this contig.
    let mut keep: HashSet<String> = HashSet::new();
    for batch in ParquetRecordBatchReaderBuilder::try_new(File::open(tx_p)?)?.build()? {
        let b = batch?;
        let (tid, ch, bt) = (
            cstr(&b, "transcript_id"),
            cstr(&b, "chrom"),
            cstr(&b, "biotype"),
        );
        let (mane, canon) = (cbool(&b, "mane_select"), cbool(&b, "canonical"));
        for i in 0..b.num_rows() {
            if ch.value(i) == chrom
                && bt.value(i) == "protein_coding"
                && (mane.value(i) || canon.value(i))
            {
                keep.insert(tid.value(i).to_string());
            }
        }
    }

    // Exons (genomic) per kept transcript.
    let mut exons: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    for batch in ParquetRecordBatchReaderBuilder::try_new(File::open(exons_p)?)?.build()? {
        let b = batch?;
        let (tid, s, e) = (
            cstr(&b, "transcript_id"),
            ci64(&b, "start"),
            ci64(&b, "end"),
        );
        for i in 0..b.num_rows() {
            let t = tid.value(i);
            if keep.contains(t) {
                exons
                    .entry(t.to_string())
                    .or_default()
                    .push((s.value(i), e.value(i)));
            }
        }
    }

    let fa = FastaReader::from_reader(File::open(fasta_p).context("open fasta")?)?;
    let mut recs: BTreeSet<Rec> = BTreeSet::new();
    for ex in exons.values() {
        let mut ex = ex.clone();
        ex.sort();
        for w in ex.windows(2) {
            let (istart, iend) = (w[0].1 + 1, w[1].0 - 1); // intron (genomic)
            if iend - istart < 40 {
                continue;
            }
            del(
                &fa,
                &chrom,
                &mut recs,
                istart - 1,
                istart + 1,
                "del_span_donor5p",
            );
            del(
                &fa,
                &chrom,
                &mut recs,
                istart + 3,
                istart + 5,
                "del_donor_region",
            );
            del(
                &fa,
                &chrom,
                &mut recs,
                istart + 4,
                istart + 4,
                "del_5th_base",
            );
            del(
                &fa,
                &chrom,
                &mut recs,
                iend - 1,
                iend + 1,
                "del_span_acceptor3p",
            );
            del(
                &fa,
                &chrom,
                &mut recs,
                iend - 12,
                iend - 8,
                "del_polypyrimidine",
            );
            del(
                &fa,
                &chrom,
                &mut recs,
                istart + 2,
                istart + 7,
                "del_splice_region",
            );
            ins(&fa, &chrom, &mut recs, istart + 1, "AT", "ins_donor");
            ins(
                &fa,
                &chrom,
                &mut recs,
                iend - 3,
                "GC",
                "ins_acceptor_region",
            );
            del(
                &fa,
                &chrom,
                &mut recs,
                istart - 1,
                istart - 1,
                "del_fs1_exonend",
            );
            del(
                &fa,
                &chrom,
                &mut recs,
                istart - 4,
                istart - 2,
                "del_inframe3_exonend",
            );
            if let Some(r) = seq(&fa, &chrom, iend + 6, iend + 8) {
                if r.len() == 3 {
                    let comp: String = r
                        .chars()
                        .map(|c| match c {
                            'A' => 'C',
                            'C' => 'G',
                            'G' => 'T',
                            'T' => 'A',
                            o => o,
                        })
                        .collect();
                    recs.insert((iend + 6, "mnv3_exon".into(), r, comp));
                }
            }
        }
    }

    let mut out = File::create(out_p)?;
    writeln!(
        out,
        "##fileformat=VCFv4.2\n##source=tools/gen-hard-variants\n##contig=<ID={chrom}>"
    )?;
    writeln!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    for (pos, tag, r, alt) in &recs {
        writeln!(out, "{chrom}\t{pos}\t{tag}\t{r}\t{alt}\t.\t.\t.")?;
    }
    eprintln!(
        "generated {} hard variants across {} transcripts",
        recs.len(),
        exons.len()
    );
    Ok(())
}
