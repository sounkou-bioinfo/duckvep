//! Prototype D (docs/kernel-algorithm.md): the codon kernel on the CDS bucket.
//! Packs all CDS sequences (457 M bases) into an nt2 pool (2 bits/base → ~114 MB),
//! then runs 256_000 codon classifications — the measured CDS-bucket size — with
//! RANDOM access (the realistic scattered-CDS-variant cache pattern): fetch the
//! reference codon, apply an alt base, translate ref+alt via a 64-entry LUT,
//! classify synonymous/missense/stop. This closes the cost picture: is the
//! sequence-dependent kernel (the only genuinely expensive part) cheap on the
//! small CDS bucket?
//!
//!   cargo run --release --example codon_bench

use arrow::array::{Array, StringArray, UInt32Array};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::fs::File;
use std::time::Instant;

// Standard genetic code, indexed by codon = b0*16 + b1*4 + b2 with A=0,C=1,G=2,T=3.
const AA: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

#[inline(always)]
fn code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

fn main() {
    const OPS: u64 = 256_000; // measured CDS-bucket size
    let t0 = Instant::now();

    // Build the nt2 packed pool (4 bases/byte) + per-CDS (offset_in_bases, len).
    let r = ParquetRecordBatchReaderBuilder::try_new(File::open("/tmp/cds_seq.parquet").unwrap())
        .unwrap()
        .build()
        .unwrap();
    let mut pool: Vec<u8> = Vec::new();
    let mut off: Vec<u64> = Vec::new();
    let mut len: Vec<u32> = Vec::new();
    let mut cur_base: u64 = 0;
    let mut acc: u8 = 0;
    let mut nb: u8 = 0;
    for b in r {
        let b = b.unwrap();
        let seqs = b.column(1).as_any().downcast_ref::<StringArray>().unwrap();
        for i in 0..seqs.len() {
            let s = seqs.value(i).as_bytes();
            off.push(cur_base);
            len.push(s.len() as u32);
            for &ch in s {
                acc |= code(ch) << (2 * nb);
                nb += 1;
                if nb == 4 {
                    pool.push(acc);
                    acc = 0;
                    nb = 0;
                }
                cur_base += 1;
            }
        }
    }
    if nb > 0 {
        pool.push(acc);
    }
    let pack_s = t0.elapsed().as_secs_f64();
    eprintln!(
        "packed {} CDS, {:.1} M bases -> {:.1} MB nt2 pool in {:.2}s",
        off.len(),
        cur_base as f64 / 1e6,
        pool.len() as f64 / 1e6,
        pack_s
    );

    #[inline(always)]
    fn base_at(pool: &[u8], idx: u64) -> u8 {
        (pool[(idx >> 2) as usize] >> (2 * (idx & 3))) & 3
    }

    // 256k codon classifications, random CDS + random in-frame codon (scattered
    // access = realistic). LCG keeps it deterministic and allocation-free.
    let n_cds = off.len() as u64;
    let mut rng: u64 = 0x9E3779B97F4A7C15;
    let mut rnd = || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };

    let t1 = Instant::now();
    let (mut syn, mut mis, mut stop) = (0u64, 0u64, 0u64);
    let mut sink: u64 = 0;
    for _ in 0..OPS {
        let ti = (rnd() % n_cds) as usize;
        let l = len[ti];
        if l < 3 {
            continue;
        }
        let n_codons = (l / 3) as u64;
        let ci = rnd() % n_codons; // codon index
        let base0 = off[ti] + ci * 3;
        let b0 = base_at(&pool, base0);
        let b1 = base_at(&pool, base0 + 1);
        let b2 = base_at(&pool, base0 + 2);
        let ref_codon = (b0 << 4) | (b1 << 2) | b2;
        // alt: flip one base of the codon to a different value
        let which = (rnd() % 3) as u8;
        let newb = ((rnd() as u8) & 3).wrapping_add(1) & 3; // any base
        let (mut a0, mut a1, mut a2) = (b0, b1, b2);
        match which {
            0 => a0 = newb,
            1 => a1 = newb,
            _ => a2 = newb,
        }
        let alt_codon = (a0 << 4) | (a1 << 2) | a2;
        let ref_aa = AA[ref_codon as usize];
        let alt_aa = AA[alt_codon as usize];
        sink ^= (ref_aa as u64) | ((alt_aa as u64) << 8);
        if ref_aa == alt_aa {
            syn += 1;
        } else if alt_aa == b'*' {
            stop += 1;
        } else {
            mis += 1;
        }
    }
    let kern_s = t1.elapsed().as_secs_f64();

    eprintln!("--- codon kernel: {OPS} CDS-bucket classifications (single thread) ---");
    eprintln!("synonymous : {syn}");
    eprintln!("missense   : {mis}");
    eprintln!("stop       : {stop}  (sink={sink})");
    eprintln!(
        "wall       : {:.4}s   ({:.0}M codons/s)",
        kern_s,
        OPS as f64 / 1e6 / kern_s
    );
    eprintln!("=> the entire sequence-dependent kernel for HG002 WGS is this bucket.");
}
