//! Region parallelization of the sweep+classify (prototype C), to answer "these
//! were all single-threaded — how does it scale?". Each chromosome is an
//! independent work unit (the haplotype-safe tile); threads pull units
//! work-stealing-style. The sweep is sequential + cache-local (active set ≤ 525,
//! L1/L2-resident), so unlike the old random-access cgranges index it should NOT
//! hit the memory-bandwidth wall that collapsed efficiency to 0.33.
//!
//!   cargo run --release --example sweep_parallel

use arrow::array::{Array, Int8Array, UInt16Array, UInt32Array};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::fs::File;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

fn col_u32(b: &arrow::record_batch::RecordBatch, i: usize) -> Vec<u32> {
    b.column(i).as_any().downcast_ref::<UInt32Array>().unwrap().values().to_vec()
}

fn read_tx(path: &str) -> (Vec<u16>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<i8>, Vec<u32>, Vec<u32>) {
    let r = ParquetRecordBatchReaderBuilder::try_new(File::open(path).unwrap()).unwrap().build().unwrap();
    let (mut c, mut s, mut e, mut x, mut st, mut cs, mut ce) =
        (vec![], vec![], vec![], vec![], vec![], vec![], vec![]);
    for b in r {
        let b = b.unwrap();
        c.extend(b.column(0).as_any().downcast_ref::<UInt16Array>().unwrap().values());
        s.extend(col_u32(&b, 1)); e.extend(col_u32(&b, 2)); x.extend(col_u32(&b, 3));
        st.extend(b.column(4).as_any().downcast_ref::<Int8Array>().unwrap().values());
        cs.extend(col_u32(&b, 5)); ce.extend(col_u32(&b, 6));
    }
    (c, s, e, x, st, cs, ce)
}
fn read3(path: &str) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let r = ParquetRecordBatchReaderBuilder::try_new(File::open(path).unwrap()).unwrap().build().unwrap();
    let (mut a, mut b2, mut c) = (vec![], vec![], vec![]);
    for b in r {
        let b = b.unwrap();
        if let Some(u) = b.column(0).as_any().downcast_ref::<UInt16Array>() {
            a.extend(u.values().iter().map(|&v| v as u32));
        } else { a.extend(col_u32(&b, 0)); }
        b2.extend(col_u32(&b, 1)); c.extend(col_u32(&b, 2));
    }
    (a, b2, c)
}

const DIST: u32 = 5000;

#[inline]
fn classify(pos: u32, k: usize, txs: &[u32], txe: &[u32], txst: &[i8], txcs: &[u32],
            txce: &[u32], txx: &[u32], eoff: &[u32], ecnt: &[u32], es: &[u32], ee: &[u32]) -> u32 {
    let (ts, te) = (txs[k], txe[k]);
    let fwd = txst[k] >= 0;
    if pos < ts { return if fwd { 1 } else { 2 }; }
    if pos > te { return if fwd { 2 } else { 1 }; }
    let ord = txx[k] as usize;
    let (off, cnt) = (eoff[ord] as usize, ecnt[ord] as usize);
    let (cs, ce) = (txcs[k], txce[k]);
    let (mut in_exon, mut splice) = (false, false);
    for e in off..off + cnt {
        in_exon |= pos >= es[e] && pos <= ee[e];
        splice |= (pos as i64 - es[e] as i64).unsigned_abs() <= 8
            || (pos as i64 - ee[e] as i64).unsigned_abs() <= 8;
    }
    let mut m = if in_exon {
        if cs != 0 && pos >= cs && pos <= ce { 16 } else if cs != 0 { 32 } else { 8 }
    } else { 4 };
    if splice { m |= 64; }
    m
}

/// Sweep one chromosome's [v_lo,v_hi) variants against [tx_lo,tx_hi) transcripts.
#[allow(clippy::too_many_arguments)]
fn sweep_chrom(v_lo: usize, v_hi: usize, tx_lo: usize, tx_hi: usize,
    vp: &[u32], txs: &[u32], txe: &[u32], txst: &[i8], txcs: &[u32], txce: &[u32],
    txx: &[u32], eoff: &[u32], ecnt: &[u32], es: &[u32], ee: &[u32]) -> u64 {
    let mut pairs = 0u64;
    let mut active: Vec<(u32, u32)> = Vec::new();
    let mut add = tx_lo;
    for vi in v_lo..v_hi {
        let pos = vp[vi];
        let at = pos.saturating_add(DIST);
        let et = pos.saturating_sub(DIST);
        while add < tx_hi && txs[add] <= at { active.push((txe[add], add as u32)); add += 1; }
        let mut k = 0;
        while k < active.len() { if active[k].0 < et { active.swap_remove(k); } else { k += 1; } }
        for &(_, kk) in &active {
            let _ = classify(pos, kk as usize, txs, txe, txst, txcs, txce, txx, eoff, ecnt, es, ee);
            pairs += 1;
        }
    }
    pairs
}

fn main() {
    let (vc, vp, _ve) = read3("/tmp/var_iv.parquet");
    let (txc, txs, txe, txx, txst, txcs, txce) = read_tx("/tmp/tx_iv.parquet");
    let (ex_tx, es, ee) = read3("/tmp/exon_iv.parquet");
    let n_ord = (*ex_tx.iter().max().unwrap() as usize) + 1;
    let (mut eoff, mut ecnt) = (vec![0u32; n_ord], vec![0u32; n_ord]);
    let mut i = 0;
    while i < ex_tx.len() {
        let t = ex_tx[i] as usize; eoff[t] = i as u32;
        let mut j = i; while j < ex_tx.len() && ex_tx[j] == ex_tx[i] { j += 1; }
        ecnt[t] = (j - i) as u32; i = j;
    }

    // per-chromosome work units (chrom, v_lo, v_hi, tx_lo, tx_hi)
    let mut units = Vec::new();
    let (n_v, n_tx) = (vp.len(), txs.len());
    let (mut vi, mut ti) = (0usize, 0usize);
    while vi < n_v {
        let ch = vc[vi] as u16;
        let v_lo = vi; while vi < n_v && (vc[vi] as u16) == ch { vi += 1; }
        while ti < n_tx && txc[ti] < ch { ti += 1; }
        let tx_lo = ti; let mut tx_hi = tx_lo;
        while tx_hi < n_tx && txc[tx_hi] == ch { tx_hi += 1; }
        units.push((v_lo, vi, tx_lo, tx_hi));
        ti = tx_hi;
    }
    eprintln!("{} chromosomes, {} variants, {} transcripts", units.len(), n_v, txs.len());

    for &threads in &[1usize, 2, 4, 8, 16] {
        let next = AtomicUsize::new(0);
        let total = AtomicUsize::new(0);
        let t = Instant::now();
        std::thread::scope(|scope| {
            for _ in 0..threads {
                let (next, total, units) = (&next, &total, &units);
                let (vp, txs, txe, txst, txcs, txce, txx, eoff, ecnt, es, ee) =
                    (&vp, &txs, &txe, &txst, &txcs, &txce, &txx, &eoff, &ecnt, &es, &ee);
                scope.spawn(move || {
                    let mut local = 0u64;
                    loop {
                        let u = next.fetch_add(1, Ordering::Relaxed);
                        if u >= units.len() { break; }
                        let (vl, vh, tl, th) = units[u];
                        local += sweep_chrom(vl, vh, tl, th, vp, txs, txe, txst, txcs, txce, txx, eoff, ecnt, es, ee);
                    }
                    total.fetch_add(local as usize, Ordering::Relaxed);
                });
            }
        });
        let s = t.elapsed().as_secs_f64();
        eprintln!("threads={:>2}  wall={:.4}s  pairs={}  speedup={:.1}x",
            threads, s, total.load(Ordering::Relaxed), 0.673 / s);
    }
}
