//! Prototype C (docs/kernel-algorithm.md): tile-sweep + a branch-light SoA region
//! classifier. For every (variant, transcript) pair the sweep produces, classify
//! the coarse structural relation (upstream/downstream/intron/exon/CDS/UTR/splice)
//! from compact integer columns + an offset/count exon pool — NO allocations, NO
//! strings, NO object graph. This measures the per-pair STRUCTURAL kernel cost (the
//! libsais "fused, branch-light pass" discipline) against the current scalar
//! `annotate_over` path's ~14 s kernel slice.
//!
//!   cargo run --release --example sweep_classify

use arrow::array::{Array, Int8Array, UInt16Array, UInt32Array};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::fs::File;
use std::time::Instant;

fn col_u32(b: &arrow::record_batch::RecordBatch, i: usize) -> impl Iterator<Item = u32> + '_ {
    b.column(i).as_any().downcast_ref::<UInt32Array>().unwrap().values().iter().copied()
}

/// transcripts: chrom_id u16, start u32, end u32, tx_idx u32, strand i8, cds_start u32, cds_end u32
fn read_tx(path: &str) -> (Vec<u16>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<i8>, Vec<u32>, Vec<u32>) {
    let r = ParquetRecordBatchReaderBuilder::try_new(File::open(path).unwrap()).unwrap().build().unwrap();
    let (mut c, mut s, mut e, mut x, mut st, mut cs, mut ce) =
        (vec![], vec![], vec![], vec![], vec![], vec![], vec![]);
    for b in r {
        let b = b.unwrap();
        c.extend(b.column(0).as_any().downcast_ref::<UInt16Array>().unwrap().values());
        s.extend(col_u32(&b, 1));
        e.extend(col_u32(&b, 2));
        x.extend(col_u32(&b, 3));
        st.extend(b.column(4).as_any().downcast_ref::<Int8Array>().unwrap().values());
        cs.extend(col_u32(&b, 5));
        ce.extend(col_u32(&b, 6));
    }
    (c, s, e, x, st, cs, ce)
}

/// 3 u32-ish columns (u16 first for var/chrom).
fn read3(path: &str) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let r = ParquetRecordBatchReaderBuilder::try_new(File::open(path).unwrap()).unwrap().build().unwrap();
    let (mut a, mut b2, mut c) = (vec![], vec![], vec![]);
    for b in r {
        let b = b.unwrap();
        // column 0 may be u16 (chrom_id) or u32 (tx_idx)
        if let Some(u16c) = b.column(0).as_any().downcast_ref::<UInt16Array>() {
            a.extend(u16c.values().iter().map(|&v| v as u32));
        } else {
            a.extend(col_u32(&b, 0));
        }
        b2.extend(col_u32(&b, 1));
        c.extend(col_u32(&b, 2));
    }
    (a, b2, c)
}

// region mask bits
const UP: u32 = 1;
const DOWN: u32 = 2;
const INTRON: u32 = 4;
const NC_EXON: u32 = 8;
const CDS: u32 = 16;
const UTR: u32 = 32;
const SPLICE: u32 = 64;

fn main() {
    const DIST: u32 = 5000;
    let t0 = Instant::now();
    let (vc, vp, _ve) = read3("/tmp/var_iv.parquet");
    let (txc, txs, txe, txx, txst, txcs, txce) = read_tx("/tmp/tx_iv.parquet");
    let (ex_tx, ex_s, ex_e) = read3("/tmp/exon_iv.parquet"); // tx_idx, start, end (sorted by tx_idx)
    eprintln!("loaded {} variants, {} tx, {} exons in {:.2}s",
        vp.len(), txs.len(), ex_s.len(), t0.elapsed().as_secs_f64());

    // exon offset/count indexed by tx_idx (ordinal)
    let n_tx_ord = (*ex_tx.iter().max().unwrap() as usize) + 1;
    let mut exon_off = vec![0u32; n_tx_ord];
    let mut exon_cnt = vec![0u32; n_tx_ord];
    let mut i = 0;
    while i < ex_tx.len() {
        let t = ex_tx[i] as usize;
        exon_off[t] = i as u32;
        let mut j = i;
        while j < ex_tx.len() && ex_tx[j] == ex_tx[i] { j += 1; }
        exon_cnt[t] = (j - i) as u32;
        i = j;
    }

    let n_tx = txs.len();
    let n_v = vp.len();
    let t1 = Instant::now();
    let mut pairs: u64 = 0;
    let mut hist = [0u64; 7];
    let mut active: Vec<(u32, u32)> = Vec::new(); // (tx_end, tx_array_index)

    let mut ti = 0usize;
    let mut vi = 0usize;
    while vi < n_v {
        let chrom = vc[vi] as u16;
        while ti < n_tx && txc[ti] < chrom { ti += 1; }
        let tx_lo = ti;
        let mut tx_hi = tx_lo;
        while tx_hi < n_tx && txc[tx_hi] == chrom { tx_hi += 1; }
        let mut add = tx_lo;
        active.clear();
        while vi < n_v && (vc[vi] as u16) == chrom {
            let pos = vp[vi];
            let add_thresh = pos.saturating_add(DIST);
            let evict_thresh = pos.saturating_sub(DIST);
            while add < tx_hi && txs[add] <= add_thresh {
                active.push((txe[add], add as u32));
                add += 1;
            }
            let mut k = 0;
            while k < active.len() {
                if active[k].0 < evict_thresh { active.swap_remove(k); } else { k += 1; }
            }
            for &(_, kk) in &active {
                let kk = kk as usize;
                let mask = classify(pos, kk, &txs, &txe, &txst, &txcs, &txce, &txx,
                                    &exon_off, &exon_cnt, &ex_s, &ex_e);
                pairs += 1;
                let mut m = mask;
                while m != 0 {
                    let b = m.trailing_zeros() as usize;
                    if b < 7 { hist[b] += 1; }
                    m &= m - 1;
                }
            }
            vi += 1;
        }
        ti = tx_hi;
    }
    let sweep_s = t1.elapsed().as_secs_f64();

    eprintln!("--- sweep + region mask (single thread) ---");
    eprintln!("pairs classified: {pairs}");
    eprintln!("wall            : {sweep_s:.3}s   ({:.1}M pairs/s)", pairs as f64 / 1e6 / sweep_s);
    let names = ["upstream", "downstream", "intron", "nc_exon", "CDS", "UTR", "splice"];
    for (n, h) in names.iter().zip(hist.iter()) {
        eprintln!("  {:<11}: {:>10}", n, h);
    }
    eprintln!("vs DuckDB join-only 31.2s/11c  &  full join+kernel 45.0s/12c");
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn classify(
    pos: u32, k: usize,
    txs: &[u32], txe: &[u32], txst: &[i8], txcs: &[u32], txce: &[u32], txx: &[u32],
    exon_off: &[u32], exon_cnt: &[u32], ex_s: &[u32], ex_e: &[u32],
) -> u32 {
    let ts = txs[k];
    let te = txe[k];
    let fwd = txst[k] >= 0;
    if pos < ts { return if fwd { UP } else { DOWN }; }
    if pos > te { return if fwd { DOWN } else { UP }; }
    // within transcript bounds: scan this transcript's exon pool
    let ord = txx[k] as usize;
    let off = exon_off[ord] as usize;
    let cnt = exon_cnt[ord] as usize;
    let (cs, ce) = (txcs[k], txce[k]);
    let mut in_exon = false;
    let mut splice = false;
    for e in off..off + cnt {
        let (es, ee) = (ex_s[e], ex_e[e]);
        in_exon |= pos >= es && pos <= ee;
        let d_s = (pos as i64 - es as i64).unsigned_abs();
        let d_e = (pos as i64 - ee as i64).unsigned_abs();
        splice |= d_s <= 8 || d_e <= 8;
    }
    let mut mask = 0u32;
    if in_exon {
        if cs != 0 && pos >= cs && pos <= ce { mask |= CDS; }
        else if cs != 0 { mask |= UTR; }
        else { mask |= NC_EXON; }
    } else {
        mask |= INTRON;
    }
    if splice { mask |= SPLICE; }
    mask
}
