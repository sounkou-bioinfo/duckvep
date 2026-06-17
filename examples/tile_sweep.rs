//! Decisive architecture experiment (docs/kernel-algorithm.md §3): does a
//! cache-local sorted-merge SWEEP over compact integer interval columns beat the
//! DuckDB range join (measured 31.2 s / ~11 cores for candidate generation)?
//!
//! Reads two integer-column parquets (exported by DuckDB):
//!   tx_iv.parquet : chrom_id u16, start u32, end_pos u32   (sorted chrom,start)
//!   var_iv.parquet: chrom_id u16, pos u32, end_pos u32     (sorted chrom,pos)
//! and counts candidate (variant, transcript) pairs under the VEP ±distance
//! overlap, SINGLE-THREADED, with a per-chromosome two-pointer + end-heap sweep:
//! O((N+T)·log D) with the active set D≤518 (measured) staying L1/L2-resident.
//!
//!   cargo run --release --example tile_sweep

use arrow::array::{Array, UInt16Array, UInt32Array};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::fs::File;
use std::time::Instant;

/// Read a 3-column (u16, u32, u32) parquet into three SoA vectors.
fn read_iv(path: &str) -> (Vec<u16>, Vec<u32>, Vec<u32>) {
    let rdr = ParquetRecordBatchReaderBuilder::try_new(File::open(path).unwrap())
        .unwrap()
        .build()
        .unwrap();
    let (mut c, mut a, mut b) = (Vec::new(), Vec::new(), Vec::new());
    for batch in rdr {
        let batch = batch.unwrap();
        let cc = batch.column(0).as_any().downcast_ref::<UInt16Array>().unwrap();
        let aa = batch.column(1).as_any().downcast_ref::<UInt32Array>().unwrap();
        let bb = batch.column(2).as_any().downcast_ref::<UInt32Array>().unwrap();
        c.extend(cc.values());
        a.extend(aa.values());
        b.extend(bb.values());
    }
    (c, a, b)
}

fn main() {
    const DIST: u32 = 5000;
    let t_read = Instant::now();
    let (tx_chrom, tx_start, tx_end) = read_iv("/tmp/tx_iv.parquet");
    let (v_chrom, v_pos, _v_end) = read_iv("/tmp/var_iv.parquet");
    let read_s = t_read.elapsed().as_secs_f64();
    eprintln!(
        "loaded {} transcripts, {} variants in {:.2}s",
        tx_start.len(),
        v_pos.len(),
        read_s
    );

    let n_tx = tx_start.len();
    let n_v = v_pos.len();

    // EMIT version: maintain the active set as (end, tx_idx) and actually visit
    // every (variant, transcript) pair — the real cost of feeding a kernel, not
    // just counting. The checksum forces the emit to do observable work.
    let t_sweep = Instant::now();
    let mut pairs: u64 = 0;
    let mut max_active: usize = 0;
    let mut checksum: u64 = 0;
    let mut active: Vec<(u32, u32)> = Vec::new(); // (end, tx_idx)

    let mut ti = 0usize;
    let mut vi = 0usize;
    while vi < n_v {
        let chrom = v_chrom[vi];
        while ti < n_tx && tx_chrom[ti] < chrom {
            ti += 1;
        }
        let tx_lo = ti;
        let mut tx_hi = tx_lo;
        while tx_hi < n_tx && tx_chrom[tx_hi] == chrom {
            tx_hi += 1;
        }
        let mut add = tx_lo;
        active.clear();
        while vi < n_v && v_chrom[vi] == chrom {
            let pos = v_pos[vi];
            let add_thresh = pos.saturating_add(DIST);
            let evict_thresh = pos.saturating_sub(DIST);
            while add < tx_hi && tx_start[add] <= add_thresh {
                active.push((tx_end[add], add as u32)); // tx ordinal = global index
                add += 1;
            }
            // drop expired (end < evict_thresh); swap_remove is fine — order of the
            // active set does not matter for emission.
            let mut k = 0;
            while k < active.len() {
                if active[k].0 < evict_thresh {
                    active.swap_remove(k);
                } else {
                    k += 1;
                }
            }
            // EMIT every (variant vi, transcript active[j].1) pair.
            for &(_, tx_idx) in &active {
                checksum ^= (vi as u64).wrapping_mul(1_000_003).wrapping_add(tx_idx as u64);
                pairs += 1;
            }
            if active.len() > max_active {
                max_active = active.len();
            }
            vi += 1;
        }
        ti = tx_hi;
    }
    let sweep_s = t_sweep.elapsed().as_secs_f64();

    eprintln!("--- tile-sweep EMIT (single thread, every pair visited) ---");
    eprintln!("candidate pairs : {pairs}");
    eprintln!("max active set  : {max_active}");
    eprintln!("checksum        : {checksum}");
    eprintln!("sweep wall      : {sweep_s:.3}s");
    eprintln!("throughput      : {:.1}M pairs/s", pairs as f64 / 1e6 / sweep_s);
    eprintln!("vs DuckDB join-only: 31.2s / ~11 cores (candidate gen, materialized)");
}
