//! Fused multi-track sweep: can ONE pass over the variant stream do consequence
//! AND supplementary annotations (point + region) together — for ~free?
//!
//! Everything in genomics is coordinate-sorted, so we never need a variant
//! key/hash. The single sweep cursor maintains, in lockstep:
//!   - active TRANSCRIPTS      -> consequence region mask         (interval set)
//!   - active REGION SAs        -> regulatory/panel/conservation   (interval sets)
//!   - merge cursors for POINT SAs -> gnomAD AF / ClinVar / dbSNP  (sorted-merge)
//! This is the duckhts/mosdepth multi-track sweep vs fastVEP-SA's M hash-joined
//! passes (var32/bloom/block). We measure: consequence-only vs +1 point SA vs
//! +1 point +3 region SAs, all in one pass.
//!
//!   cargo run --release --example sweep_multi

use arrow::array::{Array, Int8Array, UInt16Array, UInt32Array};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::fs::File;
use std::time::Instant;

fn cu32(b: &arrow::record_batch::RecordBatch, i: usize) -> Vec<u32> {
    b.column(i).as_any().downcast_ref::<UInt32Array>().unwrap().values().to_vec()
}
fn read_tx(p: &str) -> (Vec<u16>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<i8>, Vec<u32>, Vec<u32>) {
    let r = ParquetRecordBatchReaderBuilder::try_new(File::open(p).unwrap()).unwrap().build().unwrap();
    let (mut c, mut s, mut e, mut x, mut st, mut cs, mut ce) = (vec![], vec![], vec![], vec![], vec![], vec![], vec![]);
    for b in r { let b = b.unwrap();
        c.extend(b.column(0).as_any().downcast_ref::<UInt16Array>().unwrap().values());
        s.extend(cu32(&b,1)); e.extend(cu32(&b,2)); x.extend(cu32(&b,3));
        st.extend(b.column(4).as_any().downcast_ref::<Int8Array>().unwrap().values());
        cs.extend(cu32(&b,5)); ce.extend(cu32(&b,6)); }
    (c, s, e, x, st, cs, ce)
}
fn read3(p: &str) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let r = ParquetRecordBatchReaderBuilder::try_new(File::open(p).unwrap()).unwrap().build().unwrap();
    let (mut a, mut b2, mut c) = (vec![], vec![], vec![]);
    for b in r { let b = b.unwrap();
        if let Some(u) = b.column(0).as_any().downcast_ref::<UInt16Array>() { a.extend(u.values().iter().map(|&v| v as u32)); }
        else { a.extend(cu32(&b,0)); }
        b2.extend(cu32(&b,1)); c.extend(cu32(&b,2)); }
    (a, b2, c)
}
const DIST: u32 = 5000;

fn main() {
    let (vc, vp, _ve) = read3("/tmp/var_iv.parquet");
    let (txc, txs, txe, _txx, txst, _txcs, _txce) = read_tx("/tmp/tx_iv.parquet");
    let n_v = vp.len();
    let n_tx = txs.len();

    // --- POINT SA (gnomAD-like): a sorted (chrom, pos) -> af table. Use the variant
    // positions themselves as common-variant sites (high hit rate, realistic for
    // gnomAD common variants). af synthesized from pos.
    let sa_pos = vp.clone(); // already sorted within chrom (var_iv is sorted)
    let sa_chrom = vc.clone();
    let sa_af: Vec<f32> = sa_pos.iter().map(|&p| ((p.wrapping_mul(2654435761) >> 8) & 0xffff) as f32 / 65535.0).collect();

    // --- REGION SAs (regulatory/panel/conservation-like): 3 sorted interval tracks.
    // Generate ~per-chrom intervals within the transcript span; width ~500.
    let mut chrom_lo = std::collections::HashMap::new();
    let mut chrom_hi = std::collections::HashMap::new();
    for i in 0..n_tx {
        let c = txc[i];
        chrom_lo.entry(c).or_insert(txs[i]);
        let e = chrom_hi.entry(c).or_insert(0u32);
        if txe[i] > *e { *e = txe[i]; }
    }
    let make_region = |per_chrom: u32, width: u32, seed: u64| -> (Vec<u16>, Vec<u32>, Vec<u32>) {
        let (mut rc, mut rs, mut re) = (vec![], vec![], vec![]);
        let mut chroms: Vec<u16> = chrom_lo.keys().copied().collect();
        chroms.sort_unstable();
        for &c in &chroms {
            let lo = chrom_lo[&c]; let hi = chrom_hi[&c];
            if hi <= lo { continue; }
            let span = (hi - lo) as u64;
            let mut starts: Vec<u32> = (0..per_chrom).map(|i| {
                let mut z = seed ^ (c as u64).wrapping_mul(0x9E3779B97F4A7C15) ^ (i as u64).wrapping_mul(0xC2B2AE3D27D4EB4F);
                z ^= z >> 33; z = z.wrapping_mul(0xff51afd7ed558ccd); z ^= z >> 33;
                lo + (z % span) as u32
            }).collect();
            starts.sort_unstable();
            for s in starts { rc.push(c); rs.push(s); re.push(s + width); }
        }
        (rc, rs, re)
    };
    let regions: Vec<(Vec<u16>, Vec<u32>, Vec<u32>)> = vec![
        make_region(20000, 800, 1),   // regulatory
        make_region(5000, 50000, 2),  // large domains
        make_region(40000, 200, 3),   // conserved elements
    ];
    eprintln!("model {} tx; point SA {} sites; region SAs {} tracks ({} intervals)",
        n_tx, sa_pos.len(), regions.len(), regions.iter().map(|r| r.1.len()).sum::<usize>());

    // run the fused sweep with `n_point` point SAs and `n_region` region SAs.
    let run = |n_point: usize, n_region: usize| -> (f64, u64, f64) {
        let t = Instant::now();
        let mut pairs = 0u64;
        let mut reg_hits = 0u64;
        let mut af_sum = 0f64;
        let mut active: Vec<(u32, u32)> = Vec::new();              // transcripts (end, idx)
        let mut ract: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n_region]; // region active sets
        let mut radd = vec![0usize; n_region];
        let mut sa_ptr = 0usize;                                   // point-SA merge cursor
        let mut ti = 0usize;
        let mut vi = 0usize;
        while vi < n_v {
            let chrom = vc[vi];
            // position transcript cursor + region cursors to this chrom
            while ti < n_tx && (txc[ti] as u32) < chrom { ti += 1; }
            let tx_lo = ti; let mut tx_hi = tx_lo;
            while tx_hi < n_tx && txc[tx_hi] as u32 == chrom { tx_hi += 1; }
            let mut add = tx_lo;
            active.clear();
            // region cursors
            for r in 0..n_region {
                let (rc, _rs, _re) = &regions[r];
                let mut p = radd[r];
                while p < rc.len() && (rc[p] as u32) < chrom { p += 1; }
                radd[r] = p;
                ract[r].clear();
            }
            // point-SA cursor to this chrom
            while sa_ptr < sa_chrom.len() && sa_chrom[sa_ptr] < chrom { sa_ptr += 1; }

            while vi < n_v && vc[vi] == chrom {
                let pos = vp[vi];
                let at = pos.saturating_add(DIST);
                let et = pos.saturating_sub(DIST);
                // --- consequence track: active transcripts
                while add < tx_hi && txs[add] <= at { active.push((txe[add], add as u32)); add += 1; }
                let mut k = 0; while k < active.len() { if active[k].0 < et { active.swap_remove(k); } else { k += 1; } }
                pairs += active.len() as u64;
                // --- region SA tracks: active intervals (exact overlap, no halo here)
                for r in 0..n_region {
                    let (rc, rs, re) = &regions[r];
                    let av = &mut ract[r];
                    while radd[r] < rc.len() && rc[radd[r]] as u32 == chrom && rs[radd[r]] <= pos { av.push((re[radd[r]], radd[r] as u32)); radd[r] += 1; }
                    let mut k = 0; while k < av.len() { if av[k].0 < pos { av.swap_remove(k); } else { k += 1; } }
                    reg_hits += av.len() as u64;
                }
                // --- point SA track(s): sorted-merge lookup
                for _ in 0..n_point {
                    while sa_ptr < sa_pos.len() && sa_chrom[sa_ptr] == chrom && sa_pos[sa_ptr] < pos { sa_ptr += 1; }
                    if sa_ptr < sa_pos.len() && sa_chrom[sa_ptr] == chrom && sa_pos[sa_ptr] == pos {
                        af_sum += sa_af[sa_ptr] as f64;
                    }
                }
                vi += 1;
            }
            ti = tx_hi;
        }
        (t.elapsed().as_secs_f64(), pairs.wrapping_add(reg_hits), af_sum)
    };

    eprintln!("\n{:<32} {:>9} {:>13}", "config (one pass each)", "wall_s", "cons_pairs");
    let (s0, p0, _) = run(0, 0);
    eprintln!("{:<32} {:>9.3} {:>13}", "consequence only", s0, p0);
    let (s1, _, af1) = run(1, 0);
    eprintln!("{:<32} {:>9.3}   (+point SA, af_sum={:.0})", "consequence + gnomAD AF", s1, af1);
    let (s2, _, _) = run(1, 3);
    eprintln!("{:<32} {:>9.3}   (+point +3 region SAs)", "consequence + 4 annotations", s2);
    eprintln!("\nadded cost of 4 annotations in the SAME pass: {:.3}s ({:.0}% over consequence-only)",
        s2 - s0, (s2 - s0) / s0 * 100.0);
}
