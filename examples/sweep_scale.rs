//! Scaling with variant count — "what if I annotate gnomAD (150 M variants)?".
//! Replicates the real HG002 variant set k× (preserving realistic spatial
//! distribution and fan-out) to reach 4M → 16M → 64M → ~150M variants, and runs
//! the region-parallel sweep+classify at each scale. Shows the two things that
//! matter at gnomAD scale: time is LINEAR in N (P = fan-out·N dominates), and the
//! sweep's WORKING memory is CONSTANT (active set O(D), independent of N) — the
//! variant array is the only thing that grows, and in production it STREAMS.
//!
//!   cargo run --release --example sweep_scale

use arrow::array::{Array, Int8Array, UInt16Array, UInt32Array};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::fs::File;
use std::sync::atomic::{AtomicUsize, Ordering};
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
fn rss_mb() -> u64 {
    std::fs::read_to_string("/proc/self/status").ok()
        .and_then(|s| s.lines().find(|l| l.starts_with("VmHWM")).map(|l| l.to_string()))
        .and_then(|l| l.split_whitespace().nth(1).and_then(|v| v.parse::<u64>().ok()))
        .map(|kb| kb / 1024).unwrap_or(0)
}

const DIST: u32 = 5000;
#[inline]
fn classify(pos: u32, k: usize, txs: &[u32], txe: &[u32], txst: &[i8], txcs: &[u32],
            txce: &[u32], txx: &[u32], eoff: &[u32], ecnt: &[u32], es: &[u32], ee: &[u32]) -> u32 {
    let (ts, te) = (txs[k], txe[k]); let fwd = txst[k] >= 0;
    if pos < ts { return if fwd {1} else {2}; }
    if pos > te { return if fwd {2} else {1}; }
    let ord = txx[k] as usize; let (off, cnt) = (eoff[ord] as usize, ecnt[ord] as usize);
    let (cs, ce) = (txcs[k], txce[k]); let (mut ie, mut sp) = (false, false);
    for e in off..off+cnt { ie |= pos>=es[e] && pos<=ee[e];
        sp |= (pos as i64-es[e] as i64).unsigned_abs()<=8 || (pos as i64-ee[e] as i64).unsigned_abs()<=8; }
    let mut m = if ie { if cs!=0 && pos>=cs && pos<=ce {16} else if cs!=0 {32} else {8} } else {4};
    if sp { m |= 64; } m
}

fn main() {
    let (vc0, vp0, _ve) = read3("/tmp/var_iv.parquet");
    let (txc, txs, txe, txx, txst, txcs, txce) = read_tx("/tmp/tx_iv.parquet");
    let (extx, es, ee) = read3("/tmp/exon_iv.parquet");
    let n_ord = (*extx.iter().max().unwrap() as usize) + 1;
    let (mut eoff, mut ecnt) = (vec![0u32; n_ord], vec![0u32; n_ord]);
    let mut i = 0; while i < extx.len() { let t = extx[i] as usize; eoff[t] = i as u32;
        let mut j = i; while j < extx.len() && extx[j]==extx[i] { j+=1; } ecnt[t]=(j-i) as u32; i=j; }
    eprintln!("model: {} transcripts, {} exons resident\n", txs.len(), es.len());

    let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8);
    eprintln!("{:>5} {:>10} {:>9} {:>13} {:>10} {:>9}", "scale", "variants", "wall_s", "pairs", "Mpairs/s", "RSS_MB");

    // k× replication: each (chrom,pos) repeated k times keeps the array sorted by
    // (chrom,pos) and gives realistic fan-out at k× density.
    for &k in &[1usize, 4, 16, 37] {
        let nv = vp0.len() * k;
        let mut vc = Vec::with_capacity(nv);
        let mut vp = Vec::with_capacity(nv);
        for idx in 0..vp0.len() {
            for _ in 0..k { vc.push(vc0[idx]); vp.push(vp0[idx]); }
        }
        // per-chrom work units
        let mut units = Vec::new();
        let (mut vi, mut ti) = (0usize, 0usize);
        while vi < vp.len() {
            let ch = vc[vi] as u16; let vlo = vi;
            while vi < vp.len() && (vc[vi] as u16)==ch { vi+=1; }
            while ti < txs.len() && txc[ti] < ch { ti+=1; }
            let tlo = ti; let mut thi = tlo;
            while thi < txs.len() && txc[thi]==ch { thi+=1; }
            units.push((vlo, vi, tlo, thi)); ti = thi;
        }
        let next = AtomicUsize::new(0); let total = AtomicUsize::new(0);
        let t = Instant::now();
        std::thread::scope(|sc| {
            for _ in 0..threads {
                let (next, total, units, vp) = (&next, &total, &units, &vp);
                let (txs,txe,txst,txcs,txce,txx,eoff,ecnt,es,ee) =
                    (&txs,&txe,&txst,&txcs,&txce,&txx,&eoff,&ecnt,&es,&ee);
                sc.spawn(move || {
                    let mut local = 0u64; let mut active: Vec<(u32,u32)> = Vec::new();
                    loop {
                        let u = next.fetch_add(1, Ordering::Relaxed); if u >= units.len() { break; }
                        let (vl, vh, tl, th) = units[u]; active.clear(); let mut add = tl;
                        for vi in vl..vh {
                            let pos = vp[vi]; let at = pos.saturating_add(DIST); let et = pos.saturating_sub(DIST);
                            while add < th && txs[add] <= at { active.push((txe[add], add as u32)); add+=1; }
                            let mut kk=0; while kk<active.len() { if active[kk].0<et {active.swap_remove(kk);} else {kk+=1;} }
                            for &(_, k2) in &active { let _ = classify(pos, k2 as usize, txs,txe,txst,txcs,txce,txx,eoff,ecnt,es,ee); local+=1; }
                        }
                    }
                    total.fetch_add(local as usize, Ordering::Relaxed);
                });
            }
        });
        let s = t.elapsed().as_secs_f64();
        let p = total.load(Ordering::Relaxed) as u64;
        eprintln!("{:>4}x {:>10} {:>9.3} {:>13} {:>10.0} {:>9}",
            k, nv, s, p, p as f64/1e6/s, rss_mb());
    }
    eprintln!("\n(threads={threads}; variant array materialized here — in production it STREAMS from read_vcf at O(active set) memory.)");
}
