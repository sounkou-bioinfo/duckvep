#!/usr/bin/env python3
"""Stratified VEP-conformance report — the STATISTICAL tier of the conformance framework.

Reads a dated concordance dump (`data/vep_dumps/<date>/annotations.parquet`, produced by
`correctness/vep_concordance.py`: VEP-116 `--gff` + duckvep + fastVEP, one row per
(source, pos, ref, alt, transcript_id) with the SO-term set) and reports per-(variant,
transcript) conformance of duckvep vs the VEP oracle, **stratified** by the axes the
formal class-model uses observably: consequence class (the VEP/gold most-severe SO term),
variant type, and length bin.

Per stratum it reports N (shared pairs), concordant, discordant, the point discordance
rate, and the **95% upper bound** — Clopper-Pearson, which for a zero-discordance stratum
reduces to the *rule of three* (true rate < ~3/N). That upper bound is the provable
statement: "duckvep matches VEP in this stratum at < U·10^-k (95%), N=…". fastVEP is shown
alongside so duckvep's bound is read against the engine it patches.

This is the population-confidence tier; the formal tier (covering-array witness generation
on a synthetic transcript) feeds the SAME diff and reports class COVERAGE instead of CIs.
They meet at a coverage-guided differential fuzzer — see conformance/README.md.

Usage: conformance/stratified_conformance.py [annotations.parquet]  (default: newest dump)
"""
import os, sys, subprocess, math, csv

ROOT = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True).stdout.strip()
DUCKDB = f"{ROOT}/.tools/duckdb"
DUMP = sys.argv[1] if len(sys.argv) > 1 else subprocess.run(
    ["bash", "-c", f"ls -t {ROOT}/data/vep_dumps/*/annotations.parquet | head -1"],
    capture_output=True, text=True).stdout.strip()
if not DUMP or not os.path.exists(DUMP):
    sys.exit("no annotations.parquet dump found — run correctness/vep_concordance.py first")
OUT = f"{ROOT}/conformance/data"
os.makedirs(OUT, exist_ok=True)

# variant_type + length_bin are derived from (ref, alt); the gold consequence CLASS is VEP's
# most-severe SO term. The join key is (pos, ref, alt, transcript_id) — present for vep,
# duckvep AND fastvep, so the comparison is on shared (emitted-by-both) pairs.
SQL = f"""
COPY (
  WITH src AS (SELECT source, pos, ref, alt, transcript_id, consequence FROM read_parquet('{DUMP}')),
  v AS (SELECT pos,ref,alt,transcript_id, consequence AS vc FROM src WHERE source='vep'),
  d AS (SELECT pos,ref,alt,transcript_id, consequence AS dc FROM src WHERE source='duckvep'),
  f AS (SELECT pos,ref,alt,transcript_id, consequence AS fc FROM src WHERE source='fastvep'),
  shape AS (
    SELECT v.*, d.dc, f.fc,
      -- variant type from the minimal alleles (VEP writes '-' for the empty side)
      CASE
        WHEN ref='-' THEN 'ins'
        WHEN alt='-' THEN 'del'
        WHEN length(ref)=1 AND length(alt)=1 THEN 'snv'
        WHEN length(ref)=length(alt) THEN 'mnv'
        WHEN length(alt)>length(ref) THEN 'ins'
        WHEN length(alt)<length(ref) THEN 'del'
        ELSE 'delins' END AS var_type,
      -- signed length change (alt - ref), '-' counted as 0
      ( (CASE WHEN alt='-' THEN 0 ELSE length(alt) END)
      - (CASE WHEN ref='-' THEN 0 ELSE length(ref) END) ) AS net
    FROM v JOIN d USING(pos,ref,alt,transcript_id)
           JOIN f USING(pos,ref,alt,transcript_id)
  ),
  binned AS (
    SELECT *,
      -- consequence CLASS = the most-severe-ish single token (the gold VEP set's first sorted
      -- HIGH term if any, else the whole set) kept coarse: use the full VEP SO-term set string.
      vc AS consequence_class,
      CASE
        WHEN net=0 THEN '0'
        WHEN abs(net) BETWEEN 1 AND 3 THEN (CASE WHEN net<0 THEN '-' ELSE '+' END)||abs(net)::VARCHAR
        WHEN abs(net) BETWEEN 4 AND 10 THEN (CASE WHEN net<0 THEN '-' ELSE '+' END)||'4..10'
        WHEN abs(net) BETWEEN 11 AND 50 THEN (CASE WHEN net<0 THEN '-' ELSE '+' END)||'11..50'
        ELSE (CASE WHEN net<0 THEN '-' ELSE '+' END)||'>50' END AS length_bin
    FROM shape
  )
  SELECT consequence_class, var_type, length_bin,
    count(*) AS n,
    count(*) FILTER (WHERE dc = vc) AS dv_concordant,
    count(*) FILTER (WHERE dc <> vc) AS dv_discordant,
    count(*) FILTER (WHERE fc = vc) AS fv_concordant,
    count(*) FILTER (WHERE fc <> vc) AS fv_discordant
  FROM binned
  GROUP BY consequence_class, var_type, length_bin
  ORDER BY dv_discordant DESC, n DESC
) TO '{OUT}/stratified_conformance.csv' (HEADER);
"""
subprocess.run([DUCKDB, "-c", SQL], check=True)


def cp_upper(k, n, alpha=0.05):
    """Clopper-Pearson upper 95% bound on the discordance rate from k discordant of n.
    For k=0 this is 1-(alpha/2)^(1/n) ~= the rule of three (~3/n)."""
    if n == 0:
        return 1.0
    if k == 0:
        return 1.0 - (alpha / 2) ** (1.0 / n)
    # invert the Beta quantile via bisection (no scipy dependency)
    from math import lgamma
    def betacdf(x, a, b):
        # regularized incomplete beta via continued fraction (Lentz)
        if x <= 0: return 0.0
        if x >= 1: return 1.0
        lbeta = lgamma(a) + lgamma(b) - lgamma(a + b)
        front = math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) / a
        f, c, d = 1.0, 1.0, 0.0
        for i in range(0, 200):
            m = i // 2
            if i == 0: num = 1.0
            elif i % 2 == 0: num = (m * (b - m) * x) / ((a + 2*m - 1) * (a + 2*m))
            else: num = -((a + m) * (a + b + m) * x) / ((a + 2*m) * (a + 2*m + 1))
            d = 1.0 + num * d
            if abs(d) < 1e-30: d = 1e-30
            d = 1.0 / d
            c = 1.0 + num / c
            if abs(c) < 1e-30: c = 1e-30
            f *= d * c
            if abs(1.0 - d * c) < 1e-10: break
        return front * (f - 1.0)
    lo, hi = k / n, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        # P(X <= k) with p=mid should equal alpha/2 at the upper bound
        if 1.0 - betacdf(mid, k + 1, n - k) > alpha / 2:
            hi = mid
        else:
            lo = mid
    return hi


rows = list(csv.DictReader(open(f"{OUT}/stratified_conformance.csv")))
for r in rows:
    n = int(r["n"]); k = int(r["dv_discordant"])
    r["dv_upper95"] = f"{cp_upper(k, n):.2e}"
# rewrite with the CI column
with open(f"{OUT}/stratified_conformance.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

tot_n = sum(int(r["n"]) for r in rows)
tot_k = sum(int(r["dv_discordant"]) for r in rows)
tot_fk = sum(int(r["fv_discordant"]) for r in rows)
print(f"Stratified VEP-116 conformance (duckvep)  —  {len(rows)} strata, {tot_n} shared (variant,transcript) pairs")
print(f"  overall: duckvep {tot_n-tot_k}/{tot_n} concordant ({tot_k} discordant, <= {cp_upper(tot_k,tot_n):.2e} @95%)"
      f"  |  fastVEP {tot_fk} discordant")
print(f"  {'consequence_class':42} {'type':6} {'len':7} {'N':>7} {'dv_disc':>7} {'dv_<=95%':>9} {'fv_disc':>7}")
shown = [r for r in rows if int(r["dv_discordant"]) > 0] + \
        sorted([r for r in rows if int(r["dv_discordant"]) == 0], key=lambda r: -int(r["n"]))[:8]
for r in shown[:22]:
    print(f"  {r['consequence_class'][:42]:42} {r['var_type']:6} {r['length_bin']:7} "
          f"{r['n']:>7} {r['dv_discordant']:>7} {r['dv_upper95']:>9} {r['fv_discordant']:>7}")
print(f"  full table -> {OUT}/stratified_conformance.csv")
