
<!-- results.md is generated from results.Rmd — edit the .Rmd and run `make benchmarks`.
     Tables are rendered live from the recorded data via duckknit. -->

# duckvep benchmarks

Reproducible benchmark + Ensembl-VEP concordance, rendered from recorded
data. duckvep wraps the same consequence engine as fastVEP, so it is a
fair head-to-head on a real human gene model; the differences are
data-engineering (caching, streaming), not the science.

**Setup:** Ensembl GRCh38.112 GFF3 (~280k transcripts); release builds;
single run each on the dev box. Reproduce with `scripts/bench.sh` and
`scripts/vep_concordance.py` (see `benchmarks/README.md`).

## Throughput vs fastVEP

``` duckdb
SELECT dataset, variants, tool, wall_clock AS "wall", peak_rss_mb AS "rss (MB)"
FROM 'benchmarks/data/timings.csv'
ORDER BY variants DESC, tool;
```

| dataset                   | variants | tool                            | wall    | rss (MB) |
|---------------------------|---------:|---------------------------------|---------|---------:|
| GIAB HG002 (whole genome) |  4048342 | duckvep (warm cache; streaming) | 0:09.02 |     1985 |
| GIAB HG002 (whole genome) |  4048342 | fastVEP CLI                     | 0:19.46 |     1113 |
| ClinVar chr17             |   267534 | duckvep (cold / parses GFF3)    | 0:07.23 |     1370 |
| ClinVar chr17             |   267534 | duckvep (warm cache)            | 0:01.88 |     1276 |
| ClinVar chr17             |   267534 | fastVEP CLI                     | 0:06.07 |      694 |

duckvep is faster on both the full GIAB HG002 whole genome (4.0M
variants) and the ClinVar chr17 subset, because `vep_load_cache`
**caches the parsed transcript model** while fastVEP re-parses the GFF3
every run. duckvep’s higher memory on GIAB is the `read_vcf` eager-load
(4M variants materialized) — the streaming follow-up closes it; fastVEP
streams.

## Ensembl-VEP concordance

Per-(variant, transcript) consequence agreement vs the **live Ensembl
VEP** (REST API), over sampled real ClinVar variants annotated with the
reference FASTA (so coding calls — synonymous vs missense — are exact).
Both duckvep **and** fastVEP (the underlying engine) are compared
against Ensembl VEP:

``` duckdb
SELECT date, engine, n_variants AS "variants", pairs AS "shared pairs",
       agree, pct AS "concordance %"
FROM 'data/vep_dumps/concordance_log.csv'
ORDER BY date DESC, engine;
```

| date       | engine  | variants | shared pairs | agree | concordance % |
|------------|---------|---------:|-------------:|------:|--------------:|
| 2026-06-14 | duckvep |      200 |         2677 |  2677 |         100.0 |
| 2026-06-14 | fastvep |      200 |         2677 |  2677 |         100.0 |

duckvep matches Ensembl VEP because it wraps the same engine; the table
also shows fastVEP itself vs Ensembl VEP. Remaining disagreements (at
larger samples) are transcript-boundary/version edges (GFF3 release vs
VEP’s transcript set), not systematic. Dated annotation dumps accumulate
under `data/vep_dumps/<date>/annotations.parquet`.
