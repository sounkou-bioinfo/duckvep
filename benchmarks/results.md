
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

| dataset                   | variants | tool                 | wall    | rss (MB) |
|---------------------------|---------:|----------------------|---------|---------:|
| GIAB HG002 (whole genome) |  4048342 | duckvep (warm cache) | 0:09.10 |     5568 |
| GIAB HG002 (whole genome) |  4048342 | fastVEP CLI          | 0:21.58 |      589 |
| ClinVar chr17             |   267534 | duckvep (warm cache) | 0:01.4  |      737 |
| ClinVar chr17             |   267534 | fastVEP CLI          | 0:09.7  |     1114 |

duckvep is faster on both the full GIAB HG002 whole genome (4.0M
variants) and the ClinVar chr17 subset, because `vep_load_cache`
**caches the parsed transcript model** while fastVEP re-parses the GFF3
every run. duckvep’s higher memory on GIAB is the `read_vcf` eager-load
(4M variants materialized) — the streaming follow-up closes it; fastVEP
streams.

## Ensembl-VEP concordance

Per-(variant, transcript) consequence agreement between duckvep and the
**live Ensembl VEP** (REST API), over sampled real ClinVar variants —
annotated with the reference FASTA so coding calls (synonymous vs
missense) are exact.

``` duckdb
SELECT date, n_variants AS "variants", pairs AS "shared pairs",
       agree, pct AS "concordance %"
FROM 'data/vep_dumps/concordance_log.csv'
ORDER BY date DESC;
```

| date       | variants | shared pairs | agree | concordance % |
|------------|---------:|-------------:|------:|--------------:|
| 2026-06-14 |      500 |         7083 |  7067 |         99.77 |

Disagreements are transcript-boundary/version edges (GFF3 release vs
VEP’s transcript set), not systematic. Dated annotation dumps accumulate
under `data/vep_dumps/<date>/annotations.parquet`.
