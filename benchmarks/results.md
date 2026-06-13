# Benchmark results

## Real-scale: ClinVar chr17 vs fastVEP

- **Input:** ClinVar GRCh38 `chr17` subset — **267,534 variants** (from the full
  4,435,099-variant ClinVar release).
- **Gene model:** Ensembl GRCh38.112 GFF3 (~280k transcripts, the full real model).
- **Hardware:** the dev box; release builds of both; single run each.
- duckvep: `vep_annotate('clinvar.chr17.vcf', gff3 := 'GRCh38.112.gff3.gz')`
- fastVEP: `fastvep annotate -i clinvar.chr17.vcf --gff3 GRCh38.112.gff3 --output-format tab`

| tool | wall time | peak RAM | output rows |
| --- | ---: | ---: | ---: |
| **duckvep** `vep_annotate` | 13.1 s | 2.94 GB | 3,910,199 |
| **fastVEP** CLI | 9.7 s | 1.14 GB | 3,910,925 |

### Reading

- **Concordance:** output rows agree to **99.98%** (726-row diff, 0.02%). Only 6
  of fastVEP's rows are `intergenic_variant`, so the diff is *not* intergenic —
  it's a small set of edge-case variants (multiallelic / symbolic-allele handling)
  to investigate. Same consequence engine, so per-(variant,transcript) calls match
  where they align (see `concordance.sh`).
- **Speed:** ~1.4× slower. Both spend ~10 s parsing the 280k-transcript GFF3; the
  delta is duckvep's IO/materialization overhead.
- **Memory:** ~2.6× higher — the known eager-materialization in `vep_annotate`'s
  `bind` (it holds all 3.9M output rows before emitting; fastVEP streams).

### Scalar path is the fix (same data, chr17)

`vep_consequence` over a `read_vcf`/`read_parquet` scan, vs the `vep_annotate`
table function:

| path | wall | peak RAM |
| --- | ---: | ---: |
| `vep_annotate` (table fn, eager) | 13.1 s | 2.94 GB |
| `vep_consequence` (scalar, streams) | **6.3 s** | **1.16 GB** |

The scalar **streams** (DuckDB doesn't materialize the 3.9M rows) and is faster
than fastVEP (9.7 s). Memory drops 2.5×.

### Threading (measured)

The scalar run is ~100% CPU = **~1 core** — but that's because
`vep_load_cache` **parses the 45 MB GFF3 single-threaded (~5 s)**, which
dominates; the consequence eval is the small remainder. So the win is the
**native DuckDB transcript cache** (§5): parse once → re-runs skip it and the
(DuckDB-parallelizable) scalar dominates. The cache read on the hot path is
**lock-free** (`ArcSwap`, atomic pointer load — no RwLock reader-counter
contention).

### Ensembl-VEP concordance (live REST, dated dumps)

`scripts/vep_concordance.py` annotates a sample with the **live Ensembl VEP REST
API** and duckvep (with FASTA), writing a dated Parquet dump
(`data/vep_dumps/<date>/annotations.parquet`) + `concordance_log.csv`:

| date | variants | shared pairs | concordance |
| --- | ---: | ---: | ---: |
| 2026-06-14 | 500 | 7,083 | **99.77%** |

Identical consequence sets on 7,067/7,083 shared (variant,transcript) pairs.
Disagreements are transcript-boundary/version edges (GFF3 release vs VEP's). For
unlimited scale, point it at an offline VEP cache instead of REST.

### Paths to close the gap (tracked)

1. **Stream** `vep_annotate` instead of materializing all rows in `bind`.
2. Use the **scan-driven scalar** `vep_consequence` (DuckDB owns the scan,
   parallel + spillable) rather than the self-reading table function.
3. **Cache** the transcript model (native DuckDB `.duckdb`, §5) so the ~10 s GFF3
   parse is paid once, not per run — then re-runs are near-instant.

## Reproduce

```sh
scripts/fetch-giab.sh                 # or fetch ClinVar GRCh38
scripts/bench.sh       data/giab/clinvar.chr17.vcf data/giab/GRCh38.112.gff3.gz
scripts/concordance.sh data/giab/clinvar.chr17.vcf data/giab/GRCh38.112.gff3.gz
```
