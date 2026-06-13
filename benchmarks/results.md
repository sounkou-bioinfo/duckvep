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
