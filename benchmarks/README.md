# Benchmarks & concordance

Three harnesses (`../scripts/`), all runnable on the small test data and scaling
to full-size GIAB.

## Concordance — `correctness/concordance.sh <vcf> <gff3> [--fasta ref.fa] [--vep vep.tsv]`

Per-(variant, transcript) consequence agreement: duckvep vs fastVEP CLI (or an
Ensembl VEP `--vep` tsv). Writes `concordance.out/diffs.tsv`.

On the bundled test set (`DuckfastVEP/tests/test.vcf`): **every joined
(variant, transcript) pair agrees** (duckvep wraps the same engine). The only
non-joins are indels, because Ensembl VEP / fastVEP **left-align and minimize**
them (e.g. `AA>A` → Location `+1`, Allele `-`) while duckvep currently emits the
raw VCF allele. Normalizing both sides before the join (DESIGN.md §3.4) closes
this — and is required for correct annotation joins anyway.

## Benchmark — `scripts/bench.sh <vcf> <gff3> [--fasta ref.fa]`

End-to-end throughput, duckvep vs fastVEP (hyperfine if installed, else `time`).

## Full-size GIAB — `scripts/fetch-giab.sh`

Downloads HG002 GRCh38 benchmark VCF + Ensembl GRCh38 GFF3 + reference into
`data/giab/` (gitignored, ~GB). Then:

```sh
scripts/fetch-giab.sh
scripts/bench.sh       data/giab/HG002.vcf.gz data/giab/GRCh38.116.gff3.gz --fasta data/giab/GRCh38.fa
correctness/concordance.sh data/giab/HG002.vcf.gz data/giab/GRCh38.116.gff3.gz --fasta data/giab/GRCh38.fa
```

> Full GIAB runs are gated behind the manual fetch (not CI). For Ensembl-VEP
> concordance, produce a VEP tab file for the same VCF and pass `--vep`.
