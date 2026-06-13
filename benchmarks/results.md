
<!-- results.md is generated from results.Rmd — edit the .Rmd and run `make benchmarks`.
     Tables are read **directly from the recorded CSVs** and rendered with knitr. -->

# duckvep benchmarks

R Markdown report — every table below is read **directly from a recorded
CSV** under `benchmarks/data/` (regenerate the data with
`scripts/benchmark.sh` / `scripts/vep_concordance.py`, then
`make benchmarks`). duckvep wraps the same consequence engine as
fastVEP, so head-to-heads are fair; the differences are data-engineering
(caching, streaming, columnar), not the science. **Setup:** Ensembl GFF3
gene models; release builds (`opt-level=3`, `lto=fat`); steady state
(warm transcript cache); single run each.

> Scope: duckvep currently emits **consequence + amino-acids/codons**;
> HGVS strings and supplementary annotations are not yet wired, so
> throughput is “consequence”, not “full annotation”.

## Throughput (consequence)

| organism         | assembly | transcripts | variants  | source     | time | variants/s |
|:-----------------|:---------|:------------|:----------|:-----------|:-----|:-----------|
| Human (full WGS) | GRCh38   | 252,980     | 4,048,342 | GIAB HG002 | 9.0s | 449,816    |
| Human (chr17)    | GRCh38   | 252,980     | 267,534   | ClinVar    | 1.9s | 140,807    |

Add organisms by running
`scripts/benchmark.sh "<organism>" <assembly> <gff3> <vcf> [fasta]`
(yeast / fly / arabidopsis / mouse / human).

## Head-to-head vs fastVEP (full output, same gene model)

| dataset                   | variants | tool                            | wall_clock | peak_rss_mb | output | note                            |
|:--------------------------|---------:|:--------------------------------|:-----------|------------:|:-------|:--------------------------------|
| GIAB HG002 (whole genome) |  4048342 | fastVEP CLI                     | 0:19.46    |        1113 | full   | streams + re-parses GFF3        |
| GIAB HG002 (whole genome) |  4048342 | duckvep (warm cache; streaming) | 0:09.02    |        1985 | full   | cached GFF3; streaming read_vcf |
| ClinVar chr17             |   267534 | fastVEP CLI                     | 0:06.07    |         694 | full   | +fasta                          |
| ClinVar chr17             |   267534 | duckvep (cold / parses GFF3)    | 0:07.23    |        1370 | full   | +fasta; parity with fastVEP     |
| ClinVar chr17             |   267534 | duckvep (warm cache)            | 0:01.88    |        1276 | full   | +fasta; caching win             |

duckvep’s edge is **caching** the parsed transcript model (warm); a
*cold* run (parsing the GFF3) is ~parity with fastVEP — same engine.
Streaming `read_vcf` keeps memory bounded (full GIAB 2.0 GB vs 7.3 GB
eager). Offline Ensembl-VEP timing is pending its cache install;
concordance vs the live VEP is below.

## Footprint

| tool        | dependencies                      | footprint              | runtime_note                                                           |
|:------------|:----------------------------------|:-----------------------|:-----------------------------------------------------------------------|
| duckvep     | DuckDB (+ a 1.4 MB extension)     | 1.4 MB extension       | runs in any DuckDB (CLI/R/Python/WASM); transcript cache 45 MB Parquet |
| fastVEP     | none                              | ~5.5 MB binary         | standalone Rust binary                                                 |
| Ensembl VEP | Perl 5.22+, DBI, 10+ CPAN modules | ~200 MB + ~25 GB cache | interpreted Perl + species cache                                       |

## Memory & parallelism

- **Streaming `read_vcf`** bounds memory to one ~2048-row chunk: full
  GIAB **7.3 GB → 2.0 GB**.
- **Columnar Parquet transcript cache** (zstd ~45 MB,
  `read_parquet`-able); warm load ~1.1 s vs ~5.7 s GFF3 parse.
- **`Arc<str>`** shares the repeated categorical fields across output
  rows.
- **Parallel scan:** the consequence scalar is thread-safe (lock-free
  `ArcSwap` cache), so over `read_parquet` DuckDB runs it on multiple
  cores (**100% → 199% CPU**); `read_vcf` itself streams
  single-threaded.

## Ensembl-VEP concordance

Per-(variant, transcript) consequence agreement vs the **live Ensembl
VEP** (REST API), over sampled real ClinVar variants annotated with the
reference FASTA (so synonymous vs missense is exact). Both duckvep
**and** fastVEP vs Ensembl VEP:

| date       | engine  | n_variants | pairs | agree | pct |
|:-----------|:--------|-----------:|------:|------:|----:|
| 2026-06-14 | duckvep |        200 |  2677 |  2677 | 100 |
| 2026-06-14 | fastvep |        200 |  2677 |  2677 | 100 |

duckvep matches Ensembl VEP because it wraps the same engine; the table
also shows fastVEP itself vs Ensembl VEP. Dated annotation dumps
accumulate under `data/vep_dumps/<date>/annotations.parquet`.
