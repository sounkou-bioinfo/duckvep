
<!-- README.md is generated from README.Rmd — edit the .Rmd and run `make readme`.
     SQL blocks below are executed live against the built extension via duckknit. -->

# duckvep — VEP-style variant consequences in DuckDB

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

duckvep is a loadable DuckDB extension, written in Rust, that implements
Ensembl-VEP-style consequence prediction and HGVS as SQL functions. It
reads genomics formats with
[noodles](https://github.com/zaeleus/noodles) and treats annotation
databases as ordinary Parquet/DuckDB tables that the query optimizer
joins. The consequence engine is a Rust port of
[fastVEP](https://github.com/Huang-lab/fastVEP) with accuracy patches
documented in [docs/PATCHES.md](docs/PATCHES.md). It is under active
development; local builds require DuckDB v1.5.3 with unsigned extension
loading.

Concordance with Ensembl VEP is measured rather than assumed. On a
controlled `--gff` setup, where VEP, duckvep and the vendored fastVEP
read the same gene model and results are keyed to the original input
variant, duckvep diverges from VEP on 161 of roughly 1.44M
variant/transcript pairs (N=50,000 ClinVar variants): 132 with a
different consequence on a shared pair plus 29 emission differences. Of
the shared-pair disagreements, 39 are duckvep-specific (fastVEP matches
VEP where duckvep does not); fastVEP itself diverges on 5,089 pairs. The
remaining duckvep cases are boundary indels where VEP 3’-shifts the
allele before calling consequence. The reports under
[`conformance/`](conformance/) are stratified by consequence term,
variant type and length bin with exact 95% confidence intervals; the
framework is a differential falsifier and an empirical bound, not a
proof. Known gaps are tracked in
[`correctness/`](correctness/correctness.md) and
[`conformance/`](conformance/).

Beyond single-variant calls, `vep_haplotype_consequence` applies phased
variants together, in the style of bcftools-`csq` and haplosaurus, so a
silent SNV becomes missense when phased with its neighbour; this path is
experimental. Annotation sources — gnomAD, ClinVar, conservation, a
custom BED file — are Parquet or DuckDB tables joined on
`normalize_variant`’s key and composed in SQL through a `sources.sql`
manifest and an `annotate_vcf()` macro, rather than a command-line flag
per source, and duckvep can share a DuckDB session with the duckhts
extension. The transcript cache is built from Ensembl MySQL dumps and
inherits the curated MANE, `cds_start_NF` and selenocysteine flags. The
implementation is pure Rust; a WASM target exists but is experimental
and not yet covered by CI.

See [docs/DESIGN.md](docs/DESIGN.md), [NEWS.md](NEWS.md),
[docs/PATCHES.md](docs/PATCHES.md) and
[correctness/](correctness/correctness.md).

## Build

``` sh
make debug      # builds build/debug/duckvep.duckdb_extension (native)
make test       # runs the SQL test suite
```

A WASM target exists, but browser/DuckDB-WASM packaging is experimental
and not yet covered by CI.

## Load it

``` sql
LOAD 'build/debug/duckvep.duckdb_extension';
```

## `read_vcf` — VCF/BCF as a SQL table

One row per variant. `alt` and `filter` are lists; `end_pos` carries the
variant interval (`INFO/END` for SV/CNV, else `pos + len(ref) - 1`).

``` sql
SELECT chrom, pos, end_pos, ref, alt, qual, filter
FROM read_vcf('test/data/sites.vcf')
LIMIT 5;
```

| chrom |      pos |  end_pos | ref | alt   | qual | filter   |
|-------|---------:|---------:|-----|-------|-----:|----------|
| 17    | 43124090 | 43124090 | A   | \[G\] | 30.0 | \[PASS\] |
| 17    | 43106500 | 43106500 | T   | \[C\] | 30.0 | \[PASS\] |
| 17    | 43125300 | 43125300 | C   | \[T\] | 30.0 | \[PASS\] |
| 17    | 43120000 | 43120000 | C   | \[T\] | 30.0 | \[PASS\] |
| 17    | 43043000 | 43043000 | T   | \[C\] | 30.0 | \[PASS\] |

### Structural variants

Symbolic (`<DEL>`, `<CNV>`), breakend, and multiallelic alleles survive
as list elements, and `end_pos` makes interval filters work:

``` sql
SELECT pos, end_pos, ref, alt, end_pos - pos AS span
FROM read_vcf('test/data/sv.vcf')
WHERE end_pos - pos > 100
ORDER BY pos;
```

|   pos | end_pos | ref | alt       | span |
|------:|--------:|-----|-----------|-----:|
|  5000 |    8000 | N   | \[&lt;DEL&gt;\] | 3000 |
| 12000 |   20000 | A   | \[&lt;CNV&gt;\] | 8000 |

### Multi-sample, phased genotypes

`gt` is the positional per-sample genotype list (phasing preserved as
`0|1`):

``` sql
SELECT pos, ref, alt, gt
FROM read_vcf('test/data/ms.vcf');
```

|  pos | ref | alt      | gt                  |
|-----:|-----|----------|---------------------|
| 1000 | A   | \[G\]    | \[0\|1, 1/1, 0\|0\] |
| 2000 | C   | \[T, A\] | \[1\|2, 0/1, ./.\]  |

`vcf_samples()` gives the header sample order, so genotypes explode to
tidy per-sample rows on demand — no re-annotation, because annotation is
site-wise:

``` sql
SELECT v.pos, s.sample, g.gt
FROM read_vcf('test/data/ms.vcf') v,
     UNNEST(v.gt) WITH ORDINALITY AS g(gt, idx)
     JOIN vcf_samples('test/data/ms.vcf') s USING (idx)
ORDER BY v.pos, s.idx;
```

|  pos | sample | gt   |
|-----:|--------|------|
| 1000 | NA1    | 0\|1 |
| 1000 | NA2    | 1/1  |
| 1000 | NA3    | 0\|0 |
| 2000 | NA1    | 1\|2 |
| 2000 | NA2    | 0/1  |
| 2000 | NA3    | ./.  |

### Region filter

``` sql
SELECT count(*) AS n
FROM read_vcf('test/data/sv.vcf', region := 'chr1:4000-13000');
```

|   n |
|----:|
|   3 |

## Variant effect prediction

`vep_load_cache('<gene model>', '<fasta>')` builds the consequence
engine **once** into the connection’s state, so later
`vep_consequence`/`vep_annotate` calls reuse it (the result `loaded` is
just a confirmation). The gene model is either a GFF3 or duckvep’s
**columnar Parquet cache**: a GFF3 is parsed *and* written to
`<gff3>.transcripts.parquet` for next time; pass that `.parquet` to skip
the parse. The FASTA (`''` = none) is needed for exact coding calls
(synonymous vs missense).

``` sql
-- first call parses test.gff3 and writes test.gff3.transcripts.parquet;
```

``` sql
-- next time, load the cache directly: vep_load_cache('test.gff3.transcripts.parquet','')
SELECT vep_load_cache('test/data/test.gff3', '');
```

| vep_load_cache(‘test/data/test.gff3’, ’’) |
|-------------------------------------------|
| loaded                                    |

Then `vep_consequence(chrom, pos, ref, alt)` is a scan-driven scalar
returning a native `LIST<STRUCT>` — `UNNEST` it. Driven by DuckDB’s
scan, so variants come from `read_vcf`, `read_parquet`, or any relation.

Each struct also carries HGVS notation: `hgvsg` (genomic), `hgvsc`
(coding) and `hgvsp` (protein, when a FASTA is loaded). On the chr17
benchmark dataset the HGVSc and HGVSp output matches fastVEP for all
recorded comparisons (see
[`benchmarks/results.md`](benchmarks/results.md)):

``` sql
SELECT c.transcript_id, c.consequence, c.impact, c.hgvsg, c.hgvsc
FROM UNNEST(vep_consequence('17', 43124090, 'A', 'G')) AS u(c)
WHERE c.gene_symbol = 'BRCA1'
ORDER BY c.canonical DESC
LIMIT 4;
```

| transcript_id   | consequence                       | impact   | hgvsg             | hgvsc                     |
|-----------------|-----------------------------------|----------|-------------------|---------------------------|
| ENST00000357654 | \[missense_variant\]              | MODERATE | 17:g.43124090A\>G | ENST00000357654.9:c.7T\>C |
| ENST00000921914 | \[non_coding_transcript_variant\] | MODIFIER | 17:g.43124090A\>G |                           |
| ENST00000471181 | \[non_coding_transcript_variant\] | MODIFIER | 17:g.43124090A\>G |                           |
| ENST00000352993 | \[non_coding_transcript_variant\] | MODIFIER | 17:g.43124090A\>G |                           |

Or annotate a whole VCF in one call with `vep_annotate(vcf, gff3 := …)`,
one row per (variant, transcript, allele). Concordance against Ensembl
VEP is reported below:

``` sql
SELECT pos, ref, alt, gene_symbol, transcript_id, consequence, impact
FROM vep_annotate('test/data/sites.vcf', gff3 := 'test/data/test.gff3')
WHERE canonical AND impact <> 'MODIFIER'
ORDER BY pos LIMIT 5;
```

|      pos | ref | alt | gene_symbol | transcript_id   | consequence            | impact   |
|---------:|-----|-----|-------------|-----------------|------------------------|----------|
|  7675088 | C   | T   | TP53        | ENST00000269305 | \[missense_variant\]   | MODERATE |
| 43045700 | T   | C   | BRCA1       | ENST00000357654 | \[missense_variant\]   | MODERATE |
| 43106500 | T   | C   | BRCA1       | ENST00000357654 | \[missense_variant\]   | MODERATE |
| 43124090 | A   | G   | BRCA1       | ENST00000357654 | \[missense_variant\]   | MODERATE |
| 43124090 | AA  | A   | BRCA1       | ENST00000357654 | \[frameshift_variant\] | HIGH     |

### Structural variants

Symbolic SVs (`<DEL>`, `<DUP>`, `<CNV>`, `<INV>`) and breakends are
dispatched to the SV consequence predictor over their full `INFO/END`
span — `vep_annotate` reads `END` from the VCF, and the scalar takes an
END-aware 5-argument form `vep_consequence(chrom, pos, end, ref, alt)`
(drive it with `read_vcf`’s `end_pos`). A `<DEL>` covering a transcript
yields `transcript_ablation`, a partial one `feature_truncation` +
`splice_*` + `copy_number_decrease`, and a `<DUP>`
`transcript_amplification` — the Ensembl SV vocabulary:

``` sql
SELECT c.consequence, c.impact
FROM UNNEST(vep_consequence('17', 43044295, 43125483, 'N', '<DEL>')) AS u(c)
WHERE c.transcript_id = 'ENST00000357654';
```

| consequence             | impact |
|-------------------------|--------|
| \[transcript_ablation\] | HIGH   |

## Variants from any source

Because the variant table is just columns, anything DuckDB can read is a
valid variant provider — `read_vcf`, `read_parquet`, `read_csv`, an
attached DB, or a literal relation. The downstream UDFs consume the
columns, not the format:

``` sql
SELECT chrom, pos, ref, alt
FROM (VALUES ('chr1', 100, 'A', ['G']),
             ('chr2', 200, 'C', ['T', 'TA'])) AS t(chrom, pos, ref, alt);
```

| chrom | pos | ref | alt       |
|-------|----:|-----|-----------|
| chr1  | 100 | A   | \[G\]     |
| chr2  | 200 | C   | \[T, TA\] |

## Composing with duckhts

duckvep can be used with
**[duckhts](https://github.com/sounkou-bioinfo/duckhts)** (htslib-backed
format readers, a DuckDB community extension) in the same DuckDB
session: load both extensions and join their SQL outputs. duckhts brings
`read_gff` (a queryable attribute `MAP`, tabix region scans),
`read_vcf`/`read_bcf`, etc.; duckvep brings the consequence engine. The
block below is executed live in this README: gene metadata comes from
duckhts, consequence and HGVS from duckvep.

``` sql
INSTALL duckhts FROM community;
```

``` sql
LOAD duckhts;
```

``` sql
WITH csq AS (
  SELECT c.gene_id, c.consequence, c.hgvsc
  FROM UNNEST(vep_consequence('17', 43124090, 'A', 'G')) AS u(c)
  WHERE c.canonical
), genes AS (
  SELECT attributes_map['gene_id'] AS gene_id, attributes_map['Name'] AS gene_name
  FROM read_gff('test/data/test.gff3.gz', region := '17:43044295-43170245',
                attributes_map := true)
  WHERE feature = 'gene'
)
SELECT g.gene_name, csq.consequence, csq.hgvsc
FROM csq JOIN genes g USING (gene_id);
```

| gene_name | consequence          | hgvsc                     |
|-----------|----------------------|---------------------------|
| BRCA1     | \[missense_variant\] | ENST00000357654.9:c.7T\>C |

A FASTA-backed, tabix-region version (adding `hgvsp`) is in
[`scripts/duckhts_integration_demo.R`](scripts/duckhts_integration_demo.R).

## Correctness

The accuracy oracle is **Ensembl VEP**. Concordance is measured
version-matched (same Ensembl release for VEP code, VEP cache and the
gene model), keyed to the original input variant, and split by IMPACT ×
variant class so an aggregate percentage cannot hide where the misses
are (high-impact, clinically actionable sites). The full report,
rendered directly from recorded CSVs, is
**[`correctness/correctness.md`](correctness/correctness.md)**
(`make correctness`).

The table below reports error rate **per 100K pairs** vs offline Ensembl
VEP, duckvep and fastVEP side by side, generated from
`correctness/data/concordance_by_impact.csv`. SNVs have low error rates
at every impact tier; the remaining discordance is concentrated in
indels and MNVs at exon/splice boundaries.

| impact   | class | duckvep /100K | fastVEP /100K |
|:---------|:------|:--------------|:--------------|
| HIGH     | del   | 65/100K       | 8713/100K     |
| HIGH     | ins   | 383/100K      | 4930/100K     |
| HIGH     | mnv   | 0/100K        | 37179/100K    |
| HIGH     | snv   | 2/100K        | 46/100K       |
| MODERATE | del   | 0/100K        | 5042/100K     |
| MODERATE | ins   | 2304/100K     | 5707/100K     |
| MODERATE | mnv   | 0/100K        | 20957/100K    |
| MODERATE | snv   | 0/100K        | 0/100K        |
| LOW      | del   | 426/100K      | 15449/100K    |
| LOW      | ins   | 0/100K        | 11556/100K    |
| LOW      | mnv   | 0/100K        | 575/100K      |
| LOW      | snv   | 0/100K        | 17/100K       |
| MODIFIER | del   | 0/100K        | 55/100K       |
| MODIFIER | ins   | 0/100K        | 47/100K       |
| MODIFIER | mnv   | 0/100K        | 0/100K        |
| MODIFIER | snv   | 0/100K        | 5/100K        |

Error rate **per 100,000 (variant, transcript) pairs** vs Ensembl VEP
(version-matched), by impact × class, duckvep vs fastVEP — generated
from `correctness/data/concordance_by_impact.csv`.

Impact × class is too coarse, though: it hides *which* SO categories
fail and *where duckvep regresses vs fastVEP*. Splitting the
discordances by what VEP calls → what duckvep calls, flagged by whether
fastVEP matches VEP (`regression` = duckvep worse than upstream;
`shared` = inherited engine gap) — top by pair count, generated from
`correctness/data/error_transitions.csv`:

| type       | impact   | VEP calls                                                                              | duckvep calls                                                                                                                |   n |
|:-----------|:---------|:---------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------|----:|
| shared     | HIGH     | intron_variant&splice_donor_region_variant&splice_donor_variant                        | intron_variant&splice_donor_region_variant                                                                                   |  33 |
| shared     | LOW      | stop_retained_variant                                                                  | stop_lost                                                                                                                    |  28 |
| regression | MODERATE | NMD_transcript_variant&inframe_insertion&splice_region_variant                         | NMD_transcript_variant&protein_altering_variant&splice_region_variant                                                        |  20 |
| regression | MODERATE | inframe_insertion&splice_region_variant                                                | protein_altering_variant&splice_region_variant                                                                               |  19 |
| shared     | HIGH     | NMD_transcript_variant&stop_lost                                                       | NMD_transcript_variant&inframe_deletion&stop_lost                                                                            |   6 |
| shared     | HIGH     | NMD_transcript_variant&intron_variant&splice_donor_region_variant&splice_donor_variant | NMD_transcript_variant&intron_variant&splice_donor_region_variant                                                            |   6 |
| shared     | MODERATE | inframe_insertion&stop_retained_variant                                                | inframe_insertion&stop_gained                                                                                                |   3 |
| shared     | HIGH     | transcript_ablation                                                                    | intron_variant&non_coding_transcript_exon_variant&splice_acceptor_variant&splice_donor_5th_base_variant&splice_donor_variant |   2 |

The open work splits into **duckvep-specific regressions** (mainly
splice sub-term precedence on boundary indels, where VEP 3’-shifts the
allele before calling consequence) and **shared engine gaps** inherited
from the vendored engine (`frameshift&stop_gained`, `missense`→
`synonymous`). Full per-SO-term and transition tables are in
[`correctness/correctness.md`](correctness/correctness.md); the accuracy
patches over the vendored engine are in
[`docs/PATCHES.md`](docs/PATCHES.md).

## Benchmarks

Wall-clock and peak memory for duckvep and the fastVEP CLI on the
benchmark inputs, generated from `benchmarks/data/timings.csv`:

| dataset                   |  variants | tool                            | wall    | peak RSS (MB) |
|:--------------------------|----------:|:--------------------------------|:--------|--------------:|
| GIAB HG002 (whole genome) | 4,048,342 | fastVEP CLI                     | 0:19.46 |         1,113 |
| GIAB HG002 (whole genome) | 4,048,342 | duckvep (cold; parses GFF3)     | 0:26.55 |         4,810 |
| GIAB HG002 (whole genome) | 4,048,342 | duckvep (warm cache; streaming) | 0:08.70 |         1,677 |
| ClinVar chr17             |   267,534 | fastVEP CLI                     | 0:06.07 |           694 |
| ClinVar chr17             |   267,534 | duckvep (cold / parses GFF3)    | 0:07.23 |         1,370 |
| ClinVar chr17             |   267,534 | duckvep (warm cache)            | 0:01.88 |         1,276 |

Full throughput / footprint / HGVS-concordance tables, and the
methodology, are in **[`benchmarks/results.md`](benchmarks/results.md)**
(`make benchmarks`). Setup and inputs:
[`scripts/fetch-data.R`](scripts/fetch-data.R) (pinned versions).

## Roadmap

The main accuracy work is to port VEP’s `TranscriptVariationAllele`
coordinate layer (the cDNA/CDS/protein mapper, `codon()`/`peptide()`,
and the start/stop/frameshift predicates). VEP 3’-shifts an indel before
calling consequence; duckvep does not yet, which is the source of the 39
duckvep-specific boundary-indel divergences in the conformance report.
Porting that layer replaces the current windowed-peptide surrogate with
VEP’s own coordinate model and is what would turn the differential
conformance framework from a falsifier into a soundness argument.
Secondary items — multi-transcript/opposite-strand witnesses, the
haplotype path’s compound-block flush, and a native-DuckDB transcript
cache — are tracked in [docs/refinements.md](docs/refinements.md).

## License

GPL-2.0-or-later. duckvep vendors parts of
[fastVEP](https://github.com/Huang-lab/fastVEP) (Apache-2.0, compatible
with GPL ≥ 2) — see [`vendor/NOTICE.md`](vendor/NOTICE.md) and
[`docs/PATCHES.md`](docs/PATCHES.md) for the vendored crates and our
divergences.
