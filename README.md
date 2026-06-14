
<!-- README.md is generated from README.Rmd — edit the .Rmd and run `make readme`.
     SQL blocks below are executed live against the built extension via duckknit. -->

# duckvep — DuckDB-native VEP: SQL genomics readers, VEP UDFs, Parquet sources

A loadable DuckDB extension (Rust,
[`duckdb-rs`](https://github.com/duckdb/duckdb-rs)) that reads genomics
formats via [noodles](https://github.com/zaeleus/noodles), exposes the
VEP consequence / HGVS / ACMG engine as SQL functions, and treats
annotation databases as plain Parquet/DuckDB tables joined by the
optimizer — instead of hand-rolled file formats. See
[docs/DESIGN.md](docs/DESIGN.md) for the design, [NEWS.md](NEWS.md) for
the changelog, and [docs/PATCHES.md](docs/PATCHES.md) for our accuracy
patches.

> **Status (v0.3.0).** Implemented: `read_vcf`/`vcf_samples`;
> `vep_consequence` (scan-driven scalar, incl. an END-aware form for
> structural variants) and `vep_annotate`; HGVS g./c./p.;
> `normalize_variant` (canonical minimal key for valid cross-annotator
> joins); a columnar Parquet transcript cache **built directly from
> Ensembl MySQL dumps** (inheriting MANE / `cds_start_NF` /
> selenocysteine / regulatory flags). Accuracy patches make duckvep
> *more* Ensembl-VEP-concordant than fastVEP (see
> [docs/PATCHES.md](docs/PATCHES.md)); concordance is reported
> version-matched and **split by impact × class** (see
> [correctness/](correctness/correctness.md)). Next: close the
> high-impact indel/MNV engine gap; supplementary-annotation joins;
> chrM/PAR/sex-chromosome correctness. Composes with the **duckhts**
> community extension in one session.

## Build

``` sh
make debug      # builds build/debug/duckvep.duckdb_extension (native)
make test       # runs the SQL test suite
```

The extension also builds to **WASM**, so the same readers run in
DuckDB-WASM in the browser with no server.

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

### Structural variants are first-class

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
|  5000 |    8000 | N   | \[<DEL>\] | 3000 |
| 12000 |   20000 | A   | \[<CNV>\] | 8000 |

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

Each struct also carries HGVS notation — `hgvsg` (genomic), `hgvsc`
(coding) and `hgvsp` (protein, when a FASTA is loaded) — exact-matching
fastVEP’s HGVS on chr17 (HGVSc 19,828/19,828, HGVSp 9,449/9,449):

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

Or annotate a whole VCF in one call with `vep_annotate(vcf, gff3 := …)`
— one row per (variant, transcript, allele), at parity with Ensembl VEP
/ fastVEP (same engine):

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

duckvep composes with
**[duckhts](https://github.com/sounkou-bioinfo/duckhts)** (htslib-backed
format readers, a DuckDB **community extension**) by just loading both
and writing a `JOIN` — no glue code, no shared format. duckhts brings
`read_gff` (a queryable attribute `MAP`, tabix region scans),
`read_vcf`/ `read_bcf`, etc.; duckvep brings the VEP engine. They share
nothing but the DuckDB ABI and SQL. This block is **executed live** in
this README — gene metadata comes from duckhts, consequence + HGVS from
duckvep:

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
[`scripts/duckhts_integration_demo.sh`](scripts/duckhts_integration_demo.sh).

## Correctness

Our accuracy oracle is **Ensembl VEP itself**. Concordance is measured
version-matched (same Ensembl release for VEP code, VEP cache and the
gene model) and is **always split by IMPACT × variant class** — an
aggregate % hides where the misses are (high-impact, clinically
actionable sites). The full report, rendered directly from recorded
CSVs, is **[`correctness/correctness.md`](correctness/correctness.md)**
(`make correctness`).

Headline — error rate **per 100K pairs** vs offline Ensembl VEP, duckvep
and fastVEP side by side, **generated from
`correctness/data/concordance_by_impact.csv`** (not hand-written): SNVs
are near-perfect at every impact tier; the open frontier is high-impact
indels/MNVs (shared with fastVEP — an engine gap, not a duckvep bug).

| impact   | class | duckvep /100K | fastVEP /100K |
|:---------|:------|:--------------|:--------------|
| HIGH     | del   | 234/100K      | 9641/100K     |
| HIGH     | ins   | 272/100K      | 6110/100K     |
| HIGH     | mnv   | 3415/100K     | 48711/100K    |
| HIGH     | snv   | 0/100K        | 0/100K        |
| MODERATE | del   | 0/100K        | 1535/100K     |
| MODERATE | ins   | 114/100K      | 1767/100K     |
| MODERATE | mnv   | 0/100K        | 28631/100K    |
| MODERATE | snv   | 0/100K        | 5/100K        |
| LOW      | del   | 352/100K      | 16454/100K    |
| LOW      | ins   | 335/100K      | 22920/100K    |
| LOW      | mnv   | 4975/100K     | 22886/100K    |
| LOW      | snv   | 2/100K        | 15/100K       |
| MODIFIER | del   | 9/100K        | 43/100K       |
| MODIFIER | ins   | 0/100K        | 220/100K      |
| MODIFIER | mnv   | 200/100K      | 200/100K      |
| MODIFIER | snv   | 3/100K        | 3/100K        |

Error rate **per 100,000 (variant, transcript) pairs** vs Ensembl VEP
(version-matched), by impact × class, duckvep vs fastVEP — generated
from `correctness/data/concordance_by_impact.csv`.

Impact × class is too coarse, though: it hides *which* SO categories
fail and *where duckvep regresses vs fastVEP*. Splitting the
discordances by what VEP calls → what duckvep calls, flagged by whether
fastVEP matches VEP (`regression` = duckvep worse than upstream;
`shared` = inherited engine gap) — top by pair count, generated from
`correctness/data/error_transitions.csv`:

| type   | impact | VEP calls                                                                                                                           | duckvep calls                                                                                                                         |   n |
|:-------|:-------|:------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------|----:|
| shared | HIGH   | frameshift_variant&splice_donor_variant                                                                                             | frameshift_variant&splice_region_variant                                                                                              |  25 |
| shared | HIGH   | stop_gained                                                                                                                         | protein_altering_variant&stop_gained                                                                                                  |  19 |
| shared | HIGH   | 3_prime_UTR_variant&intron_variant&splice_acceptor_variant&splice_donor_5th_base_variant&splice_donor_variant&stop_retained_variant | 3_prime_UTR_variant&coding_sequence_variant&intron_variant&splice_acceptor_variant&splice_donor_5th_base_variant&splice_donor_variant |  18 |
| shared | LOW    | 3_prime_UTR_variant&stop_retained_variant                                                                                           | 3_prime_UTR_variant&stop_lost                                                                                                         |  18 |
| shared | HIGH   | intron_variant&splice_acceptor_variant&stop_retained_variant                                                                        | coding_sequence_variant&intron_variant&splice_acceptor_variant                                                                        |  14 |
| shared | HIGH   | coding_sequence_variant&intron_variant&splice_acceptor_variant                                                                      | coding_sequence_variant&splice_acceptor_variant                                                                                       |  13 |
| shared | LOW    | NMD_transcript_variant&intron_variant&splice_donor_region_variant                                                                   | NMD_transcript_variant&intron_variant&splice_donor_5th_base_variant                                                                   |  12 |
| shared | LOW    | non_coding_transcript_variant&splice_region_variant                                                                                 | non_coding_transcript_exon_variant&splice_region_variant                                                                              |  12 |

So the open work splits cleanly: **regressions** (the splice sub-term
precedence our interval rewrite broke — being fixed to match Ensembl
`VariationEffect.pm`) and **shared engine gaps**
(`frameshift&stop_gained`, `missense`→`synonymous`). Full per-SO-term
and transition tables:
[`correctness/correctness.md`](correctness/correctness.md). duckvep’s
accuracy patches over the vendored engine are in
[`docs/PATCHES.md`](docs/PATCHES.md).

## Benchmarks

Wall-clock + peak memory, duckvep vs fastVEP CLI (same consequence
engine; the difference is data-engineering — caching, streaming,
columnar), **generated from `benchmarks/data/timings.csv`**:

| dataset                   |  variants | tool                            | wall    | peak RSS (MB) |
|:--------------------------|----------:|:--------------------------------|:--------|--------------:|
| GIAB HG002 (whole genome) | 4,048,342 | fastVEP CLI                     | 0:19.46 |         1,113 |
| GIAB HG002 (whole genome) | 4,048,342 | duckvep (warm cache; streaming) | 0:09.02 |         1,985 |
| ClinVar chr17             |   267,534 | fastVEP CLI                     | 0:06.07 |           694 |
| ClinVar chr17             |   267,534 | duckvep (cold / parses GFF3)    | 0:07.23 |         1,370 |
| ClinVar chr17             |   267,534 | duckvep (warm cache)            | 0:01.88 |         1,276 |

Full throughput / footprint / HGVS-concordance tables, and the
methodology, are in **[`benchmarks/results.md`](benchmarks/results.md)**
(`make benchmarks`). Setup and inputs:
[`scripts/fetch-data.sh`](scripts/fetch-data.sh) (pinned versions).

## License

GPL-2.0-or-later. duckvep vendors parts of
[fastVEP](https://github.com/Huang-lab/fastVEP) (Apache-2.0, compatible
with GPL ≥ 2) — see [`vendor/NOTICE.md`](vendor/NOTICE.md) and
[`docs/PATCHES.md`](docs/PATCHES.md) for the vendored crates and our
divergences.
