
<!-- README.md is generated from README.Rmd — edit the .Rmd and run `make readme`.
     SQL blocks below are executed live against the built extension via duckknit. -->

# duckvep

**DuckDB-native variant effect prediction.** A loadable DuckDB extension
(Rust, [`duckdb-rs`](https://github.com/duckdb/duckdb-rs)) that reads
genomics formats via [noodles](https://github.com/zaeleus/noodles),
exposes the VEP consequence / HGVS / ACMG engine as SQL functions, and
treats annotation databases as plain Parquet/DuckDB tables joined by the
optimizer — instead of hand-rolled file formats. See
[DESIGN.md](DESIGN.md) for the full design and rationale.

> Status: `read_vcf`/`vcf_samples`, `vep_consequence` (scan-driven
> scalar) and `vep_annotate`, plus a columnar Parquet transcript cache,
> are implemented and Ensembl-VEP-concordant. HGVS (g./c./p.) is wired
> and 100%-concordant with fastVEP; supplementary-annotation joins are
> next.

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
(coding) and `hgvsp` (protein, when a FASTA is loaded) — at 100%
concordance with fastVEP:

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

duckvep is a DuckDB-native extension, so it composes with
**[duckhts](https://github.com/sounkou-bioinfo/duckhts)** (htslib-backed
format readers) by just loading both and writing a `JOIN` — no glue
code, no shared format. duckhts brings `read_gff` (a queryable attribute
`MAP`, tabix region scans), `read_vcf`/`read_bcf`, etc.; duckvep brings
the VEP engine. They share nothing but the DuckDB v1.5.3 ABI and SQL:

``` sql
LOAD 'build/release/duckvep.duckdb_extension';
LOAD 'duckhts.duckdb_extension';
SELECT vep_load_cache('GRCh38.112.gff3.gz', 'chr17.fa');

-- gene metadata from duckhts read_gff, consequence + HGVS from duckvep:
WITH csq AS (
  SELECT c.gene_id, c.consequence, c.impact, c.hgvsc, c.hgvsp
  FROM UNNEST(vep_consequence('17', 43124090, 'A', 'G')) AS u(c)
  WHERE c.canonical
), genes AS (
  SELECT attributes_map['gene_id'] AS gene_id,
         attributes_map['Name']    AS gene_name
  FROM read_gff('chr17.gff3.gz', region := '17:43044295-43125483',
                attributes_map := true)
  WHERE feature = 'gene'
)
SELECT g.gene_name, csq.consequence, csq.hgvsc, csq.hgvsp
FROM csq JOIN genes g USING (gene_id);
-- BRCA1 | [synonymous_variant] | ENST00000357654.9:c.7T>C | ENSP00000350283.1:p.Leu3=
```

A runnable version is in
[`scripts/duckhts_integration_demo.sh`](scripts/duckhts_integration_demo.sh).

## License

Apache-2.0.
