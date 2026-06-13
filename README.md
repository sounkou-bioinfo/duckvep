
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

> Status: early. `read_vcf` and `vcf_samples` are implemented; the
> consequence / HGVS / ACMG UDFs and the Parquet annotation pipeline are
> in progress.

## Build

``` sh
make debug      # builds build/debug/duckvep.duckdb_extension (native)
make test       # runs the SQL test suite
```

The extension also builds to **WASM**, so the same readers run in
DuckDB-WASM in the browser with no server.

## Load it

``` duckdb
LOAD 'build/debug/duckvep.duckdb_extension';
```

## `read_vcf` — VCF/BCF as a SQL table

One row per variant. `alt` and `filter` are lists; `end_pos` carries the
variant interval (`INFO/END` for SV/CNV, else `pos + len(ref) - 1`).

``` duckdb
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

``` duckdb
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

``` duckdb
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

``` duckdb
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

``` duckdb
SELECT count(*) AS n
FROM read_vcf('test/data/sv.vcf', region := 'chr1:4000-13000');
```

|   n |
|----:|
|   3 |

## Variants from any source

Because the variant table is just columns, anything DuckDB can read is a
valid variant provider — `read_vcf`, `read_parquet`, `read_csv`, an
attached DB, or a literal relation. The downstream UDFs consume the
columns, not the format:

``` duckdb
SELECT chrom, pos, ref, alt
FROM (VALUES ('chr1', 100, 'A', ['G']),
             ('chr2', 200, 'C', ['T', 'TA'])) AS t(chrom, pos, ref, alt);
```

| chrom | pos | ref | alt       |
|-------|----:|-----|-----------|
| chr1  | 100 | A   | \[G\]     |
| chr2  | 200 | C   | \[T, TA\] |

## License

Apache-2.0.
