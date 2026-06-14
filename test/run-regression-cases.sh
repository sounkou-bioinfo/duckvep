#!/usr/bin/env bash
# Regression corpus: the SIGNIFICANT real-data (variant, transcript) mismatches vs
# Ensembl VEP that we root-caused and fixed — re-annotated and asserted against the
# expected (VEP-concordant) consequence. Needs the full cache (gff3 + fasta), so it is
# gated like the GIAB benchmarks (not the bundled CI sqllogictests, which cover the
# positional splice/intron cases in test/sql/vep_splice.test).
#
# Cases in test/data/regression_cases.tsv are GENERATED from the concordance dump
# (rows where duckvep == Ensembl VEP), one+ per hard category — never hand-typed.
#
# Usage: test/run-regression-cases.sh [gff3] [fasta]
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DUCKDB=${DUCKVEP_DUCKDB:-$ROOT/.tools/duckdb}
EXT=${DUCKVEP_EXT:-$ROOT/build/release/duckvep.duckdb_extension}
GFF3=${1:-$ROOT/data/giab/GRCh38.116.gff3.gz}
FASTA=${2:-$ROOT/data/giab/GRCh38.primary.fa}
TSV=$ROOT/test/data/regression_cases.tsv

[ -f "$EXT" ] || { echo "extension not built: $EXT (make release)" >&2; exit 1; }
[ -f "$GFF3" ] || { echo "gff3 not found: $GFF3 — run scripts/fetch-data.sh (gated)" >&2; exit 77; }

out=$("$DUCKDB" -unsigned -noheader -list -c "
LOAD '$EXT'; SELECT vep_load_cache('$GFF3','$FASTA');
WITH got AS (
  SELECT c.category, c.chrom, c.pos, c.transcript_id, c.expected,
         list_aggregate(list_sort(x.consequence),'string_agg','&') AS actual
  FROM read_csv('$TSV', delim='\t', header=true) c,
       UNNEST(vep_consequence(c.chrom, c.pos, c.ref, c.alt)) u(x)
  WHERE x.transcript_id = c.transcript_id
)
SELECT CASE WHEN expected = actual THEN 'PASS' ELSE 'FAIL' END AS status,
       category, chrom||':'||pos AS locus, transcript_id, expected, actual
FROM got ORDER BY status DESC, category;")

echo "$out"
fails=$(printf '%s\n' "$out" | grep -c '^FAIL' || true)
total=$(printf '%s\n' "$out" | grep -cE '^(PASS|FAIL)' || true)
echo "----"
echo "regression cases: $((total - fails))/$total passed"
[ "$fails" -eq 0 ] || { echo "REGRESSION: $fails case(s) no longer match Ensembl VEP" >&2; exit 1; }
