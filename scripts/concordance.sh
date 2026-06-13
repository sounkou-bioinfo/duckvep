#!/usr/bin/env bash
# Per-(variant, transcript) consequence concordance: duckvep vs a reference
# annotator (fastVEP CLI by default; Ensembl VEP when --vep <tsv> is given).
#
# Usage:
#   scripts/concordance.sh <input.vcf> <genes.gff3> [--fasta ref.fa] [--vep vep.tsv]
#
# Prints an agreement summary and writes per-row diffs to concordance.out/.
set -euo pipefail

VCF=${1:?input.vcf}; GFF3=${2:?genes.gff3}; shift 2
FASTA=""; VEP_TSV=""
while [ $# -gt 0 ]; do case "$1" in
  --fasta) FASTA=$2; shift 2;;
  --vep)   VEP_TSV=$2; shift 2;;
  *) echo "unknown arg $1" >&2; exit 2;;
esac; done

ROOT=$(cd "$(dirname "$0")/.." && pwd)
DUCKDB=${DUCKVEP_DUCKDB:-$ROOT/.tools/duckdb}
EXT=$ROOT/build/debug/duckvep.duckdb_extension
OUT=concordance.out; mkdir -p "$OUT"

# duckvep: one row per (variant, transcript) with the most-severe SO term.
fasta_arg=""; [ -n "$FASTA" ] && fasta_arg=", fasta := '$FASTA'"
"$DUCKDB" -unsigned -c "
LOAD '$EXT';
COPY (
  SELECT chrom||':'||pos||':'||ref||'>'||alt AS variant, transcript_id,
         list_aggregate(list_sort(consequence), 'string_agg', '&') AS csq
  FROM vep_annotate('$VCF', gff3 := '$GFF3'$fasta_arg)
  ORDER BY variant, transcript_id
) TO '$OUT/duckvep.tsv' (HEADER, DELIMITER '\t');
"

# Reference: fastVEP CLI tab output -> (variant, transcript, csq), unless an
# Ensembl VEP tsv was supplied.
if [ -n "$VEP_TSV" ]; then
  cp "$VEP_TSV" "$OUT/ref.raw.tsv"
else
  ( cd "$ROOT/../DuckfastVEP" 2>/dev/null || cd "$ROOT"
    cargo run -q -p fastvep-cli -- annotate -i "$VCF" --gff3 "$GFF3" \
      ${FASTA:+--fasta "$FASTA"} --output-format tab ) > "$OUT/ref.raw.tsv" 2>/dev/null
fi

# Drop VEP '##' metadata lines, keep the '#Uploaded_variation' header.
grep -v '^##' "$OUT/ref.raw.tsv" > "$OUT/ref.clean.tsv"

# Normalize the reference into (loc, allele, transcript, csq) the same way.
"$DUCKDB" -c "
COPY (
  SELECT \"Location\" AS loc, \"Allele\" AS allele, \"Feature\" AS transcript_id,
         list_aggregate(list_sort(string_split(\"Consequence\", ',')), 'string_agg', '&') AS csq
  FROM read_csv('$OUT/ref.clean.tsv', delim='\t', header=true, ignore_errors=true)
) TO '$OUT/ref.tsv' (HEADER, DELIMITER '\t');
"

# Compare on (Location, Allele, transcript) which both sides share.
"$DUCKDB" -c "
CREATE TABLE dv AS SELECT * FROM read_csv('$OUT/duckvep.tsv', delim='\t', header=true);
CREATE TABLE rf AS SELECT * FROM read_csv('$OUT/ref.tsv', delim='\t', header=true);
-- duckvep 'variant' is chrom:pos:ref>alt; ref 'loc' is chrom:pos, 'allele' is alt.
CREATE TABLE dvn AS
  SELECT split_part(variant,':',1)||':'||split_part(variant,':',2) AS loc,
         split_part(variant,'>',2) AS allele, transcript_id, csq FROM dv;
CREATE TABLE j AS
  SELECT COALESCE(d.loc,r.loc) loc, COALESCE(d.transcript_id,r.transcript_id) tx,
         d.csq dv_csq, r.csq ref_csq
  FROM dvn d FULL OUTER JOIN rf r
    ON d.loc=r.loc AND d.allele=r.allele AND d.transcript_id=r.transcript_id;
COPY (SELECT * FROM j WHERE dv_csq IS DISTINCT FROM ref_csq) TO '$OUT/diffs.tsv' (HEADER, DELIMITER '\t');
SELECT
  count(*) AS pairs,
  count(*) FILTER (WHERE dv_csq = ref_csq) AS agree,
  count(*) FILTER (WHERE dv_csq IS NULL) AS only_ref,
  count(*) FILTER (WHERE ref_csq IS NULL) AS only_duckvep,
  round(100.0*count(*) FILTER (WHERE dv_csq = ref_csq)/nullif(count(*),0),2) AS pct_agree
FROM j;
"
echo "diffs -> $OUT/diffs.tsv"
