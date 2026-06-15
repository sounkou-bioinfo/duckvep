#!/usr/bin/env bash
# Generate test/data/regression_cases.tsv from a concordance dump: the SIGNIFICANT
# real-data (variant, transcript) cases where duckvep AGREES with Ensembl VEP, one+
# per hard category (the categories each accuracy fix targets). The expected value is
# VEP's call (= duckvep's, since concordant), so the corpus is "code as memory" — never
# hand-typed. Re-run after a fresh `correctness/vep_concordance.py` to refresh it.
#
# Usage: correctness/gen-regression-cases.sh [annotations.parquet] [sample.vcf]
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DUCKDB=${DUCKVEP_DUCKDB:-$ROOT/.tools/duckdb}
EXT=${DUCKVEP_EXT:-$ROOT/build/release/duckvep.duckdb_extension}
DUMP=${1:-$(ls -t "$ROOT"/data/vep_dumps/*/annotations.parquet | head -1)}
SAMPLE=${2:-/tmp/sample.vcf}
OUT=$ROOT/test/data/regression_cases.tsv

"$DUCKDB" -unsigned -csv -c "
LOAD '$EXT';
CREATE TEMP TABLE keyed AS
  SELECT v.chrom, v.pos AS opos, v.ref AS oref, a.alt AS oalt,
    normalize_variant(v.pos,v.ref,a.alt).pos AS npos,
    CASE WHEN normalize_variant(v.pos,v.ref,a.alt).ref='' THEN '-' ELSE normalize_variant(v.pos,v.ref,a.alt).ref END AS nref,
    CASE WHEN normalize_variant(v.pos,v.ref,a.alt).alt='' THEN '-' ELSE normalize_variant(v.pos,v.ref,a.alt).alt END AS nalt
  FROM read_vcf('$SAMPLE') v, UNNEST(v.alt) a(alt);
COPY (
  WITH dd AS (SELECT pos,ref,alt,transcript_id,consequence dc FROM read_parquet('$DUMP') WHERE source='duckvep'),
       vv AS (SELECT pos,ref,alt,transcript_id,consequence vc FROM read_parquet('$DUMP') WHERE source='vep'),
  j AS (
    SELECT k.chrom, k.opos, k.oref, k.oalt, dd.transcript_id, dd.dc,
      CASE
        WHEN dd.dc LIKE '%frameshift_variant&stop_gained%' THEN 'insertion_stop_gained'
        WHEN dd.dc='start_lost' AND k.chrom='MT' THEN 'mito_start_lost'
        WHEN dd.dc LIKE '%5_prime_UTR_variant%' AND dd.dc LIKE '%start_lost%' THEN 'utr5_start_lost_cooccur'
        WHEN dd.dc LIKE '%3_prime_UTR_variant%' AND (dd.dc LIKE '%stop_lost%' OR dd.dc LIKE '%stop_retained%') THEN 'utr3_stop_cooccur'
        WHEN (dd.dc LIKE '%splice_donor_variant%' OR dd.dc LIKE '%splice_acceptor_variant%') AND length(k.oalt) > length(k.oref) AND k.oref<>'-' THEN 'delins_alt_extends_to_splice'
        WHEN dd.dc LIKE '%protein_altering_variant%' THEN 'inframe_delins_protein_altering'
        WHEN dd.dc LIKE '%intron_variant%splice_donor_variant%' AND k.nalt='-' THEN 'boundary_del_intron_cooccur'
        WHEN dd.dc LIKE '%splice_donor_region_variant%' AND length(k.oref)=length(k.oalt) AND length(k.oref)>1 THEN 'mnv_splice_differing_regions'
        WHEN dd.dc LIKE '%splice_polypyrimidine%' AND length(k.oref)<length(k.oalt) THEN 'insertion_polypyrimidine'
        WHEN dd.dc='stop_gained' THEN 'snv_stop_gained'
        WHEN dd.dc='start_lost' THEN 'nuclear_start_lost'
        ELSE NULL END AS category
    FROM keyed k
    JOIN dd ON k.npos=dd.pos AND k.nref=dd.ref AND k.nalt=dd.alt
    JOIN vv ON dd.pos=vv.pos AND dd.ref=vv.ref AND dd.alt=vv.alt AND dd.transcript_id=vv.transcript_id
    WHERE dd.dc = vv.vc
  ),
  ranked AS (SELECT *, row_number() OVER (PARTITION BY category ORDER BY opos) rn FROM j WHERE category IS NOT NULL)
  SELECT chrom, opos AS pos, oref AS ref, oalt AS alt, transcript_id, dc AS expected, category
  FROM ranked WHERE rn<=2 ORDER BY category, opos
) TO '$OUT' (DELIMITER E'\t', HEADER);
"
echo "wrote $OUT ($(($(wc -l < "$OUT") - 1)) cases)"
