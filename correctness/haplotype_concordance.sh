#!/usr/bin/env bash
# Haplotype concordance — the hill-climb analog of correctness/vep_concordance.py, for
# the multi-edit (haplotype) path. It needs NO new oracle: a same-length MNV is exactly a
# phased set of single-base SNVs in one transcript window, so the proven single-variant
# kernel (vep_consequence, already ~VEP-116-concordant) IS the oracle.
#
#   for each coding MNV (ref/alt same length, >=2, >=2 changed bases):
#     components = the changed bases as individual SNVs
#     ASSERT  vep_haplotype_consequence(components)  ==  coding terms of vep_consequence(MNV)
#
# A mismatch means the multi-edit CodingContext combination diverges from applying the
# same change as one variant — a real haplotype-kernel bug. Coding-term subset only
# (the haplotype path emits coding consequences; splice/intron/UTR are per-variant).
# Requires the full FASTA cache (coding sequences), like test/run-regression-cases.sh.
#
# Usage: correctness/haplotype_concordance.sh [sample.vcf]
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DUCKDB=${DUCKVEP_DUCKDB:-$ROOT/.tools/duckdb}
EXT=${DUCKVEP_EXT:-$ROOT/build/release/duckvep.duckdb_extension}
GFF3=${DUCKVEP_GFF3:-$ROOT/data/giab/GRCh38.116.gff3.gz}
FASTA=${DUCKVEP_FASTA:-$ROOT/data/giab/GRCh38.primary.fa}
SAMPLE=${1:-/tmp/sample.vcf}

if [[ ! -f "$FASTA" ]]; then
  echo "SKIP: FASTA cache not present ($FASTA) — coding haplotypes need it." >&2
  exit 0
fi

THRESHOLD=${HAPLO_CONCORDANCE_MIN:-98.0}
PCT=$("$DUCKDB" -unsigned -noheader -csv -c "
LOAD '$EXT';
SELECT vep_load_cache('$GFF3','$FASTA');
-- coding SO terms the haplotype kernel emits (everything else is per-variant)
CREATE TEMP MACRO coding_only(lst) AS list_filter(lst, x -> x IN (
  'missense_variant','synonymous_variant','stop_gained','stop_lost','stop_retained_variant',
  'start_lost','start_retained_variant','inframe_insertion','inframe_deletion',
  'frameshift_variant','protein_altering_variant','coding_sequence_variant',
  'incomplete_terminal_codon_variant'));
WITH mnvs AS (
  SELECT v.chrom AS chrom, v.pos AS pos, v.ref AS r, a.alt AS al
  FROM read_vcf('$SAMPLE') v, UNNEST(v.alt) a(alt)
  WHERE length(v.ref)=length(a.alt) AND length(v.ref)>=2 AND v.ref<>a.alt
),
comps AS (
  SELECT chrom,pos,r,al,
    string_agg((pos+(i-1))||':'||substring(r,i,1)||':'||substring(al,i,1), ';' ORDER BY i) AS comp_str
  FROM mnvs, range(1,length(r)+1) t(i)
  WHERE substring(r,i,1) <> substring(al,i,1)
  GROUP BY chrom,pos,r,al
),
mnv_csq AS (
  SELECT m.chrom,m.pos,m.r,m.al, c.transcript_id,
    list_sort(coding_only(c.consequence)) AS mnv_coding
  FROM mnvs m, UNNEST(vep_consequence(m.chrom,m.pos,m.r,m.al)) u(c)
),
joined AS (
  SELECT mc.transcript_id, mc.mnv_coding,
    list_sort(string_split(vep_haplotype_consequence(mc.chrom, mc.transcript_id, cp.comp_str),'&')) AS hap
  FROM mnv_csq mc JOIN comps cp USING(chrom,pos,r,al)
  WHERE len(mc.mnv_coding) > 0
)
SELECT
  count(*) AS tested,
  count(*) FILTER (WHERE hap = mnv_coding) AS concordant,
  count(*) FILTER (WHERE hap <> mnv_coding) AS discordant,
  round(100.0*count(*) FILTER (WHERE hap = mnv_coding)/count(*), 3) AS pct
FROM joined;
" 2>/dev/null | sed 's/\x1b\[[0-9;]*m//g' | grep -E '^[0-9]+,[0-9]+,[0-9]+,[0-9.]+$' | tail -1)

read -r TESTED CONCORDANT DISCORDANT PCT <<<"$(echo "$PCT" | tr ',' ' ')"
echo "haplotype concordance (multi-edit vs proven MNV kernel): ${CONCORDANT}/${TESTED} = ${PCT}% (${DISCORDANT} divergent)"
# Gate: the multi-edit CodingContext must stay consistent with the single-variant kernel.
awk -v p="$PCT" -v t="$THRESHOLD" 'BEGIN { exit !(p+0 >= t+0) }' || {
  echo "FAIL: haplotype concordance ${PCT}% < ${THRESHOLD}% threshold" >&2
  exit 1
}
