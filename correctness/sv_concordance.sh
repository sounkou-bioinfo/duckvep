#!/usr/bin/env bash
# Structural-variant concordance vs Ensembl VEP — the SV path was never checked against
# the oracle ("Potemkin"). This runs symbolic SVs (<DEL>/<DUP>/<INV>) over a known gene
# through BOTH VEP (--gff, same gene model) and duckvep's END-aware vep_consequence, and
# diffs the consequence sets. Requires the controlled GFF + FASTA (scripts/fetch-data.sh).
#
# Usage: correctness/sv_concordance.sh
# Current state (recorded 2026-06-16, TP53 ENST00000269305) — 4/4 MATCH:
#   DEL whole-gene -> transcript_ablation         : MATCH
#   DUP whole-gene -> transcript_amplification     : MATCH
#   DEL partial    -> coding_sequence_variant&feature_truncation&intron_variant : MATCH
#                    (<DEL> classifies as Deletion not CopyNumberLoss, so no
#                     copy_number_decrease; splice heuristic replaced by intron-overlap)
#   INV partial    -> 5_prime_UTR_variant&coding_sequence_variant&intron_variant : MATCH
#                    (no feature_truncation — an inversion does not truncate; reuses the
#                     coding/intron/UTR interval predicates over the SV span)
# NOTE: only a 4-case TP53 harness — get GIAB SV Tier1 coverage before further type changes.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DUCKDB=${DUCKVEP_DUCKDB:-$ROOT/.tools/duckdb}
EXT=${DUCKVEP_EXT:-$ROOT/build/release/duckvep.duckdb_extension}
GFF3=${DUCKVEP_GFF3:-$ROOT/data/giab/GRCh38.116.gff3.gz}
CTRL=${DUCKVEP_CTRL_GFF:-$ROOT/data/giab/GRCh38.116.controlled.gff3.gz}
FASTA=${DUCKVEP_FASTA:-$ROOT/data/giab/GRCh38.primary.fa}
TR=${SV_TRANSCRIPT:-ENST00000269305}
[[ -f "$FASTA" && -f "$CTRL" ]] || { echo "SKIP: need controlled GFF + FASTA" >&2; exit 0; }

SVVCF=/tmp/sv_concordance.vcf
cat > "$SVVCF" <<'VCF'
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">
##INFO=<ID=END,Number=1,Type=Integer,Description="">
##ALT=<ID=DEL,Description="">
##ALT=<ID=DUP,Description="">
##ALT=<ID=INV,Description="">
##contig=<ID=17>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
17	7668000	del_whole	N	<DEL>	.	.	SVTYPE=DEL;END=7688000
17	7676000	del_part	N	<DEL>	.	.	SVTYPE=DEL;END=7676500
17	7668000	dup_whole	N	<DUP>	.	.	SVTYPE=DUP;END=7688000
17	7670000	inv	N	<INV>	.	.	SVTYPE=INV;END=7680000
VCF

conda run -n vep vep -i "$SVVCF" --gff "$CTRL" --fasta "$FASTA" --symbol --json \
  -o /tmp/sv_vep.json --force_overwrite --no_stats 2>/dev/null
[[ -s /tmp/sv_vep.json ]] || { echo "FAIL: VEP produced no output" >&2; exit 1; }
Rscript -e '
  suppressMessages(library(jsonlite)); tr <- commandArgs(TRUE)[1]
  for (l in readLines("/tmp/sv_vep.json")) { r <- fromJSON(l, simplifyVector = FALSE)
    for (t in r$transcript_consequences) if (!is.null(t$transcript_id) && t$transcript_id == tr)
      cat(if (!is.null(r$id)) r$id else "?",
          paste(sort(unlist(t$consequence_terms), method = "radix"), collapse = "&"),
          sep = "\t", fill = TRUE) }' "$TR" > /tmp/sv_vep.tsv

# duckvep reads the SAME controlled GFF as VEP (so the only free variable is the engine).
"$DUCKDB" -unsigned -noheader -list -c "
LOAD '$EXT'; SELECT vep_load_cache('$CTRL','$FASTA');
SELECT id||chr(9)||list_aggregate(list_sort(c.consequence),'string_agg','&')
FROM (VALUES ('del_whole',7668000,7688000,'<DEL>'),('del_part',7676000,7676500,'<DEL>'),
             ('dup_whole',7668000,7688000,'<DUP>'),('inv',7670000,7680000,'<INV>')) t(id,pos,e,alt),
     UNNEST(vep_consequence('17',pos,e,'N',alt)) u(c)
WHERE c.transcript_id='$TR';" 2>/dev/null | grep -E "^(del|dup|inv)" > /tmp/sv_dv.tsv
[[ $(wc -l < /tmp/sv_dv.tsv) -eq 4 ]] || { echo "FAIL: duckvep did not annotate all 4 SVs ($(wc -l < /tmp/sv_dv.tsv)/4)" >&2; exit 1; }

echo "SV concordance vs Ensembl VEP ($TR):"
match=0; total=0
while IFS=$'\t' read -r id vep; do
  dv=$(awk -F'\t' -v i="$id" '$1==i{print $2}' /tmp/sv_dv.tsv)
  total=$((total+1))
  if [[ "$vep" == "$dv" ]]; then match=$((match+1)); st="MATCH"; else st="DIFF "; fi
  printf "  %-5s %-10s VEP=%s  duckvep=%s\n" "$st" "$id" "$vep" "$dv"
done < /tmp/sv_vep.tsv
echo "  ---- $match/$total match"
