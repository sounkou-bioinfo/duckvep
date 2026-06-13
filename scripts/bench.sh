#!/usr/bin/env bash
# End-to-end throughput: duckvep vs fastVEP on a VCF + GFF3 (+ optional FASTA).
# Uses hyperfine if available, else `time`. For full-size GIAB:
#   scripts/fetch-giab.sh && scripts/bench.sh data/giab/HG002.vcf.gz \
#       data/giab/GRCh38.112.gff3.gz --fasta data/giab/GRCh38.fa
set -euo pipefail
VCF=${1:?input.vcf}; GFF3=${2:?genes.gff3}; shift 2
FASTA=""; [ "${1:-}" = "--fasta" ] && { FASTA=$2; shift 2; }

ROOT=$(cd "$(dirname "$0")/.." && pwd)
DUCKDB=${DUCKVEP_DUCKDB:-$ROOT/.tools/duckdb}
EXT=$ROOT/build/debug/duckvep.duckdb_extension
fasta_arg=""; [ -n "$FASTA" ] && fasta_arg=", fasta := '$FASTA'"

DUCKVEP_CMD="$DUCKDB -unsigned -c \"LOAD '$EXT'; SELECT count(*) FROM vep_annotate('$VCF', gff3 := '$GFF3'$fasta_arg);\""
FASTVEP_CMD="cargo run -q --release -p fastvep-cli --manifest-path $ROOT/../DuckfastVEP/Cargo.toml -- annotate -i $VCF --gff3 $GFF3 ${FASTA:+--fasta $FASTA} --output-format tab >/dev/null"

if command -v hyperfine >/dev/null; then
  hyperfine --warmup 1 --runs 3 -n duckvep "$DUCKVEP_CMD" -n fastvep "$FASTVEP_CMD"
else
  echo "(hyperfine not found; single timed run each)"
  echo "== duckvep =="; eval "time $DUCKVEP_CMD"
  echo "== fastvep =="; eval "time $FASTVEP_CMD"
fi
