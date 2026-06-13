#!/usr/bin/env bash
# Multi-organism throughput harness (fastVEP-paper style): annotate <vcf> against
# <gff3> [+fasta] with duckvep, append (organism, transcripts, variants, time,
# v/s) to benchmarks/data/throughput.csv. Re-run per organism for the table.
#
# Usage: scripts/benchmark.sh "<organism>" <assembly> <gff3(.gz)> <vcf(.gz)> [fasta]
set -euo pipefail
ORG=$1; ASM=$2; GFF3=$3; VCF=$4; FASTA=${5:-}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
EXT=$ROOT/build/release/duckvep.duckdb_extension
DUCKDB=${DUCKVEP_DUCKDB:-$ROOT/.tools/duckdb}

# Warm the transcript cache once (steady state).
$DUCKDB -unsigned -c "LOAD '$EXT'; SELECT vep_load_cache('$GFF3','$FASTA');" >/dev/null 2>&1
ntx=$($DUCKDB -noheader -list -c "SELECT count(*) FROM '${GFF3}.transcripts.parquet'" 2>/dev/null | tail -1)
nvar=$(zcat -f "$VCF" | grep -vc '^#')
t=$( { /usr/bin/time -f "%e" "$DUCKDB" -unsigned -c \
  "LOAD '$EXT'; SELECT vep_load_cache('$GFF3','$FASTA');
   SELECT count(*) FROM (SELECT UNNEST(vep_consequence(v.chrom,v.pos,v.ref,a.alt))
                         FROM read_vcf('$VCF') v, UNNEST(v.alt) a(alt));" >/dev/null; } 2>&1 )
vs=$(python3 -c "print(round($nvar/$t))")
echo "$ORG,$ASM,$ntx,$nvar,custom,$t,consequence" >> "$ROOT/benchmarks/data/throughput.csv"
echo "$ORG ($ASM): $nvar variants / ${t}s = $vs v/s ($ntx transcripts)"
