#!/usr/bin/env bash
# Thread-scaling / CPU-affinity benchmark for the duckvep kernel.
#
# Measures how `vep_consequence` scales with DuckDB worker threads — the
# "are we using our allotted resources?" question. For each thread count N the run
# is PINNED to N physical cores with taskset (clean affinity, no oversubscription)
# and DuckDB is told `SET threads=N`. /usr/bin/time -v gives wall clock, the
# Percent-of-CPU (direct utilization: ~N*100% means N cores saturated) and peak RSS.
#
# A serial "load-only" baseline (cache build + trivial select) is measured so the
# parallel KERNEL time can be isolated from the one-off GFF3 load (Amdahl-honest).
#
# Usage: scripts/bench-threads.sh <vcf> <gff3> [--fasta ref.fa] [--threads "1 2 4 8 16"]
# Writes: benchmarks/data/thread_scaling.csv  (recorded; rendered by results.Rmd)
set -euo pipefail
VCF=${1:?input.vcf}; GFF3=${2:?genes.gff3}; shift 2
FASTA=""
THREADS_LIST=""
while [ $# -gt 0 ]; do
  case "$1" in
    --fasta) FASTA=$2; shift 2;;
    --threads) THREADS_LIST=$2; shift 2;;
    *) echo "unknown arg $1" >&2; exit 2;;
  esac
done

ROOT=$(cd "$(dirname "$0")/.." && pwd)
DUCKDB=${DUCKVEP_DUCKDB:-$ROOT/.tools/duckdb}
EXT=${DUCKVEP_EXT:-$ROOT/build/release/duckvep.duckdb_extension}
OUT=$ROOT/benchmarks/data/thread_scaling.csv
NCPU=$(nproc)
[ -z "$THREADS_LIST" ] && THREADS_LIST="1 2 4 8 $([ "$NCPU" -ge 16 ] && echo 16) $NCPU"
THREADS_LIST=$(echo "$THREADS_LIST" | tr ' ' '\n' | awk '!seen[$0]++' | sort -n | tr '\n' ' ')

[ -f "$EXT" ] || { echo "extension not built: $EXT (run: make release)" >&2; exit 1; }
fasta_arg=""; [ -n "$FASTA" ] && fasta_arg=", '$FASTA'"
fasta_load="SELECT vep_load_cache('$GFF3'${FASTA:+, '$FASTA'});"
[ -z "$FASTA" ] && fasta_load="SELECT vep_load_cache('$GFF3', '');"

nvar=$(zcat -f "$VCF" | grep -vc '^#' || true)
echo "VCF=$VCF variants=$nvar  cores=$NCPU  threads: $THREADS_LIST" >&2

# /usr/bin/time -v field parsers
parse_wall() { awk -F': ' '/Elapsed \(wall clock\)/{n=split($2,a,":"); print (n==3)?a[1]*3600+a[2]*60+a[3]:a[1]*60+a[2]}' "$1"; }
parse_cpu()  { awk -F': ' '/Percent of CPU/{gsub("%","",$2); print $2+0}' "$1"; }
parse_rss()  { awk -F': ' '/Maximum resident set size/{print int($2/1024)}' "$1"; }

run() { # threads  query  -> echoes "wall cpu rss"
  local n=$1 q=$2 tf; tf=$(mktemp)
  local pin="0-$((n-1))"
  taskset -c "$pin" /usr/bin/time -v "$DUCKDB" -unsigned -c \
    "SET threads=$n; LOAD '$EXT'; $fasta_load $q" >/dev/null 2>"$tf" || { cat "$tf" >&2; rm -f "$tf"; exit 1; }
  echo "$(parse_wall "$tf") $(parse_cpu "$tf") $(parse_rss "$tf")"
  rm -f "$tf"
}

KERNEL_Q="SELECT count(*) FROM read_vcf('$VCF') v, UNNEST(v.alt) a(alt), UNNEST(vep_consequence(v.chrom, v.pos, v.ref, a.alt)) u(c);"
LOAD_Q="SELECT 1;"

echo "date,variants,threads,phase,wall_s,cpu_pct,rss_mb,throughput_per_s,speedup,efficiency" > "$OUT"
DATE=$(date +%F)
base_wall=""
for n in $THREADS_LIST; do
  # load-only baseline (serial GFF3 parse) at this affinity
  read -r lw lc lr <<<"$(run "$n" "$LOAD_Q")"
  # full kernel scan
  read -r fw fc fr <<<"$(run "$n" "$KERNEL_Q")"
  # isolate parallel kernel time (subtract serial load)
  kw=$(awk -v f="$fw" -v l="$lw" 'BEGIN{d=f-l; print (d>0)?d:f}')
  [ -z "$base_wall" ] && base_wall=$kw
  thr=$(awk -v v="$nvar" -v w="$kw" 'BEGIN{print (w>0)?int(v/w):0}')
  sp=$(awk -v b="$base_wall" -v w="$kw" 'BEGIN{printf "%.2f", (w>0)?b/w:0}')
  eff=$(awk -v s="$sp" -v n="$n" 'BEGIN{printf "%.2f", s/n}')
  printf "%s,%s,%s,load,%.2f,%.0f,%s,,,\n" "$DATE" "$nvar" "$n" "$lw" "$lc" "$lr" >> "$OUT"
  printf "%s,%s,%s,kernel,%.2f,%.0f,%s,%s,%s,%s\n" "$DATE" "$nvar" "$n" "$kw" "$fc" "$fr" "$thr" "$sp" "$eff" >> "$OUT"
  echo "  threads=$n  kernel_wall=${kw}s  cpu=${fc}%  rss=${fr}MB  thr=${thr}/s  speedup=${sp}x  eff=${eff}" >&2
done
echo "wrote $OUT" >&2
