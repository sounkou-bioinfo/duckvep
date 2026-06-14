#!/usr/bin/env bash
# Build duckvep's transcript + regulatory cache from Ensembl's published MySQL
# **flat-file dumps** (pub/release-N/mysql/<db>/) — NO live MySQL server, so no
# "Server has gone away" flakiness. Each Ensembl table is a headerless TSV
# (`<table>.txt.gz`, NULLs as \N); the `<db>.sql.gz` schema gives column order.
# We load them straight into a local DuckDB with read_csv, then assemble.sql does
# pure-local joins -> columnar Parquet cache.
#
# Organism- and build-agnostic (DB name built from args):
#   GRCh38 (default):  build-cache.sh homo_sapiens 116 38
#   GRCh37:            build-cache.sh homo_sapiens 113 37     # GRCh37 is frozen; use its release
#   mouse:             build-cache.sh mus_musculus 116 39
# Quick iteration on one chromosome:  CHROM=17 build-cache.sh ...
#
# Inherits Ensembl's curated knowledge as columns — MANE tags, cds_start_NF /
# cds_end_NF / selenocysteine flags, the regulatory build. The VEP cache stays the
# oracle we validate against.
set -euo pipefail
SP=${1:-homo_sapiens}; REL=${2:-116}; BUILD=${3:-38}
CORE="${SP}_core_${REL}_${BUILD}"; FUNCGEN="${SP}_funcgen_${REL}_${BUILD}"
FTP="https://ftp.ensembl.org/pub/release-${REL}/mysql"
DIR="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$DIR/data/cache/${SP}.${REL}.${BUILD}"; RAW="$OUT/raw"; DB="$OUT/ensembl.duckdb"
DUCKDB="${DUCKVEP_DUCKDB:-$DIR/.tools/duckdb}"
CHROMFILTER=""; [ -n "${CHROM:-}" ] && CHROMFILTER="AND sr.name='${CHROM}'"
mkdir -p "$RAW"; rm -f "$DB"

fetch() { [ -f "$RAW/$2" ] || curl -fsSL "$1" -o "$RAW/$2"; }       # cache downloads
schema()  { fetch "$FTP/$1/$1.sql.gz" "$1.sql.gz"; }
# Quoted, comma-separated column names (in declared order) for a table, parsed
# from the gzipped CREATE TABLE in the schema. Column lines start with `  ` + a
# backtick-quoted name; the table block ends at a line beginning with `)`.
names_of() { zcat "$RAW/$1.sql.gz" | awk -v t="$2" '
  $0 ~ ("CREATE TABLE `" t "`") { f=1; next }
  f && substr($1,1,1) == "`" { c=$1; gsub(/`/,"",c); printf "%s%c%s%c", (n++?",":""), 39, c, 39 }
  f && substr($0,1,1) == ")" { exit }'; }

# Download <db>/<table>.txt.gz and load as a local DuckDB table (read_csv).
load_table() {                       # $1=db  $2=remote table  $3=local name
  local db="$1" rt="$2" name="${3:-$2}"
  local file="${db}__${rt}.txt.gz"
  schema "$db"
  fetch "$FTP/$db/$rt.txt.gz" "$file"               # db-prefixed: no cross-DB collisions
  local names; names=$(names_of "$db" "$rt")
  "$DUCKDB" -unsigned "$DB" -c \
    "CREATE OR REPLACE TABLE $name AS
       SELECT * FROM read_csv('$RAW/$file', delim='\t', header=false,
                              nullstr='\N', names=[$names], ignore_errors=true);" >/dev/null
  echo "   + $name ($(du -h "$RAW/$file" | cut -f1))"
}

echo ">> [1/2] loading Ensembl dumps for $CORE + $FUNCGEN (server-free)"
# Funcgen regulatory features reference the CORE seq_region id space directly
# (verified: rf.seq_region_id 131550 -> core chr '1'), so no funcgen seq_region.
for t in gene transcript exon exon_transcript translation seq_region coord_system \
         attrib_type transcript_attrib translation_attrib xref; do load_table "$CORE" "$t"; done
load_table "$FUNCGEN" feature_type
load_table "$FUNCGEN" regulatory_feature

echo ">> [2/2] assembling columnar cache -> $OUT/*.parquet"
sed -e "s|@OUT@|$OUT|g" -e "s|@CHROMFILTER@|$CHROMFILTER|g" \
    "$DIR/correctness/cache-build/assemble.sql" | "$DUCKDB" -unsigned "$DB"
echo ">> done: $OUT"
