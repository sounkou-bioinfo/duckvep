#!/usr/bin/env bash
# Seamless duckvep <-> duckhts integration: BOTH DuckDB extensions in one session,
# composed by the optimizer. duckhts (a DuckDB *community extension*) provides the
# rich format readers (read_gff with a queryable attribute MAP, read_vcf/read_bcf,
# tabix region scans); duckvep provides the VEP engine (vep_consequence /
# vep_annotate / HGVS). They are independent extensions sharing nothing but the
# DuckDB ABI and SQL — the integration is just a JOIN, the point of being
# DuckDB-native.
#
# Usage: scripts/duckhts_integration_demo.sh [duckhts.duckdb_extension]
#   no arg  -> INSTALL duckhts FROM community (needs network once)
#   arg     -> LOAD a local duckhts.duckdb_extension build
set -euo pipefail
ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
DUCKDB="$ROOT/.tools/duckdb"
DUCKVEP="$ROOT/build/release/duckvep.duckdb_extension"
if [ $# -ge 1 ]; then
  LOAD_DUCKHTS="LOAD '$1';"
else
  LOAD_DUCKHTS="INSTALL duckhts FROM community; LOAD duckhts;"
fi
GFF="$ROOT/data/giab/GRCh38.112.gff3.gz"
FA="$ROOT/data/giab/chr17.fa"
# A tabix-indexed GFF for duckhts region scans (bgzip + .tbi). Built on demand.
GFF_TBX=/tmp/chr17.gff3.gz
if [ ! -f "$GFF_TBX.tbi" ]; then
  zcat "$GFF" | awk -F'\t' '/^#/ || $1=="17"' | sort -k1,1 -k4,4n | bgzip > "$GFF_TBX"
  tabix -p gff "$GFF_TBX"
fi

"$DUCKDB" -unsigned -c "
LOAD '$DUCKVEP';
$LOAD_DUCKHTS
SELECT vep_load_cache('$GFF', '$FA');

-- duckhts read_gff (rich attribute MAP, tabix region scan) JOINED to duckvep
-- vep_consequence (consequence + HGVS g./c./p.) — one session, one SQL plan.
WITH csq AS (
  SELECT c.gene_id, c.transcript_id, c.consequence, c.impact, c.hgvsc, c.hgvsp
  FROM UNNEST(vep_consequence('17', 43124090, 'A', 'G')) AS u(c)
  WHERE c.canonical
), genes AS (
  SELECT attributes_map['gene_id']     AS gene_id,
         attributes_map['Name']        AS gene_name,
         attributes_map['description'] AS description
  FROM read_gff('$GFF_TBX', region := '17:43044295-43125483', attributes_map := true)
  WHERE feature = 'gene'
)
SELECT g.gene_name, csq.consequence, csq.impact, csq.hgvsc, csq.hgvsp, g.description
FROM csq JOIN genes g USING (gene_id);
"
