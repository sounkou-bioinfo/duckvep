#!/usr/bin/env Rscript
# Seamless duckvep <-> duckhts integration: BOTH DuckDB extensions in one session, composed by
# the optimizer. duckhts (a DuckDB *community extension*) provides the rich format readers
# (read_gff with a queryable attribute MAP, read_vcf/read_bcf, tabix region scans); duckvep
# provides the VEP engine (vep_consequence / vep_annotate / HGVS). Independent extensions
# sharing nothing but the DuckDB ABI and SQL — the integration is just a JOIN.
#
# Usage: scripts/duckhts_integration_demo.R [duckhts.duckdb_extension]
#   no arg -> INSTALL duckhts FROM community (needs network once)
#   arg    -> LOAD a local duckhts.duckdb_extension build
suppressMessages(library(optparse))
pa <- parse_args(OptionParser(usage = "%prog [duckhts.duckdb_extension]"), positional_arguments = c(0, 1))

root    <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
DUCKDB  <- file.path(root, ".tools/duckdb")
DUCKVEP <- file.path(root, "build/release/duckvep.duckdb_extension")
load_duckhts <- if (length(pa$args)) sprintf("LOAD '%s';", pa$args[1]) else "INSTALL duckhts FROM community; LOAD duckhts;"
GFF <- file.path(root, "data/giab/GRCh38.116.gff3.gz")
FA  <- file.path(root, "data/giab/chr17.fa")

# A tabix-indexed GFF for duckhts region scans (bgzip + .tbi). Built on demand.
GFF_TBX <- "/tmp/chr17.gff3.gz"
if (!file.exists(paste0(GFF_TBX, ".tbi"))) {
  system2("bash", c("-c", sprintf("zcat '%s' | awk -F'\\t' '/^#/ || $1==\"17\"' | sort -k1,1 -k4,4n | bgzip > '%s' && tabix -p gff '%s'",
                    GFF, GFF_TBX, GFF_TBX)))
}

sql <- sprintf("
LOAD '%s';
%s
SELECT vep_load_cache('%s', '%s');

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
  FROM read_gff('%s', region := '17:43044295-43125483', attributes_map := true)
  WHERE feature = 'gene'
)
SELECT g.gene_name, csq.consequence, csq.impact, csq.hgvsc, csq.hgvsp, g.description
FROM csq JOIN genes g USING (gene_id);", DUCKVEP, load_duckhts, GFF, FA, GFF_TBX)

f <- tempfile(fileext = ".sql"); writeLines(sql, f); on.exit(unlink(f))
status <- system2(DUCKDB, c("-unsigned", "-f", f))
quit(status = status)
