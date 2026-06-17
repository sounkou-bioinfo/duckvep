#!/usr/bin/env Rscript
# Per-(variant, transcript) consequence concordance: duckvep vs a reference annotator
# (fastVEP CLI by default; Ensembl VEP when --vep <tsv> is given). R port (no bash slop) —
# the data plane is DuckDB SQL via the version-matched `.tools/duckdb` CLI (the unstable
# duckvep extension loads only in its own v1.5.3 build); fastVEP is an external process.
# Prints an agreement summary and writes per-row diffs to concordance.out/.
#
# Usage: correctness/concordance.R <input.vcf> <genes.gff3> [--fasta ref.fa] [--vep vep.tsv]
suppressMessages(library(optparse))

op <- OptionParser(usage = "%prog <input.vcf> <genes.gff3> [options]")
op <- add_option(op, "--fasta", default = "", help = "reference FASTA (enables coding consequences)")
op <- add_option(op, "--vep",   default = "", help = "Ensembl VEP tab output to compare against instead of fastVEP")
pa <- parse_args(op, positional_arguments = 2)
vcf <- pa$args[1]; gff3 <- pa$args[2]; FASTA <- pa$options$fasta; VEP_TSV <- pa$options$vep

root   <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
DUCKDB <- Sys.getenv("DUCKVEP_DUCKDB", file.path(root, ".tools/duckdb"))
EXT    <- file.path(root, "build/debug/duckvep.duckdb_extension")
OUT    <- "concordance.out"; dir.create(OUT, showWarnings = FALSE)

## SQL via the version-matched CLI; pass through a temp -f FILE so no shell re-parses the SQL.
sql_run <- function(sql, unsigned = TRUE) {
  f <- tempfile(fileext = ".sql"); writeLines(sql, f); on.exit(unlink(f))
  system2(DUCKDB, c(if (unsigned) "-unsigned", "-f", f), stdout = TRUE, stderr = TRUE)
}

# duckvep: one row per (variant, transcript) with the sorted SO-term set.
fasta_arg <- if (nzchar(FASTA)) sprintf(", fasta := '%s'", FASTA) else ""
sql_run(sprintf("LOAD '%s';
COPY (
  SELECT chrom||':'||pos||':'||ref||'>'||alt AS variant, transcript_id,
         list_aggregate(list_sort(consequence), 'string_agg', '&') AS csq
  FROM vep_annotate('%s', gff3 := '%s'%s)
  ORDER BY variant, transcript_id
) TO '%s/duckvep.tsv' (HEADER, DELIMITER '\t');", EXT, vcf, gff3, fasta_arg, OUT))

# Reference: fastVEP CLI tab output, unless an Ensembl VEP tsv was supplied.
if (nzchar(VEP_TSV)) {
  file.copy(VEP_TSV, file.path(OUT, "ref.raw.tsv"), overwrite = TRUE)
} else {
  dir <- if (dir.exists(file.path(root, "../DuckfastVEP"))) file.path(root, "../DuckfastVEP") else root
  raw <- system2("cargo", c("run", "-q", "-p", "fastvep-cli", "--", "annotate", "-i", vcf,
                 "--gff3", gff3, if (nzchar(FASTA)) c("--fasta", FASTA), "--output-format", "tab"),
                 stdout = TRUE, stderr = FALSE, env = c(sprintf("CARGO_MANIFEST_DIR=%s", dir)))
  writeLines(raw, file.path(OUT, "ref.raw.tsv"))
}

# Drop VEP '##' metadata lines, keep the '#Uploaded_variation' header.
ref_raw <- readLines(file.path(OUT, "ref.raw.tsv"))
writeLines(ref_raw[!grepl("^##", ref_raw)], file.path(OUT, "ref.clean.tsv"))

# Normalize the reference into (loc, allele, transcript, csq) the same way, then compare on
# (Location, Allele, transcript) which both sides share.
out <- sql_run(sprintf("
COPY (
  SELECT \"Location\" AS loc, \"Allele\" AS allele, \"Feature\" AS transcript_id,
         list_aggregate(list_sort(string_split(\"Consequence\", ',')), 'string_agg', '&') AS csq
  FROM read_csv('%s/ref.clean.tsv', delim='\t', header=true, ignore_errors=true)
) TO '%s/ref.tsv' (HEADER, DELIMITER '\t');
CREATE TABLE dv AS SELECT * FROM read_csv('%s/duckvep.tsv', delim='\t', header=true);
CREATE TABLE rf AS SELECT * FROM read_csv('%s/ref.tsv', delim='\t', header=true);
CREATE TABLE dvn AS
  SELECT split_part(variant,':',1)||':'||split_part(variant,':',2) AS loc,
         split_part(variant,'>',2) AS allele, transcript_id, csq FROM dv;
CREATE TABLE j AS
  SELECT COALESCE(d.loc,r.loc) loc, COALESCE(d.transcript_id,r.transcript_id) tx,
         d.csq dv_csq, r.csq ref_csq
  FROM dvn d FULL OUTER JOIN rf r
    ON d.loc=r.loc AND d.allele=r.allele AND d.transcript_id=r.transcript_id;
COPY (SELECT * FROM j WHERE dv_csq IS DISTINCT FROM ref_csq) TO '%s/diffs.tsv' (HEADER, DELIMITER '\t');
SELECT
  count(*) AS pairs,
  count(*) FILTER (WHERE dv_csq = ref_csq) AS agree,
  count(*) FILTER (WHERE dv_csq IS NULL) AS only_ref,
  count(*) FILTER (WHERE ref_csq IS NULL) AS only_duckvep,
  round(100.0*count(*) FILTER (WHERE dv_csq = ref_csq)/nullif(count(*),0),2) AS pct_agree
FROM j;", OUT, OUT, OUT, OUT, OUT))
cat(out, sep = "\n")
cat(sprintf("\ndiffs -> %s/diffs.tsv\n", OUT))
