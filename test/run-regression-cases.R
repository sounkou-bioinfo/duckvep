#!/usr/bin/env Rscript
# Regression corpus: the SIGNIFICANT real-data (variant, transcript) mismatches vs Ensembl VEP
# that we root-caused and fixed — re-annotated and asserted against the expected (VEP-concordant)
# consequence. Needs the full cache (gff3 + fasta), gated like the GIAB benchmarks. Cases in
# test/data/regression_cases.tsv are GENERATED from the concordance dump (rows where duckvep ==
# Ensembl VEP), one+ per hard category — never hand-typed.
#
# Native R: duckdb-r v1.5.3 (release) LOADs the duckvep extension in-process (no CLI).
suppressMessages({ library(optparse); library(duckdb); library(DBI) })

root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
op <- OptionParser(usage = "%prog [options]   (regression corpus vs Ensembl VEP)")
op <- add_option(op, "--gff3",  default = file.path(root, "data/giab/GRCh38.116.gff3.gz"), help = "gene model GFF3 [%default]")
op <- add_option(op, "--fasta", default = file.path(root, "data/giab/GRCh38.primary.fa"),  help = "reference FASTA [%default]")
op <- add_option(op, "--ext",   default = Sys.getenv("DUCKVEP_EXT", file.path(root, "build/release/duckvep.duckdb_extension")), help = "extension [%default]")
opt <- parse_args(op)
tsv <- file.path(root, "test/data/regression_cases.tsv")
if (!file.exists(opt$ext))  { message("extension not built: ", opt$ext, " (make release)"); quit(status = 1) }
if (!file.exists(opt$gff3)) { message("gff3 not found: ", opt$gff3, " — run scripts/fetch-data.R (gated)"); quit(status = 77) }

con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE))
dbExecute(con, sprintf("LOAD '%s'", opt$ext))
invisible(dbGetQuery(con, sprintf("SELECT vep_load_cache('%s','%s')", opt$gff3, opt$fasta)))

# One row per case; a LEFT JOIN keeps EVERY case even when duckvep emits no row for its
# transcript, so a dropped transcript FAILS (not silently dropped from the denominator).
res <- dbGetQuery(con, sprintf("
WITH cases AS (SELECT row_number() OVER () AS rn, * FROM read_csv('%s', delim='\t', header=true)),
ann AS (
  SELECT c.rn, list_aggregate(list_sort(x.consequence),'string_agg','&') AS actual
  FROM cases c, UNNEST(vep_consequence(c.chrom::VARCHAR, c.pos::BIGINT, c.ref::VARCHAR, c.alt::VARCHAR)) u(x)
  WHERE x.transcript_id = c.transcript_id)
SELECT CASE WHEN c.expected = coalesce(ann.actual,'(no output)') THEN 'PASS' ELSE 'FAIL' END AS status,
       c.category, c.chrom||':'||c.pos AS locus, c.transcript_id, c.expected,
       coalesce(ann.actual,'(no output)') AS actual
FROM cases c LEFT JOIN ann USING (rn) ORDER BY status DESC, c.category", tsv))

for (i in seq_len(nrow(res))) with(res[i, ],
  cat(sprintf("%s|%s|%s|%s|%s|%s\n", status, category, locus, transcript_id, expected, actual)))
fails <- sum(res$status == "FAIL"); total <- nrow(res)
cat("----\n"); cat(sprintf("regression cases: %d/%d passed\n", total - fails, total))
if (fails > 0) { message(sprintf("REGRESSION: %d case(s) no longer match Ensembl VEP", fails)); quit(status = 1) }
