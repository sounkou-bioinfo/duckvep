#!/usr/bin/env Rscript
# Multi-organism throughput harness (fastVEP-paper style): annotate <vcf> against <gff3>
# [+fasta] with duckvep, append (organism, transcripts, variants, time, v/s) to
# benchmarks/data/throughput.csv. Re-run per organism for the table.
#
# Usage: scripts/benchmark.R "<organism>" <assembly> <gff3(.gz)> <vcf(.gz)> [fasta]
suppressMessages(library(optparse))
pa <- parse_args(OptionParser(usage = "%prog <organism> <assembly> <gff3> <vcf> [fasta]"),
                 positional_arguments = c(4, 5))
ORG <- pa$args[1]; ASM <- pa$args[2]; GFF3 <- pa$args[3]; VCF <- pa$args[4]
FASTA <- if (length(pa$args) >= 5) pa$args[5] else ""

root   <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
EXT    <- file.path(root, "build/release/duckvep.duckdb_extension")
DUCKDB <- Sys.getenv("DUCKVEP_DUCKDB", file.path(root, ".tools/duckdb"))
ddb <- function(sql, ...) {
  f <- tempfile(fileext = ".sql"); writeLines(sql, f); on.exit(unlink(f))
  system2(DUCKDB, c("-unsigned", ..., "-f", f), stdout = TRUE, stderr = FALSE)
}

# Warm the transcript cache once (steady state).
ddb(sprintf("LOAD '%s'; SELECT vep_load_cache('%s','%s');", EXT, GFF3, FASTA))
ntx  <- as.integer(tail(ddb(sprintf("SELECT count(*) FROM '%s.transcripts.parquet'", GFF3), "-noheader", "-list"), 1))
nvar <- length(grep("^#", readLines(gzfile(VCF)), invert = TRUE, value = TRUE))

t <- system.time(ddb(sprintf("LOAD '%s'; SELECT vep_load_cache('%s','%s');
  SELECT count(*) FROM (SELECT UNNEST(vep_consequence(v.chrom,v.pos,v.ref,a.alt))
                        FROM read_vcf('%s') v, UNNEST(v.alt) a(alt));", EXT, GFF3, FASTA, VCF)))[["elapsed"]]
vs <- round(nvar / t)
csv <- file.path(root, "benchmarks/data/throughput.csv")
cat(sprintf("%s,%s,%d,%d,custom,%.2f,consequence\n", ORG, ASM, ntx, nvar, t), file = csv, append = TRUE)
cat(sprintf("%s (%s): %d variants / %.2fs = %d v/s (%d transcripts)\n", ORG, ASM, nvar, t, vs, ntx))
