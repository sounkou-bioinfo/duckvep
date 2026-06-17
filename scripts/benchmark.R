#!/usr/bin/env Rscript
# Multi-organism throughput harness (fastVEP-paper style): annotate <vcf> against <gff3> [+fasta]
# with duckvep and record COLD and WARM timings, appending to benchmarks/data/throughput.csv.
#
#   COLD = first ever run for this gene model: vep_load_cache must PARSE the GFF3 and build the
#          columnar `<gff3>.transcripts.parquet`, then annotate. This is the one-time cost.
#   WARM = the parquet cache already exists, so vep_load_cache just loads it; annotate only. This
#          is the steady state you get on every subsequent VCF, which is why it's the headline —
#          but the cold cost is reported too so the one-time build is not hidden.
# read_vcf streams the VCF (constant memory), so peak RSS stays bounded regardless of VCF size.
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
PARQUET <- paste0(GFF3, ".transcripts.parquet")
ddb <- function(sql, ...) {
  f <- tempfile(fileext = ".sql"); writeLines(sql, f); on.exit(unlink(f))
  system2(DUCKDB, c("-unsigned", ..., "-f", f), stdout = TRUE, stderr = FALSE)
}
annotate_sql <- sprintf("LOAD '%s'; SELECT vep_load_cache('%s','%s');
  SELECT count(*) FROM (SELECT UNNEST(vep_consequence(v.chrom,v.pos,v.ref,a.alt))
                        FROM read_vcf('%s') v, UNNEST(v.alt) a(alt));", EXT, GFF3, FASTA, VCF)
nvar <- length(grep("^#", readLines(gzfile(VCF)), invert = TRUE, value = TRUE))

# COLD: remove any prebuilt cache so vep_load_cache parses the GFF3 from scratch, then annotate.
if (file.exists(PARQUET)) file.rename(PARQUET, paste0(PARQUET, ".bak"))
t_cold <- system.time(ddb(annotate_sql))[["elapsed"]]
ntx <- as.integer(tail(ddb(sprintf("SELECT count(*) FROM '%s'", PARQUET), "-noheader", "-list"), 1))
# WARM: the parquet now exists, so this measures steady-state annotation only.
t_warm <- system.time(ddb(annotate_sql))[["elapsed"]]

csv <- file.path(root, "benchmarks/data/throughput.csv")
for (ph in list(c("cold", t_cold), c("warm", t_warm))) {
  phase <- ph[1]; t <- as.numeric(ph[2])
  cat(sprintf("%s,%s,%d,%d,%s,%s,%.2f,consequence\n", ORG, ASM, ntx, nvar, "custom", phase, t),
      file = csv, append = TRUE)
  cat(sprintf("%s (%s) %-4s: %d variants / %.2fs = %d v/s (%d transcripts)\n",
              ORG, ASM, phase, nvar, t, round(nvar / t), ntx))
}
cat(sprintf("  one-time cache-build overhead (cold - warm): %.2fs\n", t_cold - t_warm))
