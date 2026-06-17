#!/usr/bin/env Rscript
# Thread-scaling / CPU-affinity benchmark for the duckvep kernel.
#
# Measures how `vep_consequence` scales with DuckDB worker threads. For each thread count N the
# run is PINNED to N physical cores with taskset (clean affinity, no oversubscription) and
# DuckDB is told `SET threads=N`. /usr/bin/time -v gives wall clock, Percent-of-CPU (direct
# utilization: ~N*100% means N cores saturated) and peak RSS. A serial "load-only" baseline
# (cache build + trivial select) isolates the parallel KERNEL time from the one-off GFF3 load.
#
# Usage: scripts/bench-threads.R <vcf> <gff3> [--fasta ref.fa] [--threads "1 2 4 8 16"]
# Writes: benchmarks/data/thread_scaling.csv  (rendered by results.Rmd)
suppressMessages(library(optparse))
op <- OptionParser(usage = "%prog <vcf> <gff3> [options]")
op <- add_option(op, "--fasta",   default = "")
op <- add_option(op, "--threads", default = "", help = "space-separated thread counts")
pa <- parse_args(op, positional_arguments = 2)
VCF <- pa$args[1]; GFF3 <- pa$args[2]; FASTA <- pa$options$fasta

root   <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
DUCKDB <- Sys.getenv("DUCKVEP_DUCKDB", file.path(root, ".tools/duckdb"))
EXT    <- Sys.getenv("DUCKVEP_EXT",   file.path(root, "build/release/duckvep.duckdb_extension"))
OUT    <- file.path(root, "benchmarks/data/thread_scaling.csv")
NCPU   <- as.integer(system2("nproc", stdout = TRUE))
if (!file.exists(EXT)) { message(sprintf("extension not built: %s (run: make release)", EXT)); quit(status = 1) }

threads <- if (nzchar(pa$options$threads)) {
  as.integer(strsplit(trimws(pa$options$threads), "\\s+")[[1]])
} else c(1, 2, 4, 8, if (NCPU >= 16) 16, NCPU)
threads <- sort(unique(threads[threads >= 1]))

fasta_load <- sprintf("SELECT vep_load_cache('%s', '%s');", GFF3, FASTA)
nvar <- length(grep("^#", readLines(gzfile(VCF)), invert = TRUE, value = TRUE))
message(sprintf("VCF=%s variants=%d  cores=%d  threads: %s", VCF, nvar, NCPU, paste(threads, collapse = " ")))

# Run `query` pinned to N cores under /usr/bin/time -v; return c(wall_s, cpu_pct, rss_mb).
run <- function(n, query) {
  sqlf <- tempfile(fileext = ".sql"); tf <- tempfile()
  writeLines(sprintf("SET threads=%d; LOAD '%s'; %s %s", n, EXT, fasta_load, query), sqlf)
  on.exit(unlink(c(sqlf, tf)))
  rc <- system2("taskset", c("-c", sprintf("0-%d", n - 1), "/usr/bin/time", "-v",
                DUCKDB, "-unsigned", "-f", sqlf), stdout = FALSE, stderr = tf)
  ln <- readLines(tf)
  if (rc != 0) { cat(ln, sep = "\n", file = stderr()); quit(status = 1) }
  grab <- function(pat) sub(".*: ", "", grep(pat, ln, value = TRUE)[1])
  el <- strsplit(grab("Elapsed \\(wall clock\\)"), ":")[[1]]; el <- as.numeric(el)
  wall <- if (length(el) == 3) el[1]*3600 + el[2]*60 + el[3] else el[1]*60 + el[2]
  c(wall = wall, cpu = as.numeric(sub("%", "", grab("Percent of CPU"))),
    rss = floor(as.numeric(grab("Maximum resident set size")) / 1024))
}

KERNEL_Q <- sprintf("SELECT count(*) FROM read_vcf('%s') v, UNNEST(v.alt) a(alt), UNNEST(vep_consequence(v.chrom, v.pos, v.ref, a.alt)) u(c);", VCF)
LOAD_Q   <- "SELECT 1;"

rows <- list(); base_wall <- NA_real_; DATE <- format(Sys.Date())
for (n in threads) {
  l <- run(n, LOAD_Q)            # load-only baseline (serial GFF3 parse) at this affinity
  f <- run(n, KERNEL_Q)          # full kernel scan
  kw <- if (f["wall"] - l["wall"] > 0) f["wall"] - l["wall"] else f["wall"]  # isolate parallel kernel time
  if (is.na(base_wall)) base_wall <- kw
  thr <- if (kw > 0) floor(nvar / kw) else 0
  sp  <- if (kw > 0) base_wall / kw else 0
  eff <- sp / n
  rows[[length(rows)+1]] <- sprintf("%s,%d,%d,load,%.2f,%.0f,%.0f,,,", DATE, nvar, n, l["wall"], l["cpu"], l["rss"])
  rows[[length(rows)+1]] <- sprintf("%s,%d,%d,kernel,%.2f,%.0f,%.0f,%d,%.2f,%.2f", DATE, nvar, n, kw, f["cpu"], f["rss"], thr, sp, eff)
  message(sprintf("  threads=%d  kernel_wall=%.2fs  cpu=%.0f%%  rss=%.0fMB  thr=%d/s  speedup=%.2fx  eff=%.2f",
                  n, kw, f["cpu"], f["rss"], thr, sp, eff))
}
writeLines(c("date,variants,threads,phase,wall_s,cpu_pct,rss_mb,throughput_per_s,speedup,efficiency",
             unlist(rows)), OUT)
message(sprintf("wrote %s", OUT))
