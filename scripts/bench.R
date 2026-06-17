#!/usr/bin/env Rscript
# End-to-end throughput: duckvep vs fastVEP on a VCF + GFF3 (+ optional FASTA).
# Uses hyperfine if available, else a single timed run each. For full-size GIAB:
#   scripts/fetch-data.R && scripts/bench.R data/giab/HG002.vcf.gz \
#       data/giab/GRCh38.116.gff3.gz --fasta data/giab/GRCh38.fa
suppressMessages(library(optparse))
op <- add_option(OptionParser(usage = "%prog <input.vcf> <genes.gff3> [--fasta ref.fa]"),
                 "--fasta", default = "")
pa <- parse_args(op, positional_arguments = 2)
VCF <- pa$args[1]; GFF3 <- pa$args[2]; FASTA <- pa$options$fasta

root   <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
DUCKDB <- Sys.getenv("DUCKVEP_DUCKDB", file.path(root, ".tools/duckdb"))
EXT    <- file.path(root, "build/debug/duckvep.duckdb_extension")
fasta_arg <- if (nzchar(FASTA)) sprintf(", fasta := '%s'", FASTA) else ""

duckvep_cmd <- sprintf("%s -unsigned -c \"LOAD '%s'; SELECT count(*) FROM vep_annotate('%s', gff3 := '%s'%s);\"",
                       DUCKDB, EXT, VCF, GFF3, fasta_arg)
fastvep_cmd <- sprintf("cargo run -q --release -p fastvep-cli --manifest-path %s/../DuckfastVEP/Cargo.toml -- annotate -i %s --gff3 %s %s --output-format tab >/dev/null",
                       root, VCF, GFF3, if (nzchar(FASTA)) sprintf("--fasta %s", FASTA) else "")

if (nzchar(Sys.which("hyperfine"))) {
  system2("hyperfine", c("--warmup", "1", "--runs", "3", "-n", "duckvep", duckvep_cmd,
                         "-n", "fastvep", fastvep_cmd))
} else {
  cat("(hyperfine not found; single timed run each)\n")
  cat("== duckvep ==\n"); print(system.time(system(duckvep_cmd, ignore.stdout = TRUE)))
  cat("== fastvep ==\n"); print(system.time(system(fastvep_cmd, ignore.stdout = TRUE)))
}
