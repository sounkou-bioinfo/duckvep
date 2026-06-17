#!/usr/bin/env Rscript
# Regenerate the patch files capturing duckvep's divergence from the pristine vendored fastVEP
# (Huang-lab/fastVEP @ 785922e). One .patch per crate; applying them to the upstream src/
# reproduces vendor/<crate>/src/. This makes our changes reviewable and reproducible (not just
# prose in ../../PATCHES.md).
#
# Usage: vendor/patches/regen-patches.R [path-to-fastVEP-checkout-at-785922e]
#   default upstream: a sibling DuckfastVEP checkout at HEAD 785922e
suppressMessages(library(optparse))
pa <- parse_args(OptionParser(usage = "%prog [fastVEP-checkout]"), positional_arguments = c(0, 1))

HERE   <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
VENDOR <- normalizePath(file.path(HERE, ".."))
UP     <- file.path(if (length(pa$args)) pa$args[1] else "/root/DuckfastVEP", "crates")
BASE   <- "785922e"

for (c in c("fastvep-core", "fastvep-genome", "fastvep-cache", "fastvep-consequence", "fastvep-hgvs")) {
  up_src <- file.path(UP, c, "src"); v_src <- file.path(VENDOR, c, "src")
  if (!dir.exists(up_src)) { cat(sprintf("skip %s (no upstream src at %s)\n", c, file.path(UP, c))); next }
  # a/ = pristine upstream, b/ = our vendored copy; rewrite absolute paths to a/<crate>/src and
  # b/<crate>/src so `git apply -p1` works from the repo root. diff exits 1 when files differ.
  d <- system2("diff", c("-ruN", up_src, v_src), stdout = TRUE, stderr = FALSE)
  d <- gsub(up_src, file.path("a", c, "src"), d, fixed = TRUE)
  d <- gsub(v_src,  file.path("b", c, "src"), d, fixed = TRUE)
  patch <- file.path(HERE, paste0(c, ".patch")); writeLines(d, patch)
  cat(sprintf("  %s  (%d added lines)\n", basename(patch), sum(grepl("^\\+", d))))
}
cat(sprintf("regenerated patches vs fastVEP@%s\n", BASE))
