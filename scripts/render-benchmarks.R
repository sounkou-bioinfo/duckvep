#!/usr/bin/env Rscript
# Render benchmarks/results.Rmd -> benchmarks/results.md with live tables from the
# recorded data (timings.csv, concordance_log.csv) via duckknit. Run from repo
# root: `Rscript scripts/render-benchmarks.R` (or `make benchmarks`).

options(duckknit.duckdb = normalizePath("scripts/duckdb-unsigned"))
knitr::knit_engines$set(duckdb = duckknit::eng_duckdb)
rmarkdown::render("benchmarks/results.Rmd", output_file = "results.md",
                  knit_root_dir = normalizePath("."), quiet = TRUE)
duckknit::duckknit_kill_all_sessions()
cat("rendered benchmarks/results.md\n")
