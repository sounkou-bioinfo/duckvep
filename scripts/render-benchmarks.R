#!/usr/bin/env Rscript
# Render benchmarks/results.Rmd -> benchmarks/results.md. Tables are read directly
# from the recorded CSVs under benchmarks/data/ (no DuckDB needed). Run from repo
# root: `Rscript scripts/render-benchmarks.R` (or `make benchmarks`).

rmarkdown::render("benchmarks/results.Rmd", output_file = "results.md",
                  knit_root_dir = normalizePath("."), quiet = TRUE)
cat("rendered benchmarks/results.md\n")
