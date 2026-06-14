#!/usr/bin/env Rscript
# Render correctness/correctness.Rmd -> correctness/correctness.md. The table is
# read directly from the recorded CSV under correctness/data/ (no DuckDB needed).
# Run from repo root: `Rscript scripts/render-correctness.R` (or `make correctness`).

rmarkdown::render("correctness/correctness.Rmd", output_file = "correctness.md",
                  knit_root_dir = normalizePath("."), quiet = TRUE)
cat("rendered correctness/correctness.md\n")
