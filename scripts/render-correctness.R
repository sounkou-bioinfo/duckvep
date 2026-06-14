#!/usr/bin/env Rscript
# Render the correctness docs from their .Rmd sources (make correctness). Every
# table/number is read directly from a recorded CSV (no hand-typed figures, no
# drift): correctness/data/ for the concordance, correctness/cache-build/data/ for
# the cache stats. Run from repo root.

root <- normalizePath(".")
rmarkdown::render("correctness/correctness.Rmd", output_file = "correctness.md",
                  knit_root_dir = root, quiet = TRUE)
rmarkdown::render("correctness/cache-build/README.Rmd", output_file = "README.md",
                  knit_root_dir = root, quiet = TRUE)
cat("rendered correctness/correctness.md and correctness/cache-build/README.md\n")
