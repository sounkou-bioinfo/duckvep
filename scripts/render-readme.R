#!/usr/bin/env Rscript
# Render README.Rmd -> README.md with live duckvep SQL executed via duckknit
# (rundel/duckknit). Run from the repo root: `Rscript scripts/render-readme.R`
# (or `make readme`). Requires a debug build: `make debug` first.

stopifnot(file.exists("README.Rmd"))
if (!file.exists("build/debug/duckvep.duckdb_extension")) {
  stop("build the extension first: make debug")
}

# duckknit drives the DuckDB CLI; point it at our wrapper that injects
# -unsigned so the locally built (unsigned) extension can be LOADed.
options(duckknit.duckdb = normalizePath("scripts/duckdb-unsigned"))
knitr::knit_engines$set(duckdb = duckknit::eng_duckdb)

rmarkdown::render("README.Rmd", output_file = "README.md", quiet = TRUE)
duckknit::duckknit_kill_all_sessions()
cat("rendered README.md\n")
