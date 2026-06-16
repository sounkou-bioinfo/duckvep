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
knitr::knit_engines$set(sql = duckknit::eng_duckdb)

rmarkdown::render("README.Rmd", output_file = "README.md", quiet = TRUE)
duckknit::duckknit_kill_all_sessions()

# Symbolic ALTs (e.g. `<DEL>`, `<CNV>`, `<DUP:TANDEM>`) in the rendered OUTPUT tables look like
# HTML tags, so GitHub/markdown silently strips them ([<DEL>] -> []). Escape just that
# `<TOKEN>` pattern to HTML entities, on markdown TABLE-ROW lines only (lines starting with
# `|`) — leaving HGVS `\>` (already escaped by knitr), prose, and code spans untouched.
md <- readLines("README.md", warn = FALSE)
is_row <- grepl("^\\|", md)
md[is_row] <- gsub("<([A-Za-z][A-Za-z0-9:._]*)>", "&lt;\\1&gt;", md[is_row], perl = TRUE)
writeLines(md, "README.md")
cat("rendered README.md\n")
