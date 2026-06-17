#!/usr/bin/env Rscript
# Build duckvep's transcript + regulatory cache from Ensembl's published MySQL **flat-file
# dumps** (pub/release-N/mysql/<db>/) — NO live MySQL server, so no "Server has gone away"
# flakiness. Each Ensembl table is a headerless TSV (`<table>.txt.gz`, NULLs as \N); the
# `<db>.sql.gz` schema gives column order. We load them straight into a local DuckDB with
# read_csv, then assemble.sql does pure-local joins -> columnar Parquet cache.
#
# Organism- and build-agnostic (DB name built from args):
#   GRCh38 (default):  build-cache.R homo_sapiens 116 38
#   GRCh37:            build-cache.R homo_sapiens 113 37     # GRCh37 is frozen; use its release
#   mouse:             build-cache.R mus_musculus 116 39
# Quick iteration on one chromosome:  CHROM=17 build-cache.R ...
#
# Inherits Ensembl's curated knowledge as columns — MANE tags, cds_start_NF / cds_end_NF /
# selenocysteine flags, the regulatory build. The VEP cache stays the oracle we validate against.
suppressMessages(library(optparse))
pa <- parse_args(OptionParser(usage = "%prog [species] [release] [build]"), positional_arguments = c(0, 3))
a <- pa$args
SP <- if (length(a) >= 1) a[1] else "homo_sapiens"
REL <- if (length(a) >= 2) a[2] else "116"
BUILD <- if (length(a) >= 3) a[3] else "38"

CORE    <- sprintf("%s_core_%s_%s", SP, REL, BUILD)
FUNCGEN <- sprintf("%s_funcgen_%s_%s", SP, REL, BUILD)
FTP     <- sprintf("https://ftp.ensembl.org/pub/release-%s/mysql", REL)
HERE    <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
DIR     <- normalizePath(file.path(HERE, "../.."))
OUT     <- file.path(DIR, "data/cache", sprintf("%s.%s.%s", SP, REL, BUILD))
RAW     <- file.path(OUT, "raw"); DB <- file.path(OUT, "ensembl.duckdb")
DUCKDB  <- Sys.getenv("DUCKVEP_DUCKDB", file.path(DIR, ".tools/duckdb"))
CHROM   <- Sys.getenv("CHROM", "")
CHROMFILTER <- if (nzchar(CHROM)) sprintf("AND sr.name='%s'", CHROM) else ""
dir.create(RAW, recursive = TRUE, showWarnings = FALSE); unlink(DB)

ddb <- function(sql) { f <- tempfile(fileext = ".sql"); writeLines(sql, f); on.exit(unlink(f))
  system2(DUCKDB, c("-unsigned", DB, "-f", f), stdout = TRUE, stderr = TRUE) }
fetch <- function(url, name) {                                # cache downloads
  out <- file.path(RAW, name)
  if (!file.exists(out) && system2("curl", c("-fsSL", url, "-o", out)) != 0) stop(sprintf("fetch failed: %s", url))
  out
}
schema <- function(db) fetch(sprintf("%s/%s/%s.sql.gz", FTP, db, db), sprintf("%s.sql.gz", db))

# Quoted, comma-separated column names (declared order) for a table, parsed from the gzipped
# CREATE TABLE: column lines start with two spaces + a backtick-quoted name; the block ends at
# a line beginning with `)`.
names_of <- function(db, t) {
  ln <- readLines(gzfile(file.path(RAW, sprintf("%s.sql.gz", db))), warn = FALSE)
  start <- grep(sprintf("CREATE TABLE `%s`", t), ln, fixed = TRUE)[1]
  cols <- character(0)
  for (l in ln[(start + 1):length(ln)]) {
    if (startsWith(l, ")")) break
    if (grepl("^\\s+`", l)) cols <- c(cols, sub("^\\s+`([^`]+)`.*$", "\\1", l))
  }
  paste(sprintf("'%s'", cols), collapse = ",")
}

# Download <db>/<table>.txt.gz and load as a local DuckDB table (read_csv).
load_table <- function(db, rt, name = rt) {
  file <- sprintf("%s__%s.txt.gz", db, rt)              # db-prefixed: no cross-DB collisions
  schema(db)
  path <- fetch(sprintf("%s/%s/%s.txt.gz", FTP, db, rt), file)
  ddb(sprintf("CREATE OR REPLACE TABLE %s AS
       SELECT * FROM read_csv('%s', delim='\t', header=false, nullstr='\\N', names=[%s], ignore_errors=true);",
       name, path, names_of(db, rt)))
  cat(sprintf("   + %s (%s)\n", name, sub("\\s.*", "", system2("du", c("-h", path), stdout = TRUE))))
}

cat(sprintf(">> [1/2] loading Ensembl dumps for %s + %s (server-free)\n", CORE, FUNCGEN))
# Funcgen regulatory features reference the CORE seq_region id space directly, so no funcgen seq_region.
for (t in c("gene", "transcript", "exon", "exon_transcript", "translation", "seq_region",
            "coord_system", "attrib_type", "transcript_attrib", "translation_attrib", "xref"))
  load_table(CORE, t)
load_table(FUNCGEN, "feature_type")
load_table(FUNCGEN, "regulatory_feature")

cat(sprintf(">> [2/2] assembling columnar cache -> %s/*.parquet\n", OUT))
asm <- readLines(file.path(DIR, "correctness/cache-build/assemble.sql"))
asm <- gsub("@OUT@", OUT, asm, fixed = TRUE); asm <- gsub("@CHROMFILTER@", CHROMFILTER, asm, fixed = TRUE)
cat(ddb(paste(asm, collapse = "\n")), sep = "\n")
cat(sprintf(">> done: %s\n", OUT))
