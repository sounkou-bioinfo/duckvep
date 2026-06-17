#!/usr/bin/env Rscript
# Differential fuzzer over the FORMAL-tier witnesses — the shared evaluator. Runs Ensembl VEP
# (--gff, the controlled oracle) and duckvep on the witness VCF from generate_witnesses.R, joins
# on the canonical normalized key, and asserts an EXACT SO-term-SET match per (variant,transcript)
# — reported per equivalence CLASS (the witness INFO/CLASS). A class with pairs and 0 divergence
# is COVERED-and-passing; any divergence is a fuzz finding. (Statistical tier = exact CIs over a
# real corpus; this formal tier = class coverage. Same diff.) Native duckdb-r v1.5.3 loads the
# extension; VEP is a subprocess; SO-term sets sort byte-order (radix) to match DuckDB list_sort.
suppressMessages({ library(optparse); library(duckdb); library(DBI) })
options(rlang_backtrace_on_error = "none")

root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
op <- OptionParser()
op <- add_option(op, "--vcf",   default = file.path(root, "conformance/data/witnesses.vcf"))
op <- add_option(op, "--gff",   default = file.path(root, "data/giab/GRCh38.116.controlled.gff3.gz"))
op <- add_option(op, "--fasta", default = file.path(root, "data/giab/GRCh38.primary.fa"))
op <- add_option(op, "--ext",   default = Sys.getenv("DUCKVEP_EXT", file.path(root, "build/release/duckvep.duckdb_extension")))
op <- add_option(op, "--fork",  default = as.character(min(8L, parallel::detectCores())))  # VEP is the slow step; fork 8 (capped to cores)
opt <- parse_args(op)
vep_cmd <- strsplit(Sys.getenv("VEP_CMD", "conda run -n vep vep"), " ")[[1]]

# class + TARGET transcript per ORIGINAL (chrom,pos,ref,alt) from the witness INFO. The CLASS
# geometry is defined on TX only; we attach it solely to that transcript's pairs (pi P1-8).
wl <- readLines(opt$vcf); wl <- wl[!grepl("^#", wl)]
wf <- do.call(rbind, strsplit(wl, "\t"))
cls <- data.frame(opos = as.integer(wf[, 2]), oref = wf[, 4], oalt = wf[, 5],
                  class = sub(".*CLASS=([^;]*).*", "\\1", wf[, 8]),
                  tx    = sub(".*TX=([^;]*).*", "\\1", wf[, 8]), stringsAsFactors = FALSE)

con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE))
invisible(dbExecute(con, sprintf("LOAD '%s'", opt$ext)))
invisible(dbExecute(con, "INSTALL json")); invisible(dbExecute(con, "LOAD json"))  # read_json for VEP/fastVEP output
invisible(dbGetQuery(con, sprintf("SELECT vep_load_cache('%s','%s')", opt$gff, opt$fasta)))

# duckvep rows (native), keyed by ORIGINAL witness identity (v.pos,v.ref,a.alt) — NOT by
# normalize_variant. Comparing on the original input variant (the witness's ground-truth
# identity) avoids the trap (pi P0-5) where duckvep's and VEP's insertion-normalization differ
# and the SAME pair surfaces as a false "extra" on one side / "missing" on the other.
dv <- dbGetQuery(con, sprintf("
  SELECT v.pos opos, v.ref oref, a.alt oalt, c.transcript_id tx,
         list_aggregate(list_sort(c.consequence),'string_agg','&') dc
  FROM read_vcf('%s') v, UNNEST(v.alt) a(alt),
       UNNEST(vep_consequence(v.chrom,v.pos,v.ref,a.alt)) u(c)
  GROUP BY ALL", opt$vcf))

# VEP --gff on the same witnesses. Fresh temp output, hard-fail on a nonzero exit or empty
# JSON — a swallowed VEP failure would otherwise fake a clean run (pi P0-6).
vep_json <- tempfile(fileext = ".json"); unlink(vep_json)
vrc <- system2(vep_cmd[1], c(vep_cmd[-1], "-i", opt$vcf, "--gff", opt$gff, "--fasta", opt$fasta,
        "--distance", "5000", "--json", "-o", vep_json, "--fork", opt$fork,
        "--force_overwrite", "--no_stats"), stdout = FALSE, stderr = FALSE)
if (vrc != 0 || !file.exists(vep_json) || file.info(vep_json)$size == 0)
  stop(sprintf("VEP --gff failed (exit %s) or produced no JSON at %s", vrc, vep_json))
# Parse VEP's NDJSON with duckdb-r's read_json (not an R per-row loop) and UNNEST in SQL. Each
# record echoes the ORIGINAL VCF line in `input` (fields pos=2, ref=4, alt=5); key VEP rows to
# that original identity. vmap (orig -> VEP-normalized start/allele_string) bridges fastVEP (no
# `input` field, but it shares VEP's normalization) back to original.
vep <- dbGetQuery(con, sprintf("
  SELECT CAST(split_part(input,chr(9),2) AS BIGINT) opos, split_part(input,chr(9),4) oref,
         split_part(input,chr(9),5) oalt, tc.transcript_id tx,
         list_aggregate(list_sort(tc.consequence_terms),'string_agg','&') vc
  FROM read_json('%s', format='newline_delimited', sample_size=-1), UNNEST(transcript_consequences) u(tc)", vep_json))
vmap <- dbGetQuery(con, sprintf("
  SELECT DISTINCT CAST(split_part(input,chr(9),2) AS BIGINT) opos, split_part(input,chr(9),4) oref,
         split_part(input,chr(9),5) oalt, start||' '||allele_string vkey
  FROM read_json('%s', format='newline_delimited', sample_size=-1)", vep_json))

# fastVEP (the vendored engine) on the same witnesses -> term-fair "duckvep-specific" flag.
# fastvep_ok gates ALL "shared with fastVEP" claims: if fastVEP is absent or its JSON fails to
# parse, we must NOT silently fold its silence into "shared" (pi P0-4) — divergences then count
# as fastvep_unknown, not as duckvep-non-specific.
FASTVEP <- Sys.getenv("FASTVEP", file.path(root, "../DuckfastVEP/target/release/fastvep"))
fv <- NULL; fastvep_ok <- FALSE
if (file.exists(FASTVEP)) {
  fv_json <- tempfile(fileext = ".json"); unlink(fv_json)
  frc <- system2(FASTVEP, c("annotate", "-i", opt$vcf, "--gff3", opt$gff, "--fasta", opt$fasta,
                            "--output-format", "json"), stdout = fv_json, stderr = FALSE)
  if (frc == 0 && file.exists(fv_json) && file.info(fv_json)$size > 0) {
    fvn <- tryCatch(dbGetQuery(con, sprintf("       -- read + UNNEST fastVEP JSON in DuckDB
      SELECT f.start||' '||f.allele_string vkey, tc.transcript_id tx,   -- bridge key = VEP-style 'start allele_string'
             list_aggregate(list_sort(tc.consequence_terms),'string_agg','&') fc
      FROM read_json('%s', sample_size=-1) f, UNNEST(f.transcript_consequences) u(tc)", fv_json)),
      error = function(e) { message("WARN: fastVEP JSON parse failed: ", conditionMessage(e)); NULL })
    if (!is.null(fvn) && nrow(fvn)) {
      # join fastVEP's normalized vkey to vmap to recover the original variant identity
      fv <- merge(fvn, vmap, by = "vkey")[, c("opos","oref","oalt","tx","fc")]
      fastvep_ok <- TRUE
    }
  }
} else message("WARN: fastVEP binary not found (", FASTVEP, ") — 'duckvep-specific' undeterminable.")

# FULL-OUTER union of (original-variant, transcript) keys — denominator is every pair EITHER
# engine emits, so emission misses/extras count as divergence, not as silently-dropped inner-join
# rows (pi P0-3 / P1-11). All three engines are keyed to the ORIGINAL witness identity (pi P0-5).
# NA in a column = that engine did not emit the pair.
key <- function(d) paste(d$opos, d$oref, d$oalt, d$tx)
m <- unique(rbind(vep[, c("opos","oref","oalt","tx")], dv[, c("opos","oref","oalt","tx")]))
mk <- key(m)
m$vc <- vep$vc[match(mk, key(vep))]
m$dc <- dv$dc[match(mk, key(dv))]
m$fc <- if (fastvep_ok) fv$fc[match(mk, key(fv))] else NA_character_
# Join on (opos,oref,oalt,TX) so the class attaches ONLY to the target transcript's pair; every
# other overlapping transcript at the locus has different geometry and becomes '(other)' (pi P1-8).
m <- merge(m, cls, by = c("opos","oref","oalt","tx"), all.x = TRUE)
m$class[is.na(m$class)] <- "(other)"

eq <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)   # NA = "not emitted"
m$ok <- !is.na(m$vc) & !is.na(m$dc) & m$vc == m$dc
m$status <- ifelse(m$ok, "match",
             ifelse(is.na(m$dc), "duckvep_missing",      # VEP emitted, duckvep did not
             ifelse(is.na(m$vc), "duckvep_extra",        # duckvep emitted, VEP did not
                                  "term_mismatch")))      # both emitted, term sets differ
# term-fair: a divergence is duckvep-SPECIFIC only when fastVEP ran AND made VEP's exact call
# (same term set, or same emission/non-emission). Without fastVEP it is UNDETERMINABLE, never "shared".
m$dv_specific <- fastvep_ok & !m$ok & eq(m$fc, m$vc)
m$fastvep_unknown <- !fastvep_ok & !m$ok

agg <- function(d) data.frame(n = nrow(d), discordant = sum(!d$ok), dv_spec = sum(d$dv_specific))
rep <- do.call(rbind, by(m, m$class, agg)); rep$class <- rownames(rep)
rep <- rep[order(-rep$discordant, -rep$n), ]
cat(sprintf("Differential fuzz over %d witnesses -> %d (variant,transcript) pairs (union of both engines), %d classes\n",
            nrow(cls), nrow(m), length(unique(m$class))))
cat(sprintf("  duckvep ≡ VEP on %d/%d pairs; %d divergences across %d classes\n",
            sum(m$ok), nrow(m), sum(!m$ok), sum(rep$discordant > 0)))
cat(sprintf("    by status: term_mismatch=%d  duckvep_missing=%d  duckvep_extra=%d\n",
            sum(m$status == "term_mismatch"), sum(m$status == "duckvep_missing"), sum(m$status == "duckvep_extra")))
cat(if (fastvep_ok) sprintf("    %d duckvep-SPECIFIC (fastVEP≡VEP≠duckvep); rest shared-or-both-wrong with fastVEP\n", sum(m$dv_specific))
    else sprintf("    duckvep-specific UNDETERMINABLE (fastVEP unavailable) for all %d divergences\n", sum(!m$ok)))
cat(sprintf("  %-26s %6s %6s %9s\n", "class", "n", "disc", "dv_spec"))
for (i in seq_len(nrow(rep))) cat(sprintf("  %-26s %6d %6d %9d\n", rep$class[i], rep$n[i], rep$discordant[i], rep$dv_spec[i]))
d <- m[!m$ok, c("class","status","opos","oref","oalt","tx","vc","dc")]
if (nrow(d)) { cat("--- divergences (fuzz findings) ---\n")
  for (i in seq_len(min(20, nrow(d)))) with(d[i, ],
    cat(sprintf("  [%s/%s] %s %s>%s %s  VEP=%s  duckvep=%s\n", class, status, opos, oref, oalt, tx,
        ifelse(is.na(vc), "(none)", vc), ifelse(is.na(dc), "(none)", dc)))) }
write.csv(m, file.path(root, "conformance/data/fuzz_results.csv"), row.names = FALSE)
quit(status = if (sum(!m$ok) > 0) 1 else 0)
