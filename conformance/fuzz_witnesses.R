#!/usr/bin/env Rscript
# Differential fuzzer over the FORMAL-tier witnesses — the shared evaluator. Runs Ensembl VEP
# (--gff, the controlled oracle) and duckvep on the witness VCF from generate_witnesses.R, joins
# on the canonical normalized key, and asserts an EXACT SO-term-SET match per (variant,transcript)
# — reported per equivalence CLASS (the witness INFO/CLASS). A class with pairs and 0 divergence
# is COVERED-and-passing; any divergence is a fuzz finding. (Statistical tier = exact CIs over a
# real corpus; this formal tier = class coverage. Same diff.) Native duckdb-r v1.5.3 loads the
# extension; VEP is a subprocess; SO-term sets sort byte-order (radix) to match DuckDB list_sort.
suppressMessages({ library(optparse); library(jsonlite); library(duckdb); library(DBI) })
options(rlang_backtrace_on_error = "none")

root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
op <- OptionParser()
op <- add_option(op, "--vcf",   default = file.path(root, "conformance/data/witnesses.vcf"))
op <- add_option(op, "--gff",   default = file.path(root, "data/giab/GRCh38.116.controlled.gff3.gz"))
op <- add_option(op, "--fasta", default = file.path(root, "data/giab/GRCh38.primary.fa"))
op <- add_option(op, "--ext",   default = Sys.getenv("DUCKVEP_EXT", file.path(root, "build/release/duckvep.duckdb_extension")))
op <- add_option(op, "--fork",  default = as.character(max(1, parallel::detectCores() - 2)))
opt <- parse_args(op)
vep_cmd <- strsplit(Sys.getenv("VEP_CMD", "conda run -n vep vep"), " ")[[1]]
setlist <- function(x) paste(sort(x, method = "radix"), collapse = "&")   # byte order, matches DuckDB

# class per ORIGINAL (chrom,pos,ref,alt) from the witness INFO/CLASS
wl <- readLines(opt$vcf); wl <- wl[!grepl("^#", wl)]
wf <- do.call(rbind, strsplit(wl, "\t"));
cls <- data.frame(opos = as.integer(wf[, 2]), oref = wf[, 4], oalt = wf[, 5],
                  class = sub(".*CLASS=", "", wf[, 8]), stringsAsFactors = FALSE)

con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE))
invisible(dbExecute(con, sprintf("LOAD '%s'", opt$ext)))
invisible(dbGetQuery(con, sprintf("SELECT vep_load_cache('%s','%s')", opt$gff, opt$fasta)))

# normalized key per original witness -> attach class to (npos,nref,nalt)
nk <- dbGetQuery(con, sprintf("
  SELECT v.pos opos, v.ref oref, a.alt oalt,
    normalize_variant(v.pos,v.ref,a.alt).pos npos,
    CASE WHEN normalize_variant(v.pos,v.ref,a.alt).ref='' THEN '-' ELSE normalize_variant(v.pos,v.ref,a.alt).ref END nref,
    CASE WHEN normalize_variant(v.pos,v.ref,a.alt).alt='' THEN '-' ELSE normalize_variant(v.pos,v.ref,a.alt).alt END nalt
  FROM read_vcf('%s') v, UNNEST(v.alt) a(alt)", opt$vcf))
nk <- merge(nk, cls, by = c("opos","oref","oalt"))

# duckvep rows (native), canonical key
dv <- dbGetQuery(con, sprintf("
  SELECT nv.pos pos, CASE WHEN nv.ref='' THEN '-' ELSE nv.ref END AS \"ref\",
         CASE WHEN nv.alt='' THEN '-' ELSE nv.alt END AS \"alt\", c.transcript_id tx,
         list_aggregate(list_sort(c.consequence),'string_agg','&') dc
  FROM read_vcf('%s') v, UNNEST(v.alt) a(alt),
       UNNEST(vep_consequence(v.chrom,v.pos,v.ref,a.alt)) u(c),
       LATERAL (SELECT normalize_variant(v.pos,v.ref,a.alt) nv)
  GROUP BY ALL", opt$vcf))

# VEP --gff on the same witnesses
system2(vep_cmd[1], c(vep_cmd[-1], "-i", opt$vcf, "--gff", opt$gff, "--fasta", opt$fasta,
        "--distance", "5000", "--json", "-o", "/tmp/wit_vep.json", "--fork", opt$fork,
        "--force_overwrite", "--no_stats"), stdout = FALSE, stderr = FALSE)
recs <- stream_in(file("/tmp/wit_vep.json"), verbose = FALSE)
V <- list()
for (i in seq_len(nrow(recs))) {
  r <- recs[i, ]; tcs <- r$transcript_consequences[[1]]; if (is.null(tcs) || !length(tcs)) next
  nref <- strsplit(r$allele_string, "/")[[1]][1]
  for (j in seq_len(nrow(tcs))) { tc <- tcs[j, ]
    V[[length(V)+1]] <- data.frame(pos = r$start, ref = nref,
      alt = if (!is.null(tc$variant_allele)) tc$variant_allele else "",
      tx = tc$transcript_id, vc = setlist(tc$consequence_terms[[1]]), stringsAsFactors = FALSE) }
}
vep <- do.call(rbind, V)

# join VEP <-> duckvep on canonical key, attach class, compare sets
m <- merge(vep, dv, by = c("pos","ref","alt","tx"))
m <- merge(m, nk[, c("npos","nref","nalt","class")], by.x = c("pos","ref","alt"),
           by.y = c("npos","nref","nalt"), all.x = TRUE)
m$class[is.na(m$class)] <- "(other)"
m$ok <- m$vc == m$dc

by <- aggregate(ok ~ class, m, function(x) c(n = length(x), disc = sum(!x)))
rep <- data.frame(class = by$class, n = by$ok[, "n"], discordant = by$ok[, "disc"])
rep <- rep[order(-rep$discordant, -rep$n), ]
cat(sprintf("Differential fuzz over %d witnesses -> %d (variant,transcript) pairs, %d classes covered\n",
            nrow(cls), nrow(m), length(unique(m$class))))
cat(sprintf("  duckvep ≡ VEP on %d/%d pairs; %d divergences across %d classes\n",
            sum(m$ok), nrow(m), sum(!m$ok), sum(rep$discordant > 0)))
for (i in seq_len(nrow(rep))) cat(sprintf("  %-26s n=%-5d disc=%d\n", rep$class[i], rep$n[i], rep$discordant[i]))
d <- m[!m$ok, c("class","pos","ref","alt","tx","vc","dc")]
if (nrow(d)) { cat("--- divergences (fuzz findings) ---\n")
  for (i in seq_len(min(20, nrow(d)))) with(d[i, ],
    cat(sprintf("  [%s] %s:%s %s>%s %s  VEP=%s  duckvep=%s\n", class, "", pos, ref, alt, tx, vc, dc))) }
write.csv(m, file.path(root, "conformance/data/fuzz_results.csv"), row.names = FALSE)
quit(status = if (sum(!m$ok) > 0) 1 else 0)
