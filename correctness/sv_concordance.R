#!/usr/bin/env Rscript
# Structural-variant concordance vs Ensembl VEP — the SV path was never checked against the
# oracle ("Potemkin"). Runs symbolic SVs (<DEL>/<DUP>/<INV>) over a known gene through BOTH
# VEP (--gff, same gene model) and duckvep's END-aware vep_consequence, and diffs the
# consequence sets. Requires the controlled GFF + FASTA (scripts/fetch-data.R).
#
# Usage: correctness/sv_concordance.R
# Current state (recorded 2026-06-16, TP53 ENST00000269305) — 4/4 MATCH:
#   DEL whole-gene -> transcript_ablation; DUP whole-gene -> transcript_amplification;
#   DEL partial -> coding_sequence_variant&feature_truncation&intron_variant;
#   INV partial -> 5_prime_UTR_variant&coding_sequence_variant&intron_variant.
# NOTE: only a 4-case TP53 harness — get GIAB SV Tier1 coverage before further type changes.
suppressMessages(library(jsonlite))

root   <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
DUCKDB <- Sys.getenv("DUCKVEP_DUCKDB", file.path(root, ".tools/duckdb"))
EXT    <- Sys.getenv("DUCKVEP_EXT",   file.path(root, "build/release/duckvep.duckdb_extension"))
CTRL   <- Sys.getenv("DUCKVEP_CTRL_GFF", file.path(root, "data/giab/GRCh38.116.controlled.gff3.gz"))
FASTA  <- Sys.getenv("DUCKVEP_FASTA", file.path(root, "data/giab/GRCh38.primary.fa"))
TR     <- Sys.getenv("SV_TRANSCRIPT", "ENST00000269305")
vep_cmd <- strsplit(Sys.getenv("VEP_CMD", "conda run -n vep vep"), " ")[[1]]
if (!file.exists(FASTA) || !file.exists(CTRL)) { message("SKIP: need controlled GFF + FASTA"); quit(status = 0) }
setlist <- function(x) paste(sort(x, method = "radix"), collapse = "&")

# The four symbolic SVs (id, pos, end, alt) — one VCF for VEP, the same tuples inline for duckvep.
svs <- data.frame(
  id  = c("del_whole", "del_part", "dup_whole", "inv"),
  pos = c(7668000L,    7676000L,   7668000L,    7670000L),
  end = c(7688000L,    7676500L,   7688000L,    7680000L),
  alt = c("<DEL>",     "<DEL>",    "<DUP>",     "<INV>"),
  stringsAsFactors = FALSE)

SVVCF <- "/tmp/sv_concordance.vcf"
writeLines(c("##fileformat=VCFv4.2",
  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">',
  '##INFO=<ID=END,Number=1,Type=Integer,Description="">',
  '##ALT=<ID=DEL,Description="">', '##ALT=<ID=DUP,Description="">', '##ALT=<ID=INV,Description="">',
  "##contig=<ID=17>", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
  sprintf("17\t%d\t%s\tN\t%s\t.\t.\tSVTYPE=%s;END=%d", svs$pos, svs$id, svs$alt,
          gsub("[<>]", "", svs$alt), svs$end)), SVVCF)

system2(vep_cmd[1], c(vep_cmd[-1], "-i", SVVCF, "--gff", CTRL, "--fasta", FASTA, "--symbol",
        "--json", "-o", "/tmp/sv_vep.json", "--force_overwrite", "--no_stats"), stdout = FALSE, stderr = FALSE)
if (!file.exists("/tmp/sv_vep.json") || file.info("/tmp/sv_vep.json")$size == 0) {
  message("FAIL: VEP produced no output"); quit(status = 1) }
vep <- list()
for (l in readLines("/tmp/sv_vep.json")) { r <- fromJSON(l, simplifyVector = FALSE)
  for (t in r$transcript_consequences) if (!is.null(t$transcript_id) && t$transcript_id == TR)
    vep[[if (!is.null(r$id)) r$id else "?"]] <- setlist(unlist(t$consequence_terms)) }

# duckvep reads the SAME controlled GFF as VEP (so the only free variable is the engine).
f <- tempfile(fileext = ".sql")
writeLines(sprintf("LOAD '%s'; SELECT vep_load_cache('%s','%s');
SELECT id||chr(9)||list_aggregate(list_sort(c.consequence),'string_agg','&')
FROM (VALUES %s) t(id,pos,e,alt), UNNEST(vep_consequence('17',pos,e,'N',alt)) u(c)
WHERE c.transcript_id='%s';", EXT, CTRL, FASTA,
  paste(sprintf("('%s',%d,%d,'%s')", svs$id, svs$pos, svs$end, svs$alt), collapse = ","), TR), f)
dv_raw <- system2(DUCKDB, c("-unsigned", "-noheader", "-list", "-f", f), stdout = TRUE, stderr = FALSE)
unlink(f)
dv_raw <- dv_raw[grepl("^(del|dup|inv)", dv_raw)]
if (length(dv_raw) != 4) { message(sprintf("FAIL: duckvep did not annotate all 4 SVs (%d/4)", length(dv_raw))); quit(status = 1) }
dv <- setNames(sub("^[^\t]*\t", "", dv_raw), sub("\t.*$", "", dv_raw))

cat(sprintf("SV concordance vs Ensembl VEP (%s):\n", TR))
match <- 0
for (id in names(vep)) {
  ok <- identical(vep[[id]], unname(dv[id]))
  match <- match + ok
  cat(sprintf("  %-5s %-10s VEP=%s  duckvep=%s\n", if (ok) "MATCH" else "DIFF ", id, vep[[id]], dv[id]))
}
cat(sprintf("  ---- %d/%d match\n", match, length(vep)))
quit(status = if (match == length(vep)) 0 else 1)
