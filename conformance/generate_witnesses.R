#!/usr/bin/env Rscript
# Formal-tier witness generator (R-native, via Rduckhts/duckhts — stable C API, loads in any
# DuckDB build). Reads a transcript's structure from a GFF and the reference sequence from a
# FASTA *range*, then tiles the equivalence classes where consequence bugs live — each exon/intron
# boundary (splice donor/acceptor/5th-base/region), the start codon, the stop codon, exon interior
# and UTR — crossed with allele shapes: SNV ×3, 1 bp ins + 1/2 bp del (frameshift), 3 bp del
# (in-frame), and a 2 bp MNV. Emits a witness VCF; the differential fuzzer (VEP --gff ⟂ duckvep) asserts an exact
# SO-term-set match on it. This is the FORMAL (covering) tier; conformance/stratified_conformance.R
# is the STATISTICAL tier — they share the one fuzzer. See conformance/README.md.
suppressMessages({ library(optparse); library(Rduckhts); library(duckdb); library(DBI) })
options(rlang_backtrace_on_error = "none")

root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
op <- OptionParser(usage = "%prog --tx ENST... [options]")
op <- add_option(op, "--tx",    default = "ENST00000269305", help = "transcript id (default: TP53) [%default]")
op <- add_option(op, "--gff",   default = file.path(root, "data/giab/GRCh38.116.controlled.gff3.gz"))
op <- add_option(op, "--fasta", default = file.path(root, "data/giab/GRCh38.primary.fa"))
op <- add_option(op, "--out",   default = file.path(root, "conformance/data/witnesses.vcf"))
op <- add_option(op, "--duckhts", default = Sys.getenv("DUCKHTS_EXT", "/root/Rtinycc/Rducks/duckhts.duckdb_extension"))
opt <- parse_args(op)

con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE))
invisible(rduckhts_load(con, extension_path = opt$duckhts))

# 1. transcript exons + CDS from the GFF (attributes carry Parent=transcript:ENST...).
invisible(dbExecute(con, sprintf("CREATE TABLE feats AS SELECT seqname AS chrom, feature, start, \"end\", strand,
  regexp_extract(attributes,'(ENST[0-9]+)',1) AS tx FROM read_gff('%s')
  WHERE feature IN ('exon','CDS','five_prime_UTR','three_prime_UTR')", opt$gff)))
tx <- dbGetQuery(con, sprintf("SELECT chrom,feature,start,\"end\",strand FROM feats WHERE tx='%s' ORDER BY start", opt$tx))
if (nrow(tx) == 0) stop("no features for ", opt$tx, " — is it in the GFF region?")
chrom <- tx$chrom[1]; strand <- tx$strand[1]
exons <- tx[tx$feature == "exon", c("start", "end")]; exons <- exons[order(exons$start), ]
cds   <- tx[tx$feature == "CDS",  c("start", "end")]; cds   <- cds[order(cds$start), ]
cds_lo <- if (nrow(cds)) min(cds$start) else NA; cds_hi <- if (nrow(cds)) max(cds$end) else NA

# 2. class-relevant genomic POSITIONS (the covering-array offsets).
pos <- integer(0); lab <- character(0)
add <- function(p, l) { pos <<- c(pos, p); lab <<- c(lab, rep(l, length(p))) }
for (i in seq_len(nrow(exons))) {
  s <- exons$start[i]; e <- exons$end[i]
  add(c(s, s + 1, e - 1, e), "exon_edge")                 # exon boundary interior
  if (i < nrow(exons)) add(c(e + 1, e + 2, e + 5), "donor")     # intron 5' (donor/5th/region)
  if (i > 1)           add(c(s - 1, s - 2, s - 8), "acceptor")  # intron 3' (acceptor/ppt)
  add(round((s + e) / 2), "exon_mid")                     # deep exon interior
}
if (!is.na(cds_lo)) {                                     # start & stop codons (strand-aware)
  start_codon <- if (strand == "+") cds_lo:(cds_lo + 2) else (cds_hi - 2):cds_hi
  stop_codon  <- if (strand == "+") (cds_hi - 2):cds_hi else cds_lo:(cds_lo + 2)
  add(start_codon, "start_codon"); add(stop_codon, "stop_codon")
}
# Dedup positions keeping the MOST SPECIFIC label (a position that is both an exon edge and the
# start codon must keep `start_codon`, not whichever was added first — pi P0/P1-8). Lower rank
# = more specific = wins.
pos <- pmax(pos, 1L)
rank <- c(start_codon = 1, stop_codon = 2, donor = 3, acceptor = 4, exon_edge = 5, exon_mid = 6)
ord <- order(pos, rank[lab]); pos <- pos[ord]; lab <- lab[ord]
keep <- !duplicated(pos); pos <- pos[keep]; lab <- lab[keep]

# 3. reference base at each position (FASTA range read via duckhts -> the range's SEQUENCE).
p0 <- min(pos); span <- sprintf("%s:%d-%d", chrom, p0, max(pos) + 4)
invisible(rduckhts_fasta(con, "refseq", opt$fasta, region = span, overwrite = TRUE))
seqstr <- toupper(dbGetQuery(con, "SELECT SEQUENCE FROM refseq LIMIT 1")$SEQUENCE[1])
nuc <- function(p, n) substr(seqstr, p - p0 + 1, p - p0 + n)
base_at <- function(p) { s <- nuc(p, 1); if (s %in% c("A","C","G","T")) s else NA_character_ }
twobase <- function(p) nuc(p, 2)

# 4. tile allele shapes per position -> witness rows.
W <- list()
for (k in seq_along(pos)) {
  p <- pos[k]; rb <- base_at(p); if (is.na(rb)) next
  for (alt in setdiff(c("A","C","G","T"), rb)) W[[length(W)+1]] <- c(p, rb, alt, paste0(lab[k], "_snv"))   # SNV ×3
  tb <- twobase(p); if (grepl("^[ACGT]{2}$", tb)) {
    W[[length(W)+1]] <- c(p, tb, substr(tb,1,1), paste0(lab[k], "_del1"))          # 1bp del (frameshift)
    W[[length(W)+1]] <- c(p, rb, paste0(rb,"T"), paste0(lab[k], "_ins1"))          # 1bp ins (frameshift)
    W[[length(W)+1]] <- c(p, tb, chartr("ACGT","TGCA",tb), paste0(lab[k], "_mnv2"))# 2bp MNV (both bases flip; complement guarantees both differ)
  }
  fb <- nuc(p, 3)
  if (grepl("^[ACGT]{3}$", fb)) W[[length(W)+1]] <- c(p, fb, substr(fb,1,1), paste0(lab[k], "_del2"))  # 2bp del (frameshift)
  qb <- nuc(p, 4)
  if (grepl("^[ACGT]{4}$", qb)) W[[length(W)+1]] <- c(p, qb, substr(qb,1,1), paste0(lab[k], "_del3"))  # 3bp del (in-frame)
}
wdf <- as.data.frame(do.call(rbind, W), stringsAsFactors = FALSE)
names(wdf) <- c("pos","ref","alt","class"); wdf$pos <- as.integer(wdf$pos)
wdf <- wdf[order(wdf$pos), ]

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
vc <- file(opt$out, "w")
# Tag the TARGET transcript: each CLASS describes geometry on THIS transcript only. The fuzzer
# applies the class label solely to the target transcript's pairs (other overlapping transcripts
# at the same locus have different geometry, so they are '(other)', not mislabeled — pi P1-8).
writeLines(c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s>", chrom),
  '##INFO=<ID=CLASS,Number=1,Type=String,Description="witness equivalence class (on the target transcript)">',
  '##INFO=<ID=TX,Number=1,Type=String,Description="target transcript the CLASS geometry refers to">',
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), vc)
writeLines(sprintf("%s\t%d\t.\t%s\t%s\t.\t.\tCLASS=%s;TX=%s", chrom, wdf$pos, wdf$ref, wdf$alt, wdf$class, opt$tx), vc)
close(vc)
cat(sprintf("wrote %d witnesses for %s (%s, %d classes) -> %s\n",
            nrow(wdf), opt$tx, chrom, length(unique(wdf$class)), opt$out))
