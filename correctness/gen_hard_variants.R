#!/usr/bin/env Rscript
# Generate a synthetic corpus of PATHOLOGICAL / hard variants for the consequence
# classes we target — placed at *exact* offsets relative to transcript features
# (which random sampling almost never hits), to stress-test duckvep vs Ensembl VEP.
#
# Classes: splice-boundary indels that SPAN a donor/acceptor/5th-base/donor-region/
# polypyrimidine/splice-region boundary (the interval-vs-point case), CDS frameshift
# (1 bp) / inframe (3 bp) indels, and in-codon MNVs.
#
# Stack (no Python): DuckDB SQL for coordinates, samtools faidx for ref bases,
# vcfppR's writer for output — the same tooling as ~/duckhts.
#
# Usage: gen_hard_variants.R <exons.parquet> <fasta> <transcripts.parquet> <chrom> <out.vcf>
suppressMessages(library(vcfppR))
a <- commandArgs(trailingOnly = TRUE)
EXONS <- a[1]; FASTA <- a[2]; TX <- a[3]; CHROM <- a[4]; OUT <- a[5]
DUCKDB <- Sys.getenv("DUCKVEP_DUCKDB", ".tools/duckdb")

duck <- function(sql) {
  txt <- system2(DUCKDB, c("-csv", "-noheader", "-c", sql), stdout = TRUE)
  if (length(txt) == 0) return(data.frame())
  read.csv(text = paste(txt, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE)
}

# Coding, multi-exon transcripts on this contig (both strands), MANE/canonical.
tx <- duck(sprintf("SELECT transcript_id FROM read_parquet('%s')
  WHERE chrom='%s' AND biotype='protein_coding' AND (mane_select OR canonical)
  ORDER BY gene_symbol LIMIT 40", TX, CHROM))[[1]]
ex <- duck(sprintf("SELECT transcript_id, start, \"end\" FROM read_parquet('%s')
  WHERE transcript_id IN ('%s') ORDER BY transcript_id, start",
  EXONS, paste(tx, collapse = "','")))
names(ex) <- c("tid", "start", "end")

# Introns = gaps between position-sorted exons; keep windows around each boundary.
introns <- do.call(rbind, lapply(split(ex, ex$tid), function(d) {
  d <- d[order(d$start), ]
  if (nrow(d) < 2) return(NULL)
  data.frame(istart = head(d$end, -1) + 1, iend = tail(d$start, -1) - 1)
}))
introns <- introns[introns$iend - introns$istart >= 40, ]

# One samtools faidx call for all donor/acceptor windows; index by "chr:a-b".
regions <- unique(c(sprintf("%s:%d-%d", CHROM, introns$istart - 6, introns$istart + 8),
                    sprintf("%s:%d-%d", CHROM, introns$iend - 18, introns$iend + 8)))
fa <- system2("samtools", c("faidx", FASTA, regions), stdout = TRUE)
h <- grepl("^>", fa); ids <- sub("^>", "", fa[h])
seqs <- toupper(tapply(fa[!h], cumsum(h)[!h], paste, collapse = ""))
win <- setNames(as.character(seqs), ids)
# Genomic [a,b] inclusive from whichever window contains it (1-based slice).
getseq <- function(a, b) {
  for (nm in names(win)) {
    p <- as.integer(strsplit(sub(".*:", "", nm), "-")[[1]])
    if (a >= p[1] && b <= p[2]) return(substr(win[[nm]], a - p[1] + 1, b - p[2] + nchar(win[[nm]])))
  }
  NA_character_
}

recs <- list(); seen <- new.env()
add <- function(pos, ref, alt, tag) {
  if (is.na(ref) || is.na(alt) || pos < 2 || grepl("N", ref) || grepl("N", alt)) return()
  k <- paste(pos, ref, alt); if (!is.null(seen[[k]])) return(); seen[[k]] <- TRUE
  recs[[length(recs) + 1]] <<- sprintf("%s\t%d\t%s\t%s\t%s\t.\t.\t.", CHROM, pos, tag, ref, alt)
}
del <- function(s, e, tag) add(s - 1, getseq(s - 1, e), getseq(s - 1, s - 1), tag)  # VCF anchor-base del
ins <- function(p, b, tag) { a <- getseq(p, p); add(p, a, paste0(a, b), tag) }
comp <- c(A = "T", T = "A", C = "G", G = "C")

for (i in seq_len(nrow(introns))) {
  is <- introns$istart[i]; ie <- introns$iend[i]
  del(is - 1, is + 1, "del_span_donor5p")      # exon|GT (spans 5' boundary)
  del(is + 3, is + 5, "del_donor_region")
  del(is + 4, is + 4, "del_5th_base")
  del(ie - 1, ie + 1, "del_span_acceptor3p")   # AG|exon (spans 3' boundary)
  del(ie - 12, ie - 8, "del_polypyrimidine")
  del(is + 2, is + 7, "del_splice_region")
  ins(is + 1, "AT", "ins_donor")
  ins(ie - 3, "GC", "ins_acceptor_region")
  del(is - 1, is - 1, "del_fs1_exonend")        # 1 bp -> frameshift
  del(is - 4, is - 2, "del_inframe3_exonend")   # 3 bp -> inframe
  r <- getseq(ie + 6, ie + 8)                   # in-codon MNV in next exon
  if (!is.na(r) && nchar(r) == 3) add(ie + 6, r, paste(comp[strsplit(r, "")[[1]]], collapse = ""), "mnv3_exon")
}

w <- methods::new(vcfppR::vcfwriter, OUT, "VCFv4.2")
w$addLine("##source=correctness/gen_hard_variants.R")
w$addContig(CHROM)
ord <- order(as.integer(sub(".*\t(\\d+)\t.*", "\\1", sub("^[^\t]*\t", "", recs))))
for (r in recs[ord]) w$writeline(r)
w$close()
cat(sprintf("generated %d hard variants across %d transcripts\n", length(recs), length(tx)), file = stderr())
