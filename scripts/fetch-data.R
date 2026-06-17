#!/usr/bin/env Rscript
# Data receipt — fetch EVERY input the benchmarks + correctness pipelines use, with pinned
# sources/versions so anyone reproduces the same data. Everything lands in data/giab/ +
# data/vep_cache/ (both gitignored; large, ~GB). Run once. Override any URL via env. Provenance
# is the point: no mystery data.
#
# Pinned versions (match what the committed reports were generated against):
#   Ensembl release : 116 (GRCh38)        -> GFF3 + primary-assembly FASTA + VEP cache
#   ClinVar         : fileDate 2026-06-06 -> GRCh38 VCF (NCBI, Ensembl-style numeric contigs)
#   GIAB            : HG002 v4.2.1 benchmark (AshkenazimTrio)
#   Ensembl VEP     : install code matching the cache (116) via micromamba/condathis
root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
setwd(root)
REL  <- Sys.getenv("ENSEMBL_RELEASE", "116")
dir.create("data/giab", recursive = TRUE, showWarnings = FALSE)
dir.create("data/vep_cache", recursive = TRUE, showWarnings = FALSE)
setwd("data/giab")

EFTP    <- sprintf("https://ftp.ensembl.org/pub/release-%s", REL)
REF     <- Sys.getenv("REF",     sprintf("%s/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", EFTP))
GFF3    <- Sys.getenv("GFF3",    sprintf("%s/gff3/homo_sapiens/Homo_sapiens.GRCh38.%s.gff3.gz", EFTP, REL))
VEP_CACHE_URL <- Sys.getenv("VEP_CACHE_URL", sprintf("%s/variation/indexed_vep_cache/homo_sapiens_vep_%s_GRCh38.tar.gz", EFTP, REL))
CLINVAR  <- Sys.getenv("CLINVAR",  "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20260606.vcf.gz")
GIAB_VCF <- Sys.getenv("GIAB_VCF", "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")

fetch <- function(url, out) {
  if (file.exists(out)) { cat(sprintf("have %s\n", out)); return(invisible()) }
  cat(sprintf("↓ %s  (%s)\n", out, url))
  if (system2("curl", c("-fSL", "--retry", "3", "-o", out, url)) != 0) stop(sprintf("download failed: %s", url))
}

fetch(GFF3,     sprintf("GRCh38.%s.gff3.gz", REL))
fetch(REF,      "GRCh38.primary.fa.gz")
fetch(CLINVAR,  "clinvar.vcf.gz")
fetch(GIAB_VCF, "HG002.vcf.gz")

# ── Diverse real-data validation cohort (gated by DIVERSE=1). Broadens concordance vs VEP across
#    sex (M/F), ancestry (Ashkenazi + Han Chinese + 1000G all superpopulations incl. African),
#    structural variants, and PHASED haplotypes. Pinned: GIAB v4.2.1, HG002 SV Tier1 v0.6, 1000G
#    3202 high-cov phased panel (2020-08-05). chr21 slice keeps the phased panel tractable.
if (Sys.getenv("DIVERSE", "1") == "1") {
  G <- Sys.getenv("GIABREL", "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab")
  fetch(sprintf("%s/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz", G), "HG003.vcf.gz")  # AJ father (male)
  fetch(sprintf("%s/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz", G), "HG004.vcf.gz")  # AJ mother (female)
  fetch(sprintf("%s/release/ChineseTrio/HG005_NA24631_son/latest/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz", G), "HG005.vcf.gz")        # Han Chinese son (male)
  fetch(sprintf("%s/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz", G), "HG002.SV.Tier1.vcf.gz")          # structural variants
  TGP <- Sys.getenv("TGP", "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased")
  fetch(sprintf("%s/CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz", TGP), "1000G.3202.phased.chr21.vcf.gz")   # phased haplotypes, all superpops
}

# Primary-assembly FASTA: decompress + faidx for random access; carve chr17 subset.
if (!file.exists("GRCh38.primary.fa")) { cat("decompressing reference…\n"); system2("gunzip", c("-k", "GRCh38.primary.fa.gz")) }
if (nzchar(Sys.which("samtools"))) {
  if (!file.exists("GRCh38.primary.fa.fai")) system2("samtools", c("faidx", "GRCh38.primary.fa"))
  if (!file.exists("chr17.fa")) {
    writeLines(system2("samtools", c("faidx", "GRCh38.primary.fa", "17"), stdout = TRUE), "chr17.fa")
    system2("samtools", c("faidx", "chr17.fa"))
  }
} else cat("note: install samtools, then 'samtools faidx GRCh38.primary.fa' and carve chr17.fa\n")

# Offline Ensembl VEP cache (release-matched) -> data/vep_cache/homo_sapiens/${REL}_GRCh38
if (!dir.exists(sprintf("../vep_cache/homo_sapiens/%s_GRCh38", REL))) {
  cat(sprintf("↓ VEP cache %s (~27GB)…\n", REL))
  if (system2("curl", c("-fSL", "--retry", "3", VEP_CACHE_URL, "-o", "/tmp/vep_cache.tar.gz")) == 0) {
    system2("tar", c("-xzf", "/tmp/vep_cache.tar.gz", "-C", "../vep_cache")); unlink("/tmp/vep_cache.tar.gz")
  }
}

cat(sprintf("inputs ready under data/ :\n  giab/GRCh38.%s.gff3.gz  giab/GRCh38.primary.fa(.fai)  giab/chr17.fa  giab/clinvar.vcf.gz  giab/HG002.vcf.gz\n  vep_cache/homo_sapiens/%s_GRCh38/\n\nEnsembl VEP code must match the cache (%s). Preferred (condathis/micromamba):\n  micromamba create -y -n vep -c bioconda -c conda-forge ensembl-vep=%s\n", REL, REL, REL, REL))
