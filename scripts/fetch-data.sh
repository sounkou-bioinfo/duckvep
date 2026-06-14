#!/usr/bin/env bash
# Data receipt — fetch EVERY input the benchmarks + correctness pipelines use, with
# pinned sources/versions so anyone reproduces the same data. Everything lands in
# data/giab/ + data/vep_cache/ (both gitignored; large, ~GB). Run once. Override
# any URL via env. Provenance is the point: no mystery data.
#
# Pinned versions (match what the committed reports were generated against):
#   Ensembl release : 116 (GRCh38)         -> GFF3 + primary-assembly FASTA + VEP cache
#   ClinVar         : fileDate 2026-06-06  -> GRCh38 VCF (NCBI, Ensembl-style numeric contigs)
#   GIAB            : HG002 v4.2.1 benchmark (AshkenazimTrio)
#   Ensembl VEP     : install code matching the cache (116) via micromamba/condathis
set -euo pipefail
cd "$(dirname "$0")/.."
REL=${ENSEMBL_RELEASE:-116}
mkdir -p data/giab data/vep_cache && cd data/giab

EFTP="https://ftp.ensembl.org/pub/release-${REL}"
REF=${REF:-$EFTP/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}
GFF3=${GFF3:-$EFTP/gff3/homo_sapiens/Homo_sapiens.GRCh38.${REL}.gff3.gz}
VEP_CACHE_URL=${VEP_CACHE_URL:-$EFTP/variation/indexed_vep_cache/homo_sapiens_vep_${REL}_GRCh38.tar.gz}
CLINVAR=${CLINVAR:-https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20260606.vcf.gz}
GIAB_VCF=${GIAB_VCF:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}

fetch() { local url=$1 out=$2; [ -f "$out" ] && { echo "have $out"; return; }; echo "↓ $out  ($url)"; curl -fSL --retry 3 -o "$out" "$url"; }

fetch "$GFF3"     "GRCh38.${REL}.gff3.gz"
fetch "$REF"      GRCh38.primary.fa.gz
fetch "$CLINVAR"  clinvar.vcf.gz
fetch "$GIAB_VCF" HG002.vcf.gz

# Primary-assembly FASTA: decompress + faidx for random access; carve chr17 subset.
[ -f GRCh38.primary.fa ] || { echo "decompressing reference…"; gunzip -k GRCh38.primary.fa.gz; }
if command -v samtools >/dev/null; then
  [ -f GRCh38.primary.fa.fai ] || samtools faidx GRCh38.primary.fa
  [ -f chr17.fa ] || { samtools faidx GRCh38.primary.fa 17 > chr17.fa; samtools faidx chr17.fa; }
else
  echo "note: install samtools, then 'samtools faidx GRCh38.primary.fa' and carve chr17.fa"
fi

# Offline Ensembl VEP cache (release-matched) -> data/vep_cache/homo_sapiens/${REL}_GRCh38
if [ ! -d "../vep_cache/homo_sapiens/${REL}_GRCh38" ]; then
  echo "↓ VEP cache ${REL} (~27GB)…"; curl -fSL --retry 3 "$VEP_CACHE_URL" -o /tmp/vep_cache.tar.gz
  tar -xzf /tmp/vep_cache.tar.gz -C ../vep_cache && rm -f /tmp/vep_cache.tar.gz
fi

cat <<EOF
inputs ready under data/ :
  giab/GRCh38.${REL}.gff3.gz  giab/GRCh38.primary.fa(.fai)  giab/chr17.fa  giab/clinvar.vcf.gz  giab/HG002.vcf.gz
  vep_cache/homo_sapiens/${REL}_GRCh38/

Ensembl VEP code must match the cache (${REL}). Preferred (condathis/micromamba):
  micromamba create -y -n vep -c bioconda -c conda-forge ensembl-vep=${REL}
EOF
