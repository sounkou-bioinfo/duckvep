#!/usr/bin/env bash
# Fetch full-size GIAB + reference inputs for benchmarking/concordance into
# data/giab/ (gitignored). Large downloads (~GB); run once. Override URLs via env.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/giab && cd data/giab

GIAB_VCF=${GIAB_VCF:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}
REF=${REF:-https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}
GFF3=${GFF3:-https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.gff3.gz}

fetch() { local url=$1 out=$2; [ -f "$out" ] && { echo "have $out"; return; }; echo "↓ $out"; curl -fSL --retry 3 -o "$out" "$url"; }

fetch "$GIAB_VCF" HG002.vcf.gz
fetch "$GFF3"     GRCh38.112.gff3.gz
fetch "$REF"      GRCh38.fa.gz

# Reference needs to be bgzipped + faidx for random access (see DESIGN.md §5).
if [ ! -f GRCh38.fa ]; then echo "decompressing reference…"; gunzip -k GRCh38.fa.gz; fi
command -v samtools >/dev/null && [ ! -f GRCh38.fa.fai ] && samtools faidx GRCh38.fa || \
  echo "note: install samtools and run 'samtools faidx GRCh38.fa' for FASTA random access"
echo "GIAB inputs ready in $(pwd)"
