# Refinements & open items

Tracked refinements and future work, kept **out of** [`DESIGN.md`](DESIGN.md) (which
holds the stable design rationale) so the design doc doesn't accumulate drift.
Current, measured correctness figures are in the generated report
[`../correctness/correctness.md`](../correctness/correctness.md); divergences from
upstream fastVEP (with patch files) are in [`PATCHES.md`](PATCHES.md).

## Open correctness gaps (paper-relevant)

- **High-impact indel / MNV engine gap** — the measured frontier vs Ensembl VEP,
  *shared with fastVEP* (so it's engine accuracy, not a duckvep-specific bug).
  Concrete next pieces: a frameshift that introduces a premature stop should add
  `stop_gained` (we emit only `frameshift_variant`); MNV codon handling.
- **Mitochondrial start codons** — the chrM codon table (NCBI table 2) is wired
  (`predictor.rs::ct`), but non-ATG mitochondrial *start* codons (ATA/ATT/GTG) are
  not yet recognised by the `start_lost` gate — the residual chrM SNV discordance.
- **chrX/Y haploid + pseudoautosomal regions (PAR)** — PAR genes appear on both X
  and Y; hemizygous calling and PAR boundaries are untested. HG002 is male; add a
  **female sample (e.g. NA12878)** so chrX-diploid + PAR are exercised. Needs
  regression tests in the `correctness/` suite (genome-wide cache + primary FASTA).
- **Haplotype-aware consequences (haplosaurus)** — the kernel is per-variant; a
  multi-variant apply over `translateable_seq` (phased GT is already first-class in
  `read_vcf`) generalises it. See AGENTS.md.

## Cache builder

- **RefSeq xrefs** + the `otherfeatures` RefSeq gene models (for `--merged`).
- **MySQL-dump text escaping** — `read_csv` does not de-escape MySQL's `\t`/`\n`/`\\`
  in free-text columns; we only consume escape-free id/coord/code/stable_id columns.
  Revisit if importing description-heavy tables.
- **Engine loader** for the columnar Ensembl cache — currently the engine reads its
  own GFF3-derived Parquet; a loader for `transcripts/exons/translations.parquet`
  would let the curated flags (codon_table, cds_start_NF, …) feed the predictor
  directly.

## Supplementary annotations (roadmap → joins, not formats)

ClinVar / gnomAD / dbSNP / COSMIC / dbNSFP scores / OMIM / constraint — each a
Parquet/DuckDB table joined on the **`normalize_variant`** canonical key (the
load-bearing piece is already in place). See AGENTS.md feature roadmap.
