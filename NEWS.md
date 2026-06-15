# duckvep NEWS

Changelog, most recent first. (R-package style.)

## duckvep 0.3.0 (in development)

### Accuracy (oracle = Ensembl VEP, not fastVEP)

* **Concordance is always split by IMPACT × variant class** and per SO term, with
  fastVEP shown alongside — never a single aggregate. duckvep beats fastVEP on
  **every HIGH-impact SO term**. Percentages use enough precision that a non-zero
  discordance never reads as 100%.
* **The consequence engine was rebuilt to mirror Ensembl's own structure**, taking
  N=50000 ClinVar discordance vs offline Ensembl VEP 116 from **~3,876 → 78 (a 98%
  reduction) with ZERO duckvep-specific regressions** at every step (every remaining
  discordance is a *shared* gap fastVEP has too — fastVEP itself is discordant on
  thousands of the same calls). Three VEP-faithful abstractions:
  * **`CodingContext` (haplotype-ready):** coding consequences are a predicate SET
    over a peptide/codon context built from `CdsEdit`s applied to the reference CDS.
    One variant = one edit; a phased haplotype = many edits on the same CDS before
    translation (the haplosaurus / `bcftools csq` model — a capability fastVEP lacks).
    New terms are new predicates, not new branches.
  * **Declarative `OverlapConsequence` includes:** term gates are data (a
    `FeatureOverlap{exon,intron,intron_boundary}` flag set + an `include_satisfied`
    table from Ensembl `Constants.pm`), not scattered `if` guards — e.g.
    `splice_polypyrimidine` needs `exon=0,intron=1`.
  * **`_get_differing_regions` (genomic analog of `CdsEdit`):** a same-length MNV is
    split at its internal *matching* bases, so every splice/intron predicate is
    evaluated over only the actually-changed sub-intervals — e.g. `AACTC/GACCA` hits
    `splice_donor_region` but not the donor 5th base. SNVs/indels stay one region.
* **Root-caused fixes** (each verified against VEP — several by *instrumenting* VEP),
  all locked by a generated regression corpus: mitochondrial codon table + non-ATG
  `start_lost`; splice precedence and insertion handling (`(min,max)` swap, exact
  `codons` window); `intron_variant` / 5′ & 3′ UTR / `stop_lost` co-occurrence as
  unions; `protein_altering_variant` vs clean inframe ins/del; the
  `splice_polypyrimidine` exon gate; CDS-boundary straddle suppression; the
  `frameshift` stop-codon guards (a stop-deleting indel is `stop_lost` alone, not
  `frameshift`); the `coding_unknown` `X` guard (an N-padded `cds_start_NF` first codon
  is `coding_sequence_variant`, not synonymous); the **frameshift-intron 12 bp exon
  stretch** (a transcript with a ≤12 bp intron treats near-boundary intronic variants
  as exonic, suppressing a spurious polypyrimidine-tract call).
* **Regression corpus is mandatory:** every fixed divergence is captured — unit
  tests + `test/sql/vep_splice.test` + `test/data/regression_cases.tsv` (generated
  from the concordance dump by `correctness/gen-regression-cases.sh`, run by
  `test/run-regression-cases.sh`). See `docs/PATCHES.md`.
* **Open frontier (tracked):** `_get_differing_regions` for the alt-longer *delins*
  case (split into minimal del+ins sub-regions, not just the same-length MNV split
  already shipped); stop_retained-at-essential-splice (acceptor vs donor); and
  3′-shifting.

### Haplotype-aware consequence (experimental)

* **`vep_haplotype_consequence(chrom, transcript_id, 'pos:ref:alt;…') → VARCHAR`** —
  combines a set of PHASED variants on one transcript into a single consequence by
  applying them together to the reference CDS before translating once (the
  bcftools-`csq` / Haplosaurus model). Co-located variants merge — an in-codon SNV pair
  becomes one MNV, so e.g. a *silent* SNV flips to *missense* when phased with its
  neighbour — the capability fastVEP lacks. The kernel is the existing multi-edit
  `CodingContext`; merging is in transcript/CDS coordinates (strand-aware,
  intron-collapsed), and grouping by `(sample, haplotype, transcript)` stays in SQL.
* **Hill-climb harness:** `correctness/haplotype_concordance.sh` validates the multi-edit
  path against the *proven* single-variant kernel with no new oracle — a same-length MNV
  is exactly a phased SNV set, so the haplotype of its split components must equal the
  whole MNV's coding terms. Currently **98.8% (1740/1761)** on the ClinVar sample; the
  21 divergences are the `coding_unknown`/X-guard interaction with the multi-edit window.
* **Marked experimental** — known gaps tracked in code: no bcftools-`hap_finalize`
  compound-block *flush* (independent codon edits far apart are over-merged); a net-zero
  indel haplotype (deletion + restoring insertion) reads as non-indel; the string input
  and `VARCHAR` output will become a typed `LIST<STRUCT>` / structured row.

### New SQL functions

* `normalize_variant(pos, ref, alt) → STRUCT(pos, ref, alt)` — canonical minimal
  variant form (right-trim + left-trim), the load-bearing join key matching
  Ensembl VEP's representation. Makes cross-annotator comparison valid for every
  class (incl. indels) and underpins future supplementary-annotation joins.
* HGVS `g./c./p.` on `vep_consequence` and `vep_annotate` (100% concordant with
  fastVEP).
* Structural variants: `vep_consequence` gains an END-aware
  `(chrom, pos, end, ref, alt)` form; `<DEL>/<DUP>/<CNV>/<INV>/<BND>/<CN*>` →
  Ensembl SV consequence vocabulary.

### Tooling & data

* **Server-free Ensembl cache builder** (`correctness/cache-build/`): loads Ensembl's
  published MySQL flat-file dumps and assembles a columnar Parquet cache inheriting
  the curated flags (MANE Select / Plus Clinical, `cds_start_NF`/`cds_end_NF`,
  selenocysteine, TSL, APPRIS, CCDS, regulatory build). Organism/build-agnostic.
* Composes with **duckhts** (community extension) in one DuckDB session.
* Reproducible, structured layout: `benchmarks/` (perf) and `correctness/`
  (concordance + synthetic hard-variant corpus), each rendered from recorded CSVs
  and linked from the README.

### Project

* License changed to **GPL-2.0-or-later** (vendored fastVEP crates stay Apache-2.0).
* Crate/extension version set to **0.3.0**.

## duckvep 0.2.0

* Initial DuckDB-native foundation: `read_vcf`/`vcf_samples`, `vep_consequence`
  (scan-driven scalar) and `vep_annotate`, columnar Parquet transcript cache;
  vendored fastVEP consequence engine. Ensembl-VEP-concordant on SNVs.
