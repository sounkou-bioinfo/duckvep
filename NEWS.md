# duckvep NEWS

Changelog, most recent first. (R-package style.)

## duckvep 0.3.0 (in development)

### Accuracy (oracle = Ensembl VEP, not fastVEP)

* **Concordance is always split by IMPACT × variant class** — never a single
  aggregate. An average hides high-impact discordances (the clinically actionable
  sites); the harness and the rendered report surface them. Percentages use enough
  precision that a non-zero discordance never reads as 100%.
* **Accuracy patches over the vendored engine** make duckvep *more* concordant
  with Ensembl VEP than fastVEP (which carries the same bugs):
  * `start_lost` requires a real ATG start codon (was firing on non-ATG / `cds_start_NF`
    annotated starts).
  * incomplete 5′ CDS (`cds_start_NF`) → `coding_sequence_variant`, not `start_lost`.
  * incomplete terminal codon (`cds_end_NF`) → `incomplete_terminal_codon_variant`.
  * splice predicates are **interval-aware** (variant `[start,end]` vs the splice
    bands, Ensembl's `overlap()`), so indels/MNVs/haplotypes spanning a boundary
    are classified correctly.
  * Version-matched SNV concordance: 13 discordant in 1,315,496 (99.999%); the open
    frontier is high-impact indels and MNVs. See `PATCHES.md`.

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
