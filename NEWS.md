# duckvep NEWS

Changelog, most recent first. (R-package style.)

## duckvep 0.3.0 (in development)

### Accuracy (oracle = Ensembl VEP, not fastVEP)

* **Concordance is always split by IMPACT × variant class** and per SO term, with
  fastVEP shown alongside — never a single aggregate. duckvep beats fastVEP on
  **every HIGH-impact SO term**. Percentages use enough precision that a non-zero
  discordance never reads as 100%.
* **The consequence engine was rebuilt to better match Ensembl's output structure**, taking
  N=50000 ClinVar divergence vs **controlled** Ensembl VEP 116 (VEP run with `--gff` on the
  *same* gene model the engines read, keyed to the original input variant — so only the engine
  differs) to **239 total divergence** (210 discordant on shared pairs + 29 emission), vs
  fastVEP's **5,063** — duckvep has ~21× fewer divergences than the vendored engine. Of those
  239, **92 are duckvep-specific** (fastVEP matches VEP there): boundary indels where VEP
  3'-shifts the allele before calling consequence and duckvep does not. (An earlier "39 total /
  0 duckvep-specific" headline was a **normalized-key measurement artifact** — VEP and duckvep
  align indels differently, so a normalized-key join compared mismatched pairs and hid these
  regressions as emission. Keying both engines to the original input variant
  (`correctness/vep_concordance.R`) reveals the true counts; see `correctness/correctness.md`.)
  The VEP-faithful abstractions:
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
  * **CDS-boundary straddle `coding_unknown` fallback:** a variant straddling the CDS
    boundary still overlaps the coding sequence, so when the codon-window peptide is
    undeterminable and no boundary term survives, Ensembl emits the generic
    `coding_sequence_variant` (not nothing) — e.g. a 5′UTR-into-CDS frameshift deletion is
    `5_prime_UTR_variant&coding_sequence_variant`.
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
  as exonic, suppressing a spurious polypyrimidine-tract call); and the **`consider_ins_len`
  genomic stop model** (`_overlaps_stop_codon_cil` + `_ins_del_stop_altered_cil` — a
  deletion reaching an essential splice site keeps a determinable `stop_retained`/
  `stop_lost` from the genomic stop-codon overlap, not the blanket `coding_sequence_variant`);
  the **phase-correct non-ATG `start_lost`** (VEP's `_overlaps_start_codon` AND
  `translation_start==1` as separate facts — a non-ATG / non-zero-phase annotated start
  still gets `start_lost`); the **`_ins_del_start_altered` genomic start model** (the start
  analog of `cil_stop` — a length-changing variant overlapping the start codon is `start_lost`
  iff a 5′UTR+CDS reconstruction shows the start codon altered and it is not an inframe
  ins/del, so a 5′UTR-into-start deletion that preserves the start no longer over-fires and a
  frameshift-start deletion is no longer missed); and the **insertion-aware non-coding exon
  membership** (Ensembl
  `non_coding_exon_variant` checks the RAW unsorted variant bounds, so an insertion at an
  exon 5′ edge is `non_coding_transcript_variant`, not `…_exon_variant`).
* **Regression corpus is mandatory:** every fixed divergence is captured — unit
  tests + `test/sql/vep_splice.test` + `test/data/regression_cases.tsv` (generated
  from the concordance dump by `correctness/gen-regression-cases.R`, run by
  `Rscript test/run-regression-cases.R`). See `docs/PATCHES.md`.
* **Boundary predicates are migrating to a `VepAlleleContext`** — a phase-correct
  coordinate/protein projection (cDNA/CDS/protein coords, start/stop-codon windows,
  circular/codon-table aware) that the start/stop predicates query, instead of the
  CDS-codon-index `CodingContext`. This is why the non-ATG `start_lost` landed cleanly
  where three codon-index patches had regressed.
* **Open frontier (tracked, all SHARED with fastVEP — none duckvep-specific):**
  `mature_miRNA_variant` (a feature region not yet in the cache); frameshift / 3′UTR-straddle
  deletions at the stop codon needing VEP's exact `_get_peptide_alleles` window (313154,
  73419638); large multi-exon deletions spanning the start codon (a `start_lost` miss); and
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
* **Hill-climb harness:** `correctness/haplotype_concordance.R` validates the multi-edit
  path against the current single-variant kernel with no new oracle — a same-length MNV
  is exactly a phased SNV set, so the haplotype of its split components must equal the
  whole MNV's coding terms. 100% (1761/1761) **on this MNV-split ClinVar harness**.
* **Marked experimental** — known gaps tracked in code: no bcftools-`hap_finalize`
  compound-block *flush* (independent codon edits far apart are over-merged); a net-zero
  indel haplotype (deletion + restoring insertion) reads as non-indel; the string input
  and `VARCHAR` output will become a typed `LIST<STRUCT>` / structured row.

### New SQL functions

* `normalize_variant(pos, ref, alt) → STRUCT(pos, ref, alt)` — duckvep's minimal
  (anchor-trimmed) variant form, the join key for annotation sources normalized by
  this same function. Note it is **not** a universal cross-engine key: VEP and duckvep
  3'-shift indels differently, so differential VEP conformance keys on the original
  input variant, not the normalized form (see `correctness/vep_concordance.R`).
* HGVS `g./c./p.` on `vep_consequence` and `vep_annotate` (matches fastVEP's recorded
  HGVSc/HGVSp on the chr17 benchmark; see `benchmarks/results.md`).
* Structural variants: `vep_consequence` gains an END-aware
  `(chrom, pos, end, ref, alt)` form; `<DEL>/<DUP>/<CNV>/<INV>/<BND>/<CN*>` →
  Ensembl SV consequence vocabulary.

### Tooling & data

* **Server-free Ensembl cache builder** (`correctness/cache-build/`): loads Ensembl's
  published MySQL flat-file dumps and assembles a columnar Parquet cache inheriting
  the curated flags (MANE Select / Plus Clinical, `cds_start_NF`/`cds_end_NF`,
  selenocysteine, TSL, APPRIS, CCDS, regulatory build). Organism/build-agnostic.
* Can coexist with **duckhts** (community extension) in one DuckDB session; duckvep is
  self-contained and does not depend on it.
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
