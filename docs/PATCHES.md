# Divergences from upstream fastVEP

> **The patches as files (reproducible):** the complete divergence is captured as
> per-crate diffs against the pristine base **Huang-lab/fastVEP@785922e** in
> [`vendor/patches/`](vendor/patches/) (`*.patch`, `git apply -p1` from the repo
> root). Regenerate with [`vendor/patches/regen-patches.sh`](vendor/patches/regen-patches.sh).
> This prose is the *rationale*; the patch files are the *source of truth*.

duckvep vendors the lean compute closure of [fastVEP](https://github.com/Huang-lab/fastVEP)
(`fastvep-core`, `-genome`, `-cache`, `-consequence`, `-hgvs`; see
[`vendor/NOTICE.md`](vendor/NOTICE.md)) and **patches it**. Two principles drive
the patches:

1. **Our accuracy oracle is Ensembl VEP itself, not fastVEP.** Where fastVEP
   disagrees with Ensembl VEP, we fix the engine toward VEP — so duckvep is
   *more* concordant with Ensembl VEP than upstream fastVEP is.
2. **duckvep has advantages fastVEP does not** — a columnar SQL data plane
   (variants from `read_vcf`/`read_parquet`/any relation; annotation DBs as
   Parquet/DuckDB tables joined by the optimizer), composition with the
   [duckhts](https://github.com/sounkou-bioinfo/duckhts) community extension
   (rich GFF attribute `MAP`, BCF/CRAM), a WASM target, and set-based QC in SQL.
   These raise the accuracy ceiling: annotations fastVEP's fixed pipeline can't
   reach (e.g. miRBase mature-miRNA regions) are a join away.

## Architecture: peptide-level `CodingContext` (haplotype-ready)

The largest structural change. Upstream `predict_coding_consequence` returned a
**single** consequence chosen by `if/else`, so co-occurring Ensembl terms
(`frameshift_variant&stop_gained`, the generic `coding_sequence_variant`, the
intron co-occurrence on boundary indels) and MNV codon handling all required
special cases — and computing a frame off the *un-normalized* allele kept
re-introducing anchor-base bugs. Replaced with the shape Ensembl actually uses:

- A **`CdsEdit`** = `{cds_idx, ref_bases, alt_bases}` in transcript orientation,
  derived from the **normalized** (anchor-trimmed) allele — so the reading frame
  is correct by construction.
- **`CodingContext`** applies a *set* of edits to the reference CDS, then
  translates → `ref_pep`/`alt_pep` over the affected window. Built **once**.
- **`coding_consequence_terms`** is a flat predicate set that **collects every
  applicable** term (Ensembl's `@OverlapConsequences`), not one.

**Why a set of edits, not `(ref, alt)`:** one variant is one edit; a *phased
haplotype* is several edits applied to the **same** CDS before translation — the
Ensembl haplosaurus / `bcftools csq` model. Co-located variants on one haplotype
therefore combine into one correct protein consequence (an in-codon MNV is the
degenerate "local haplotype"). fastVEP cannot do haplotype-aware consequence;
here it falls out of the abstraction for free (`test_haplotype_two_edits_one_codon`).
Result: MNV `missense↔synonymous` discordances **325 → 0**, with **0**
duckvep-specific regressions vs Ensembl VEP. Two boundary rules are kept exact: a
coding variant reaching an essential splice site straddles the exon/intron edge →
generic `coding_sequence_variant` (VEP can't determine the peptide); insertion
`stop_gained` is deferred (needs VEP's exact insertion window).

## Accuracy-driven patches (toward Ensembl VEP)

These close systematic disagreements found by genome-wide concordance against
offline Ensembl VEP 112 (the residual was ~0.01%, and it was *not* noise — it
was concentrated in incomplete-CDS transcripts).

### Incomplete 5′ CDS — `start_lost` over-calling (`cds_start_NF`)
`vendor/fastvep-consequence/src/predictor.rs`, `predict_coding_consequence`.
Ensembl flags transcripts whose annotated start is not a real ATG
(`cds_start_NF`, encoded as a non-zero phase on the first CDS). Upstream fastVEP
unconditionally calls `start_lost` for any first-codon change. Ensembl VEP never
does on these transcripts. Patch: when `codon_table_start_phase != 0`, codon 0 is
the incomplete leading codon → emit `coding_sequence_variant`; when only the
`cds_start_NF` flag is set, suppress `start_lost` and fall through to the normal
missense/synonymous decision.

### Incomplete terminal codon (`cds_end_NF`)
Same file. When the CDS length is not a multiple of 3, the last codon is
incomplete and cannot be translated. Upstream translated the partial codon into a
spurious `missense_variant`. Patch: detect the codon that runs past the
translateable sequence and emit `incomplete_terminal_codon_variant`, paired with
`coding_sequence_variant` at the call site (matching VEP's two-term output).

### `start_lost` follows the Met initiator, not a canonical-ATG test
Same file. Upstream emitted `start_lost` whenever the *alt* first codon was not a
start. Ensembl's `start_lost` (`VariationEffect.pm`) instead compares peptides: the
annotated first codon of a complete CDS *is* the initiator, which always encodes
Met — independent of the triplet (in the vertebrate-mito table `ATT`/`ATC` are
valid starts that *translate* to Ile). So `start_lost` fires when the alt codon no
longer yields Met, and `start_retained_variant` when it still does. Patch: at codon 0 of a complete
CDS (`!cds_start_NF`, start phase 0) whose reference codon is a real initiator *in
the transcript's codon table* (`CodonTable::is_start` is now table-aware: ATG for the
standard code, plus ATT/ATC/ATA/GTG for mito), emit `start_lost` iff `alt_aa != 'M'`,
else `start_retained_variant`. The table-aware start set both recognises non-ATG mito
starts the prior standard-table-only `is_start(ref) && !is_start(alt)` gate dropped
(`ATT→ATC` at MT-ND2 chrM:4472 is start_lost even though both are mito starts) **and**
keeps a real gate so a first codon that is not a genuine start stays
missense/synonymous. Standard ATG starts are unchanged (no SNV keeps ATG → always
`start_lost`). Regression: `predictor.rs::test_mito_start_lost_non_atg` /
`…_start_retained_to_met`.

### Interval-aware splice predicates (indels/MNVs/haplotypes)
`splice.rs`. Upstream tested a single `genomic_pos` against each splice band;
Ensembl tests the *variant interval* `[r_start, r_end]` (`overlap(...)` in
`BaseTranscriptVariationAllele.pm`). Patch: every splice predicate
(`is_splice_donor`/`acceptor`/`donor_region`/`5th_base`/`polypyrimidine_tract`/
`region`) takes `(start, end)` and uses interval overlap, so an indel spanning a
boundary is classified correctly. (SNVs unchanged: `start == end`.)

### Insertion handling at the polypyrimidine tract / splice region
`splice.rs`. An insertion normalizes to a zero-width interval `start = end + 1`
(VEP's convention). Two predicates need insertion-specific handling to match
`_intron_effects`: (1) **polypyrimidine tract** — VEP swaps the interval to
`(min, max)` before the overlap (`first loop`), widening it by one base each end, so
an insertion touching the window edge IS in the tract; (2) **splice region** — VEP's
`_intron_overlap` adds explicit insertion boundary-touch cases
(`start==intron_start || end==intron_end || start==intron_start+2 || end==intron_end-2`).
Without these, insertions at intronic splice landmarks were under-called to bare
`intron_variant` — the dominant duckvep-specific regression (283 pairs) the
error-transition analysis surfaced. Now eliminated; regressions in
`test/sql/vep_splice.test` (insertion PPT-edge cases).

### Mitochondrial codon table (chrM / NCBI table 2)
`predictor.rs`. Upstream translates everything with `CodonTable::standard()`, so
**every chrM coding variant was mis-called** (TGA→stop instead of Trp, ATA→Ile
instead of Met, AGA/AGG not recognised as stop). Patch: a per-transcript codon
table — `ct(transcript)` returns the vertebrate-mitochondrial table
(`from_ncbi_table(2)`) for chrM/MT transcripts, the standard table otherwise.
Validated on ClinVar chrM vs offline Ensembl VEP (the non-ATG mitochondrial *start*
codons are now handled too — see the `start_lost` patch above).

**Impact — current, valid, version-matched figures are in the generated report
[`correctness/correctness.md`](correctness/correctness.md)** (split by impact ×
class, per-100K error rates, all read from `correctness/data/concordance_by_impact.csv`
so they never drift). Headline: SNVs near-perfect at every impact tier and ahead of
fastVEP; the open frontier is high-impact indels/MNVs (shared with fastVEP — a real
engine gap vs Ensembl VEP, not a measurement artifact).

The remaining engine frontier (tracked, paper-relevant): a frameshift that
introduces a premature stop should add `stop_gained` (duckvep emits only
`frameshift_variant`); and MNV codon handling. These are **shared with fastVEP** vs
Ensembl VEP (both differ from Ensembl) — engine accuracy work, measured per-100K in
the generated report. After the splice-insertion and start-codon patches there are
**zero duckvep-specific regressions** vs the upstream engine in the concordance run:
every remaining discordance is a shared engine gap, not something our patches broke.

## Feature patches (parity with `fastvep annotate`, kept lean)

- **HGVS g./c./p.** assembled in the DuckDB-free engine (`engine.rs::build_hgvs`)
  directly from `fastvep-hgvs`, bypassing the `fastvep-annotate` god-object.
  Exact match with fastVEP's own HGVS on chr17 (HGVSc 19,828/19,828, HGVSp
  9,449/9,449 — agree==total).
- **Structural-variant consequences** wired via `sv_predictor::predict_sv_consequences`
  over the full `INFO/END` span (`engine.rs::annotate_variant_spanned`), with an
  END-aware `vep_consequence(chrom, pos, end, ref, alt)` scalar.
- Exposed `predictor` / `AlleleConsequenceResult` (made `pub`) so the engine can
  route around `AnnotationContext`.

## Performance / engineering patches

- **Columnar Parquet transcript cache** (`vep/tcache.rs`) replacing fastVEP's
  bincode+zstd cache — portable and itself a queryable table.
- **Lock-free engine cache** (`ArcSwapOption`) for concurrent DuckDB worker-thread
  reads (no `RwLock` contention on the hot scalar path).
- **`Arc<str>`** shared categorical fields in `AnnotatedRow` (refcount bump, not
  re-allocation, across millions of rows).
- **Streaming `read_vcf`** (bounded memory) and **magic-byte format detection**
  (content, not file extension).
- Direct-write of the nested `LIST<VARCHAR>` consequence entry array to work
  around duckdb-rs's per-vector (2048) wrapper bound (`vep/consequence.rs`).

All modifications are also recorded in this repository's git history.
