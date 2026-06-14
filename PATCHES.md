# Divergences from upstream fastVEP

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

**Impact:** on a genome-wide ClinVar sample these two patches resolved 42 of the
67 (~63%) discordant transcript pairs — every incomplete-CDS case — lifting
duckvep above fastVEP on Ensembl-VEP concordance. The remaining residual is a
subtler start-codon case, the splice polypyrimidine-tract boundary, and
`mature_miRNA_variant` (which needs miRBase data, a join duckvep can add and
fastVEP cannot).

## Feature patches (parity with `fastvep annotate`, kept lean)

- **HGVS g./c./p.** assembled in the DuckDB-free engine (`engine.rs::build_hgvs`)
  directly from `fastvep-hgvs`, bypassing the `fastvep-annotate` god-object.
  100% concordant with fastVEP's own HGVS.
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
