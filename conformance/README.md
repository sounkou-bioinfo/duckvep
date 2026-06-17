# Conformance framework — proving duckvep ≡ Ensembl VEP

Two tiers that answer different questions and **meet at a coverage-guided differential
fuzzer**. The oracle in both is Ensembl VEP 116 run `--gff` on the *same* gene model the
engine reads (so only the engine differs). The unit is the **SO-term SET** per
`(variant, transcript)` — both *emission* (is the pair produced at all) and *consequence*
(which terms).

```
            class model (from VEP source)
                     │ defines the axes / equivalence classes
        ┌────────────┴─────────────┐
  FORMAL tier                 STATISTICAL tier
  enumerative generation      corpus sampling
  (covering array over the    (Ensembl MySQL gold +
   classes, synthetic tx)      TOPMed / gnomAD / ClinVar)
        └────────────┬─────────────┘
                     ▼
        ONE differential fuzzer:  VEP(--gff)  ⟂  duckvep   →  exact set match?
                     │   coverage signal = predicate-equivalence-class vector
        ┌────────────┴─────────────┐
  coverage report             stratified CI report
  (every reachable class hit, (per class × type × length bin:
   0 divergence ⇒ falsified?)   N, discordant, 95% upper bound)
```

> **Status: this is a differential FALSIFIER + empirical bound, NOT a proof.** The "class
> coverage ⇒ identity" step is sound only once the engine evaluates VEP's *predicates* (the
> [[port-faithfulness-pivot-to-TVA]]); today the class is VEP's *output* SO-term set, so the
> windowed-peptide surrogate can pass a class yet diverge within it. The oracle is VEP **`--gff`
> controlled mode** — it does not exercise full cache-mode curated attributes (NMD/NF/seleno/
> codon-table/regulatory). Treat every number here as "0 found, bounded," never "proven 0."

## Why they meet at fuzzing

A *differential fuzzer* feeds inputs to two implementations and flags any disagreement. Its
**coverage signal** is what distinguishes the two tiers — and here that signal is the
**predicate-equivalence-class vector** (region/allele/peptide/transcript-context predicates,
derived from VEP's `Constants.pm` + `VariationEffect.pm`). The FORMAL tier *generates to
maximize* class coverage (a covering array — one witness per reachable cell, t-wise for
interactions); the STATISTICAL tier *samples a real distribution* and *measures* which
classes it hits. Same loop, same oracle, same coverage metric — only the generator and the
stopping criterion differ. Coverage-guided mutation (steer toward unhit cells) unifies them.

## Soundness (the crux)

"Pass on the witnesses ⇒ identical on the class" is a **proof only if the engine genuinely
computes VEP's predicates** — a *surrogate* (the current windowed-peptide model) can match
the witnesses yet diverge *within* a cell (e.g. the KDM5A 313154 frameshift-stop was a
within-cell divergence the categories didn't witness). So the formal proof **requires the
TVA-coordinate-layer port** (predicates ported 1:1 from VEP — see `../docs/refinements.md`).
Until then this framework is a strong *falsifier + bound*, not yet a proof.

## Components & status

| component | what | status |
|---|---|---|
| **class model** | parse `Constants.pm` OverlapConsequence table + `VariationEffect.pm` predicate signatures → the equivalence-class axes | **partial** — the 41-term spec (SO term, accession, IMPACT, rank, tier, predicate) is extracted from VEP 116's `Constants.pm` into [`data/so_consequences.tsv`](data/so_consequences.tsv) by [`extract_so_spec.pl`](extract_so_spec.pl); deriving the predicate-equivalence classes from the `VariationEffect.pm` sub signatures is the remaining step |
| **formal generator** | tile equivalence classes on a real transcript (exon/intron boundaries, splice, start/stop) × allele shapes (SNV×3, 1bp ins, 1/2bp del frameshift, 3bp del in-frame, 2bp MNV) | **`generate_witnesses.R`** ✅ R-native via **Rduckhts** (GFF structure + FASTA ranges; duckhts stable C API loads in any DuckDB build). TP53 → 968 witnesses / 36 classes. Still single-transcript: extend to multi-transcript/opposite-strand + a full covering array. |
| **statistical corpus** | gold labels from Ensembl `transcript_variation` MySQL for the release; variant distributions from TOPMed/gnomAD/ClinVar + the pinned GIAB diverse cohort | ClinVar in place; Ensembl-MySQL puller TODO |
| **differential oracle** | VEP `--gff` ⟂ duckvep, exact SO-term-set match per pair | `../correctness/vep_concordance.R` (dumps `annotations.parquet`) |
| **stratified CI report** | per (consequence class × variant type × length bin): N, discordant, Clopper-Pearson 95% upper bound (rule-of-three at k=0) | **`stratified_conformance.R`** ✅ runnable now (DuckDB SQL + `binom.test`) |
| **differential fuzzer** | VEP --gff ⟂ duckvep on the witness VCF; all engines keyed to the **original witness identity** (VEP `input` / fastVEP bridged via VEP's orig→norm map), **union** denominator (emission misses/extras count), exact set match per pair, per CLASS, fastVEP-gated `duckvep_specific` | **`fuzz_witnesses.R`** ✅ R-native. TP53 968 witnesses → 39,048 union pairs: **71 divergences, 0 emission miss/extra, 0 duckvep-specific** — all `start_codon_*` + `stop_codon_del1`, and fastVEP diverges identically on each (the genuine start/stop frontier). The earlier "1,854 boundary-indel" tail was a normalized-key artifact, fixed by original-identity keying. |

## The statistical statement

For a stratum with **0 discordances out of N**, the **two-sided** 95% Clopper-Pearson upper
bound (what `binom.test(0,N)$conf.int[2]` returns) is `1-(α/2)^(1/N) ≈ 3.7/N` (the one-sided
rule-of-three is the looser `≈3/N` from `1-α^(1/N)`; we report the two-sided number the code
computes). So certifying a stratum at `< 10⁻⁵` needs `~3.7·10⁵` examples *in that stratum* —
the "minimal number," made precise and **per class**. This bounds the discordance rate **on
the sampled corpus's distribution only**: the samples are deterministic hash draws from
filtered ClinVar/GIAB and `(variant,transcript)` pairs cluster, so the CI is an empirical
corpus bound, not a guarantee over all deployable variants. Errors concentrate in rare cells (the class model says which), so the corpus
**over-samples** the rare large-indel / boundary / SV strata where uniform sampling
under-powers. Current run (`make`-able from the N=50000 ClinVar dump): the common SNV classes
are certified `≤ ~1–9·10⁻⁵ @95%`; the discordances are all isolated in rare large-indel/SV
strata (N=1–19, CI≈1.0) — honestly *under-powered*, the targets for over-sampling and for the
formal witnesses.

## Run

```sh
# statistical tier (reads the newest correctness dump). DuckDB SQL for the data, R for the
# stats/report (exact Clopper-Pearson via binom.test) — no hand-rolled CIs.
Rscript conformance/stratified_conformance.R
# -> conformance/data/stratified_conformance.csv  (class × type × length-bin, N, discordant, 95% bound)
```
