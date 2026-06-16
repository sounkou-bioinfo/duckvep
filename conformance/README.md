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
   0 divergence  ⇒ PROOF)      N, discordant, 95% upper bound)
```

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
| **class model** | parse `Constants.pm` OverlapConsequence table + `VariationEffect.pm` predicate signatures → the equivalence-class axes | TODO (derive from VEP source) |
| **formal generator** | covering array over the classes on a synthetic transcript (known exons/introns/UTRs/start+stop/codon-table/NF flags) | `../tools/gen_hard_variants.rs` (splice/exon cells; extend to full tiling) |
| **statistical corpus** | gold labels from Ensembl `transcript_variation` MySQL for the release; variant distributions from TOPMed/gnomAD/ClinVar + the pinned GIAB diverse cohort | ClinVar in place; Ensembl-MySQL puller TODO |
| **differential oracle** | VEP `--gff` ⟂ duckvep, exact SO-term-set match per pair | `../correctness/vep_concordance.py` (dumps `annotations.parquet`) |
| **stratified CI report** | per (consequence class × variant type × length bin): N, discordant, Clopper-Pearson 95% upper bound (rule-of-three at k=0) | **`stratified_conformance.R`** ✅ runnable now (DuckDB SQL + `binom.test`) |
| **coverage report** | every reachable class hit? (formal completeness) | TODO (needs the class model + labeler) |

## The statistical statement

For a stratum with **0 discordances out of N**, the 95% Clopper-Pearson upper bound on the
true discordance rate is `1-(α/2)^(1/N) ≈ 3/N` (the *rule of three*). So certifying a stratum
at `< 10⁻⁵` needs `~3·10⁵` examples *in that stratum* — the "minimal number," made precise and
**per class**. Errors concentrate in rare cells (the class model says which), so the corpus
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
