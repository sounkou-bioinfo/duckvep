# Vendored fastVEP crates

These crates are hard-copied (vendored) from **fastVEP**
(<https://github.com/Huang-lab/fastVEP>), Apache-2.0, at commit
`785922ebcaacd3f646d5f1edf374f40f1a39efe5` (v0.2.0). Upstream license:
[`LICENSE-fastVEP.md`](LICENSE-fastVEP.md).

We vendor rather than depend so we can adapt decisions that are suboptimal for a
DuckDB-native, columnar use case (see `../DESIGN.md` §1, §6):

- **Lean closure only.** We vendor `fastvep-core`, `fastvep-genome`,
  `fastvep-cache`, `fastvep-consequence`, `fastvep-hgvs` and deliberately route
  around `fastvep-annotate` (a god-object that also pulls `fastvep-sa`,
  `fastvep-classification`, `fastvep-io`). The engine builds the transcript
  provider + reference + predictor directly and assembles HGVS from the lean
  `fastvep-hgvs` functions itself (see `src/vep/engine.rs::build_hgvs`).
- **Planned divergence:** replace the bincode+zstd transcript cache with Parquet;
  drop the hand-rolled `fastvep-sa` formats in favour of Parquet joins; unify on
  noodles 0.87. Tracked in `DESIGN.md`.

Local modifications, if any, are recorded in this repo's git history.
