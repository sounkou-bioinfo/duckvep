# duckvep — design

DuckDB-native variant effect prediction. A loadable DuckDB extension (Rust,
`duckdb-rs`) that reads genomics formats via **noodles**, exposes the VEP
consequence / HGVS / ACMG engine as SQL functions, and treats annotation
databases as plain Parquet/DuckDB tables joined by the optimizer.

This is the successor to fastVEP. It deletes the hand-rolled file formats and
operations and lets DuckDB do the data engineering it is already good at.

> One doc. No ADR sprawl. Decisions live here and in commit messages.

---

## 1. Thesis: stop re-inventing what a columnar engine already does

The old `fastvep-sa` (~7.7k LOC) is, structurally, a hand-rolled query engine:

| fastSA mechanism | What it actually is | DuckDB / Parquet equivalent |
| --- | --- | --- |
| `kmer16`, `var32`, zigzag | A **packed variant key** for exact `(chrom,pos,ref,alt)` joins | **No key needed** — hash/merge join on the natural columns |
| `block` / `chunk` binning, `bloom` | Hand-rolled **zone maps** to skip data | Parquet row-group **min/max statistics** → predicate pushdown |
| `.osi` interval index | Hand-rolled **range index** | Range join / `ASOF` join / `range_overlap` |
| `.osa`/`.osa2` LRU block cache | Manual decompression cache | Engine buffer pool + columnar projection pushdown |
| per-source `sources/*.rs` builders | One encoder per database | One generic `VCF → Parquet` importer + metadata |

Every one of these is something DuckDB's optimizer chooses automatically once
the data is sorted Parquet. So we delete the lot and let SQL plan the joins.

**The variant key is a pruning device — keep the pruning, drop the format.**
A naive equi-join on all four columns (two of them variable-length `ref`/`alt`
strings) is genuinely not free, and that is exactly what echtvar's packed key
optimizes: it makes pruning cheap. The point is that the pruning belongs in the
**physical layout and query plan**, all expressible in SQL/Parquet, not in a
bespoke on-disk encoding:

- Sort each source by `(chrom,pos)` → Parquet row-group min/max act as zone
  maps; a position predicate skips most row groups before any string compare.
- Join *position-first* (cheap integer/zone-mapped), resolve `ref`/`alt`
  equality only on the small surviving set — the optimizer does this given the
  sorted layout; partitioning by chrom helps further.
- For hot paths, optionally **materialize a generic integer key column**
  (chrom-dictionary + pos, or a hash) so the equi-join is single-column — but as
  an ordinary computed Parquet column, not a human-hardcoded bit layout.

So we keep echtvar's *benefit* (pruning) via sort + zone maps + an optional
generic key column, and delete its *format*. No bespoke variant/region key is
part of the default.

The genuinely irreplaceable part — consequence prediction, HGVS, ACMG — stays
in Rust and is exposed as UDFs.

**Organism-agnostic by construction.** The hand-rolled variant key (kmer16/
Var32, echtvar-style) bakes in *human* assumptions — enumerated chromosomes and
position bit-budgets. We instead join on the natural `(chrom,pos,ref,alt)`
columns, which works for any organism (fastVEP already advertises any-organism
GFF3 and is validated on mouse). If a packed key is ever added purely for
speed, it must be **generic** (chrom-dictionary or hash based), never a
human-hardcoded bit layout.

---

## 1b. Layering: core library → extension → CLI → C API

One DuckDB-free **engine core** (`src/vep/engine.rs` + the vendored fastVEP
crates) with thin frontends over it:

- **Engine core** — consequence/HGVS/ACMG + IO. Zero `duckdb` dependency.
- **DuckDB extension** (`cdylib`) — the only target that touches DuckDB, via
  `duckdb` with the `loadable-extension` feature, which resolves DuckDB symbols
  from the *host* process at load (it doesn't even bundle libduckdb).
- **CLI binary** — links the core only, so **no libduckdb**. A small
  `fastvep annotate`-style compute tool. SQL access is "run `duckdb` and
  `LOAD duckvep`", so the CLI never pulls DuckDB.
- **C API** — another `extern "C"` shim over the core; also no libduckdb.

This is why the engine must be a DuckDB-free module: it keeps libduckdb linkage
isolated to the extension target alone.

**Kernel vs data plane (so DuckDB-free does *not* mean slow scans).** "Engine
core" is really two things, and only the first is DuckDB-free:

- **Compute kernel** — the consequence/HGVS algorithm: `(variant, overlapping
  transcripts, ref_seq) → rows`. Pure CPU (interval lookup, coordinate math,
  codon translation). DuckDB cannot speed this up; it does no scan/join/filter.
  It must take data *in* and do **no IO** — DuckDB calls it, never the reverse.
- **Data plane** — scanning variants, reading the transcript model + reference,
  filtering, joining, parallelism, spilling. This is **DuckDB's job** when we run
  as the extension (vectorized, parallel, pushdown, any source), and noodles/
  Arrow in the standalone CLI.

Keeping the kernel DuckDB-free is what *enables* the good scans: the kernel is a
scalar over DuckDB-scanned rows. The idiomatic path is therefore the scalar
`vep_consequence(chrom,pos,ref,alt)` driven by `read_vcf`/`read_parquet` (§3.2),
**not** a table function that re-reads the VCF itself.

**Implemented:** `vep_load_cache(gff3, fasta)` builds the engine once into a
global; `vep_consequence(chrom,pos,ref,alt)` is a raw `VScalar` over
DuckDB-scanned rows (any source) returning a native **`LIST<STRUCT>`** (typed:
`consequence VARCHAR[]`, `canonical BOOLEAN`, `protein_pos BIGINT`, …), directly
`UNNEST(...) AS u(c)` then `c.consequence`/`c.impact` — verified at parity with
`vep_annotate`/fastVEP over a multi-row scan. `vep_annotate(vcf, …)` remains a
one-call convenience; the scalar is the primary, scan-driven path.

> **No Arrow on the output path.** The scalar writes DuckDB vectors directly via
> the core `ListVector`/`StructVector` API (same as the VTabs). `VArrowScalar`'s
> Arrow→vector conversion mangles nested structs for multi-row batches (struct
> fields come back NULL); we use Arrow only to *read* the input batch.

**Vendored, not depended.** The fastVEP crates are hard-copied under `vendor/`
(Apache-2.0, `Huang-lab/fastVEP@785922e`) so we can fix decisions that are
suboptimal for a columnar/DuckDB use case — chiefly: route around the
`AnnotationContext` god-object (which drags in `fastvep-sa`/`-classification`/
`-io`/`-hgvs`) and build the lean transcript+reference+predictor closure
directly. See `vendor/NOTICE.md`.

## 1c. Ensembl is the source of truth; persist the slices we need into local DuckDB

**Scope decision:** duckvep targets the organisms Ensembl's VEP MySQL (or a
compatible service / local mirror) offers — we do **not** chase general-organism
GFF3. Riding on Ensembl's canonical core schema buys **exact `ensembl-vep` and
`haplosaurus` compatibility by construction**: same transcript set, same
attributes, same HGVS exceptions — not a re-derived approximation. (VEP's own
cache is Perl `Storable` *derived from this same DB*; we skip that format
entirely.)

**Runtime artifact = a local DuckDB database** holding the slices we require,
synced once from Ensembl MySQL. At annotation time duckvep reads the local
slices (`ATTACH … (READ_ONLY)`) — **no MySQL, no network, no duckhts**. This
local DB *is* the cache; it replaces the bincode `transcript_cache.rs` and is
strictly richer.

**Slices to persist** (`correctness/cache-build/`):

- **Coordinates & names:** `coord_system`, `seq_region`, `seq_region_synonym`
  (contig-name aliasing), `seq_region_attrib`.
- **Model:** `gene`, `transcript`, `exon`, `exon_transcript`, `translation`.
- **HGVS exceptions VEP applies:** `transcript_attrib`, `translation_attrib`,
  `attrib_type` (decode) — RNA edits `_rna_edit`, selenocysteine, incomplete CDS
  `cds_start_NF`/`cds_end_NF`, ribosomal slippage.
- **Aliases:** `xref`, `external_synonym`, `object_xref`, `external_db`.

Materialized as DuckDB tables (or nested Parquet), sorted by `(chrom, pos)` for
zone-map pruning. Reference **sequence** stays a bgzip+faidx FASTA (Ensembl's
per-release FASTA) — sequence is the one thing not worth holding in the DB.

**GFF3+FASTA** remains only a degraded fallback for non-Ensembl inputs (core
consequence, no HGVS-exception fidelity).

**Self-contained — duckhts is not a dependency or an integration target.**
duckvep ships its own readers, interval overlap, and variant keys. We may
*borrow patterns* (cgranges-style interval search, key encoding) and reimplement
them here, but nothing requires loading another extension at build or run time.

## 2. Architecture (decided)

- **Engine:** DuckDB loadable extension, `duckvep.duckdb_extension`, built on
  `extension-template-rs` + `duckdb-rs` (scalar UDFs + `VTab` table functions).
  Single optimizer, single vectorized runtime.
- **IO:** Rust port of the **duckhts** reader surface, backed by **noodles**
  (not htslib, not rust-htslib). Pure Rust → also compiles to **WASM**.
- **VEP compute:** keep `fastvep-consequence`, `fastvep-hgvs`,
  `fastvep-classification`, `fastvep-genome`; wrap as UDFs.
- **Annotation sources:** Parquet-primary (range-readable, CDN-friendly), with
  optional `ATTACH 'duckvep.duckdb'` to bundle everything as one artifact.
- **Semi-structured values:** DuckDB **`VARIANT`** (1.5.3+) for VCF `INFO` and
  for annotation payloads — no more JSON strings. `STRUCT`/JSON fallback while
  `VARIANT` is new.
- **Front-ends, both first-class:** the loadable extension *and* a thin
  `fastvep` CLI that `LOAD duckvep` and runs SQL.
- **WASM / no server:** `duckvep.wasm` + DuckDB-WASM + Parquet caches over HTTP
  range requests makes fastVEP.org a static page. Deletes `fastvep-web` in the
  general case.
- **`exon` (wheretrue):** design reference only — its UDTF pattern and
  `chrom_optimizer_rule` region pushdown — **not** a dependency. It is a
  DataFusion engine; we picked DuckDB's optimizer, and we don't run two.

### Why DuckDB over the DataFusion (exon) path
DuckDB gives the optimizer plus the whole ecosystem (Python / CLI / **WASM** /
`httpfs` / Parquet / joins) for zero binding work, which is exactly the stated
goal ("the optimizer works for us"). DataFusion/exon is more ergonomic for
custom Rust operators but makes us own more of the stack. We keep an Arrow
bridge open so a DataFusion component could feed in later, but DuckDB is the
one engine.

---

## 3. API surface

### 3.0 The variant table is the interface (providers conform), genotypes are tidy

`read_vcf` is just one **variant provider**. Two complementary shapes, both
plain DuckDB relations, so any source (VCF, BCF, Parquet/Arrow variant table,
phased or multi-sample provider) conforms and every downstream UDF
(`vep_consequence`, `normalize_variant`, `acmg_classify`) consumes them
source-agnostically.

**(a) Sites** — one row per variant (the contract that annotation joins to):

| column | type | notes |
| --- | --- | --- |
| `chrom` | VARCHAR | |
| `pos` | BIGINT | 1-based start |
| `end_pos` | BIGINT | interval end; `INFO/END` for SV/CNV, else `pos+len(ref)-1` |
| `id` | VARCHAR | |
| `ref` | VARCHAR | |
| `alt` | LIST\<VARCHAR\> | multiallelic / symbolic `<DEL>` / breakend |
| `qual` | DOUBLE | nullable |
| `filter` | LIST\<VARCHAR\> | |
| `info` | VARIANT | (VARCHAR fallback for now) |

**(b) Genotypes are a presentation choice, because annotation is site-wise.**
Consequence/SA depend only on the variant, never on which samples carry it, so
we annotate each *unique site once* and never re-annotate per sample. That frees
genotypes to be carried in whatever layout is cheapest for the query, with no
correctness impact:

- **Default — compact `gt LIST<VARCHAR>` on the site row** (positional, in
  sample order; phasing kept in the `0|1` text). Sample names come from a
  metadata function `vcf_samples(path) -> LIST<VARCHAR>`. No row multiplication,
  and it rides along with the site you already annotated.
- **Tidy long on demand** — `UNNEST` the `gt` list zipped with `vcf_samples(...)`
  to get one row per variant × sample (duckhts-style) for per-sample filtering,
  aggregation, or phased haplotype-aware consequence (§3.5). So tidy is a *view*
  over the compact form, not a separate ingest.
- **Rich per-sample fields** (DP/GQ/AD/PS) — `LIST<STRUCT>` when needed.

Net: **trivially annotate multi-sample BCF/VCF** — annotate sites once, keep
`gt` as a list beside them, explode to tidy only where a query needs it.

### 3.1 IO table functions (noodles, region-aware — duckhts parity)

```sql
read_vcf(path, region := 'chr1:1-1000', samples := true)   read_bcf(...)
read_bam(path, region := ...)   read_cram(...)             read_fastq(path)
read_fasta(path, region := ...) read_gff(path) read_gtf(path) read_bed(path)
read_tabix(path, region := ...)
-- utilities: faidx(), bgzip()/bgunzip(), *_index()
```

- Region syntax `CHROM:START-END`, comma-separated multi-region (duckhts-compatible).
- `INFO` and per-sample `FORMAT` materialize as `VARIANT` (or `STRUCT`).
- Region is an explicit named parameter (not reliant on filter pushdown) in v1;
  pushdown is a later optimization.

### 3.2 VEP compute — scalar UDF + `UNNEST` (not a relation-arg table function)

A variant yields *many* rows (one per transcript/allele), and DuckDB table
functions can't cleanly take a relation as an argument. So consequence is a
**scalar UDF returning `LIST<STRUCT>`**, composed with any row source via
`UNNEST`. The transcript model + reference are loaded once via a setup call and
held in extension state (see §5), so the per-row UDF takes only the variant:

```sql
SELECT 'duckvep.transcripts.parquet' AS cache,
       vep_load_cache('grch38.transcripts.parquet', 'grch38.fa');  -- once

SELECT v.chrom, v.pos, v.ref, alt.alt, c.*
FROM   read_vcf('in.vcf.gz') v,
       UNNEST(v.alt) AS alt(alt),
       UNNEST(vep_consequence(v.chrom, v.pos, v.ref, alt.alt)) AS c;
--   c => STRUCT(gene, transcript, source, consequence LIST<VARCHAR>,
--               impact, hgvsc, hgvsp, canonical, ...)
```

Because consequence is a scalar over plain columns, **variants can come from any
source DuckDB can read** — `read_vcf(...)`, but equally `read_parquet`,
`read_csv`, JSON, Arrow, or an attached DB. `read_vcf` is a convenience, not a
requirement; `SELECT chrom,pos,ref,alt FROM 'variants.parquet'` feeds
`vep_consequence` just as well. That falls out of the design for free.

Scalar helpers compose too: `vep_most_severe(list)`, `hgvs_c(...)`,
`hgvs_p(...)`, `revcomp(seq)`, `region_overlap(chrom, pos, end, 'chr1:..')`.

### 3.3 Annotation = joins (the "trivial to extend" goal)

No per-source Rust. Sources are tables; the optimizer matches them:

```sql
SELECT a.*, clinvar.clnsig, gnomad.af_grpmax
FROM   annotated a
LEFT JOIN clinvar USING (chrom, pos, ref, alt)          -- exact join
LEFT JOIN gnomad  USING (chrom, pos, ref, alt)
LEFT JOIN sv ON range_overlap(a.chrom, a.pos, a.end, sv.chrom, sv.start, sv.end);  -- region join
```

Adding a source = produce a Parquet + one line in the SQL manifest (§4).

### 3.4 Variant normalization (required, not optional)

A join-based annotation design is only correct if both sides are normalized the
**same** way — otherwise `(chrom,pos,ref,alt)` equi-joins silently miss
(`AT/A` vs left-aligned `T/`). So normalization is a first-class, reference-aware
kernel (like the duckhts/bcftools-norm functions), applied at **both**
`fastvep build` time (sources) and query time (variants):

```sql
normalize_variant(chrom, pos, ref, alt, fasta)  -- left-align + parsimony/trim → STRUCT
split_multiallelic(ref, alt)                     -- → LIST of biallelic records
```

Canonicalization rule (left-align + minimal representation) is fixed and
documented so sources built today match queries run later.

### 3.5 Haplotype-aware consequence (bcftools csq-style) — roadmap

Per-variant prediction (VEP default, what fastVEP does today) is the v1 path.
A later phase adds **haplotype-aware** prediction: phase co-located variants
within a transcript (needs phased `GT`) and predict the *combined* protein
effect — adjacent SNVs sharing a codon, or frameshifts that compensate. This is
a capability upgrade to the consequence engine, exposed as a mode flag on
`vep_consequence`, not a new format.

### 3.6 ACMG

A table function over the joined relation: `acmg_classify(relation, config)`,
backed by `fastvep-classification`. Trio/compound-het via extra params.

---

## 4. Source registration: Parquet-KV + declarative SQL manifest (option C)

No YAML. Two layers with clear roles:

- **Parquet footer KV metadata = intrinsic, travels with the file:** which
  columns are `chrom/pos/ref/alt`, source name/version/assembly,
  allele-vs-positional, output key. A source file is self-describing; DuckDB
  reads it via `parquet_kv_metadata()`.
- **Declarative SQL manifest = deployment choices:** which sources are active,
  exact vs. range-overlap join, output aliases. `sources.sql` is *code*
  (`CREATE MACRO` / `CREATE VIEW` / `ATTACH`), versioned and diffable.

`fastvep build` converts a source VCF → sorted Parquet (by `chrom,pos`) with KV
metadata, or into the bundled `duckdb` file. This replaces `sa-build` and every
`sources/*.rs`.

---

## 5. The cache is a native DuckDB file (core-only runtime)

The runtime cache is a **native DuckDB database** (`.duckdb`) holding relational
tables, read with `ATTACH … (READ_ONLY)`. **Runtime requires DuckDB core only —
no Parquet, no MySQL, no extensions** (verified: a GFF-built cache reads back in
plain `duckdb` with no `-unsigned` and nothing loaded). This replaces the bincode
`transcript_cache.rs`.

**Schema (relational, mirrors Ensembl):** `transcripts(transcript_id, chrom,
start, end_pos, strand, biotype, gene_id, gene_symbol, canonical, coding, tsl,
appris, flags[])`, `exons(transcript_id, rank, start, end_pos, phase,
end_phase)`, `translations(...)`, `chrom_alias(alias, chrom, source)`,
plus HGVS-exception attribs on the Ensembl path. Sorted by `(chrom, start)` for
zone-map pruning.

**Two importers, one schema** (the cache *is* the contract):

- **GFF importer** — `read_gff_transcripts(gff3)` (+ `read_gff_exons`) table
  functions reuse `parse_gff3`; build with `CREATE TABLE transcripts AS SELECT *
  FROM read_gff_transcripts('x.gff3')`. Portable/offline, any organism. ✅
- **Ensembl cache builder** — `correctness/cache-build/` (`build-cache.R` +
  `assemble.sql`) loads Ensembl's published flat-file MySQL **dumps** (no live
  server, no flakiness) and assembles a columnar Parquet cache inheriting the
  curated flags (MANE, cds_start_NF/cds_end_NF, selenocysteine, regulatory build).
  Organism/build-agnostic. Exact ensembl-vep fidelity. ✅ (see correctness/cache-build/README.md)

MySQL/Parquet are *build-time only*; either way the runtime artifact is the same
core-only `.duckdb`.

- **Chromosome aliases are first-class.** `chrom_alias` maps input contig names
  (`1`/`chr17`/`NC_000017.11`/`CM000679.2`) to the cache's `chrom` — the naming
  reconciliation ensembl-vep does. Ensembl fills it from `seq_region_synonym`;
  the GFF path takes an optional synonyms file (else chr-prefix normalization).
- **FASTA → bgzip + faidx.** Sequence is the one thing not in the DB: O(1) region
  access via noodles for coding consequences. Also `read_fasta(path, region)`.
- **Runtime:** the kernel loads transcripts from the attached cache once into the
  in-memory interval index (`get_transcripts` overlap), held in extension state
  keyed by `path+mtime`, reused across rows/queries.

Kept verbatim: `fastvep_genome::Transcript` coordinate logic. Only the cache
*substrate* changes (bincode → native DuckDB tables).

**Runtime pattern = ATTACH the cache `.duckdb`** (DuckDB-core ATTACH, no
extension):

```sql
ATTACH 'grch38.cache.duckdb' AS cache (READ_ONLY);   -- transcripts/exons/chrom_alias
ATTACH 'clinvar.duckdb'      AS clinvar (READ_ONLY);  -- annotation sources, same way
-- annotate read_vcf variants against cache.transcripts (kernel loads from the
-- attached cache), then LEFT JOIN clinvar.* — all in one query.
```

Attachments are per-session (re-attach each run) and may be local **or remote**
(`ATTACH 'https://…/grch38.cache.duckdb'` is read-only) — so caches can be
CDN-hosted, which is also what makes the WASM/no-server path work (§2).

---

## 5b. GFF3 conformance (honest)

The GFF importer uses fastVEP's `parse_gff3`, which is **permissive**: it parses
well-formed Ensembl/RefSeq GFF3 correctly but warns-and-skips malformed lines and
defaults a missing `phase` to a sentinel — it is **not** a strict spec validator
(unlike duckhts's `GFFBase`). For our inputs (Ensembl GFF3) this is fine and
matches fastVEP/VEP; strict validation (borrowing duckhts's approach) is a
possible hardening. Benchmarks are `.Rmd`-driven (`benchmarks/results.Rmd`,
duckknit) in the same spirit as duckhts's conformance reports.

## 6. Migration: wrap-then-delete (decided)

Never lose a working product mid-flight.

1. Scaffold `duckvep` from `extension-template-rs`; CI builds native + WASM. ✅
2. **IO:** `read_vcf` / `read_gff` / `read_fasta` on noodles to input parity;
   then `read_bam/bcf/cram/fastq/bed/gtf/tabix`.
3. **Wrap compute:** `vep_consequence`, HGVS scalars, `acmg_classify` as UDFs.
   **Validate output parity** vs current `fastvep annotate` on `tests/` +
   `validation/` + GIAB.
4. **Sources → Parquet:** `fastvep build`; replace `sa-build` + `sources/*.rs`.
5. **Delete only after parity:** `var32/zigzag/bloom/kmer16/chunk/block/`
   `writer*/reader*/index`, `.osa/.osa2/.osi/.oga`, all `sources/*.rs`, custom
   VCF parsing in `fastvep-io`, `transcript_cache.rs`, and `fastvep-web`.
6. **CLI rewire:** `fastvep annotate` becomes a SQL driver over `duckvep`.

**Kept:** consequence, hgvs, classification, genome, GFF3 parsing.
**Deleted (~7k+ LOC):** the hand-rolled formats, per-source builders, custom IO,
web server.

---

## 7. Testing & benchmarking

- **Unit:** Rust tests per reader/UDF.
- **SQL behavior:** DuckDB sqllogictest under `test/sql/`.
- **Living docs:** `README.Rmd` rendered with **duckknit** (rundel/duckknit) —
  example SQL runs against the built extension, outputs are real.
- **Parity:** harness diffing `duckvep` output vs current `fastvep annotate`
  (golden VCF/tab/JSON) on `tests/`, `validation/human`, `validation/mouse`.
- **Full-size GIAB:** fetch scripts for HG001/HG002 (GRCh38) truth VCFs; a
  bench harness (hyperfine for end-to-end, criterion for hot Rust paths)
  measuring throughput vs fastVEP and VEP. Full-size runs are gated (env flag /
  nightly), not on every CI run.
- **CI:** build + test native and WASM targets.

---

## 8. Design-review findings (resolved)

1. **Relation-arg table function won't work.** DuckDB table functions can't take
   a subquery/relation as an argument. → consequence is a **scalar UDF →
   `LIST<STRUCT>` + `UNNEST`**; cache passed via `vep_load_cache()` state. (§3.2)
2. **Per-query cache reload is wasteful.** → load-once into extension state,
   keyed by `path+mtime`. (§5)
3. **Region predicate pushdown is hard in v1.** → explicit `region :=` param
   first (duckhts/exon style); optimizer pushdown later. (§3.1)
4. **Async noodles would drag tokio into the extension.** DuckDB calls are sync.
   → use noodles **sync** readers (bgzf/vcf/...); no embedded runtime.
5. **`VARIANT` is new (1.5.3).** → use it for INFO/annotations, keep a
   `STRUCT`/JSON fallback; confirm WASM support before relying on it.
6. **WASM IO has no local FS.** → caches must be range-readable Parquet over
   `httpfs`; that is already the chosen format, so it aligns.
7. **Unstable C ABI is version-pinned.** Built artifact loads only on its
   `TARGET_DUCKDB_VERSION` (currently v1.5.3). → pin CI's DuckDB to the crate's
   version; document it.
8. **Parity before deletion is load-bearing.** Step 5 must not start until the
   step-3 parity harness is green on GIAB. (§6)
