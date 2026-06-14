# Ensembl cache builder

Builds duckvep's transcript + regulatory cache **directly from Ensembl**, so we
inherit Ensembl's curated knowledge (MANE tags, `cds_start_NF`/`cds_end_NF` flags,
the regulatory build) as queryable columns instead of re-deriving it — lossily —
from GFF3. Ensembl VEP stays the *oracle we validate against*; this is the
*source we build from*.

## Why flat-file dumps, not live MySQL

Ensembl exposes a public MySQL server (`ensembldb.ensembl.org`), and DuckDB's
`mysql` extension can `ATTACH` it — but bulk pulls hit "Server has gone away"
flakiness and the scanner can lose database scoping. Ensembl **also publishes the
same data as flat-file dumps** (`pub/release-N/mysql/<db>/`): one headerless TSV
per table (`<table>.txt.gz`, NULLs as `\N`) plus the `<db>.sql.gz` schema for
column order. We download the handful of tables we need (~120 MB for human, vs the
27 GB VEP cache) and load them with `read_csv` — **no server, fully reproducible,
offline after download.** This is also the exact pattern for importing any
supplementary-annotation source (ClinVar/gnomAD/dbSNP/COSMIC, or the Ensembl
`variation` DB's `phenotype`/`phenotype_feature`) into Parquet.

## Usage

```sh
correctness/cache-build/build-cache.sh [species] [release] [build] [host] [port]
```

Organism- and build-agnostic — the database name is built from the args, so this
works for **every species and assembly in Ensembl**:

| target | command |
|--------|---------|
| GRCh38 (default) | `build-cache.sh homo_sapiens 116 38` |
| GRCh37 (frozen)  | `build-cache.sh homo_sapiens 113 37` |
| mouse            | `build-cache.sh mus_musculus 116 39` |
| one chromosome (fast iteration) | `CHROM=17 build-cache.sh homo_sapiens 116 38` |

Output (gitignored) lands in `data/cache/<species>.<release>.<build>/`:
`transcripts.parquet`, `exons.parquet`, `translations.parquet`, `regulatory.parquet`
(plus `raw/` dumps and an `ensembl.duckdb` staging db).

## Files

- `build-cache.sh` — downloads dumps + loads them into a local DuckDB (stage 1),
  then runs `assemble.sql` (stage 2).
- `assemble.sql` — pure-local joins turning the raw Ensembl tables into the
  columnar cache. This is where curated flags become boolean columns.

## What you get for free

`transcripts.parquet` carries, per transcript, Ensembl's curated knowledge as
columns:

- **Selection / reporting:** `canonical`, `mane_select`, `mane_plus_clinical`,
  `gencode_basic`, `gencode_primary`, `ccds`, `tsl` (transcript support level),
  `appris` (principal-isoform tier).
- **Incomplete CDS:** `cds_start_nf`, `cds_end_nf`, `mrna_start_nf`, `mrna_end_nf`,
  `upstream_atg`, `readthrough` (the flags the consequence engine should consult
  instead of inferring from phase).
- **Correctness-critical (translated sequence diverges from naive genomic
  translation — a naive engine mis-calls these):** `selenocysteine` (UGA→Sec),
  `stop_codon_readthrough`, `rna_edit`, `amino_acid_sub`.

(GRCh38 r116: 19,299 MANE Select, 14,507 `cds_start_NF`, 115k CCDS, 104
selenocysteine, 14 stop-codon-readthrough, 5 RNA-edited.) `regulatory.parquet`
carries the regulatory build (promoter/enhancer/CTCF/open chromatin/TFBS) with SO
terms.

## Known refinements (TODO)

- **RefSeq xrefs** and the `otherfeatures` RefSeq gene models (for `--merged`).
- **MySQL-dump text escaping**: `read_csv` does not de-escape MySQL's `\t`/`\n`/`\\`
  in free-text columns; we only consume id/coord/code/stable_id columns, which are
  escape-free. Revisit if importing description-heavy tables.
- Engine integration: a loader that builds the in-memory transcript provider from
  this columnar cache (currently the engine reads its own GFF3-derived Parquet).
