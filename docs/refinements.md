# Refinements & open items

Tracked refinements and future work, kept **out of** [`DESIGN.md`](DESIGN.md) (which
holds the stable design rationale) so the design doc doesn't accumulate drift.
Current, measured correctness figures are in the generated report
[`../correctness/correctness.md`](../correctness/correctness.md); divergences from
upstream fastVEP (with patch files) are in [`PATCHES.md`](PATCHES.md).

## Port faithfulness — the architectural pivot (paper-relevant headline)

A pi port-faithfulness review (vs VEP-116 Perl) found that the engine reached **0
duckvep-specific divergences / 15 consequence / 39 total** vs controlled VEP-116 by
matching VEP's **output** through a *windowed-peptide surrogate + boolean proxies*
(`peptide_defined`, `pep_determinable`, `at_start`, `vep_start_overlap`, `cil_stop_term`),
**not** VEP's actual machinery. The corpus number is a lower bound on faithfulness, not
proof of it. The real road to 100% (and to closing the residual tail at the root) is to
**port VEP's coordinate/codon/peptide layer**, in priority order:

1. **`TranscriptVariationAllele` coordinate + codon/peptide layer** (`BaseTranscriptVariation.pm`,
   `TranscriptVariationAllele.pm`): mapper `Coordinate`/`Gap` cDNA/CDS/protein arrays (not scalar
   `Option` endpoints), allele/feature seq, shift hash, `codon()`, `peptide()`,
   `_get_alternate_cds()`. Structurally fixes stop/frameshift/inframe/partial-codon/HGVS at once.
2. **Start/stop predicates verbatim** (`VariationEffect.pm:851-1541`) → kills the residual tail
   (KDM5A 313154 frameshift-stop on the reverse strand; multi-exon start-codon deletions).
3. **Generate the `OverlapConsequence` include/tier table** from `Constants.pm`; run predicates by
   include/tier, not manual push/filter.
4. **Intron/exon overlap trees + `shift_hash`** before more `splice.rs` geometry patches.
5. **SVs through the same engine** (classify from SVTYPE/class/CN, not ALT strings).
6. Cache the missing biological attributes: codon_table from slice attrs (not chrom name),
   cds_start/end_NF, `_rna_edit`, selenocysteine, mature_miRNA.

Until then: claim **"VEP-concordant on N=50000 ClinVar," not "faithful VEP port."**

## Open correctness gaps (residual tail — all SHARED with fastVEP)

The formal-tier fuzzer (`conformance/fuzz_witnesses.R`, original-identity keyed) pins this tail
precisely on TP53: **71 divergences / 39,048 union pairs, all `start_codon_*` + `stop_codon_del1`,
0 duckvep-specific** (fastVEP diverges from VEP identically). An earlier "1,854 boundary-indel"
figure was a normalized-key measurement artifact (the two engines align indels differently),
removed by keying VEP/duckvep/fastVEP to the original witness identity — so the real frontier is
~26× smaller and is start/stop-codon engine accuracy, not splice-boundary indel bugs.

- **High-impact indel / MNV tail** — a frameshift / 3′UTR-straddle deletion at the stop
  (KDM5A 313154: VEP `stop_lost`, duckvep `stop_retained` — the windowed peptide is too short
  to see the read-through; needs the full `peptide()` from item 1 above), and multi-exon
  start-codon deletions. Shared with fastVEP (engine accuracy, not duckvep-specific).
- **`mature_miRNA_variant`** — a feature region not yet in the cache (a join away).
- **chrX/Y haploid + PAR** — now have a **female sample (HG004) + HG003/HG005 + 1000G** in the
  pinned diverse cohort (`scripts/fetch-data.R`, `DIVERSE=1`); add chrX-diploid + PAR
  regression coverage to the `correctness/` suite.
- **Haplotype experimental edges** — the multi-edit path is 100% on the MNV-split harness but
  over-merges as a single block (no bcftools-`csq` `hap_finalize` flush) and can't yet yield the
  stop exception; needs a real haplosaurus/bcftools-csq oracle before block-scoping.
- **SV `<INS>` mis-routes** (classifies as a small `Insertion`) — give it a structural route;
  validate against the 36.6K real GIAB HG002 SV Tier1 insertions now fetched.

## Cache builder

- **RefSeq xrefs** + the `otherfeatures` RefSeq gene models (for `--merged`).
- **MySQL-dump text escaping** — `read_csv` does not de-escape MySQL's `\t`/`\n`/`\\`
  in free-text columns; we only consume escape-free id/coord/code/stable_id columns.
  Revisit if importing description-heavy tables.
- **Engine loader** for the columnar Ensembl cache — currently the engine reads its
  own GFF3-derived Parquet; a loader for `transcripts/exons/translations.parquet`
  would let the curated flags (codon_table, cds_start_NF, …) feed the predictor
  directly.

## Supplementary annotations — SQL-native plugin system (designed)

ClinVar / gnomAD / dbSNP / COSMIC / dbNSFP / OMIM / conservation — each a Parquet/DuckDB
table joined on the **`normalize_variant`** canonical key (the load-bearing piece is in
place). A pi design consultation settled the **ergonomic API** (the SQL-native replacement
for fastVEP's per-source CLI flags — and VEP's `--plugin`/`--custom`). Three layers:

1. **Import** → normalized typed tables (Parquet with footer `KV_METADATA`: name, kind,
   assembly, normalization version, key). Source classes: `exact_allele`
   (`USING(chrom,pos,ref,alt)`), `interval` (range join), `point_track` (BETWEEN join),
   `gene/transcript` (join to consequence rows). Keys + hot fields are typed columns;
   `STRUCT`/`VARIANT` only for payload tails. (Avoid a universal `(…,field,value)` long table.)
2. **Register** → `CALL vep_register_source(name, relation := …, kind := …)` PRAGMA: validates
   assembly/coords/normalization/columns into a session registry; does not annotate.
3. **Manifest** → a `sources.sql` (`ATTACH` + `CREATE VIEW` + `CALL vep_register_source`) =
   the reproducible plugin config.

Query ergonomics: **hybrid** — a VEP-one-call `annotate_vcf('x.vcf')` TABLE MACRO defined by
the manifest, plus the canonical, most-idiomatic **CTE/join** ("tidyverse-pipe") form. A literal
`|> add_frequencies()` fluent pipe is not worth faking in SQL (no grammar to add; `duckdb-rs` has
no `ParserExtension`) — leave that to R/Python wrappers. A `FROM 'x.vcf'` **replacement scan**
(FFI `duckdb_add_replacement_scan`) is later UX polish, only for `.vcf/.bcf/.bed/.gff`.
Prereqs to build first: `vcf_info_float/int/string(info,key,allele_no)` accessors +
`read_bed`/`read_bigwig` readers. The big-source read side streams bounded-memory via duckhts.
