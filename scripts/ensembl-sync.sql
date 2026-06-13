-- Sync the Ensembl core slices duckvep needs into a LOCAL DuckDB database, once.
-- At runtime duckvep reads this local DB (ATTACH … READ_ONLY) — no MySQL, no
-- network. Riding Ensembl's canonical schema gives exact ensembl-vep /
-- haplosaurus compatibility (DESIGN.md §1c). This local DB *is* the cache.
--
-- Run with the target file as the default database (so CREATE TABLE persists),
-- ens.* being the attached MySQL source:
--   duckdb data/ensembl.duckdb -unsigned -c ".read scripts/ensembl-sync.sql"
--
-- The public ensembldb server rate-limits heavy pulls; for a whole-genome sync
-- use a local Ensembl MySQL mirror, or sync per chromosome by adding a
-- `WHERE sr.name = '<chrom>'` to the `transcript` slice below.

INSTALL mysql; LOAD mysql;
SET VARIABLE ens_db = coalesce(try(getvariable('ens_db')), 'homo_sapiens_core_112_38');
ATTACH 'host=ensembldb.ensembl.org port=5306 user=anonymous database=' ||
       getvariable('ens_db') AS ens (TYPE mysql, READ_ONLY);

-- 1. Model (transcripts first; dependent slices filter by these ids).
CREATE OR REPLACE TABLE transcript AS SELECT * FROM ens.transcript;
CREATE OR REPLACE TABLE gene       AS SELECT * FROM ens.gene
  WHERE gene_id IN (SELECT gene_id FROM transcript);
CREATE OR REPLACE TABLE exon_transcript AS SELECT * FROM ens.exon_transcript
  WHERE transcript_id IN (SELECT transcript_id FROM transcript);
CREATE OR REPLACE TABLE exon AS SELECT * FROM ens.exon
  WHERE exon_id IN (SELECT exon_id FROM exon_transcript);
CREATE OR REPLACE TABLE translation AS SELECT * FROM ens.translation
  WHERE transcript_id IN (SELECT transcript_id FROM transcript);

-- 2. Coordinates & contig-name aliasing.
CREATE OR REPLACE TABLE seq_region AS SELECT * FROM ens.seq_region;
CREATE OR REPLACE TABLE seq_region_synonym AS SELECT * FROM ens.seq_region_synonym;
CREATE OR REPLACE TABLE coord_system AS SELECT * FROM ens.coord_system;

-- 3. HGVS exceptions VEP applies (rna_edit, selenocysteine, incomplete CDS, …).
CREATE OR REPLACE TABLE attrib_type AS SELECT * FROM ens.attrib_type;
CREATE OR REPLACE TABLE transcript_attrib AS SELECT * FROM ens.transcript_attrib
  WHERE transcript_id IN (SELECT transcript_id FROM transcript);
CREATE OR REPLACE TABLE translation_attrib AS SELECT * FROM ens.translation_attrib
  WHERE translation_id IN (SELECT translation_id FROM translation);

-- 4. Gene/transcript aliases & synonyms.
CREATE OR REPLACE TABLE xref AS SELECT * FROM ens.xref;
CREATE OR REPLACE TABLE external_synonym AS SELECT * FROM ens.external_synonym;
CREATE OR REPLACE TABLE external_db AS SELECT * FROM ens.external_db;

DETACH ens;
