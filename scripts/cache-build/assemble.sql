-- Stage 2 — ASSEMBLE. Pure LOCAL joins over the loaded Ensembl dump tables (no
-- network). Produces the columnar, Ensembl-sourced cache duckvep loads. @OUT@ and
-- @CHROMFILTER@ are filled by build-cache.sh. The curated knowledge becomes
-- queryable columns — the flags we kept rediscovering (cds_start_NF, MANE, …) are
-- inherited from the source of truth, not re-derived (lossily) from GFF3.

-- Per-transcript attrib codes collapsed to a list (incomplete-CDS flags, MANE, …).
CREATE OR REPLACE TEMP VIEW tx_flags AS
  SELECT ta.transcript_id, list(aty.code) AS codes
  FROM transcript_attrib ta JOIN attrib_type aty ON ta.attrib_type_id = aty.attrib_type_id
  GROUP BY ta.transcript_id;

-- One row per transcript on a chromosome: coords, gene, symbol, canonical, MANE,
-- and the incomplete-CDS / selenocysteine flags as booleans.
COPY (
  SELECT t.stable_id AS transcript_id, t.version, t.biotype, t.transcript_id AS internal_id,
         sr.name AS chrom, t.seq_region_start AS start, t.seq_region_end AS "end",
         t.seq_region_strand AS strand,
         g.stable_id AS gene_id, g.biotype AS gene_biotype, x.display_label AS gene_symbol,
         (g.canonical_transcript_id = t.transcript_id) AS canonical,
         coalesce(list_contains(f.codes, 'MANE_Select'), false)        AS mane_select,
         coalesce(list_contains(f.codes, 'MANE_Plus_Clinical'), false) AS mane_plus_clinical,
         coalesce(list_contains(f.codes, 'gencode_basic'), false)      AS gencode_basic,
         coalesce(list_contains(f.codes, 'cds_start_NF'), false)        AS cds_start_nf,
         coalesce(list_contains(f.codes, 'cds_end_NF'), false)          AS cds_end_nf,
         coalesce(list_contains(f.codes, 'mRNA_start_NF'), false)       AS mrna_start_nf,
         coalesce(list_contains(f.codes, 'mRNA_end_NF'), false)         AS mrna_end_nf,
         coalesce(list_contains(f.codes, '_selenocysteine'), false)     AS selenocysteine
  FROM transcript t
  JOIN gene g ON t.gene_id = g.gene_id
  JOIN seq_region sr ON t.seq_region_id = sr.seq_region_id
  JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id AND cs.name = 'chromosome'
  LEFT JOIN xref x ON g.display_xref_id = x.xref_id
  LEFT JOIN tx_flags f ON f.transcript_id = t.transcript_id
  WHERE 1=1 @CHROMFILTER@
  ORDER BY chrom, start
) TO '@OUT@/transcripts.parquet' (FORMAT parquet);

-- Exons ordered within transcript by rank (the structure the engine splices).
COPY (
  SELECT t.stable_id AS transcript_id, et.rank,
         e.seq_region_start AS start, e.seq_region_end AS "end",
         e.phase, e.end_phase, et.exon_id
  FROM exon_transcript et
  JOIN exon e ON et.exon_id = e.exon_id
  JOIN transcript t ON et.transcript_id = t.transcript_id
  JOIN seq_region sr ON t.seq_region_id = sr.seq_region_id
  JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id AND cs.name = 'chromosome'
  WHERE 1=1 @CHROMFILTER@
  ORDER BY transcript_id, et.rank
) TO '@OUT@/exons.parquet' (FORMAT parquet);

-- Translation: CDS start/end exon + offsets, protein id/version (for HGVSp).
COPY (
  SELECT t.stable_id AS transcript_id, tl.stable_id AS protein_id, tl.version AS protein_version,
         tl.start_exon_id, tl.end_exon_id, tl.seq_start, tl.seq_end
  FROM translation tl JOIN transcript t ON tl.transcript_id = t.transcript_id
) TO '@OUT@/translations.parquet' (FORMAT parquet);

-- Regulatory build (promoter / enhancer / CTCF / open chromatin / TFBS).
COPY (
  SELECT rf.stable_id, sr.name AS chrom, rf.seq_region_start AS start, rf.seq_region_end AS "end",
         ft.name AS feature_type, ft.so_term, ft.so_accession
  FROM regulatory_feature rf
  JOIN feature_type ft ON rf.feature_type_id = ft.feature_type_id
  JOIN seq_region sr ON rf.seq_region_id = sr.seq_region_id
  ORDER BY chrom, start
) TO '@OUT@/regulatory.parquet' (FORMAT parquet);

SELECT 'transcripts' AS cache, count(*) AS rows FROM read_parquet('@OUT@/transcripts.parquet')
UNION ALL SELECT 'exons',        count(*) FROM read_parquet('@OUT@/exons.parquet')
UNION ALL SELECT 'translations', count(*) FROM read_parquet('@OUT@/translations.parquet')
UNION ALL SELECT 'regulatory',   count(*) FROM read_parquet('@OUT@/regulatory.parquet');
