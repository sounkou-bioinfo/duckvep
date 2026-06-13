-- Build a nested transcript-model Parquet straight from Ensembl's public MySQL
-- core DB — the SQL-native replacement for GFF3 parsing + the bincode transcript
-- cache (DESIGN.md §1c, §5). No Perl Storable cache parsing.
--
-- Usage:
--   duckdb -unsigned \
--     -c ".read scripts/ensembl-transcripts.sql"
-- Override DB/region with: -c "SET VARIABLE ens_db='homo_sapiens_core_112_38'; SET VARIABLE region_chrom='17'; .read ..."
--
-- Output: data/transcripts.parquet, one row per transcript with a nested exon
-- list. Sorted by (chrom, start) so Parquet row-group stats act as zone maps.
--
-- NOTE: the public ensembldb server rate-limits heavy queries ("Server has gone
-- away") — a whole-genome exon aggregation can drop. For the full build prefer a
-- local Ensembl MySQL mirror; for testing, filter a region first (proven to pull
-- the BRCA1 region → 83 transcripts / 860 exons cleanly).

INSTALL mysql; LOAD mysql;

-- Parameters (defaults; override with SET VARIABLE before .read).
SET VARIABLE ens_db = coalesce(try(getvariable('ens_db')), 'homo_sapiens_core_112_38');

ATTACH 'host=ensembldb.ensembl.org port=5306 user=anonymous database=' ||
       getvariable('ens_db') AS ens (TYPE mysql, READ_ONLY);

COPY (
  WITH exons AS (
    SELECT et.transcript_id,
           list(struct_pack(
                  rank        := et.rank,
                  start       := e.seq_region_start,
                  "end"       := e.seq_region_end,
                  phase       := e.phase,
                  end_phase   := e.end_phase
                ) ORDER BY et.rank) AS exons,
           min(e.seq_region_start) AS tx_start,
           max(e.seq_region_end)   AS tx_end
    FROM ens.exon_transcript et
    JOIN ens.exon e ON e.exon_id = et.exon_id
    GROUP BY et.transcript_id
  )
  SELECT
    t.stable_id                                   AS transcript_id,
    sr.name                                       AS chrom,
    x.tx_start                                    AS start,
    x.tx_end                                      AS "end",
    t.seq_region_strand                           AS strand,
    t.biotype                                     AS biotype,
    g.stable_id                                   AS gene_id,
    xr.display_label                              AS gene_symbol,
    (g.canonical_transcript_id = t.transcript_id) AS canonical,
    tl.stable_id                                  AS protein_id,
    x.exons                                       AS exons
  FROM ens.transcript t
  JOIN ens.seq_region sr ON sr.seq_region_id = t.seq_region_id
  JOIN ens.gene g        ON g.gene_id = t.gene_id
  JOIN exons x           ON x.transcript_id = t.transcript_id
  LEFT JOIN ens.translation tl ON tl.transcript_id = t.transcript_id
  LEFT JOIN ens.xref xr        ON xr.xref_id = g.display_xref_id
  -- Primary assembly only (skip haplotype/patch seq_regions): coord_system rank 1.
  WHERE sr.coord_system_id IN (
    SELECT coord_system_id FROM ens.coord_system WHERE rank = 1
  )
  ORDER BY sr.name, x.tx_start
) TO 'data/transcripts.parquet' (FORMAT parquet);
