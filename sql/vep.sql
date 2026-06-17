-- duckvep SQL surface: annotation as joins the optimizer plans, not a black-box
-- table function over a path. See docs/kernel-algorithm.md §6/§8.
--
-- Contract — `src_sql` is ANY query yielding ONE ROW PER ALT with:
--   chrom    VARCHAR
--   pos      UBIGINT   -- 1-based
--   end_pos  UBIGINT   -- full ref span; INFO/END for symbolic SVs
--   ref      VARCHAR
--   alt      VARCHAR   -- a single ALT
-- so the surface is source-agnostic (read_vcf / read_bcf / read_parquet / a
-- staged table / any query), hence `vep_annotate` over a relation, not
-- `annotate_vcf` over a path. Requires a prior vep_load_cache(gff3, fasta).

-- Normalize a VCF into the contract (one row per ALT). Use read_vcf's OWN
-- end_pos — it already resolves INFO/END → SVLEN → precise interval, so symbolic
-- SVs (<DEL>/<CNV>) span correctly. (Recomputing pos+len(ref)-1 here was a bug:
-- it collapses every symbolic SV to a 1 bp anchor.)
CREATE OR REPLACE MACRO vep_variants_from_vcf(path) AS TABLE
  SELECT v.chrom::VARCHAR   AS chrom,
         v.pos::UBIGINT     AS pos,
         v.end_pos::UBIGINT AS end_pos,
         v.ref::VARCHAR     AS ref,
         a.alt::VARCHAR     AS alt
  FROM read_vcf(path) v, UNNEST(v.alt) AS a(alt);

-- Candidate (variant, transcript) pairs via a full-span interval range join.
-- DuckDB plans the inequality predicate as an IEJoin/range join; `read_parquet`
-- of the cache is used (native scan = zone-maps + parallel) rather than the
-- vep_transcripts() table function, which has no statistics (see §8 negative result).
CREATE OR REPLACE MACRO vep_candidate_pairs(src_sql, tx_parquet, dist := 5000) AS TABLE
  SELECT v.chrom, v.pos, v.end_pos, v.ref, v.alt,
         t.transcript_id, t.gene_id, t.gene_symbol, t.biotype, t.canonical
  FROM query(src_sql) v
  JOIN read_parquet(tx_parquet) t
    ON t.chrom     =  v.chrom
   AND t.start     <= v.end_pos + dist
   AND t.end_pos   >= CASE WHEN v.pos > dist THEN v.pos - dist ELSE 1 END;

-- Full annotation. `src_table` is the NAME of a (normalized, ideally materialized
-- + ANALYZEd) variant relation — see vep_input below — so the source is scanned
-- once and DuckDB has statistics. The join is LEAN: only transcript_id flows
-- through the 47M-row join; the STRUCT scalar returns gene/biotype/HGVS, so those
-- columns are not redundantly carried (that redundancy cost ~30s at WGS scale).
-- The scalar returns a nullable STRUCT (not a 1-element LIST), so there is no
-- LATERAL UNNEST and `WHERE c IS NOT NULL` pushes the filter down.
CREATE OR REPLACE MACRO vep_annotate(src_table, tx_parquet, dist := 5000) AS TABLE
  WITH scored AS (
    SELECT v.chrom, v.pos, v.end_pos, v.ref, v.alt,
           -- positions cast to BIGINT to match the scalar signature (coords are
           -- always positive, so the UBIGINT contract narrows safely).
           vep_consequence_pair(v.chrom, v.pos::BIGINT, v.end_pos::BIGINT,
                                v.ref, v.alt, t.transcript_id) AS c
    FROM query_table(src_table::VARCHAR) v
    JOIN read_parquet(tx_parquet) t
      ON t.chrom     =  v.chrom
     AND t.start     <= v.end_pos + dist
     AND t.end_pos   >= CASE WHEN v.pos > dist THEN v.pos - dist ELSE 1 END
  )
  SELECT chrom, pos, end_pos, ref, alt,
         c.transcript_id, c.gene_id, c.gene_symbol, c.biotype, c.canonical,
         c.consequence, c.impact, c.amino_acids, c.codons, c.protein_pos,
         c.hgvsc, c.hgvsp, c.hgvsg
  FROM scored
  WHERE c IS NOT NULL;

-- Haplotype-aware: NOT a per-pair join. Phased variants on the same
-- (sample, phase_set, haplotype, transcript) interact, so they must be reduced
-- together — a partitioned GROUP BY aggregate over the multi-edit kernel. DuckDB's
-- hash-partitioned group-by IS the haplotype-safe parallel tiling (§8).
CREATE OR REPLACE MACRO vep_haplotype_annotate(src_sql, tx_parquet, dist := 5000) AS TABLE
  SELECT v.sample_id, v.phase_set, v.haplotype, v.chrom, t.transcript_id,
         vep_haplotype_consequence(
             v.chrom, t.transcript_id,
             list(struct_pack(pos := v.pos, ref := v.ref, alt := v.alt) ORDER BY v.pos)
         ) AS haplotype_consequence
  FROM query(src_sql) v
  JOIN read_parquet(tx_parquet) t
    ON t.chrom   =  v.chrom
   AND t.start   <= v.end_pos + dist
   AND t.end_pos >= CASE WHEN v.pos > dist THEN v.pos - dist ELSE 1 END
  WHERE v.haplotype IS NOT NULL
  GROUP BY v.sample_id, v.phase_set, v.haplotype, v.chrom, t.transcript_id;

-- Supplementary annotation — POINT: an exact equi-join (the DuckDB-native
-- equivalent of fastVEP's var32 + bloom + block lookup; Parquet zone-maps +
-- bloom filters prune the scan).
CREATE OR REPLACE MACRO vep_join_point_sa(src_sql, sa_sql) AS TABLE
  SELECT v.*, sa.* EXCLUDE (chrom, pos, ref, alt)
  FROM query(src_sql) v
  LEFT JOIN query(sa_sql) sa USING (chrom, pos, ref, alt);

-- Supplementary annotation — REGION: an interval range join (the DuckDB-native
-- equivalent of fastVEP's interval index / cgranges SA), collected per variant.
CREATE OR REPLACE MACRO vep_join_region_sa(src_sql, region_sql) AS TABLE
  SELECT v.chrom, v.pos, v.end_pos, v.ref, v.alt,
         list(struct_pack(name := r.name, value := r.value, start := r.start, end_pos := r.end_pos))
             FILTER (WHERE r.name IS NOT NULL) AS region_annotations
  FROM query(src_sql) v
  LEFT JOIN query(region_sql) r
    ON r.chrom     =  v.chrom
   AND r.start     <= v.end_pos
   AND r.end_pos   >= v.pos
  GROUP BY v.chrom, v.pos, v.end_pos, v.ref, v.alt;
