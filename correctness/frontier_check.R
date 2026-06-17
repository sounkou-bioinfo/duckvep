#!/usr/bin/env Rscript
# Fast frontier check — the 5 ClinVar insertions (≈92 variant/transcript pairs) that are
# duckvep-SPECIFIC vs controlled VEP-116 (fastVEP gets them right; duckvep over-calls splice/
# intron terms and inframe_insertion->protein_altering because it does not trim/3'-shift the
# allele before computing the consequence interval). Runs VEP --gff + duckvep on the tiny fixture
# (seconds, not the 10-min full run) and prints the per-pair divergence count, so the 3'-shift fix
# can be iterated tightly. Lower is better; 0 = the frontier is closed.
#
# Usage: correctness/frontier_check.R   (needs the controlled GFF + FASTA + built release extension)
suppressMessages(library(optparse))
pa <- parse_args(OptionParser(usage = "%prog [options]"), positional_arguments = 0)
root  <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
DUCKDB <- file.path(root, ".tools/duckdb")
EXT   <- Sys.getenv("DUCKVEP_EXT", file.path(root, "build/release/duckvep.duckdb_extension"))
CTRL  <- Sys.getenv("DUCKVEP_CTRL_GFF", file.path(root, "data/giab/GRCh38.116.controlled.gff3.gz"))
FASTA <- Sys.getenv("DUCKVEP_FASTA", file.path(root, "data/giab/GRCh38.primary.fa"))
VCF   <- file.path(root, "correctness/data/frontier_duckvep_specific.vcf")
vep_cmd <- strsplit(Sys.getenv("VEP_CMD", "conda run -n vep vep"), " ")[[1]]
if (!file.exists(FASTA) || !file.exists(CTRL) || !file.exists(EXT)) { message("SKIP: need controlled GFF + FASTA + release extension"); quit(status = 0) }

sql_csv <- function(sql) { f <- tempfile(fileext = ".sql"); writeLines(sql, f); on.exit(unlink(f))
  read.csv(text = paste(system2(DUCKDB, c("-unsigned", "-csv", "-f", f), stdout = TRUE, stderr = FALSE), collapse = "\n")) }

# VEP --gff -> NDJSON, parsed by DuckDB (keyed to the original input variant).
vj <- tempfile(fileext = ".json"); unlink(vj)
vrc <- system2(vep_cmd[1], c(vep_cmd[-1], "-i", VCF, "--gff", CTRL, "--fasta", FASTA, "--distance", "5000",
        "--json", "-o", vj, "--force_overwrite", "--no_stats"), stdout = FALSE, stderr = FALSE)
if (vrc != 0 || !file.exists(vj) || file.info(vj)$size == 0) stop(sprintf("VEP failed (exit %s)", vrc))

# duckvep + VEP in one DuckDB session, compare per (original variant, transcript).
cmp <- sql_csv(sprintf("LOAD '%s'; CREATE TEMP TABLE _l AS SELECT vep_load_cache('%s','%s') AS x;
  WITH vep AS (
    SELECT CAST(split_part(input,chr(9),2) AS BIGINT) opos, split_part(input,chr(9),4) oref,
           split_part(input,chr(9),5) oalt, tc.transcript_id tx,
           list_aggregate(list_sort(tc.consequence_terms),'string_agg','&') vc
    FROM read_json('%s', format='newline_delimited', sample_size=-1), UNNEST(transcript_consequences) u(tc)),
  dv AS (
    SELECT v.pos opos, v.ref oref, a.alt oalt, c.transcript_id tx,
           list_aggregate(list_sort(c.consequence),'string_agg','&') dc
    FROM read_vcf('%s') v, UNNEST(v.alt) a(alt), UNNEST(vep_consequence(v.chrom,v.pos,v.ref,a.alt)) u(c) GROUP BY ALL)
  SELECT count(*) AS pairs, count(*) FILTER (WHERE vc=dc) AS agree, count(*) FILTER (WHERE vc<>dc) AS divergent
  FROM vep JOIN dv USING(opos,oref,oalt,tx);", EXT, CTRL, FASTA, vj, VCF))

cat(sprintf("frontier check (5 duckvep-specific insertions): %d/%d pairs match VEP, %d divergent\n",
            cmp$agree, cmp$pairs, cmp$divergent))
# the per-pair detail when there are divergences, to guide the fix
if (cmp$divergent > 0) {
  d <- sql_csv(sprintf("LOAD '%s'; CREATE TEMP TABLE _l AS SELECT vep_load_cache('%s','%s') AS x;
    WITH vep AS (SELECT CAST(split_part(input,chr(9),2) AS BIGINT) opos, split_part(input,chr(9),5) oalt, tc.transcript_id tx,
           list_aggregate(list_sort(tc.consequence_terms),'string_agg','&') vc
       FROM read_json('%s', format='newline_delimited', sample_size=-1), UNNEST(transcript_consequences) u(tc)),
      dv AS (SELECT v.pos opos, a.alt oalt, c.transcript_id tx, list_aggregate(list_sort(c.consequence),'string_agg','&') dc
       FROM read_vcf('%s') v, UNNEST(v.alt) a(alt), UNNEST(vep_consequence(v.chrom,v.pos,v.ref,a.alt)) u(c) GROUP BY ALL)
    SELECT opos, tx, vc AS vep, dc AS duckvep FROM vep JOIN dv USING(opos,oalt,tx) WHERE vc<>dc ORDER BY opos LIMIT 8;",
    EXT, CTRL, FASTA, vj, VCF))
  for (i in seq_len(nrow(d))) cat(sprintf("  %s %s  VEP=%s  duckvep=%s\n", d$opos[i], d$tx[i], d$vep[i], d$duckvep[i]))
}
quit(status = if (cmp$divergent > 0) 1 else 0)
