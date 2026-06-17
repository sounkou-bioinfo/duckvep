#!/usr/bin/env Rscript
# Ensembl-VEP concordance with dated Parquet annotation dumps — R port (no Python).
#
# Annotates a deterministic sample of biallelic variants with (a) controlled Ensembl VEP 116
# (--gff on the SAME gene model the engines read), (b) duckvep, and (c) fastVEP; writes all
# three as a dated Parquet dump and the recorded concordance CSVs. Data plane = DuckDB SQL via
# the version-matched `.tools/duckdb` CLI (the unstable duckvep extension loads only in its own
# v1.5.3 build — the dev duckdb R pkg reports a git-hash version and is rejected). VEP/fastVEP
# are external processes (system2); their JSON is parsed with jsonlite. Stats/reports are R.
#
# Usage: correctness/vep_concordance.R <vcf> <gff3(.gz)> <fasta> [-n N]
# Dumps: data/vep_dumps/<YYYY-MM-DD>/annotations.parquet  (+ the correctness/data/*.csv reports)
suppressMessages({ library(optparse); library(jsonlite) })

op <- OptionParser(usage = "%prog <vcf> <gff3> <fasta> [options]")
op <- add_option(op, c("-n", "--n"), type = "integer", default = 100L, help = "sample size [%default]")
op <- add_option(op, "--oracle", default = Sys.getenv("VEP_ORACLE", "gff"), help = "VEP oracle: gff|cache [%default]")
op <- add_option(op, "--fork", default = Sys.getenv("VEP_FORK", ""), help = "VEP --fork workers")
op <- add_option(op, "--buffer", default = Sys.getenv("VEP_BUFFER", "25000"), help = "VEP --buffer_size")
pa <- parse_args(op, positional_arguments = 3)
vcf <- pa$args[1]; gff3 <- normalizePath(pa$args[2]); fasta <- normalizePath(pa$args[3]); N <- pa$options$n
oracle <- pa$options$oracle

root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
DUCKDB <- file.path(root, ".tools/duckdb"); EXT <- file.path(root, "build/release/duckvep.duckdb_extension")
DATE <- as.character(Sys.Date()); OUTDIR <- file.path(root, "data/vep_dumps", DATE)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
SAMPLE_VCF <- "/tmp/sample.vcf"
rel <- Sys.getenv("VEP_CACHE_VERSION", "116")
vep_cmd <- strsplit(Sys.getenv("VEP_CMD", "conda run -n vep vep"), " ")[[1]]
msg <- function(...) cat(..., "\n", file = stderr())

## SQL via the version-matched CLI. Pass via a temp -f FILE (not -c) so no shell re-parses the
## SQL's parens/quotes. `unsigned` allows the locally built extension.
sql_run  <- function(sql, unsigned = TRUE, csv = FALSE) {
  f <- tempfile(fileext = ".sql"); writeLines(sql, f); on.exit(unlink(f))
  a <- c(if (unsigned) "-unsigned", if (csv) "-csv", "-f", f)
  system2(DUCKDB, a, stdout = TRUE, stderr = if (csv) FALSE else TRUE)
}
sql_csv  <- function(sql) read.csv(text = paste(sql_run(sql, csv = TRUE), collapse = "\n"))

## 1. DETERMINISTIC sample (no RNG): biallelic ACGT variants from read_vcf, ordered by a hash
## of (chrom,pos,ref,alt) so any engine/run reproduces the same set, then write a sorted VCF.
msg(sprintf("sampling %d biallelic variants from %s", N, vcf))
samp_csv <- "/tmp/sample_rows.csv"
sql_run(sprintf("LOAD '%s';
  COPY (
    SELECT chrom, pos, ref, a.alt AS alt
    FROM read_vcf('%s') v, UNNEST(v.alt) AS a(alt)
    WHERE len(v.alt)=1 AND regexp_full_match(v.ref,'[ACGT]+') AND regexp_full_match(a.alt,'[ACGT]+')
          AND v.ref <> a.alt
    ORDER BY hash(chrom||':'||pos||':'||ref||'>'||a.alt) LIMIT %d
  ) TO '%s' (HEADER);", EXT, vcf, N, samp_csv))
samp <- read.csv(samp_csv, colClasses = c(chrom = "character"))
n_indel <- sum(nchar(samp$ref) != nchar(samp$alt))
msg(sprintf("sampled %d variants (%d indels)", nrow(samp), n_indel))
samp <- samp[order(samp$chrom, samp$pos), ]
contigs <- sort(unique(samp$chrom))
con <- file(SAMPLE_VCF, "w")
writeLines(c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s>", contigs),
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), con)
writeLines(sprintf("%s\t%d\t.\t%s\t%s\t.\t.\t.", samp$chrom, samp$pos, samp$ref, samp$alt), con)
close(con)

## 2. Controlled VEP --gff input: rewrite ncRNA_gene->gene (VEP's include-list omits it), sort,
## bgzip, tabix — the SINGLE gene model all three tools read.
build_controlled_gff <- function() {
  out <- file.path(root, "data/giab", sprintf("GRCh38.%s.controlled.gff3.gz", rel))
  if (file.exists(paste0(out, ".tbi")) && file.info(out)$mtime >= file.info(gff3)$mtime) return(out)
  msg(sprintf("building controlled VEP --gff input -> %s", out))
  sh <- sprintf("set -euo pipefail; zcat '%s' | grep -v '^#' | awk -F'\\t' 'BEGIN{OFS=\"\\t\"} $3==\"ncRNA_gene\"{$3=\"gene\"} {print}' | sort -k1,1 -k4,4n -S1G | bgzip > '%s' && tabix -p gff '%s'",
                gff3, out, out)
  if (system2("bash", c("-c", sh)) != 0) stop("controlled GFF build failed")
  out
}
MODEL_GFF <- if (oracle == "gff") build_controlled_gff() else gff3

## 3. Ensembl VEP (offline, controlled --gff, forked) -> JSON -> canonical rows.
fork <- if (nzchar(pa$options$fork)) pa$options$fork else as.character(max(1, parallel::detectCores() - 4))
gene_model <- if (oracle == "gff") c("--gff", MODEL_GFF, "--fasta", fasta) else
  c("--offline", "--cache", "--dir_cache", file.path(root, "data/vep_cache"),
    "--cache_version", rel, "--species", "homo_sapiens", "--assembly", "GRCh38", "--fasta", fasta)
vrc <- system2(vep_cmd[1], c(vep_cmd[-1], "-i", SAMPLE_VCF, gene_model, "--distance", "5000",
        "--symbol", "--json", "-o", "/tmp/vep_off.json", "--fork", fork, "--buffer_size",
        pa$options$buffer, "--force_overwrite", "--no_stats"), stdout = FALSE, stderr = FALSE)
if (vrc != 0 || !file.exists("/tmp/vep_off.json") || file.info("/tmp/vep_off.json")$size == 0)
  stop(sprintf("VEP failed (exit %s) or produced empty JSON — refusing to write a partial dump", vrc))
# Stream the JSON-lines; one row per (variant, transcript). Key on the ORIGINAL input variant
# (VEP echoes the VCF line in `input`: pos=2, ref=4, alt=5) — NOT VEP's normalized output —
# so a deletion/insertion that VEP and duckvep left/right-align differently still compares as
# the same pair (pi P0-5). vmap (orig -> VEP normalized start/allele_string) bridges fastVEP.
vmap <<- list()
vep_rows <- function() {
  out <- list(); recs <- stream_in(file("/tmp/vep_off.json"), verbose = FALSE)
  for (i in seq_len(nrow(recs))) {
    rec <- recs[i, ]; inf <- strsplit(rec$input, "\t")[[1]]
    opos <- as.integer(inf[2]); oref <- inf[4]; oalt <- inf[5]
    vmap[[length(vmap) + 1]] <<- data.frame(opos = opos, oref = oref, oalt = oalt,
      vkey = paste(rec$start, rec$allele_string), stringsAsFactors = FALSE)
    tcs <- rec$transcript_consequences[[1]]; if (is.null(tcs) || !length(tcs)) next
    nref <- strsplit(rec$allele_string, "/")[[1]][1]
    for (j in seq_len(nrow(tcs))) {
      tc <- tcs[j, ]
      out[[length(out) + 1]] <- data.frame(source = "vep", date = DATE, opos = opos, oref = oref, oalt = oalt,
        pos = rec$start, ref = nref, alt = if (!is.null(tc$variant_allele)) tc$variant_allele else "",
        transcript_id = tc$transcript_id, gene_symbol = if (!is.null(tc$gene_symbol)) tc$gene_symbol else "",
        consequence = paste(sort(tc$consequence_terms[[1]], method = "radix"), collapse = "&"),
        impact = if (!is.null(tc$impact)) tc$impact else "", stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}
vr <- vep_rows(); vmap <- unique(do.call(rbind, vmap))
msg(sprintf("VEP %s (%s GRCh38): %d (variant,transcript) rows", oracle, rel, nrow(vr)))
stream_out(vr, file("/tmp/vep_raw.json"), verbose = FALSE)

## 4. duckvep (CLI, one session): annotate the sample. Carry the ORIGINAL (v.pos,v.ref,a.alt)
## as the comparison key (opos/oref/oalt) AND the engine-normalized alleles (pos/ref/alt, for
## shape + the haplotype oracle). -> JSON.
sql_run(sprintf("LOAD '%s';
SELECT vep_load_cache('%s', '%s');
COPY (
  WITH dv AS (
    SELECT v.pos AS opos, v.ref AS oref, a.alt AS oalt,
           normalize_variant(v.pos, v.ref, a.alt) AS nv, c.transcript_id, c.gene_symbol, c.consequence, c.impact
    FROM read_vcf('%s') v, UNNEST(v.alt) AS a(alt), UNNEST(vep_consequence(v.chrom, v.pos, v.ref, a.alt)) AS u(c))
  SELECT 'duckvep' AS source, '%s' AS date, opos, oref, oalt, nv.pos AS pos,
         CASE WHEN nv.ref='' THEN '-' ELSE nv.ref END AS ref,
         CASE WHEN nv.alt='' THEN '-' ELSE nv.alt END AS alt,
         transcript_id, gene_symbol,
         list_aggregate(list_sort(consequence),'string_agg','&') AS consequence, impact
  FROM dv
) TO '/tmp/dv.json' (FORMAT json);", EXT, MODEL_GFF, fasta, SAMPLE_VCF, DATE))

## 5. fastVEP (the underlying engine) on the same sample -> JSON -> rows.
FASTVEP <- Sys.getenv("FASTVEP", file.path(root, "../DuckfastVEP/target/release/fastvep"))
fv_present <- file.exists(FASTVEP)
if (fv_present) {
  raw <- system2(FASTVEP, c("annotate", "-i", SAMPLE_VCF, "--gff3", MODEL_GFF, "--fasta", fasta,
                            "--output-format", "json"), stdout = TRUE, stderr = FALSE)
  arr <- tryCatch(fromJSON(paste(raw, collapse = "\n"), simplifyVector = FALSE),
                  error = function(e) stop("fastVEP JSON parse failed (not silently skipped): ", conditionMessage(e)))
  # fastVEP has no `input` field but mirrors VEP's normalization; bridge its (start,allele_string)
  # back to the original identity via VEP's vmap so all three engines share the original key.
  vlk <- setNames(seq_len(nrow(vmap)), vmap$vkey)
  out <- list()
  for (rec in arr) for (tc in rec$transcript_consequences) {
    parts <- strsplit(rec$allele_string, "/")[[1]]
    mi <- vlk[[paste(rec$start, rec$allele_string)]]; if (is.null(mi)) next   # no VEP anchor -> can't key to original
    out[[length(out) + 1]] <- data.frame(source = "fastvep", date = DATE,
      opos = vmap$opos[mi], oref = vmap$oref[mi], oalt = vmap$oalt[mi],
      pos = rec$start, ref = parts[1], alt = parts[length(parts)], transcript_id = tc$transcript_id,
      gene_symbol = if (!is.null(tc$gene_symbol)) tc$gene_symbol else "",
      consequence = paste(sort(unlist(tc$consequence_terms), method = "radix"), collapse = "&"),
      impact = if (!is.null(tc$impact)) tc$impact else "", stringsAsFactors = FALSE)
  }
  fv <- if (length(out)) do.call(rbind, out) else NULL
  if (!is.null(fv)) { stream_out(fv, file("/tmp/fv.json"), verbose = FALSE); msg(sprintf("fastVEP: %d rows", nrow(fv))) }
}

## 6. Dated Parquet dump + the recorded concordance CSVs (DuckDB SQL, unchanged from the port).
fv_union <- if (fv_present && file.exists("/tmp/fv.json")) "UNION ALL BY NAME SELECT * FROM read_json('/tmp/fv.json')" else ""
dpath <- file.path(root, "correctness/data")
# All concordance/emission joins key on the ORIGINAL identity (opos,oref,oalt,transcript_id),
# not each engine's normalized alleles (pi P0-5) — so indels the engines align differently
# still compare as the same pair, and the EXCEPT-based emission audit counts true misses/extras.
summary_sql <- sprintf("
CREATE TABLE ann AS
  SELECT * FROM read_json('/tmp/vep_raw.json', columns={source:'VARCHAR',date:'VARCHAR',opos:'BIGINT',oref:'VARCHAR',oalt:'VARCHAR',pos:'BIGINT',ref:'VARCHAR',alt:'VARCHAR',transcript_id:'VARCHAR',gene_symbol:'VARCHAR',consequence:'VARCHAR',impact:'VARCHAR'})
  UNION ALL BY NAME SELECT * FROM read_json('/tmp/dv.json') %s;
COPY (SELECT * FROM ann ORDER BY opos, transcript_id, source) TO '%s/annotations.parquet' (FORMAT parquet);
CREATE TABLE vv AS SELECT opos,oref,oalt,transcript_id,consequence,impact,
  CASE WHEN length(oref)=1 AND length(oalt)=1 THEN 'snv' WHEN length(oalt)>length(oref) THEN 'ins'
       WHEN length(oref)>length(oalt) THEN 'del' ELSE 'mnv' END AS class
  FROM ann WHERE source='vep';
CREATE TABLE pairs AS SELECT e.source AS engine, vv.impact, vv.class, vv.consequence AS vep_csq, e.consequence AS eng_csq
  FROM (SELECT * FROM ann WHERE source<>'vep') e JOIN vv USING (opos,oref,oalt,transcript_id);
COPY (SELECT '%s' AS date, engine, impact, class, %d AS n_variants, count(*) AS pairs,
        count(*) FILTER (WHERE vep_csq=eng_csq) AS agree,
        round(100.0*count(*) FILTER (WHERE vep_csq=eng_csq)/nullif(count(*),0),4) AS pct
      FROM pairs WHERE engine IN ('duckvep','fastvep') GROUP BY ALL ORDER BY engine, impact, class
) TO '%s/concordance_by_impact.csv' (HEADER, FORMAT csv);
COPY (WITH t AS (SELECT engine, impact, unnest(string_split(vep_csq,'&')) AS so_term, (vep_csq<>eng_csq) AS disc
        FROM pairs WHERE engine IN ('duckvep','fastvep'))
      SELECT '%s' AS date, engine, so_term, impact, count(*) AS pairs, count(*) FILTER (WHERE disc) AS discordant,
        round(1e5*count(*) FILTER (WHERE disc)/count(*)) AS per100k
      FROM t GROUP BY ALL HAVING count(*) >= 20 ORDER BY engine, discordant DESC
) TO '%s/discordance_by_consequence.csv' (HEADER, FORMAT csv);
COPY (WITH v AS (SELECT opos,oref,oalt,transcript_id,consequence vc,impact FROM ann WHERE source='vep'),
        dd AS (SELECT opos,oref,oalt,transcript_id,consequence dc FROM ann WHERE source='duckvep'),
        ff AS (SELECT opos,oref,oalt,transcript_id,consequence fc FROM ann WHERE source='fastvep')
      SELECT '%s' AS date, v.impact, v.vc AS vep_calls, dd.dc AS duckvep_calls,
        (ff.fc IS NOT DISTINCT FROM v.vc) AS duckvep_specific_regression, count(*) AS n
      FROM v JOIN dd USING(opos,oref,oalt,transcript_id) LEFT JOIN ff USING(opos,oref,oalt,transcript_id)
      WHERE v.vc <> dd.dc GROUP BY ALL ORDER BY n DESC LIMIT 60
) TO '%s/error_transitions.csv' (HEADER, FORMAT csv);
COPY (WITH v AS (SELECT opos,oref,oalt,transcript_id,consequence vc,impact FROM ann WHERE source='vep'),
        dd AS (SELECT opos,oref,oalt,transcript_id,consequence dc FROM ann WHERE source='duckvep'),
        ff AS (SELECT opos,oref,oalt,transcript_id,consequence fc FROM ann WHERE source='fastvep'),
      j AS (SELECT v.impact, string_split(v.vc,'&') AS vt, string_split(dd.dc,'&') AS dt, coalesce(string_split(ff.fc,'&'),[]) AS ft
        FROM v JOIN dd USING(opos,oref,oalt,transcript_id) LEFT JOIN ff USING(opos,oref,oalt,transcript_id) WHERE v.vc<>dd.dc),
      t AS (SELECT impact, list_filter(vt, x -> NOT list_contains(dt,x)) AS vep_only,
              list_filter(dt, x -> NOT list_contains(vt,x)) AS dv_only, ft FROM j)
      SELECT '%s' AS date, impact,
        coalesce(list_aggregate(list_sort(vep_only),'string_agg','&'),'(none)') AS vep_terms,
        coalesce(list_aggregate(list_sort(dv_only),'string_agg','&'),'(none)') AS duckvep_terms,
        (len(list_filter(vep_only, x -> NOT list_contains(ft,x)))=0 AND len(list_filter(dv_only, x -> list_contains(ft,x)))=0) AS duckvep_specific_regression,
        count(*) AS n FROM t GROUP BY ALL ORDER BY n DESC LIMIT 60
) TO '%s/so_term_transitions.csv' (HEADER, FORMAT csv);
COPY (WITH v AS (SELECT opos,oref,oalt,transcript_id,consequence vc FROM ann WHERE source='vep'),
        f AS (SELECT opos,oref,oalt,transcript_id,consequence fc FROM ann WHERE source='fastvep'),
        dd AS (SELECT opos,oref,oalt,transcript_id,consequence dc FROM ann WHERE source='duckvep'),
        vk AS (SELECT DISTINCT opos,oref,oalt,transcript_id FROM v), fk AS (SELECT DISTINCT opos,oref,oalt,transcript_id FROM f),
        dk AS (SELECT DISTINCT opos,oref,oalt,transcript_id FROM dd)
      SELECT * FROM (VALUES
        ('duckvep_discordant_on_shared',(SELECT count(*) FROM v JOIN dd USING(opos,oref,oalt,transcript_id) WHERE vc<>dc)),
        ('duckvep_only_pairs_emission',(SELECT count(*) FROM (SELECT * FROM dk EXCEPT SELECT * FROM vk))),
        ('duckvep_missing_pairs_emission',(SELECT count(*) FROM (SELECT * FROM vk EXCEPT SELECT * FROM dk))),
        ('duckvep_total_divergence',(SELECT count(*) FROM v JOIN dd USING(opos,oref,oalt,transcript_id) WHERE vc<>dc)
           + (SELECT count(*) FROM (SELECT * FROM dk EXCEPT SELECT * FROM vk)) + (SELECT count(*) FROM (SELECT * FROM vk EXCEPT SELECT * FROM dk))),
        ('fastvep_discordant_on_shared',(SELECT count(*) FROM v JOIN f USING(opos,oref,oalt,transcript_id) WHERE vc<>fc)),
        ('fastvep_only_pairs_emission',(SELECT count(*) FROM (SELECT * FROM fk EXCEPT SELECT * FROM vk))),
        ('fastvep_missing_pairs_emission',(SELECT count(*) FROM (SELECT * FROM vk EXCEPT SELECT * FROM fk))),
        ('fastvep_total_divergence',(SELECT count(*) FROM v JOIN f USING(opos,oref,oalt,transcript_id) WHERE vc<>fc)
           + (SELECT count(*) FROM (SELECT * FROM fk EXCEPT SELECT * FROM vk)) + (SELECT count(*) FROM (SELECT * FROM vk EXCEPT SELECT * FROM fk)))
      ) t(metric,value)
) TO '%s/methodology_audit.csv' (HEADER, FORMAT csv);",
  fv_union, OUTDIR, DATE, nrow(samp), dpath, DATE, dpath, DATE, dpath, DATE, dpath, dpath)
sql_run(summary_sql)
msg(sprintf("dump: %s/annotations.parquet  |  reports: %s/*.csv", OUTDIR, dpath))
