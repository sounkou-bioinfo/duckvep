#!/usr/bin/env Rscript
# Generate test/data/regression_cases.tsv from a concordance dump: the SIGNIFICANT real-data
# (variant, transcript) cases where duckvep AGREES with Ensembl VEP, one+ per hard category
# (the categories each accuracy fix targets). The expected value is VEP's call (= duckvep's,
# since concordant), so the corpus is "code as memory" — never hand-typed. Re-run after a fresh
# correctness/vep_concordance.R. Native R: duckdb-r v1.5.3 LOADs the extension (read_vcf +
# normalize_variant) in-process.
suppressMessages({ library(optparse); library(duckdb); library(DBI) })

root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
op <- OptionParser(usage = "%prog [options]")
op <- add_option(op, "--dump",   default = "", help = "annotations.parquet [newest]")
op <- add_option(op, "--sample", default = "/tmp/sample.vcf", help = "sample VCF [%default]")
op <- add_option(op, "--ext",    default = Sys.getenv("DUCKVEP_EXT", file.path(root, "build/release/duckvep.duckdb_extension")))
opt <- parse_args(op)
dump <- if (nzchar(opt$dump)) opt$dump else
  system2("bash", c("-c", shQuote(sprintf("ls -t %s/data/vep_dumps/*/annotations.parquet | head -1", root))), stdout = TRUE)
out <- file.path(root, "test/data/regression_cases.tsv")

con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true")))
on.exit(dbDisconnect(con, shutdown = TRUE))
dbExecute(con, sprintf("LOAD '%s'", opt$ext))
dbExecute(con, sprintf("
  CREATE TEMP TABLE keyed AS
    SELECT v.chrom, v.pos AS opos, v.ref AS oref, a.alt AS oalt,
      normalize_variant(v.pos,v.ref,a.alt).pos AS npos,
      CASE WHEN normalize_variant(v.pos,v.ref,a.alt).ref='' THEN '-' ELSE normalize_variant(v.pos,v.ref,a.alt).ref END AS nref,
      CASE WHEN normalize_variant(v.pos,v.ref,a.alt).alt='' THEN '-' ELSE normalize_variant(v.pos,v.ref,a.alt).alt END AS nalt
    FROM read_vcf('%s') v, UNNEST(v.alt) a(alt)", opt$sample))
dbExecute(con, sprintf("
  COPY (
    WITH dd AS (SELECT pos,ref,alt,transcript_id,consequence dc FROM read_parquet('%s') WHERE source='duckvep'),
         vv AS (SELECT pos,ref,alt,transcript_id,consequence vc FROM read_parquet('%s') WHERE source='vep'),
    j AS (
      SELECT k.chrom, k.opos, k.oref, k.oalt, dd.transcript_id, dd.dc,
        CASE
          WHEN dd.dc LIKE '%%frameshift_variant&stop_gained%%' THEN 'insertion_stop_gained'
          WHEN dd.dc='start_lost' AND k.chrom='MT' THEN 'mito_start_lost'
          WHEN dd.dc LIKE '%%5_prime_UTR_variant%%' AND dd.dc LIKE '%%start_lost%%' THEN 'utr5_start_lost_cooccur'
          WHEN dd.dc LIKE '%%5_prime_UTR_variant%%' AND dd.dc LIKE '%%coding_sequence_variant%%' AND dd.dc NOT LIKE '%%start_lost%%' AND length(k.oref)<>length(k.oalt) THEN 'utr5_straddle_coding_unknown'
          WHEN dd.dc LIKE '%%3_prime_UTR_variant%%' AND (dd.dc LIKE '%%stop_lost%%' OR dd.dc LIKE '%%stop_retained%%') THEN 'utr3_stop_cooccur'
          WHEN (dd.dc LIKE '%%splice_donor_variant%%' OR dd.dc LIKE '%%splice_acceptor_variant%%') AND length(k.oalt) > length(k.oref) AND k.oref<>'-' THEN 'delins_alt_extends_to_splice'
          WHEN dd.dc LIKE '%%protein_altering_variant%%' THEN 'inframe_delins_protein_altering'
          WHEN dd.dc LIKE '%%intron_variant%%splice_donor_variant%%' AND k.nalt='-' THEN 'boundary_del_intron_cooccur'
          WHEN dd.dc LIKE '%%splice_donor_region_variant%%' AND length(k.oref)=length(k.oalt) AND length(k.oref)>1 THEN 'mnv_splice_differing_regions'
          WHEN dd.dc LIKE '%%stop_lost%%' AND dd.dc NOT LIKE '%%frameshift%%' AND k.nalt='-' THEN 'stop_del_frameshift_suppressed'
          WHEN dd.dc LIKE '%%stop_retained_variant%%' AND (dd.dc LIKE '%%splice_acceptor_variant%%' OR dd.dc LIKE '%%splice_donor_variant%%') AND k.nalt='-' THEN 'stop_retained_essential_splice_cil'
          WHEN dd.dc LIKE '%%start_lost%%' AND length(k.nref)=1 AND length(k.nalt)=1 AND length(k.oref)=1 THEN 'snv_start_lost_incl_non_atg'
          WHEN dd.dc LIKE '%%frameshift_variant%%' AND dd.dc LIKE '%%start_lost%%' THEN 'indel_start_lost_reconstruction'
          WHEN dd.dc LIKE '%%non_coding_transcript_variant%%' AND dd.dc NOT LIKE '%%exon%%' AND k.nref='-' THEN 'boundary_insertion_noncoding'
          WHEN dd.dc LIKE '%%coding_sequence_variant%%' AND dd.dc NOT LIKE '%%synonym%%' AND dd.dc NOT LIKE '%%missense%%' AND dd.dc NOT LIKE '%%splice%%' AND dd.dc NOT LIKE '%%intron%%' AND length(k.nref)=1 AND length(k.nalt)=1 THEN 'cds_start_nf_coding_unknown'
          WHEN dd.dc = 'intron_variant&splice_region_variant' AND length(k.nref)=1 AND length(k.nalt)=1 THEN 'splice_region_intronic_no_ppt'
          WHEN dd.dc LIKE '%%splice_polypyrimidine%%' AND length(k.oref)<length(k.oalt) THEN 'insertion_polypyrimidine'
          WHEN dd.dc='stop_gained' THEN 'snv_stop_gained'
          WHEN dd.dc='start_lost' THEN 'nuclear_start_lost'
          ELSE NULL END AS category
      FROM keyed k
      JOIN dd ON k.npos=dd.pos AND k.nref=dd.ref AND k.nalt=dd.alt
      JOIN vv ON dd.pos=vv.pos AND dd.ref=vv.ref AND dd.alt=vv.alt AND dd.transcript_id=vv.transcript_id
      WHERE dd.dc = vv.vc),
    ranked AS (SELECT *, row_number() OVER (PARTITION BY category ORDER BY opos) rn FROM j WHERE category IS NOT NULL)
    SELECT chrom, opos AS pos, oref AS ref, oalt AS alt, transcript_id, dc AS expected, category
    FROM ranked WHERE rn<=2 ORDER BY category, opos
  ) TO '%s' (DELIMITER E'\t', HEADER)", dump, dump, out))
cat(sprintf("wrote %s (%d cases)\n", out, length(readLines(out)) - 1L))
