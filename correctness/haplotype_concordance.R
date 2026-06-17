#!/usr/bin/env Rscript
# Haplotype concordance — the hill-climb analog of correctness/vep_concordance.R, for the
# multi-edit (haplotype) path. It needs NO new oracle: a same-length MNV is exactly a phased
# set of single-base SNVs in one transcript window, so the proven single-variant kernel
# (vep_consequence, already ~VEP-116-concordant) IS the oracle.
#
#   for each coding MNV (ref/alt same length, >=2, >=2 changed bases):
#     components = the changed bases as individual SNVs
#     ASSERT  vep_haplotype_consequence(components) == coding terms of vep_consequence(MNV)
#
# A mismatch means the multi-edit CodingContext combination diverges from applying the same
# change as one variant — a real haplotype-kernel bug. Coding-term subset only. Requires the
# full FASTA cache (coding sequences), like Rscript test/run-regression-cases.R.
#
# Usage: correctness/haplotype_concordance.R [sample.vcf]
suppressMessages(library(optparse))
op <- OptionParser(usage = "%prog [sample.vcf]")
pa <- parse_args(op, positional_arguments = c(0, 1))
root   <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
DUCKDB <- Sys.getenv("DUCKVEP_DUCKDB", file.path(root, ".tools/duckdb"))
EXT    <- Sys.getenv("DUCKVEP_EXT",   file.path(root, "build/release/duckvep.duckdb_extension"))
GFF3   <- Sys.getenv("DUCKVEP_GFF3",  file.path(root, "data/giab/GRCh38.116.gff3.gz"))
FASTA  <- Sys.getenv("DUCKVEP_FASTA", file.path(root, "data/giab/GRCh38.primary.fa"))
SAMPLE <- if (length(pa$args)) pa$args[1] else "/tmp/sample.vcf"
THRESHOLD <- as.numeric(Sys.getenv("HAPLO_CONCORDANCE_MIN", "98.0"))

if (!file.exists(FASTA)) {
  message(sprintf("SKIP: FASTA cache not present (%s) — coding haplotypes need it.", FASTA))
  quit(status = 0)
}

## SQL via the version-matched CLI, returning the result rows parsed as a data frame.
sql_csv <- function(sql) {
  f <- tempfile(fileext = ".sql"); writeLines(sql, f); on.exit(unlink(f))
  txt <- system2(DUCKDB, c("-unsigned", "-csv", "-f", f), stdout = TRUE, stderr = FALSE)
  read.csv(text = paste(txt, collapse = "\n"))
}

CODING_MACRO <- "CREATE TEMP MACRO coding_only(lst) AS list_filter(lst, x -> x IN (
  'missense_variant','synonymous_variant','stop_gained','stop_lost','stop_retained_variant',
  'start_lost','start_retained_variant','inframe_insertion','inframe_deletion',
  'frameshift_variant','protein_altering_variant','coding_sequence_variant',
  'incomplete_terminal_codon_variant'));"

r <- sql_csv(sprintf("LOAD '%s';
SELECT vep_load_cache('%s','%s');
%s
WITH mnvs AS (
  SELECT v.chrom AS chrom, v.pos AS pos, v.ref AS r, a.alt AS al
  FROM read_vcf('%s') v, UNNEST(v.alt) a(alt)
  WHERE length(v.ref)=length(a.alt) AND length(v.ref)>=2 AND v.ref<>a.alt
),
comps AS (
  SELECT chrom,pos,r,al,
    string_agg((pos+(i-1))||':'||substring(r,i,1)||':'||substring(al,i,1), ';' ORDER BY i) AS comp_str
  FROM mnvs, range(1,length(r)+1) t(i)
  WHERE substring(r,i,1) <> substring(al,i,1)
  GROUP BY chrom,pos,r,al
),
mnv_csq AS (
  SELECT m.chrom,m.pos,m.r,m.al, c.transcript_id,
    list_sort(coding_only(c.consequence)) AS mnv_coding
  FROM mnvs m, UNNEST(vep_consequence(m.chrom,m.pos,m.r,m.al)) u(c)
),
joined AS (
  SELECT mc.transcript_id, mc.mnv_coding,
    list_sort(string_split(vep_haplotype_consequence(mc.chrom, mc.transcript_id, cp.comp_str),'&')) AS hap
  FROM mnv_csq mc JOIN comps cp USING(chrom,pos,r,al)
  WHERE len(mc.mnv_coding) > 0
)
SELECT count(*) AS tested,
  count(*) FILTER (WHERE hap = mnv_coding) AS concordant,
  count(*) FILTER (WHERE hap <> mnv_coding) AS discordant,
  round(100.0*count(*) FILTER (WHERE hap = mnv_coding)/count(*), 3) AS pct
FROM joined;", EXT, GFF3, FASTA, CODING_MACRO, SAMPLE))

cat(sprintf("haplotype concordance (multi-edit vs proven MNV kernel): %d/%d = %s%% (%d divergent)\n",
            r$concordant, r$tested, format(r$pct), r$discordant))
if (is.na(r$pct) || r$pct < THRESHOLD) {
  message(sprintf("FAIL: haplotype concordance %s%% < %s%% threshold", format(r$pct), THRESHOLD))
  quit(status = 1)
}

# ── Stronger check: validate against the REAL Ensembl VEP oracle (the controlled dump), not
# just our own MNV kernel. A same-length MNV's split components, applied as one haplotype, must
# equal VEP's coding consequence for the whole MNV. Skipped if no dump.
dumps <- Sys.glob(file.path(root, "data/vep_dumps/*/annotations.parquet"))
if (length(dumps)) {
  DUMP <- dumps[order(file.info(dumps)$mtime, decreasing = TRUE)][1]
  o <- sql_csv(sprintf("LOAD '%s'; SELECT vep_load_cache('%s','%s');
  %s
  WITH mnvs AS (SELECT v.chrom, v.pos, v.ref r, a.alt al FROM read_vcf('%s') v, UNNEST(v.alt) a(alt)
      WHERE length(v.ref)=length(a.alt) AND length(v.ref)>=2 AND v.ref<>a.alt),
  comps AS (SELECT chrom,pos,r,al, string_agg((pos+(i-1))||':'||substring(r,i,1)||':'||substring(al,i,1),';' ORDER BY i) cs
      FROM mnvs, range(1,length(r)+1) t(i) WHERE substring(r,i,1)<>substring(al,i,1) GROUP BY chrom,pos,r,al),
  -- the dump now carries the ORIGINAL identity (opos,oref,oalt) for VEP rows, so join directly on it
  -- (no normalize_variant round-trip / '-' alignment guesswork — the pi P0-5 keying fix).
  vep AS (SELECT opos,oref,oalt,transcript_id, list_sort(coding_only(string_split(consequence,'&'))) vc
          FROM read_parquet('%s') WHERE source='vep'),
  joined AS (SELECT
      list_sort(string_split(vep_haplotype_consequence(c.chrom, vep.transcript_id, c.cs),'&')) hap, vep.vc vep_terms
    FROM vep JOIN comps c ON c.pos=vep.opos AND c.r=vep.oref AND c.al=vep.oalt WHERE len(vep.vc)>0)
  SELECT count(*) AS tested, count(*) FILTER (WHERE hap=vep_terms) AS concordant,
         round(100.0*count(*) FILTER (WHERE hap=vep_terms)/count(*),3) AS pct FROM joined;",
    EXT, GFF3, FASTA, CODING_MACRO, SAMPLE, DUMP))
  cat(sprintf("haplotype concordance (multi-edit vs REAL Ensembl VEP oracle): %d/%d = %s%%\n",
              o$concordant, o$tested, format(o$pct)))
}
