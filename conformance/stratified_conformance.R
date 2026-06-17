#!/usr/bin/env Rscript
# Stratified VEP-conformance report — the STATISTICAL tier of the conformance framework.
# Data wrangling in DuckDB SQL (over the dated `annotations.parquet` dump from
# correctness/vep_concordance.R: VEP-116 --gff + duckvep + fastVEP, one row per
# (source, pos, ref, alt, transcript_id) with the SO-term set); statistics + report in R.
#
# Per (consequence class x variant type x length bin) stratum: N shared (variant,transcript)
# pairs, duckvep/fastVEP discordances vs the VEP oracle, and the exact Clopper-Pearson 95%
# upper bound on the duckvep discordance rate (binom.test) — the provable statement
# ("matches VEP at < U @95%, N=..."). At 0 discordances this is the rule of three (~3/N), so
# certifying a stratum at <1e-5 needs ~3e5 examples IN that stratum.
#
# Usage: conformance/stratified_conformance.R [annotations.parquet]   (default: newest dump)
suppressMessages({ library(duckdb); library(DBI) })

args <- commandArgs(trailingOnly = TRUE)
root <- system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE)
dump <- if (length(args) >= 1) args[1] else
  system2("bash", c("-c", shQuote(sprintf("ls -t %s/data/vep_dumps/*/annotations.parquet | head -1", root))),
          stdout = TRUE)
if (length(dump) == 0 || !file.exists(dump))
  stop("no annotations.parquet dump found — run correctness/vep_concordance.R first")
outdir <- file.path(root, "conformance", "data"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

con <- dbConnect(duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE))

# Stratify in SQL: variant type + signed length bin from the minimal alleles (VEP writes '-'
# for the empty side), consequence class = VEP's SO-term set; compare on shared pairs.
sql <- sprintf("
  WITH src AS (SELECT source,pos,ref,alt,transcript_id,consequence FROM read_parquet('%s')),
  v AS (SELECT pos,ref,alt,transcript_id,consequence vc FROM src WHERE source='vep'),
  d AS (SELECT pos,ref,alt,transcript_id,consequence dc FROM src WHERE source='duckvep'),
  f AS (SELECT pos,ref,alt,transcript_id,consequence fc FROM src WHERE source='fastvep'),
  shape AS (
    SELECT v.vc, d.dc, f.fc,
      CASE WHEN ref='-' THEN 'ins' WHEN alt='-' THEN 'del'
           WHEN length(ref)=1 AND length(alt)=1 THEN 'snv'
           WHEN length(ref)=length(alt) THEN 'mnv'
           WHEN length(alt)>length(ref) THEN 'ins'
           WHEN length(alt)<length(ref) THEN 'del' ELSE 'delins' END AS var_type,
      ((CASE WHEN alt='-' THEN 0 ELSE length(alt) END) -
       (CASE WHEN ref='-' THEN 0 ELSE length(ref) END)) AS net
    FROM v JOIN d USING(pos,ref,alt,transcript_id) JOIN f USING(pos,ref,alt,transcript_id))
  SELECT vc AS consequence_class, var_type,
    CASE WHEN net=0 THEN '0'
         WHEN abs(net) BETWEEN 1 AND 3 THEN (CASE WHEN net<0 THEN '-' ELSE '+' END)||abs(net)::VARCHAR
         WHEN abs(net) BETWEEN 4 AND 10 THEN (CASE WHEN net<0 THEN '-' ELSE '+' END)||'4..10'
         WHEN abs(net) BETWEEN 11 AND 50 THEN (CASE WHEN net<0 THEN '-' ELSE '+' END)||'11..50'
         ELSE (CASE WHEN net<0 THEN '-' ELSE '+' END)||'>50' END AS length_bin,
    count(*) n,
    count(*) FILTER (WHERE dc<>vc) dv_discordant,
    count(*) FILTER (WHERE fc<>vc) fv_discordant
  FROM shape GROUP BY 1,2,3 ORDER BY dv_discordant DESC, n DESC", dump)
df <- dbGetQuery(con, sql)

# Exact Clopper-Pearson 95% upper bound on the discordance rate (native binom.test).
upper95 <- function(k, n) if (n == 0) 1 else binom.test(k, n, conf.level = 0.95)$conf.int[2]
df$dv_upper95 <- mapply(upper95, df$dv_discordant, df$n)
write.csv(df, file.path(outdir, "stratified_conformance.csv"), row.names = FALSE)

N <- sum(df$n); K <- sum(df$dv_discordant); FK <- sum(df$fv_discordant)
cat(sprintf("Stratified VEP-116 conformance (duckvep) — %d strata, %d shared (variant,transcript) pairs\n",
            nrow(df), N))
cat(sprintf("  overall: duckvep %d/%d concordant (%d discordant, <= %.2e @95%%)  |  fastVEP %d discordant\n",
            N - K, N, K, upper95(K, N), FK))
# show the discordant strata, then the highest-N certified-clean strata
show <- rbind(df[df$dv_discordant > 0, ],
              head(df[df$dv_discordant == 0, ][order(-df$n[df$dv_discordant == 0]), ], 8))
fmt <- "  %-42s %-6s %-7s %7s %7s %9s %7s\n"
cat(sprintf(fmt, "consequence_class", "type", "len", "N", "dv_disc", "dv_<=95%", "fv_disc"))
for (i in seq_len(min(22, nrow(show)))) with(show[i, ],
  cat(sprintf(fmt, substr(consequence_class, 1, 42), var_type, length_bin,
              n, dv_discordant, sprintf("%.2e", dv_upper95), fv_discordant)))
cat(sprintf("  full table -> %s/stratified_conformance.csv\n", outdir))
