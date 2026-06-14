#!/usr/bin/env python3
"""Ensembl-VEP concordance with **dated Parquet annotation dumps**.

Annotates a sample of SNVs with (a) Ensembl VEP and (b) duckvep, writes both as a
dated Parquet dump (accumulating dataset), and appends a concordance summary.

VEP source: the live REST API by default (paced for rate limits), or an offline
VEP cache via `--vep-tsv <file>` (run VEP yourself for unlimited scale).

Usage:
  scripts/vep_concordance.py <vcf> <gff3(.gz)> <fasta> [N] [--vep-tsv f]
Dumps: data/vep_dumps/<YYYY-MM-DD>/annotations.parquet  (+ concordance_log.csv)
"""
import json, subprocess, sys, time, urllib.request, urllib.error, random, datetime, os

args = [a for a in sys.argv[1:] if not a.startswith("--")]
VCF = args[0]
GFF3, FASTA = os.path.abspath(args[1]), os.path.abspath(args[2])  # VEP/duckvep want absolute
N = int(args[3]) if len(args) > 3 else 100
ROOT = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True).stdout.strip()
DUCKDB, EXT = f"{ROOT}/.tools/duckdb", f"{ROOT}/build/release/duckvep.duckdb_extension"
REST = "https://rest.ensembl.org/vep/human/region"
DATE = datetime.date.today().isoformat()
OUTDIR = f"{ROOT}/data/vep_dumps/{DATE}"
os.makedirs(OUTDIR, exist_ok=True)

# 1. Sample biallelic variants — SNVs AND indels/MNVs (the indel paths exercise
# frameshift / inframe / HGVS-3'-shift, the most divergence-prone logic).
# Transparently handles bgzip/gzip input.
import gzip, re
ACGT = re.compile(r"[ACGT]+\Z")
opener = gzip.open if VCF.endswith(".gz") else open
variants = []
with opener(VCF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        ref, alt = f[3], f[4]
        if "," in alt:  # biallelic only (clean join key)
            continue
        if ACGT.match(ref) and ACGT.match(alt) and not (ref == alt):
            variants.append((f[0], int(f[1]), ref, alt))
random.seed(42)
sample = random.sample(variants, min(N, len(variants)))
n_indel = sum(1 for (_, _, r, a) in sample if len(r) != len(a))
print(f"sampled {len(sample)} biallelic variants ({n_indel} indels)", file=sys.stderr)

# Shared sorted sample VCF (read by duckvep, offline VEP, and fastVEP). Declare
# every contig present so genome-wide samples parse, and sort grouped-by-contig
# then by position (VEP requires sorted input).
SAMPLE_VCF = "/tmp/sample.vcf"
contigs = sorted({c for (c, _, _, _) in sample})
with open(SAMPLE_VCF, "w") as fh:
    fh.write("##fileformat=VCFv4.2\n")
    for c in contigs:
        fh.write(f"##contig=<ID={c}>\n")
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for (c, p, r, a) in sorted(sample, key=lambda x: (x[0], x[1])):
        fh.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t.\n")

# 2. Ensembl VEP REST, PACED (batch 100, sleep between, backoff on 429).
def vep_rest(batch, tries=4):
    body = json.dumps({"variants": [f"{c} {p} . {r} {a} . . ." for (c, p, r, a) in batch]}).encode()
    req = urllib.request.Request(REST, data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    for t in range(tries):
        try:
            return json.load(urllib.request.urlopen(req, timeout=120))
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait = int(e.headers.get("Retry-After", 2 ** t)); time.sleep(wait); continue
            raise
    raise RuntimeError("rate-limited")

# 2b. OFFLINE Ensembl VEP (preferred): no REST rate limit -> high throughput.
VEP_CACHE = os.environ.get("VEP_CACHE", f"{ROOT}/data/vep_cache")
VEP_REL = os.environ.get("VEP_CACHE_VERSION", "116")
VEP_RUN = os.environ.get("VEP_CMD", "conda run -n vep vep").split()

def vep_offline(sample):
    # JSON output (not tab) so every record carries `input` = the original VCF
    # line: we join on the ORIGINAL VCF (pos,ref,alt), which matches duckvep
    # exactly. (Tab output reports a trimmed/anchor-stripped Allele that does NOT
    # match the VCF alt for indels — a silent join miss.)
    subprocess.run(VEP_RUN + ["-i", SAMPLE_VCF, "--offline", "--cache",
        "--dir_cache", VEP_CACHE, "--cache_version", VEP_REL, "--species", "homo_sapiens",
        "--assembly", "GRCh38", "--fasta", FASTA, "--symbol", "--json", "-o", "/tmp/vep_off.json",
        "--force_overwrite", "--no_stats"], check=True, capture_output=True, text=True)
    rows = []
    for line in open("/tmp/vep_off.json"):
        rec = json.loads(line)
        inp = rec["input"].split("\t")
        pos, ref, alt = int(inp[1]), inp[3], inp[4]
        for tc in rec.get("transcript_consequences", []):
            rows.append(dict(source="vep", date=DATE, pos=pos, ref=ref, alt=alt,
                transcript_id=tc["transcript_id"], gene_symbol=tc.get("gene_symbol", ""),
                consequence="&".join(sorted(tc["consequence_terms"])), impact=tc.get("impact", "")))
    return rows

USE_OFFLINE = os.path.isdir(f"{VEP_CACHE}/homo_sapiens/{VEP_REL}_GRCh38")
if USE_OFFLINE:
    vep_rows = vep_offline(sample)
    print(f"VEP offline ({VEP_REL} GRCh38): {len(vep_rows)} (variant,transcript) rows", file=sys.stderr)
else:
    vep_rows = []
    for i in range(0, len(sample), 100):
        try:
            recs = vep_rest(sample[i:i+100])
        except Exception as e:
            print(f"  batch {i}: skipped ({e})", file=sys.stderr); continue
        for rec in recs:
            inp = rec["input"].split()
            for tc in rec.get("transcript_consequences", []):
                vep_rows.append(dict(source="vep", date=DATE, pos=int(inp[1]), ref=inp[3], alt=inp[4],
                    transcript_id=tc["transcript_id"], gene_symbol=tc.get("gene_symbol", ""),
                    consequence="&".join(sorted(tc["consequence_terms"])), impact=tc.get("impact", "")))
        time.sleep(0.5)  # pace
    print(f"VEP REST: {len(vep_rows)} (variant,transcript) rows", file=sys.stderr)
with open("/tmp/vep_raw.json", "w") as fh:
    for r in vep_rows:
        fh.write(json.dumps(r) + "\n")

# 3. duckvep (one session).
# duckvep reads the sample VCF (chunked streaming path), not a giant VALUES.
sql = f"""LOAD '{EXT}';
SELECT vep_load_cache('{GFF3}', '{FASTA}');
COPY (
  SELECT 'duckvep' AS source, '{DATE}' AS date, v.pos, v.ref, a.alt, c.transcript_id,
         c.gene_symbol, list_aggregate(list_sort(c.consequence),'string_agg','&') AS consequence, c.impact
  FROM read_vcf('{SAMPLE_VCF}') v, UNNEST(v.alt) AS a(alt),
       UNNEST(vep_consequence(v.chrom, v.pos, v.ref, a.alt)) AS u(c)
) TO '/tmp/dv.json' (FORMAT json);"""
subprocess.run([DUCKDB, "-unsigned", "-c", sql], check=True)

# 3b. fastVEP (the underlying engine) on the same sample — validates engine vs VEP.
FASTVEP = os.environ.get("FASTVEP", f"{ROOT}/../DuckfastVEP/target/release/fastvep")
fv_rows = []
if os.path.exists(FASTVEP):
    fv = subprocess.run([FASTVEP, "annotate", "-i", SAMPLE_VCF, "--gff3", GFF3,
                         "--fasta", FASTA, "--output-format", "tab"], capture_output=True, text=True)
    cols = None
    for line in fv.stdout.splitlines():
        if line.startswith("##"):
            continue
        f = line.split("\t")
        if line.startswith("#"):
            cols = {n: i for i, n in enumerate(f)}; continue
        if not cols:
            continue
        loc, allele, tid, csq = f[cols["Location"]], f[cols["Allele"]], f[cols["Feature"]], f[cols["Consequence"]]
        pos = int(loc.split(":")[1].split("-")[0])
        fv_rows.append(dict(source="fastvep", date=DATE, pos=pos, ref="", alt=allele,
            transcript_id=tid, gene_symbol="", consequence="&".join(sorted(csq.split(","))), impact=""))
    with open("/tmp/fv.json", "w") as fh:
        for r in fv_rows:
            fh.write(json.dumps(r) + "\n")
    print(f"fastVEP: {len(fv_rows)} rows", file=sys.stderr)

# 4. Dated Parquet dump (all sources) + THREE-WAY concordance vs Ensembl VEP.
fv_union = "UNION ALL BY NAME SELECT * FROM read_json('/tmp/fv.json')" if fv_rows else ""
# Join on (pos, alt, transcript_id): SNVs, and fastVEP's tab output omits ref.
summary_sql = f"""
CREATE TABLE ann AS
  SELECT * FROM read_json('/tmp/vep_raw.json', columns={{source:'VARCHAR',date:'VARCHAR',pos:'BIGINT',ref:'VARCHAR',alt:'VARCHAR',transcript_id:'VARCHAR',gene_symbol:'VARCHAR',consequence:'VARCHAR',impact:'VARCHAR'}})
  UNION ALL BY NAME SELECT * FROM read_json('/tmp/dv.json')
  {fv_union};
COPY (SELECT * FROM ann ORDER BY pos, transcript_id, source) TO '{OUTDIR}/annotations.parquet' (FORMAT parquet);
WITH v AS (
  SELECT pos,alt,transcript_id,consequence,
         CASE WHEN length(ref)=1 AND length(alt)=1 THEN 'snv'
              WHEN length(ref)>length(alt) THEN 'del'
              WHEN length(ref)<length(alt) THEN 'ins'
              ELSE 'mnv' END AS class
  FROM ann WHERE source='vep')
SELECT e.source AS engine, v.class, count(*) AS pairs,
       count(*) FILTER (WHERE v.consequence=e.consequence) AS agree,
       round(100.0*count(*) FILTER (WHERE v.consequence=e.consequence)/nullif(count(*),0),2) AS pct
FROM (SELECT * FROM ann WHERE source<>'vep') e
JOIN v USING (pos,alt,transcript_id)
GROUP BY e.source, v.class ORDER BY e.source, v.class;
"""
out = subprocess.run([DUCKDB, "-csv", "-c", summary_sql], capture_output=True, text=True).stdout.strip().splitlines()
log = f"{ROOT}/data/vep_dumps/concordance_log.csv"
new = not os.path.exists(log)
print(f"\n=== Concordance vs Ensembl VEP ({DATE}) — by variant class ===")
with open(log, "a") as fh:
    if new:
        fh.write("date,engine,class,n_variants,pairs,agree,pct\n")
    for row in out[1:]:  # skip header
        engine, vclass, pairs, agree, pct = row.split(",")
        fh.write(f"{DATE},{engine},{vclass},{len(sample)},{pairs},{agree},{pct}\n")
        print(f"  {engine:8s} {vclass:4s} vs VEP: pairs={pairs} agree={agree} concordance={pct}%")
print(f"  dump: {OUTDIR}/annotations.parquet")
