//! The annotation engine core — **DuckDB-free**.
//!
//! This module wraps fastVEP's consequence engine and produces plain owned
//! rows. It has no `duckdb` dependency, which is the point: the DuckDB
//! extension (`super::annotate`), a future CLI binary, and a future C API are
//! all thin frontends over this one function. Only the extension links the
//! DuckDB C API; the CLI/C-API do not link libduckdb at all.

use crate::vep::tcache;
use fastvep_cache::fasta::{FastaReader, MmapFastaReader};
use fastvep_cache::gff::parse_gff3;
use fastvep_cache::providers::{
    FastaSequenceProvider, IndexedTranscriptProvider, MmapFastaSequenceProvider, SequenceProvider,
    TranscriptProvider,
};
use fastvep_consequence::sv_predictor::predict_sv_consequences;
use fastvep_consequence::AlleleConsequenceResult;
use fastvep_consequence::{ConsequencePredictor, TranscriptConsequence};
use fastvep_core::{Allele, Consequence, GenomicPosition, Strand, VariantType};
use fastvep_genome::Transcript;
use noodles::vcf;
use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::Arc;

/// Detected gene-model input format (by content, not file extension).
enum ModelFormat {
    /// Our columnar Parquet transcript cache (magic `PAR1`).
    ParquetCache,
    /// gzip/bgzf-compressed GFF3 (magic `1f 8b`).
    GzippedGff3,
    /// Plain-text GFF3.
    Gff3,
}

fn detect_model_format(path: &str) -> std::io::Result<ModelFormat> {
    let mut head = [0u8; 4];
    let n = File::open(path)?.read(&mut head)?;
    Ok(if n >= 4 && &head == b"PAR1" {
        ModelFormat::ParquetCache
    } else if n >= 2 && head[0] == 0x1f && head[1] == 0x8b {
        ModelFormat::GzippedGff3
    } else {
        ModelFormat::Gff3
    })
}

pub const DEFAULT_DISTANCE: u64 = 5000;

type SeqProvider = Box<dyn SequenceProvider + Send + Sync>;

/// The lean consequence context: transcript index + reference + predictor.
/// Built directly from the cache/consequence crates, deliberately bypassing
/// fastVEP's `AnnotationContext` (which also drags in SA/ACMG/HGVS/IO).
///
/// `Send + Sync` so it can live in the extension's load-once cache and be shared
/// across DuckDB worker threads (`super::consequence`).
pub(crate) struct EngineContext {
    transcripts: IndexedTranscriptProvider,
    predictor: ConsequencePredictor,
    distance: u64,
}

pub(crate) fn build_context(
    gff3: &str,
    fasta: Option<&str>,
    distance: u64,
    zstd_level: i32,
) -> Result<EngineContext, Box<dyn Error>> {
    // The `model` argument is detected by content (not extension):
    //  - a Parquet transcript cache -> load it directly (the fast path), or
    //  - a GFF3 (plain or gzip) -> parse it, and write a columnar Parquet cache
    //    next to it (`<gff3>.transcripts.parquet`) so later loads take the fast
    //    path. Parsing ~280k transcripts dominates the cold build; the warm cache
    //    load (bincode `model` column) is ~1.5 s. `build_sequences` is cheap and
    //    FASTA-specific, so we rebuild it fresh below rather than caching it.
    //    (docs/DESIGN.md §5; cache format in vep::tcache.)
    let mut transcripts = match detect_model_format(gff3)? {
        ModelFormat::ParquetCache => tcache::load(Path::new(gff3))?,
        format => {
            let cache_path = tcache::cache_path(gff3);
            if tcache::is_fresh(&cache_path, Path::new(gff3)) {
                tcache::load(&cache_path)?
            } else {
                let gff_file = File::open(gff3)?;
                let t = match format {
                    ModelFormat::GzippedGff3 => {
                        parse_gff3(flate2::read::MultiGzDecoder::new(gff_file))?
                    }
                    _ => parse_gff3(gff_file)?,
                };
                let _ = tcache::save(&t, &cache_path, zstd_level);
                t
            }
        }
    };

    let seq: Option<SeqProvider> = match fasta {
        Some(path) if Path::new(&format!("{path}.fai")).exists() => Some(Box::new(
            MmapFastaSequenceProvider::new(MmapFastaReader::open(Path::new(path))?),
        )),
        Some(path) => Some(Box::new(FastaSequenceProvider::new(
            FastaReader::from_reader(File::open(path)?)?,
        ))),
        None => None,
    };

    // Build spliced cDNA / protein sequences so coding consequences (missense,
    // synonymous, amino acids) are exact. Requires a reference. The provider is used
    // ONLY here (sequences are baked into each transcript's `translateable_seq`); it is
    // not retained on the context, since nothing reads the reference per variant.
    if let Some(sp) = &seq {
        for t in &mut transcripts {
            let _ = t.build_sequences(|chrom, start, end| {
                sp.fetch_sequence(chrom, start, end)
                    .map_err(|e| e.to_string())
            });
        }
    }
    drop(seq);

    Ok(EngineContext {
        transcripts: IndexedTranscriptProvider::new(transcripts),
        predictor: ConsequencePredictor::new(distance, distance),
        distance,
    })
}

impl EngineContext {
    /// END-aware kernel. `end` is the variant interval end (`INFO/END` for SVs);
    /// for SNVs/indels it is `pos + len(ref) - 1`. Structural alt alleles
    /// (symbolic `<DEL>`/`<DUP>`/`<CNV>`/`<INV>`, breakends) are dispatched to the
    /// SV consequence predictor over the spanned interval; sequence alleles go
    /// through the regular predictor. A mixed multiallelic site does both.
    pub(crate) fn annotate_variant_spanned(
        &self,
        chrom: &str,
        pos: u64,
        end: u64,
        ref_str: &str,
        alt_raw: &str,
    ) -> Vec<AnnotatedRow> {
        let mut rows = Vec::new();
        if alt_raw.is_empty() || alt_raw == "." {
            return rows;
        }
        let end = end.max(pos + (ref_str.len() as u64).saturating_sub(1));
        let ref_allele = Allele::from_str(ref_str);
        let alt_alleles: Vec<Allele> = alt_raw.split(',').map(Allele::from_str).collect();
        // Partition alts: structural (symbolic/breakend) vs ordinary sequence.
        let (sv_alleles, seq_alleles): (Vec<Allele>, Vec<Allele>) = alt_alleles
            .iter()
            .cloned()
            .partition(|a| classify_sv_type(a).is_structural());

        let position = GenomicPosition::new(chrom.to_string(), pos, end, Strand::Forward);
        let query_start = pos.saturating_sub(self.distance).max(1);
        let query_end = end + self.distance;
        // JUSTIFIED default: a `get_transcripts` Err (or a contig absent from the index)
        // means there are no candidate transcripts for THIS variant — annotation is
        // best-effort per variant, so we yield no rows for it rather than aborting the whole
        // scan. The empty result is then handled by the `is_empty()` early return below.
        let overlapping = self
            .transcripts
            .get_transcripts(chrom, query_start, query_end)
            .unwrap_or_default();
        if overlapping.is_empty() {
            return rows;
        }
        // NOTE: we deliberately do NOT fetch the reference region here. Coding
        // consequences read each transcript's cached `translateable_seq` (built once at
        // cache load), and `predict_allele` ignores the `ref_seq` argument. Fetching
        // `pos ± distance` (~10 KB) per variant was a pure allocation hotspot that
        // capped thread scaling (allocator contention across DuckDB workers) for zero
        // benefit. If a future reference-validating path needs it, fetch conditionally.
        let mut transcript_consequences: Vec<TranscriptConsequence> = Vec::new();
        if !seq_alleles.is_empty() {
            transcript_consequences.extend(
                self.predictor
                    .predict(&position, &ref_allele, &seq_alleles, &overlapping, None)
                    .transcript_consequences,
            );
        }
        // Each SV allele may carry its own type (e.g. <DEL> vs <DUP>), so call the
        // SV predictor per allele with that allele's classified type.
        for sv in &sv_alleles {
            transcript_consequences.extend(predict_sv_consequences(
                chrom,
                pos,
                end,
                classify_sv_type(sv),
                std::slice::from_ref(sv),
                &overlapping,
                self.distance,
                self.distance,
            ));
        }

        // Index the overlapping transcripts by stable id so HGVS can reach the
        // full transcript model (coding bounds, spliced sequence, protein id).
        let tmap: std::collections::HashMap<&str, &Transcript> =
            overlapping.iter().map(|t| (&*t.stable_id, *t)).collect();

        let chrom_arc: Arc<str> = Arc::from(chrom);
        let empty: Arc<str> = Arc::from("");
        for tc in &transcript_consequences {
            let tr = tmap.get(&*tc.transcript_id).copied();
            for ac in &tc.allele_consequences {
                // Ensembl omits a (variant, transcript) pair that yields no consequence
                // (e.g. a candidate transcript beyond the up/downstream distance) — emit no
                // row rather than an empty one.
                if ac.consequences.is_empty() {
                    continue;
                }
                let (hgvsg, hgvsc, hgvsp) = match tr {
                    // SV/breakend HGVS (del/dup/ins ranges) is its own spec and not
                    // yet validated, so we leave it empty rather than emit guesses.
                    Some(t) if !ac.allele.is_symbolic() => {
                        build_hgvs(t, chrom, pos, end, &ref_allele, ac)
                    }
                    _ => (String::new(), String::new(), String::new()),
                };
                rows.push(AnnotatedRow {
                    chrom: chrom_arc.clone(),
                    pos: pos as i64,
                    reference: ref_str.to_string(),
                    alt: ac.allele.to_string(),
                    gene_id: tc.gene_id.clone(),
                    gene_symbol: tc.gene_symbol.clone().unwrap_or_else(|| empty.clone()),
                    transcript_id: tc.transcript_id.clone(),
                    biotype: tc.biotype.clone(),
                    canonical: tc.canonical,
                    consequence: ac
                        .consequences
                        .iter()
                        .map(|c| c.so_term().to_string())
                        .collect(),
                    impact: ac.impact.as_str().to_string(),
                    amino_acids: pair_to_string(&ac.amino_acids),
                    codons: pair_to_string(&ac.codons),
                    protein_pos: ac.protein_start.map(|p| p as i64),
                    hgvsg,
                    hgvsc,
                    hgvsp,
                });
            }
        }
        rows
    }

    /// Combine a set of PHASED variants on ONE transcript into a single haplotype coding
    /// consequence (the bcftools-csq / haplosaurus model — co-located variants merged in
    /// transcript/CDS coordinates and translated once). `variants` are `(pos, ref, alt)`;
    /// returns the sorted SO terms, or empty if the transcript isn't found in range or no
    /// variant is coding. The grouping by (sample, haplotype, transcript) is left to SQL.
    pub(crate) fn haplotype_consequence(
        &self,
        chrom: &str,
        transcript_id: &str,
        variants: &[(u64, String, String)],
    ) -> Vec<String> {
        if variants.is_empty() {
            return Vec::new();
        }
        let lo = variants.iter().map(|(p, _, _)| *p).min().unwrap_or(1);
        let hi = variants
            .iter()
            .map(|(p, r, _)| p + (r.len() as u64).saturating_sub(1))
            .max()
            .unwrap_or(lo);
        let q_start = lo.saturating_sub(self.distance).max(1);
        let q_end = hi + self.distance;
        // JUSTIFIED default (same as annotate_variant_spanned): a lookup Err / absent contig
        // means no candidate transcripts, so the haplotype has no coding context here — we
        // return no terms rather than aborting. The `find` below then yields None.
        let overlapping = self
            .transcripts
            .get_transcripts(chrom, q_start, q_end)
            .unwrap_or_default();
        let Some(tr) = overlapping.iter().find(|t| &*t.stable_id == transcript_id) else {
            return Vec::new();
        };
        let edits: Vec<(u64, u64, Allele, Allele)> = variants
            .iter()
            .map(|(p, r, a)| {
                let end = p + (r.len() as u64).saturating_sub(1);
                (*p, end, Allele::from_str(r), Allele::from_str(a))
            })
            .collect();
        match self.predictor.haplotype_coding_terms(tr, &edits) {
            Some(mut terms) => {
                terms.sort_by_key(|c| c.rank());
                terms.dedup();
                terms.iter().map(|c| c.so_term().to_string()).collect()
            }
            None => Vec::new(),
        }
    }
}

/// One annotated (variant, transcript, allele) row. Plain owned data so it can
/// cross any frontend boundary (DuckDB vectors, CLI output, C API).
pub(crate) struct AnnotatedRow {
    // Categorical/repeated fields are `Arc<str>` shared from the `Transcript`
    // model — cloning is a refcount bump, not a re-allocation, across the
    // millions of (variant, transcript) rows.
    pub chrom: Arc<str>,
    pub pos: i64,
    pub reference: String,
    pub alt: String,
    pub gene_id: Arc<str>,
    pub gene_symbol: Arc<str>,
    pub transcript_id: Arc<str>,
    pub biotype: Arc<str>,
    pub canonical: bool,
    pub consequence: Vec<String>,
    pub impact: String,
    pub amino_acids: String,
    pub codons: String,
    pub protein_pos: Option<i64>,
    // HGVS nomenclature (empty string where a notation does not apply).
    pub hgvsg: String,
    pub hgvsc: String,
    pub hgvsp: String,
}

/// Annotate every variant in `vcf_path` against the `gff3` gene model (optional
/// `fasta` for sequence-dependent consequences). The single library entry point.
pub(crate) fn annotate(
    vcf_path: &str,
    gff3: &str,
    fasta: Option<&str>,
    distance: u64,
) -> Result<Vec<AnnotatedRow>, Box<dyn Error>> {
    // The library entry point builds the cache (if stale) at the default zstd
    // level; the SQL surface exposes the level via `vep_load_cache(..., zstd)`.
    let ctx = build_context(gff3, fasta, distance, tcache::DEFAULT_ZSTD)?;
    annotate_all(&ctx, vcf_path, distance)
}

fn pair_to_string(p: &Option<(String, String)>) -> String {
    match p {
        Some((a, b)) => format!("{a}/{b}"),
        None => String::new(),
    }
}

/// Classify an alt allele's structural type from its symbolic tag or breakend
/// syntax, mirroring fastVEP's `classify_sv_type`. Non-structural alleles (plain
/// sequence, `<INS>`) return a non-structural type so they route to the regular
/// predictor. `<DEL>` is a copy-number loss (a structural type), not a small del.
fn classify_sv_type(alt: &Allele) -> VariantType {
    let s = match alt {
        Allele::Symbolic(s) => s.trim_matches(|c| c == '<' || c == '>').to_uppercase(),
        // Breakend ALTs carry mate brackets, e.g. `G]chr2:200]` or `[chr2:200[G`.
        Allele::Sequence(b) if b.contains(&b'[') || b.contains(&b']') => {
            return VariantType::TranslocationBreakend
        }
        _ => return VariantType::Unknown,
    };
    match s.as_str() {
        // A symbolic `<DEL>` is a sequence deletion (feature_truncation / transcript_ablation),
        // NOT a copy-number-loss allele — only `<CN0>`/`<CNn<2>` get the copy_number_* terms.
        "DEL" => VariantType::Deletion,
        "DUP" | "DUP:TANDEM" => VariantType::TandemDuplication,
        "INV" => VariantType::Inversion,
        "BND" => VariantType::TranslocationBreakend,
        "INS" => VariantType::Insertion,
        "CNV" => VariantType::CopyNumberVariation,
        "STR" => VariantType::ShortTandemRepeatVariation,
        s if s.starts_with("CN") => match s.trim_start_matches("CN").parse::<u32>() {
            Ok(cn) if cn < 2 => VariantType::CopyNumberLoss,
            Ok(cn) if cn > 2 => VariantType::CopyNumberGain,
            _ => VariantType::CopyNumberVariation,
        },
        _ => VariantType::Unknown,
    }
}

fn complement_allele(allele: &Allele) -> Allele {
    match allele {
        Allele::Sequence(bases) => Allele::Sequence(
            bases
                .iter()
                .map(|&b| match b {
                    b'A' | b'a' => b'T',
                    b'T' | b't' => b'A',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    other => other,
                })
                .collect(),
        ),
        other => other.clone(),
    }
}

/// Build `(HGVSg, HGVSc, HGVSp)` for one (transcript, allele) consequence,
/// mirroring fastVEP's annotate layer (`fastvep-annotate`). Empty strings where a
/// notation does not apply (e.g. HGVSp for a non-coding consequence). This is the
/// last feature-parity gap with `fastvep annotate`'s tab/JSON output.
fn build_hgvs(
    tr: &Transcript,
    chrom: &str,
    pos_start: u64,
    pos_end: u64,
    ref_allele: &Allele,
    ac: &AlleleConsequenceResult,
) -> (String, String, String) {
    let hgvsg = fastvep_hgvs::hgvsg(chrom, pos_start, pos_end, ref_allele, &ac.allele);

    let versioned_tid = match tr.version {
        Some(v) => format!("{}.{}", tr.stable_id, v),
        None => tr.stable_id.to_string(),
    };
    // HGVS bases are written in the transcript's orientation.
    let (hgvs_ref, hgvs_alt) = if tr.strand == Strand::Reverse {
        (complement_allele(ref_allele), complement_allele(&ac.allele))
    } else {
        (ref_allele.clone(), ac.allele.clone())
    };
    let intronic = |start: u64, end: u64| {
        tr.genomic_to_intronic_cdna(start)
            .map(|(cdna_pos, offset)| {
                let (end_cdna, end_offset) = if start != end {
                    tr.genomic_to_intronic_cdna(end)
                        .map(|(c, o)| (Some(c), Some(o)))
                        .unwrap_or((None, None))
                } else {
                    (None, None)
                };
                (cdna_pos, offset, end_cdna, end_offset)
            })
    };

    let mut hgvsc = None;
    if let Some(coding_start) = tr.cdna_coding_start {
        if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
            let (cs, ce) = (cs.min(ce), cs.max(ce));
            hgvsc = fastvep_hgvs::hgvsc_with_seq(
                &versioned_tid,
                cs,
                ce,
                &hgvs_ref,
                &hgvs_alt,
                coding_start,
                tr.cdna_coding_end,
                tr.spliced_seq.as_deref(),
                tr.codon_table_start_phase,
            );
        } else if ac.intron.is_some() {
            if let Some((cdna_pos, offset, end_cdna, end_offset)) = intronic(pos_start, pos_end) {
                hgvsc = fastvep_hgvs::hgvsc_intronic_range(
                    &versioned_tid,
                    cdna_pos,
                    offset,
                    end_cdna,
                    end_offset,
                    &hgvs_ref,
                    &hgvs_alt,
                    coding_start,
                    tr.cdna_coding_end,
                );
            }
        }
    } else if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
        hgvsc = fastvep_hgvs::hgvsc_noncoding(&versioned_tid, cs, ce, &hgvs_ref, &hgvs_alt);
    } else if ac.intron.is_some() {
        if let Some((cdna_pos, offset, end_cdna, end_offset)) = intronic(pos_start, pos_end) {
            hgvsc = fastvep_hgvs::hgvsc_noncoding_intronic_range(
                &versioned_tid,
                cdna_pos,
                offset,
                end_cdna,
                end_offset,
                &hgvs_ref,
                &hgvs_alt,
            );
        }
    }

    let mut hgvsp = None;
    if let (Some(aa), Some(ps)) = (&ac.amino_acids, ac.protein_start) {
        if let Some(pid) = &tr.protein_id {
            let versioned_pid = match tr.protein_version {
                Some(v) => {
                    let suffix = format!(".{v}");
                    if pid.ends_with(suffix.as_str()) {
                        pid.clone()
                    } else {
                        format!("{pid}.{v}")
                    }
                }
                None => pid.clone(),
            };
            if ac.consequences.contains(&Consequence::FrameshiftVariant) {
                if let (Some(spliced), Some(coding_start), Some(cds_s)) =
                    (&tr.spliced_seq, tr.cdna_coding_start, ac.cds_start)
                {
                    let ref_from_cds = &spliced.as_bytes()[(coding_start - 1) as usize..];
                    let cds_idx = (cds_s - 1) as usize;
                    let mut alt_from_cds = ref_from_cds.to_vec();
                    if ac.allele == Allele::Deletion {
                        let end = (cds_idx + ref_allele.len()).min(alt_from_cds.len());
                        alt_from_cds.drain(cds_idx..end);
                    } else if let Allele::Sequence(ins_bases) = &ac.allele {
                        let bases: Vec<u8> = if tr.strand == Strand::Reverse {
                            match &complement_allele(&ac.allele) {
                                Allele::Sequence(b) => b.clone(),
                                _ => ins_bases.clone(),
                            }
                        } else {
                            ins_bases.clone()
                        };
                        for (j, &b) in bases.iter().enumerate() {
                            if cds_idx + j <= alt_from_cds.len() {
                                alt_from_cds.insert(cds_idx + j, b);
                            }
                        }
                    }
                    hgvsp = fastvep_hgvs::hgvsp_frameshift(
                        &versioned_pid,
                        ref_from_cds,
                        &alt_from_cds,
                        cds_idx / 3,
                    );
                }
            } else {
                let ref_aa = aa.0.as_bytes().first().copied().unwrap_or(b'X');
                let alt_aa = aa.1.as_bytes().first().copied().unwrap_or(b'X');
                hgvsp = fastvep_hgvs::hgvsp(&versioned_pid, ps, ref_aa, alt_aa, false);
            }
        }
    }

    (hgvsg, hgvsc.unwrap_or_default(), hgvsp.unwrap_or_default())
}

fn annotate_all(
    ctx: &EngineContext,
    vcf_path: &str,
    _distance: u64,
) -> Result<Vec<AnnotatedRow>, Box<dyn Error>> {
    let mut reader = vcf::io::reader::Builder::default().build_from_path(vcf_path)?;
    let _header = reader.read_header()?;
    let mut record = vcf::Record::default();
    let mut rows = Vec::new();

    while reader.read_record(&mut record)? != 0 {
        let chrom = record.reference_sequence_name().to_string();
        // A record with no POS is malformed; skip it rather than annotating at the invalid
        // coordinate 0 (`unwrap_or(0)` would have masked it). The `?` propagates a parse Err.
        let pos = match record.variant_start().transpose()? {
            Some(p) => usize::from(p) as u64,
            None => continue,
        };
        let ref_str = record.reference_bases().to_string();
        let alt_raw = record.alternate_bases().as_ref().to_string();
        // INFO/END spans symbolic SVs (<DEL>/<CNV>/…) past their 1 bp ref anchor.
        let info = record.info();
        let end = info_end(info.as_ref()).unwrap_or(pos + (ref_str.len() as u64).saturating_sub(1));
        rows.extend(ctx.annotate_variant_spanned(&chrom, pos, end, &ref_str, &alt_raw));
    }
    Ok(rows)
}

/// `INFO/END` (the SV/CNV interval end), if present and parseable.
fn info_end(info_raw: &str) -> Option<u64> {
    info_raw
        .split(';')
        .find_map(|f| f.strip_prefix("END="))
        .and_then(|v| v.parse::<u64>().ok())
}
