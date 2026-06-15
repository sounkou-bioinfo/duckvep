use fastvep_core::{Allele, Consequence, GenomicPosition, Impact, Strand};
use fastvep_genome::codon::format_codon_change;
use fastvep_genome::{CodonTable, Transcript};
use std::sync::Arc;

use crate::splice;

/// Result of consequence prediction for a variant against a transcript.
#[derive(Debug, Clone)]
pub struct TranscriptConsequence {
    pub transcript_id: Arc<str>,
    pub gene_id: Arc<str>,
    pub gene_symbol: Option<Arc<str>>,
    pub biotype: Arc<str>,
    pub allele_consequences: Vec<AlleleConsequenceResult>,
    pub canonical: bool,
    pub strand: Strand,
}

/// Consequence result for a single allele against a transcript.
#[derive(Debug, Clone)]
pub struct AlleleConsequenceResult {
    pub allele: Allele,
    pub consequences: Vec<Consequence>,
    pub impact: Impact,
    pub cdna_start: Option<u64>,
    pub cdna_end: Option<u64>,
    pub cds_start: Option<u64>,
    pub cds_end: Option<u64>,
    pub protein_start: Option<u64>,
    pub protein_end: Option<u64>,
    pub amino_acids: Option<(String, String)>,
    pub codons: Option<(String, String)>,
    pub exon: Option<(u32, u32)>,
    pub intron: Option<(u32, u32)>,
    pub distance: Option<i64>,
}

/// Full prediction result for a variant.
#[derive(Debug, Clone)]
pub struct PredictionResult {
    pub transcript_consequences: Vec<TranscriptConsequence>,
    pub most_severe: Option<Consequence>,
}

/// The genomic interval actually changed by a variant — the VCF anchor-inclusive
/// `[start, end]` minus the common prefix/suffix of ref/alt (matching Ensembl
/// VEP's normalized `r_start/r_end`). SNVs/MNVs are unchanged; a deletion trims its
/// leading anchor base; an insertion collapses to a point (start may exceed end).
fn normalized_interval(start: u64, end: u64, r: &Allele, a: &Allele) -> (u64, u64) {
    let (rb, ab) = match (r, a) {
        (Allele::Sequence(rb), Allele::Sequence(ab)) => (rb.as_slice(), ab.as_slice()),
        _ => return (start, end),
    };
    let mut p = 0;
    while p < rb.len() && p < ab.len() && rb[p] == ab[p] {
        p += 1;
    }
    let mut s = 0;
    while s < rb.len().saturating_sub(p)
        && s < ab.len().saturating_sub(p)
        && rb[rb.len() - 1 - s] == ab[ab.len() - 1 - s]
    {
        s += 1;
    }
    (start + p as u64, end.saturating_sub(s as u64))
}

/// The feature-overlap flags a variant has against a transcript, computed ONCE
/// (Ensembl `_overlapped_introns`/`_boundary`/`_exons`). Consequence terms are gated
/// on these via their declarative `include` conditions (see `include_satisfied`),
/// rather than scattered `if` guards.
#[derive(Debug, Clone, Copy)]
struct FeatureOverlap {
    exon: bool,
    intron: bool,
    intron_boundary: bool,
}

/// Ensembl OverlapConsequence `include` conditions (`Constants.pm`): a term is kept
/// only if the variant's feature-overlap flags match. This is the declarative gate that
/// replaces hand-wired guards — e.g. splice_polypyrimidine requires `exon=0, intron=1`
/// (so a deletion spanning an exon into the tract is correctly NOT called PPT), the
/// splice donor/acceptor/region terms require `intron_boundary=1`, and the UTR /
/// non-coding-exon terms require `exon=1`. Terms with no `include` (the protein-level
/// consequences) always pass.
fn include_satisfied(c: Consequence, f: &FeatureOverlap) -> bool {
    use Consequence::*;
    match c {
        SplicePolypyrimidineTractVariant => !f.exon && f.intron,
        IntronVariant => f.intron,
        SpliceDonorVariant
        | SpliceAcceptorVariant
        | SpliceDonorFifthBaseVariant
        | SpliceDonorRegionVariant
        | SpliceRegionVariant => f.intron_boundary,
        FivePrimeUtrVariant | ThreePrimeUtrVariant | NonCodingTranscriptExonVariant => f.exon,
        _ => true,
    }
}

/// One edit applied to the reference CDS in transcript (5'->3') orientation: replace
/// `ref_bases` starting at `cds_idx` (0-based) with `alt_bases`. A single variant
/// produces one edit; a phased haplotype (Ensembl haplosaurus / `bcftools csq`)
/// produces several edits applied to the SAME CDS before translation — which is how
/// co-located variants on one haplotype combine into one correct protein
/// consequence (an in-codon multi-SNP is just a degenerate "local haplotype").
#[derive(Debug, Clone)]
pub(crate) struct CdsEdit {
    pub cds_idx: usize,
    pub ref_bases: Vec<u8>,
    pub alt_bases: Vec<u8>,
}

/// Peptide/codon view of one or more `CdsEdit`s, built once and queried by the flat
/// consequence-predicate set. Equivalent to Ensembl's TranscriptVariationAllele
/// peptide layer, but keyed on a set of edits so it serves single variants and whole
/// phased haplotypes identically.
struct CodingContext {
    /// 0-based codon number of the first affected codon.
    first_codon: usize,
    /// the affected codon span runs past the end of the translateable CDS.
    incomplete_terminal: bool,
    /// reference / alternate codon sequence spanning the change (whole codons).
    ref_codons: Vec<u8>,
    alt_codons: Vec<u8>,
    /// translated reference / alternate peptides over the affected window.
    ref_pep: Vec<u8>,
    alt_pep: Vec<u8>,
    /// total reference / alternate bases changed (summed over edits).
    ref_len: usize,
    alt_len: usize,
    /// the edits overlap the CDS stop codon (the last 3 CDS bases).
    overlaps_stop: bool,
    /// for an indel overlapping the stop: Ensembl `_ins_del_stop_altered` — applying the
    /// edit to CDS+3'UTR leaves NO stop at the old stop position (or shortens the CDS),
    /// so the stop is altered (-> stop_lost). `false` with `overlaps_stop` -> stop_retained.
    stop_altered: bool,
}

/// Trim the shared prefix/suffix of two alleles to their minimal changed bytes
/// (genomic orientation). Matches VEP's minimal representation; e.g. GAATTT/G -> (AATTT, "").
fn minimal_alleles(r: &Allele, a: &Allele) -> (Vec<u8>, Vec<u8>) {
    let rb: &[u8] = match r {
        Allele::Sequence(b) => b,
        _ => &[],
    };
    let ab: &[u8] = match a {
        Allele::Sequence(b) => b,
        _ => &[],
    };
    let mut p = 0;
    while p < rb.len() && p < ab.len() && rb[p] == ab[p] {
        p += 1;
    }
    let mut s = 0;
    while s < rb.len().saturating_sub(p)
        && s < ab.len().saturating_sub(p)
        && rb[rb.len() - 1 - s] == ab[ab.len() - 1 - s]
    {
        s += 1;
    }
    (rb[p..rb.len() - s].to_vec(), ab[p..ab.len() - s].to_vec())
}

/// Trim the shared prefix and suffix of two byte slices (Ensembl `trim_sequences`),
/// returning the differing middles. Used to test whether an inframe deletion's alt
/// codon matches a contiguous part of the ref codon (clean deletion vs delins).
fn trim_common<'a>(r: &'a [u8], a: &'a [u8]) -> (&'a [u8], &'a [u8]) {
    let mut p = 0;
    while p < r.len() && p < a.len() && r[p] == a[p] {
        p += 1;
    }
    let mut s = 0;
    while s < r.len().saturating_sub(p) && s < a.len().saturating_sub(p) && r[r.len() - 1 - s] == a[a.len() - 1 - s] {
        s += 1;
    }
    (&r[p..r.len() - s], &a[p..a.len() - s])
}

/// Convert one variant into a CDS edit in transcript orientation. `cds_start`/`cds_end`
/// are the already-normalized (anchor-trimmed) CDS coordinates from `predict_allele`,
/// so the minimal allele and the CDS index stay consistent (the anchor-base frame bug
/// is impossible by construction). Reverse-strand bases are reverse-complemented.
fn variant_to_cds_edit(
    ref_allele: &Allele,
    alt_allele: &Allele,
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    strand: Strand,
) -> Option<CdsEdit> {
    let (mref_g, malt_g) = minimal_alleles(ref_allele, alt_allele);
    let cds_lo = match (cds_start, cds_end) {
        (Some(a), Some(b)) => a.min(b),
        (Some(a), None) => a,
        (None, Some(b)) => b,
        (None, None) => return None,
    };
    if cds_lo == 0 {
        return None;
    }
    let (ref_bases, alt_bases) = match strand {
        Strand::Forward => (mref_g, malt_g),
        Strand::Reverse => (
            fastvep_genome::codon::reverse_complement(&mref_g),
            fastvep_genome::codon::reverse_complement(&malt_g),
        ),
    };
    // A pure insertion (no ref bases) goes AFTER cds_lo; a substitution/deletion
    // starts AT cds_lo. cds_lo is 1-based.
    let cds_idx = if ref_bases.is_empty() {
        cds_lo as usize
    } else {
        (cds_lo - 1) as usize
    };
    Some(CdsEdit {
        cds_idx,
        ref_bases,
        alt_bases,
    })
}

/// The consequence prediction engine.
pub struct ConsequencePredictor {
    pub upstream_distance: u64,
    pub downstream_distance: u64,
    codon_table: CodonTable,
    mito_codon_table: CodonTable,
}

impl ConsequencePredictor {
    pub fn new(upstream_distance: u64, downstream_distance: u64) -> Self {
        Self {
            upstream_distance,
            downstream_distance,
            codon_table: CodonTable::standard(),
            // Vertebrate mitochondrial code (NCBI table 2): TGA=Trp, ATA=Met,
            // AGA/AGG=stop. chrM transcripts must use this, not the standard table.
            mito_codon_table: CodonTable::from_ncbi_table(2),
        }
    }

    /// The codon table to translate a given transcript with — mitochondrial
    /// (NCBI table 2) for chrM/MT transcripts, the standard table otherwise.
    fn ct(&self, transcript: &Transcript) -> &CodonTable {
        if matches!(&*transcript.chromosome, "MT" | "chrM" | "M" | "chrMT") {
            &self.mito_codon_table
        } else {
            &self.codon_table
        }
    }

    /// Predict consequences of a variant against a set of transcripts.
    pub fn predict(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_alleles: &[Allele],
        transcripts: &[&Transcript],
        ref_seq: Option<&[u8]>,
    ) -> PredictionResult {
        let mut transcript_consequences = Vec::new();

        for transcript in transcripts {
            let tc =
                self.predict_transcript(position, ref_allele, alt_alleles, transcript, ref_seq);
            transcript_consequences.push(tc);
        }

        let all_consequences: Vec<Consequence> = transcript_consequences
            .iter()
            .flat_map(|tc| {
                tc.allele_consequences
                    .iter()
                    .flat_map(|ac| ac.consequences.iter().copied())
            })
            .collect();

        let most_severe = Consequence::most_severe(&all_consequences);

        PredictionResult {
            transcript_consequences,
            most_severe,
        }
    }

    fn predict_transcript(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_alleles: &[Allele],
        transcript: &Transcript,
        ref_seq: Option<&[u8]>,
    ) -> TranscriptConsequence {
        let allele_consequences: Vec<AlleleConsequenceResult> = alt_alleles
            .iter()
            .map(|alt| self.predict_allele(position, ref_allele, alt, transcript, ref_seq))
            .collect();

        TranscriptConsequence {
            transcript_id: transcript.stable_id.clone(),
            gene_id: transcript.gene.stable_id.clone(),
            gene_symbol: transcript.gene.symbol.clone(),
            biotype: transcript.biotype.clone(),
            allele_consequences,
            canonical: transcript.canonical,
            strand: transcript.strand,
        }
    }

    fn predict_allele(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        _ref_seq: Option<&[u8]>,
    ) -> AlleleConsequenceResult {
        // VEP computes consequences on the NORMALIZED affected interval (the
        // anchor-trimmed changed bases), not the VCF anchor-inclusive interval.
        // For indels at a splice site, the unchanged anchor base wrongly overlaps
        // the donor/acceptor and mis-classifies — match VEP by trimming it.
        let (var_start, var_end) =
            normalized_interval(position.start, position.end, ref_allele, alt_allele);
        // Splice/intron predicates use the genomic interval Ensembl's
        // `_get_differing_regions` produces: a delins whose ALT is longer than its REF
        // spans the alt length (the inserted bases extend the affected region toward a
        // splice site), so `r_end = var_start + alt_len - 1`. SNVs/MNVs/deletions and
        // pure insertions are unchanged (alt not longer, or the anchor-trimmed form).
        let splice_end = {
            let (mref, malt) = minimal_alleles(ref_allele, alt_allele);
            if !mref.is_empty() && malt.len() > mref.len() {
                var_start + malt.len() as u64 - 1
            } else {
                var_end
            }
        };

        // Ensembl `_get_differing_regions`: a same-length MNV is split at its internal
        // MATCHING bases, so splice/intron predicates see only the actually-changed
        // sub-intervals (e.g. AACTC/GACCA differs at positions 0,3,4 -> regions
        // [start] and [start+3,start+4], NOT the whole span). SNVs/indels -> one region.
        let diff_regions: Vec<(u64, u64)> = match (ref_allele, alt_allele) {
            (Allele::Sequence(r), Allele::Sequence(a)) if r.len() == a.len() && r.len() > 1 => {
                let mut regions = Vec::new();
                let mut i = 0;
                while i < r.len() {
                    if r[i] != a[i] {
                        let s = i;
                        while i < r.len() && r[i] != a[i] {
                            i += 1;
                        }
                        regions.push((position.start + s as u64, position.start + (i - 1) as u64));
                    } else {
                        i += 1;
                    }
                }
                regions
            }
            _ => Vec::new(),
        };
        let splice_regions: Vec<(u64, u64)> = if diff_regions.is_empty() {
            vec![(var_start, splice_end)]
        } else {
            diff_regions.clone()
        };
        let ref_regions: Vec<(u64, u64)> = if diff_regions.is_empty() {
            vec![(var_start, var_end)]
        } else {
            diff_regions
        };
        let splice_any = |f: fn(&Transcript, u64, u64) -> bool| {
            splice_regions.iter().any(|&(s, e)| f(transcript, s, e))
        };
        let ref_any =
            |f: fn(&Transcript, u64, u64) -> bool| ref_regions.iter().any(|&(s, e)| f(transcript, s, e));

        let tr_start = transcript.start;
        let tr_end = transcript.end;

        let mut consequences = Vec::new();
        let mut cds_start = None;
        let mut cds_end = None;
        let mut protein_start = None;
        let mut protein_end = None;
        let mut amino_acids = None;
        let mut codons = None;
        let mut distance = None;

        // 1. Check if variant overlaps the transcript at all
        let overlaps = var_start <= tr_end && var_end >= tr_start;

        if !overlaps {
            // Check upstream/downstream
            let dist = self.distance_to_transcript(var_start, var_end, transcript);
            if let Some(d) = dist {
                distance = Some(d);
                let abs_dist = d.unsigned_abs();
                if abs_dist <= self.upstream_distance {
                    match transcript.strand {
                        Strand::Forward => {
                            if var_end < tr_start {
                                consequences.push(Consequence::UpstreamGeneVariant);
                            } else {
                                consequences.push(Consequence::DownstreamGeneVariant);
                            }
                        }
                        Strand::Reverse => {
                            if var_start > tr_end {
                                consequences.push(Consequence::UpstreamGeneVariant);
                            } else {
                                consequences.push(Consequence::DownstreamGeneVariant);
                            }
                        }
                    }
                }
            }

            if consequences.is_empty() {
                consequences.push(Consequence::IntergenicVariant);
            }

            let impact = Consequence::worst_impact(&consequences).unwrap_or(Impact::Modifier);
            return AlleleConsequenceResult {
                allele: alt_allele.clone(),
                consequences,
                impact,
                cdna_start: None,
                cdna_end: None,
                cds_start: None,
                cds_end: None,
                protein_start: None,
                protein_end: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance,
            };
        }

        // 2. Map to cDNA coordinates
        let cdna_start = transcript.genomic_to_cdna(var_start);
        let cdna_end = transcript.genomic_to_cdna(var_end);

        // 3. Determine exon/intron location
        // Use range-based overlap for exon detection to handle large indels
        let exon_info = transcript
            .exon_at(var_start)
            .or_else(|| transcript.exon_overlapping(var_start, var_end))
            .map(|(i, t)| (i as u32 + 1, t as u32));
        let intron_info = transcript
            .intron_at(var_start)
            .map(|(i, t)| (i as u32 + 1, t as u32));

        let in_exon = exon_info.is_some();
        let in_intron = intron_info.is_some();

        // 4. Check splice sites (always check regardless of coding status)
        if splice_any(splice::is_splice_donor) {
            consequences.push(Consequence::SpliceDonorVariant);
        }
        if splice_any(splice::is_splice_acceptor) {
            consequences.push(Consequence::SpliceAcceptorVariant);
        }

        // Extended splice terms, with Ensembl `VariationEffect.pm` precedence (NOT
        // an essential-splice blanket suppression — that over-suppressed):
        //   5th_base: no suppression;  donor_region: suppressed by 5th_base;
        //   polypyrimidine: no suppression;
        //   splice_region: suppressed by donor/acceptor/donor_region/5th_base (NOT
        //   by polypyrimidine, so polypyrimidine + splice_region co-occur).
        let is_donor = consequences.contains(&Consequence::SpliceDonorVariant);
        let is_acceptor = consequences.contains(&Consequence::SpliceAcceptorVariant);
        let is_essential_splice = is_donor || is_acceptor;
        // Whether the REFERENCE interval (not the alt-extended splice interval) reaches an
        // essential splice site. Only then does the variant straddle the exon/intron
        // boundary so the peptide is undeterminable (-> coding_sequence_variant). A delins
        // whose alt merely extends to the donor still has a determinable coding effect
        // (the ref is in-exon -> keep frameshift/missense).
        let ref_reaches_essential = ref_any(splice::is_splice_donor)
            || ref_any(splice::is_splice_acceptor);
        let is_donor_5th = splice_any(splice::is_splice_donor_5th_base);
        let is_donor_region =
            !is_donor_5th && splice_any(splice::is_splice_donor_region);
        if is_donor_5th {
            consequences.push(Consequence::SpliceDonorFifthBaseVariant);
        }
        if is_donor_region {
            consequences.push(Consequence::SpliceDonorRegionVariant);
        }
        // polypyrimidine is a candidate from the positional predicate; its `exon=0,
        // intron=1` include gate is applied declaratively below (`include_satisfied`).
        if splice_any(splice::is_splice_polypyrimidine_tract) {
            consequences.push(Consequence::SplicePolypyrimidineTractVariant);
        }
        if !is_essential_splice
            && !is_donor_region
            && !is_donor_5th
            && splice_any(splice::is_splice_region)
        {
            consequences.push(Consequence::SpliceRegionVariant);
        }

        // 5. Coding vs non-coding transcript
        if transcript.is_coding() {
            let coding_start = transcript.coding_region_start.unwrap_or(0);
            let coding_end = transcript.coding_region_end.unwrap_or(0);

            // Map to CDS coordinates
            if let Some(cs) = cdna_start {
                cds_start = transcript.cdna_to_cds(cs);
            }
            if let Some(ce) = cdna_end {
                cds_end = transcript.cdna_to_cds(ce);
            }
            if let Some(cds_s) = cds_start {
                protein_start = Some(Transcript::cds_to_protein(cds_s));
            }
            if let Some(cds_e) = cds_end {
                protein_end = Some(Transcript::cds_to_protein(cds_e));
            }

            let in_coding_region =
                self.is_in_coding_region(var_start, var_end, coding_start, coding_end);
            // Interval-overlap UTR flags (for the co-occurrence union). A variant
            // spanning a UTR into the CDS is BOTH the UTR term AND the coding term.
            let ov_5utr = self.overlaps_5utr(var_start, var_end, transcript);
            let ov_3utr = self.overlaps_3utr(var_start, var_end, transcript);
            // Straddles the CDS boundary: part coding, part UTR. Ensembl then cannot
            // call the length-based terms (frameshift/inframe/missense: cds bounds
            // undefined / partial codon) — only the start/stop boundary terms survive.
            let straddles_cds = in_coding_region && (ov_5utr || ov_3utr);

            // 5'/3' UTR terms as a UNION (exon=1 include gate applied later).
            if ov_5utr && in_exon {
                consequences.push(Consequence::FivePrimeUtrVariant);
            }
            if ov_3utr && in_exon {
                consequences.push(Consequence::ThreePrimeUtrVariant);
            }

            if in_coding_region && in_exon {
                // Coding exonic variant. TERMS = peptide predicate set + the CDS-boundary
                // straddle rules (resolved flatly in `resolve_coding_terms`); the
                // amino-acid/codon DISPLAY comes from the HGVS-validated computation.
                let display = self.predict_coding_consequence(
                    ref_allele, alt_allele, transcript, cds_start, cds_end,
                );
                let terms = self.coding_terms_for_variant(
                    ref_allele, alt_allele, transcript, cds_start, cds_end,
                );
                consequences.extend(Self::resolve_coding_terms(
                    terms,
                    display.as_ref().map(|(c, _, _)| *c),
                    straddles_cds,
                    ref_reaches_essential,
                ));
                if let Some((_, aa, cdn)) = display {
                    amino_acids = aa;
                    codons = cdn;
                }
            }
            // intron_variant is added below as a UNION (it co-occurs with the coding
            // term for boundary-spanning indels), not in this exclusive chain.
        } else {
            // Non-coding transcript
            if in_exon {
                consequences.push(Consequence::NonCodingTranscriptExonVariant);
            } else if in_intron {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
        }

        // intron_variant (SO:0001627) co-occurs whenever the variant overlaps an intron
        // INTERIOR — Ensembl's `within_intron`/`intronic` flag. A deletion spanning the
        // exon/intron boundary is therefore BOTH the coding/splice terms AND
        // intron_variant (a union, not an exclusive branch). The essential-splice
        // dinucleotides are excluded by `is_intronic`, so a variant only at the
        // donor/acceptor is not called intron_variant.
        if splice_any(splice::is_intronic)
            && !consequences.contains(&Consequence::IntronVariant)
        {
            consequences.push(Consequence::IntronVariant);
        }

        // Declarative OverlapConsequence `include` gate: keep only terms whose feature
        // requirements match the variant's overlap flags (computed once). This replaces
        // the scattered hand-wired guards with Ensembl's actual model.
        let flags = FeatureOverlap {
            exon: in_exon,
            intron: splice_any(splice::overlaps_intron),
            intron_boundary: splice_any(splice::overlaps_intron_boundary),
        };
        consequences.retain(|c| include_satisfied(*c, &flags));

        // If still no consequences, add catch-all
        if consequences.is_empty() {
            if transcript.is_coding() {
                consequences.push(Consequence::CodingTranscriptVariant);
            } else {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
        }

        // Add NMD_transcript_variant modifier for nonsense_mediated_decay transcripts
        if &*transcript.biotype == "nonsense_mediated_decay" {
            consequences.push(Consequence::NmdTranscriptVariant);
        }

        // Deduplicate
        consequences.sort_by_key(|c| c.rank());
        consequences.dedup();

        let impact = Consequence::worst_impact(&consequences).unwrap_or(Impact::Modifier);

        AlleleConsequenceResult {
            allele: alt_allele.clone(),
            consequences,
            impact,
            cdna_start,
            cdna_end,
            cds_start,
            cds_end,
            protein_start,
            protein_end,
            amino_acids,
            codons,
            exon: exon_info,
            intron: intron_info,
            distance,
        }
    }

    /// Translate every coding consequence TERM for a variant from the peptide-level
    /// `CodingContext` (the haplotype-ready abstraction). Returns the full term SET
    /// (Ensembl collects overlap-consequences; a variant can be e.g.
    /// frameshift_variant&stop_gained), or None if no coding context could be built
    /// (caller falls back). Single-variant path: one `CdsEdit`.
    fn coding_terms_for_variant(
        &self,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
    ) -> Option<Vec<Consequence>> {
        let edit = variant_to_cds_edit(ref_allele, alt_allele, cds_start, cds_end, transcript.strand)?;
        let ctx = self.build_coding_context(transcript, &[edit])?;
        Some(self.coding_consequence_terms(&ctx, transcript))
    }

    /// Resolve the coding TERMS for an exonic-coding variant, applying the CDS-boundary
    /// straddle rules flatly (no nesting in the caller). `raw` is the peptide predicate
    /// set; `legacy` is the single-term fallback when no coding context could be built.
    ///   - reaching an essential splice site -> the peptide is undeterminable (part is
    ///     intronic) -> the generic coding_sequence_variant;
    ///   - straddling a UTR -> only the start/stop boundary terms survive;
    ///   - otherwise -> all terms, pairing incomplete_terminal with coding_sequence.
    fn resolve_coding_terms(
        raw: Option<Vec<Consequence>>,
        legacy: Option<Consequence>,
        straddles_cds: bool,
        ref_reaches_essential: bool,
    ) -> Vec<Consequence> {
        let mut out = Vec::new();
        if ref_reaches_essential {
            out.push(Consequence::CodingSequenceVariant);
            return out;
        }
        let terms = match raw.or_else(|| legacy.map(|c| vec![c])) {
            Some(ts) => ts,
            None => {
                out.push(Consequence::CodingSequenceVariant);
                return out;
            }
        };
        let is_boundary = |t: Consequence| {
            matches!(
                t,
                Consequence::StartLost
                    | Consequence::StartRetainedVariant
                    | Consequence::StopGained
                    | Consequence::StopLost
                    | Consequence::StopRetainedVariant
            )
        };
        for t in terms {
            if straddles_cds && !is_boundary(t) {
                continue;
            }
            // VEP pairs incomplete_terminal_codon_variant with coding_sequence_variant.
            if t == Consequence::IncompleteTerminalCodonVariant {
                out.push(Consequence::CodingSequenceVariant);
            }
            out.push(t);
        }
        out
    }

    /// Build the peptide/codon view of a set of edits applied to the transcript CDS.
    /// One edit = one variant; many edits = a phased haplotype (haplosaurus /
    /// `bcftools csq`): they are applied to the SAME reference CDS before translation,
    /// so co-located variants combine into one correct protein consequence.
    fn build_coding_context(
        &self,
        transcript: &Transcript,
        edits: &[CdsEdit],
    ) -> Option<CodingContext> {
        let seq = transcript.translateable_seq.as_ref()?.as_bytes();
        if edits.is_empty() {
            return None;
        }
        let mut edits = edits.to_vec();
        edits.sort_by_key(|e| e.cds_idx);
        let first_idx = edits.first().unwrap().cds_idx;
        let last = edits.last().unwrap();
        let last_ref_end = last.cds_idx + last.ref_bases.len();
        if first_idx > seq.len() {
            return None;
        }
        let first_codon = first_idx / 3;
        // Number of REFERENCE codons the edits overlap. A pure insertion between two
        // codons (boundary, first_idx % 3 == 0, no ref bases) overlaps ZERO reference
        // codons — Ensembl's `codon` is then `-` and the alt window is ONLY the inserted
        // bases (no downstream base pulled in, so no spurious stop). A mid-codon
        // insertion or any ref-spanning edit overlaps the containing codon(s).
        let total_ref_bases: usize = edits.iter().map(|e| e.ref_bases.len()).sum();
        let (codon_start, codon_len) = if total_ref_bases == 0 && first_idx % 3 == 0 {
            (first_idx, 0)
        } else {
            let cs = first_codon * 3;
            let last_codon = last_ref_end.saturating_sub(1).max(first_idx) / 3;
            (cs, (last_codon - first_codon + 1) * 3)
        };
        // incomplete_terminal_codon_variant is a GENUINE cds_end_NF case: the CDS length
        // is not a multiple of 3, so its last codon is partial, and the variant lands in
        // that partial codon. A normal CDS whose variant merely EXTENDS past the stop
        // codon into the 3'UTR is NOT incomplete-terminal — the stop codon is complete,
        // so it translates (-> stop_lost/stop_retained, kept by the straddle filter).
        let incomplete_terminal =
            seq.len() % 3 != 0 && codon_start >= (seq.len() / 3) * 3 && codon_start < seq.len();
        let ref_codons = seq[codon_start..(codon_start + codon_len).min(seq.len())].to_vec();

        // VEP bounds the alt peptide to the affected window: the ref codon span plus the
        // net length change — NOT the whole downstream frame. Sum the edit lengths up
        // front so the reconstruction below can STOP at `window_len` instead of copying
        // the entire downstream CDS per variant (the per-row allocation hotspot — an SNV
        // only needs ~3 bytes, not the whole gene).
        let (total_ref, total_alt) = edits.iter().fold((0usize, 0usize), |(r, a), e| {
            (r + e.ref_bases.len(), a + e.alt_bases.len())
        });
        let window_len =
            (codon_len as i64 + total_alt as i64 - total_ref as i64).max(0) as usize;

        // Reconstruct only the first `window_len` bases of the alternate sequence from
        // codon_start, applying every edit; `push` truncates each slice to what is still
        // needed so the downstream tail is never fully copied.
        let mut alt: Vec<u8> = Vec::with_capacity(window_len);
        let push = |alt: &mut Vec<u8>, s: &[u8]| {
            if alt.len() < window_len {
                let take = (window_len - alt.len()).min(s.len());
                alt.extend_from_slice(&s[..take]);
            }
        };
        push(&mut alt, &seq[codon_start..first_idx.min(seq.len())]);
        let mut cursor = first_idx;
        for e in &edits {
            if e.cds_idx > cursor {
                let lo = cursor.min(seq.len());
                let hi = e.cds_idx.min(seq.len());
                push(&mut alt, &seq[lo..hi]);
            }
            push(&mut alt, &e.alt_bases);
            cursor = e.cds_idx + e.ref_bases.len();
        }
        if cursor < seq.len() {
            push(&mut alt, &seq[cursor..]);
        }
        let alt_codons = alt;

        let ct = self.ct(transcript);
        let translate = |s: &[u8]| -> Vec<u8> {
            s.chunks_exact(3)
                .map(|c| ct.translate(&[c[0], c[1], c[2]]))
                .collect()
        };

        // Ensembl `_overlaps_stop_codon` / `_ins_del_stop_altered` over CDS+3'UTR — the
        // primitive behind indel stop_lost vs stop_retained. The stop codon is the last
        // 3 CDS bases [seq.len()-3, seq.len()).
        let stop_lo = seq.len().saturating_sub(3);
        let span_lo = first_idx;
        let span_hi = last_ref_end.max(first_idx + 1); // half-open; insertion = a point
        // There must BE a terminal stop codon to overlap. A cds_end_NF transcript (the GFF3
        // may not tag it) ends in a non-stop codon, so Ensembl's `_overlaps_stop_codon`
        // returns 0 — no stop_lost/retained. Checking the last codon is the robust proxy.
        let has_terminal_stop = seq.len() >= 3
            && ct.translate(&[seq[stop_lo], seq[stop_lo + 1], seq[stop_lo + 2]]) == b'*';
        let overlaps_stop = has_terminal_stop && span_lo < seq.len() && span_hi > stop_lo;
        let stop_altered = if overlaps_stop && total_ref != total_alt {
            // CDS + 3'UTR (utr = the spliced cDNA after the CDS end).
            let mut cds_utr = seq.to_vec();
            if let (Some(spliced), Some(cend)) =
                (transcript.spliced_seq.as_ref(), transcript.cdna_coding_end)
            {
                let sb = spliced.as_bytes();
                let u = (cend as usize).min(sb.len());
                cds_utr.extend_from_slice(&sb[u..]);
            }
            // Apply every edit to the CDS+3'UTR (clipped to the available sequence).
            for e in edits.iter().rev() {
                if e.cds_idx <= cds_utr.len() {
                    let end = (e.cds_idx + e.ref_bases.len()).min(cds_utr.len());
                    cds_utr.splice(e.cds_idx..end, e.alt_bases.iter().copied());
                }
            }
            if cds_utr.len() < seq.len() {
                true // shortened past the CDS -> stop altered
            } else {
                // codon at the OLD stop position in the edited sequence — still a stop?
                let c = &cds_utr[stop_lo..stop_lo + 3];
                ct.translate(&[c[0], c[1], c[2]]) != b'*'
            }
        } else {
            false
        };

        Some(CodingContext {
            first_codon,
            incomplete_terminal,
            ref_pep: translate(&ref_codons),
            alt_pep: translate(&alt_codons),
            ref_codons,
            alt_codons,
            ref_len: total_ref,
            alt_len: total_alt,
            overlaps_stop,
            stop_altered,
        })
    }

    /// The flat predicate set over a `CodingContext`. Every applicable Ensembl
    /// overlap-consequence is collected (not a single hand-picked term), so
    /// co-occurring terms like frameshift_variant&stop_gained fall out naturally.
    fn coding_consequence_terms(
        &self,
        ctx: &CodingContext,
        transcript: &Transcript,
    ) -> Vec<Consequence> {
        // --- shared facts about the peptide change (computed once) ---
        let net = ctx.alt_len as i64 - ctx.ref_len as i64;
        let is_indel = net != 0;
        let frameshift_len = is_indel && (net % 3 != 0);
        let ref_stop = ctx.ref_pep.contains(&b'*');
        let alt_stop = ctx.alt_pep.contains(&b'*');
        let incomplete = ctx.incomplete_terminal;

        // The reference first codon is a real initiator (table-aware) of a complete CDS.
        let at_start = ctx.first_codon == 0
            && transcript.codon_table_start_phase == 0
            && !transcript.flags.iter().any(|f| f == "cds_start_NF")
            && ctx.ref_codons.len() >= 3
            && self
                .ct(transcript)
                .is_start(&[ctx.ref_codons[0], ctx.ref_codons[1], ctx.ref_codons[2]]);
        let alt_first = ctx.alt_pep.first().copied().unwrap_or(b'X');

        // Codon-level peptide terms are unavailable for an incomplete terminal codon,
        // and the start-codon terms pre-empt the rest of the peptide change (Ensembl runs
        // them, then the other peptide predicates `return 0 if start_lost`). `frameshift`
        // is the one length term that co-occurs with start_lost.
        let coding_ok = !incomplete && !at_start;

        // --- protein_altering_variant exclusions (VariationEffect.pm), shared ---
        let stop_starts =
            ctx.alt_pep.first() == Some(&b'*') || ctx.ref_pep.first() == Some(&b'*');
        // inframe-like check on the UNTRIMMED alt peptide (everything past the stop).
        let pa_excluded =
            stop_starts || ctx.alt_pep.starts_with(&ctx.ref_pep) || ctx.alt_pep.ends_with(&ctx.ref_pep);
        // inframe_insertion trims at the first stop before its prefix/suffix check.
        let altp_trimmed = {
            let s = ctx.alt_pep.as_slice();
            match s.iter().position(|&b| b == b'*') {
                Some(i) => &s[..=i],
                None => s,
            }
        };
        let inframe_ins = altp_trimmed.starts_with(&ctx.ref_pep) || altp_trimmed.ends_with(&ctx.ref_pep);
        let inframe_del = {
            let (r, a) = (ctx.ref_codons.as_slice(), ctx.alt_codons.as_slice());
            r.starts_with(a) || r.ends_with(a) || {
                let (rt, at) = trim_common(r, a);
                at.is_empty() && rt.len() % 3 == 0
            }
        };

        // --- independent predicates: each carries its own exclusions, none nested ---
        let p_start_lost = at_start && alt_first != b'M';
        let p_start_retained = at_start && alt_first == b'M';
        let p_stop_gained = coding_ok && alt_stop && !ref_stop;
        // DELETION at the stop codon: Ensembl uses `_ins_del_stop_altered` over CDS+3'UTR
        // (stop_altered -> stop_lost, !stop_altered -> stop_retained), NOT the windowed
        // peptide. Insertions use a genomic `consider_ins_len` overlap (not yet ported),
        // so they keep the windowed peptide behavior; substitutions use the peptide.
        let del = is_indel && net < 0;
        let p_stop_lost = coding_ok
            && if del {
                ctx.overlaps_stop && ctx.stop_altered
            } else {
                ref_stop && !alt_stop
            };
        let p_stop_retained = coding_ok
            && if del {
                ctx.overlaps_stop && !ctx.stop_altered
            } else {
                !is_indel && ref_stop && alt_stop && ctx.ref_pep == ctx.alt_pep
            };
        // VEP `frameshift` (VariationEffect.pm) suppresses the call when the affected
        // reference peptide starts with the stop (`$ref_pep =~ /^\*/` — the stop codon is
        // the first residue hit, so the frame past it is moot) or when the stop is
        // retained. Then it's stop_lost/stop_retained alone, never frameshift.
        let p_frameshift = !incomplete
            && frameshift_len
            && ctx.ref_pep.first() != Some(&b'*')
            && !p_stop_retained;
        let inframe = is_indel && !frameshift_len;
        let p_inframe_insertion = coding_ok && inframe && net > 0 && inframe_ins;
        let p_inframe_deletion = coding_ok && inframe && net < 0 && inframe_del;
        let p_protein_altering = coding_ok
            && inframe
            && !(net > 0 && inframe_ins)
            && !(net < 0 && inframe_del)
            && !pa_excluded;
        let p_missense =
            coding_ok && !is_indel && !alt_stop && !ref_stop && ctx.ref_pep != ctx.alt_pep;
        let p_synonymous =
            coding_ok && !is_indel && !alt_stop && !ref_stop && ctx.ref_pep == ctx.alt_pep;

        // --- flat collection (filter-style); order is cosmetic, terms re-ranked later ---
        let mut out = Vec::new();
        for (fired, term) in [
            (incomplete, Consequence::IncompleteTerminalCodonVariant),
            (p_start_lost, Consequence::StartLost),
            (p_start_retained, Consequence::StartRetainedVariant),
            (p_stop_gained, Consequence::StopGained),
            (p_stop_lost, Consequence::StopLost),
            (p_stop_retained, Consequence::StopRetainedVariant),
            (p_frameshift, Consequence::FrameshiftVariant),
            (p_inframe_insertion, Consequence::InframeInsertion),
            (p_inframe_deletion, Consequence::InframeDeletion),
            (p_protein_altering, Consequence::ProteinAlteringVariant),
            (p_missense, Consequence::MissenseVariant),
            (p_synonymous, Consequence::SynonymousVariant),
        ] {
            if fired {
                out.push(term);
            }
        }
        if out.is_empty() {
            out.push(Consequence::CodingSequenceVariant);
        }
        out
    }

    /// Predict the coding consequence (missense, synonymous, frameshift, etc.)
    fn predict_coding_consequence(
        &self,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
    ) -> Option<(
        Consequence,
        Option<(String, String)>,
        Option<(String, String)>,
    )> {
        let cds_pos_start = cds_start?;

        let ref_len = ref_allele.len();
        let alt_len = alt_allele.len();

        // Check if this is a frameshift or in-frame indel
        let is_deletion = *ref_allele != Allele::Deletion && *alt_allele == Allele::Deletion;
        let is_insertion = *ref_allele == Allele::Deletion && *alt_allele != Allele::Missing;
        let is_indel = is_deletion || is_insertion || ref_len != alt_len;

        if is_indel {
            let (consequence, is_frameshift) = if is_deletion || is_insertion {
                let indel_len = if is_deletion { ref_len } else { alt_len };
                if indel_len % 3 != 0 {
                    (Consequence::FrameshiftVariant, true)
                } else if is_insertion {
                    (Consequence::InframeInsertion, false)
                } else {
                    (Consequence::InframeDeletion, false)
                }
            } else {
                let len_diff = (ref_len as i64 - alt_len as i64).unsigned_abs() as usize;
                if len_diff % 3 != 0 {
                    (Consequence::FrameshiftVariant, true)
                } else if ref_len > alt_len {
                    (Consequence::InframeDeletion, false)
                } else {
                    (Consequence::InframeInsertion, false)
                }
            };

            // For deletions on reverse strand, cds_start maps to the end of the
            // deletion in CDS space. Use the lower CDS position as the start.
            let cds_pos = if is_deletion {
                match cds_end {
                    Some(ce) => cds_pos_start.min(ce),
                    None => cds_pos_start,
                }
            } else {
                cds_pos_start
            };

            // Try to compute amino acids and codons from translateable_seq
            let (aa_pair, codon_pair) = self.compute_indel_amino_acids(
                transcript,
                cds_pos,
                ref_allele,
                alt_allele,
                is_frameshift,
            );

            return Some((consequence, aa_pair, codon_pair));
        }

        // Same length substitution (SNV or MNV)
        if let Some(ref translateable_seq) = transcript.translateable_seq {
            let seq_bytes = translateable_seq.as_bytes();

            // Get the codon containing this CDS position
            let codon_number = ((cds_pos_start - 1) / 3) as usize;
            let codon_offset = ((cds_pos_start - 1) % 3) as usize;
            let codon_start = codon_number * 3;

            // Incomplete terminal codon (Ensembl cds_end_NF / CDS length not a
            // multiple of 3): the last codon runs past the translateable sequence
            // and cannot be translated. VEP reports incomplete_terminal_codon_variant.
            if codon_start < seq_bytes.len() && codon_start + 3 > seq_bytes.len() {
                return Some((Consequence::IncompleteTerminalCodonVariant, None, None));
            }

            if codon_start + 3 <= seq_bytes.len() {
                let ref_codon = [
                    seq_bytes[codon_start],
                    seq_bytes[codon_start + 1],
                    seq_bytes[codon_start + 2],
                ];
                let mut alt_codon = ref_codon;

                // Apply the substitution
                if let Allele::Sequence(alt_bases) = alt_allele {
                    for (i, &base) in alt_bases.iter().enumerate() {
                        let pos = codon_offset + i;
                        if pos < 3 {
                            alt_codon[pos] = match transcript.strand {
                                Strand::Forward => base,
                                Strand::Reverse => complement(base),
                            };
                        }
                    }
                }

                let ref_aa = self.ct(transcript).translate(&ref_codon);
                let alt_aa = self.ct(transcript).translate(&alt_codon);

                let (ref_codon_str, alt_codon_str) = format_codon_change(&ref_codon, &alt_codon);

                let ref_aa_str = String::from(ref_aa as char);
                let alt_aa_str = String::from(alt_aa as char);

                let codon_pair = Some((ref_codon_str, alt_codon_str));
                let aa_pair = Some((ref_aa_str.clone(), alt_aa_str.clone()));

                // Start codon. `start_lost` is, per Ensembl, a change to the
                // *canonical start codon* — so it requires the REFERENCE first codon
                // to actually be a start (ATG) that the variant destroys. Checking
                // only the alt codon (as upstream did) mis-fires on transcripts whose
                // annotated CDS does not begin with ATG: incomplete 5' CDS
                // (`cds_start_NF`, non-zero start phase), or any non-ATG annotated
                // start. A non-zero start phase additionally means codon 0 is the
                // incomplete leading codon (untranslatable -> coding_sequence_variant).
                if codon_number == 0 {
                    if transcript.codon_table_start_phase != 0 {
                        return Some((Consequence::CodingSequenceVariant, None, None));
                    }
                    let incomplete_start = transcript.flags.iter().any(|f| f == "cds_start_NF");
                    // Gate on the reference first codon being a real initiator *in this
                    // transcript's codon table* (ATG, or the mito alt-starts ATT/ATC/…).
                    // This both (a) skips first codons that are not genuine starts — e.g.
                    // a transcript whose modelled CDS start does not align with VEP's
                    // initiator, which must stay missense/synonymous — and (b) recognises
                    // non-ATG mito starts the old standard-table-only check dropped.
                    if !incomplete_start && self.ct(transcript).is_start(&ref_codon) {
                        // The initiator always encodes Methionine ('M') regardless of the
                        // triplet (mito ATT/ATC translate to Ile). Ensembl's start_lost
                        // therefore tests whether the alt codon still yields Met, NOT
                        // whether the alt is itself a start codon (the old
                        // `!is_start(alt)` gate, which mis-fired on start->start changes
                        // like ATT->ATC). If the alt still encodes Met it is
                        // start_retained, not lost.
                        let start_pair = Some(("M".to_string(), alt_aa_str.clone()));
                        if alt_aa != b'M' {
                            return Some((Consequence::StartLost, start_pair, codon_pair));
                        }
                        return Some((
                            Consequence::StartRetainedVariant,
                            Some(("M".to_string(), "M".to_string())),
                            codon_pair,
                        ));
                    }
                }

                // Determine consequence type
                if ref_aa == alt_aa {
                    if ref_aa == b'*' {
                        return Some((Consequence::StopRetainedVariant, aa_pair, codon_pair));
                    }
                    return Some((Consequence::SynonymousVariant, aa_pair, codon_pair));
                }

                if alt_aa == b'*' {
                    return Some((Consequence::StopGained, aa_pair, codon_pair));
                }
                if ref_aa == b'*' {
                    return Some((Consequence::StopLost, aa_pair, codon_pair));
                }

                return Some((Consequence::MissenseVariant, aa_pair, codon_pair));
            }
        }

        // Fallback: if we can't determine the exact consequence,
        // classify based on whether it's an in-frame or frameshift change
        if ref_len == alt_len {
            Some((Consequence::MissenseVariant, None, None))
        } else {
            None
        }
    }

    /// Compute amino acids and codons affected by an indel variant.
    /// Returns (amino_acids, codons) tuples.
    /// For frameshifts: ref codon with VEP-style case formatting, truncated alt codon.
    fn compute_indel_amino_acids(
        &self,
        transcript: &Transcript,
        cds_pos: u64,
        ref_allele: &Allele,
        alt_allele: &Allele,
        is_frameshift: bool,
    ) -> (Option<(String, String)>, Option<(String, String)>) {
        let translateable_seq = match transcript.translateable_seq.as_ref() {
            Some(s) => s,
            None => return (None, None),
        };
        let seq_bytes = translateable_seq.as_bytes();
        let cds_idx = (cds_pos - 1) as usize;

        if cds_idx >= seq_bytes.len() {
            return (None, None);
        }

        // Get the codon at the affected position
        let codon_number = cds_idx / 3;
        let codon_offset = cds_idx % 3;
        let codon_start = codon_number * 3;

        if codon_start + 3 > seq_bytes.len() {
            return (None, None);
        }

        let ref_codon = [
            seq_bytes[codon_start],
            seq_bytes[codon_start + 1],
            seq_bytes[codon_start + 2],
        ];
        let ref_aa = self.ct(transcript).translate(&ref_codon);
        let ref_aa_str = String::from(ref_aa as char);

        if is_frameshift {
            // Build the alt sequence by applying the indel
            let mut alt_seq: Vec<u8> = seq_bytes.to_vec();

            match (ref_allele, alt_allele) {
                (Allele::Sequence(_), Allele::Deletion) => {
                    let del_len = ref_allele.len();
                    let end = (cds_idx + del_len).min(alt_seq.len());
                    alt_seq.drain(cds_idx..end);
                }
                (Allele::Deletion, Allele::Sequence(ins_bases)) => {
                    let mut bases: Vec<u8> = ins_bases.clone();
                    if transcript.strand == Strand::Reverse {
                        bases = bases.iter().map(|&b| complement(b)).collect();
                    }
                    for (i, &b) in bases.iter().enumerate() {
                        alt_seq.insert(cds_idx + i, b);
                    }
                }
                (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases)) => {
                    let end = (cds_idx + ref_bases.len()).min(alt_seq.len());
                    let mut replacement = alt_bases.clone();
                    if transcript.strand == Strand::Reverse {
                        replacement = replacement.iter().map(|&b| complement(b)).collect();
                    }
                    alt_seq.splice(cds_idx..end, replacement);
                }
                _ => return (None, None),
            }

            // Build codon display: VEP style with deleted base uppercase
            // ref codon: lowercase bases, uppercase at the deleted position(s)
            let mut ref_codon_display = String::with_capacity(3);
            for i in 0..3 {
                if i == codon_offset {
                    ref_codon_display.push((ref_codon[i] as char).to_ascii_uppercase());
                } else {
                    ref_codon_display.push((ref_codon[i] as char).to_ascii_lowercase());
                }
            }

            // alt codon: show only the remaining bases of the original codon after the indel
            // For a deletion at offset 2 in a 3-base codon: show only the 2 remaining bases
            let alt_codon_display: String = {
                let mut original_codon: Vec<u8> = ref_codon.to_vec();
                match (ref_allele, alt_allele) {
                    (Allele::Sequence(_), Allele::Deletion) => {
                        // Remove the deleted base(s) from the codon
                        let del_len = ref_allele.len().min(3 - codon_offset);
                        let end = (codon_offset + del_len).min(original_codon.len());
                        original_codon.drain(codon_offset..end);
                    }
                    (Allele::Deletion, Allele::Sequence(ins_bases)) => {
                        // Insert bases into the codon at the offset
                        let mut bases = ins_bases.clone();
                        if transcript.strand == Strand::Reverse {
                            bases = bases.iter().map(|&b| complement(b)).collect();
                        }
                        for (j, &b) in bases.iter().enumerate() {
                            original_codon.insert(codon_offset + j, b);
                        }
                    }
                    _ => {}
                }
                original_codon
                    .iter()
                    .map(|&b| (b as char).to_ascii_lowercase())
                    .collect()
            };

            // For frameshifts, alt amino acid is always X (unknown/frameshift)
            // For pure insertions, VEP uses "-" for ref amino acid/codon
            // and only the inserted bases for alt codon
            let (fs_ref_aa, fs_ref_codon, fs_alt_codon) = if *ref_allele == Allele::Deletion {
                let ins_codon = if let Allele::Sequence(ins_bases) = alt_allele {
                    let mut bases = ins_bases.clone();
                    if transcript.strand == Strand::Reverse {
                        bases = bases.iter().map(|&b| complement(b)).collect();
                    }
                    bases
                        .iter()
                        .map(|&b| (b as char).to_ascii_uppercase())
                        .collect::<String>()
                } else {
                    alt_codon_display
                };
                ("-".to_string(), "-".to_string(), ins_codon)
            } else {
                (ref_aa_str, ref_codon_display, alt_codon_display)
            };
            let aa_pair = Some((fs_ref_aa, "X".to_string()));
            let codon_pair = Some((fs_ref_codon, fs_alt_codon));
            (aa_pair, codon_pair)
        } else {
            // In-frame indel: build alt sequence and translate affected codons
            let mut alt_seq: Vec<u8> = seq_bytes.to_vec();
            match (ref_allele, alt_allele) {
                (Allele::Sequence(_), Allele::Deletion) => {
                    // In-frame deletion: remove bases and compare ref/alt amino acids
                    let del_len = ref_allele.len();
                    let end = (cds_idx + del_len).min(alt_seq.len());
                    alt_seq.drain(cds_idx..end);

                    // Number of complete codons deleted
                    let del_codons = del_len / 3;

                    if codon_offset == 0 {
                        // Deletion starts at codon boundary: VEP shows deleted AAs vs "-"
                        let ref_end = (codon_start + del_codons * 3).min(seq_bytes.len());
                        let ref_region = &seq_bytes[codon_start..ref_end];
                        let ref_aas: String = ref_region
                            .chunks(3)
                            .filter(|c| c.len() == 3)
                            .map(|c| self.ct(transcript).translate(&[c[0], c[1], c[2]]) as char)
                            .collect();
                        let ref_codons: String = ref_region
                            .iter()
                            .map(|&b| (b as char).to_uppercase().next().unwrap())
                            .collect();
                        let aa_pair = Some((ref_aas, "-".to_string()));
                        let codon_pair = Some((ref_codons, "-".to_string()));
                        return (aa_pair, codon_pair);
                    } else {
                        // Deletion within a codon: show affected codons ref and alt
                        let n_ref_codons = del_codons + 1;
                        let ref_end = (codon_start + n_ref_codons * 3).min(seq_bytes.len());
                        let ref_region = &seq_bytes[codon_start..ref_end];
                        let ref_aas: String = ref_region
                            .chunks(3)
                            .filter(|c| c.len() == 3)
                            .map(|c| self.ct(transcript).translate(&[c[0], c[1], c[2]]) as char)
                            .collect();
                        let alt_codon_end = (codon_start + 3).min(alt_seq.len());
                        let alt_region = &alt_seq[codon_start..alt_codon_end];
                        let alt_aas: String = if alt_region.len() == 3 {
                            String::from(self.ct(transcript).translate(&[
                                alt_region[0],
                                alt_region[1],
                                alt_region[2],
                            ]) as char)
                        } else {
                            "-".to_string()
                        };
                        let ref_codons: String = ref_region
                            .iter()
                            .map(|&b| (b as char).to_uppercase().next().unwrap())
                            .collect();
                        let alt_codons: String = if alt_aas == "-" {
                            "-".to_string()
                        } else {
                            alt_region
                                .iter()
                                .map(|&b| (b as char).to_uppercase().next().unwrap())
                                .collect()
                        };
                        let aa_pair = Some((ref_aas, alt_aas));
                        let codon_pair = Some((ref_codons, alt_codons));
                        return (aa_pair, codon_pair);
                    }
                }
                (Allele::Deletion, Allele::Sequence(ins_bases)) => {
                    // In-frame insertion: reverse-complement for reverse strand
                    let mut bases: Vec<u8> = ins_bases.clone();
                    if transcript.strand == Strand::Reverse {
                        bases = bases.iter().rev().map(|&b| complement(b)).collect();
                    }
                    // For reverse strand, the VCF insertion point maps to one base
                    // earlier in CDS space, so shift the insertion index by 1
                    let ins_idx = if transcript.strand == Strand::Reverse {
                        cds_idx + 1
                    } else {
                        cds_idx
                    };
                    for (i, &b) in bases.iter().enumerate() {
                        if ins_idx + i <= alt_seq.len() {
                            alt_seq.insert(ins_idx + i, b);
                        }
                    }

                    // Ref: the single codon at the insertion point
                    let ref_codon_str: String = ref_codon
                        .iter()
                        .map(|&b| (b as char).to_lowercase().next().unwrap())
                        .collect();

                    // Alt: translate codons spanning the insertion
                    let ins_codons = (bases.len() / 3) + 1;
                    let alt_end = (codon_start + ins_codons * 3).min(alt_seq.len());
                    let alt_region = &alt_seq[codon_start..alt_end];
                    let alt_aas: String = alt_region
                        .chunks(3)
                        .filter(|c| c.len() == 3)
                        .map(|c| self.ct(transcript).translate(&[c[0], c[1], c[2]]) as char)
                        .collect();

                    // Build alt codon string: original bases lowercase, inserted uppercase
                    let ins_offset_in_codon = if transcript.strand == Strand::Reverse {
                        codon_offset + 1
                    } else {
                        codon_offset
                    };
                    let mut alt_codon_display = String::new();
                    for (i, &b) in alt_region.iter().enumerate() {
                        let is_original = if i < ins_offset_in_codon {
                            true
                        } else if i >= ins_offset_in_codon + bases.len() {
                            true
                        } else {
                            false
                        };
                        if is_original {
                            alt_codon_display.push((b as char).to_lowercase().next().unwrap());
                        } else {
                            alt_codon_display.push((b as char).to_uppercase().next().unwrap());
                        }
                    }

                    let aa_pair = Some((ref_aa_str, alt_aas));
                    let codon_pair = Some((ref_codon_str, alt_codon_display));
                    return (aa_pair, codon_pair);
                }
                _ => {}
            }
            (Some((ref_aa_str, "X".to_string())), None)
        }
    }

    fn distance_to_transcript(
        &self,
        var_start: u64,
        var_end: u64,
        transcript: &Transcript,
    ) -> Option<i64> {
        // For insertions (end < start), use start for distance calculation
        // since start represents the actual insertion position
        let effective_start = var_start.min(var_end);
        let effective_end = var_start.max(var_end);
        if effective_end < transcript.start {
            Some((transcript.start - effective_end) as i64)
        } else if effective_start > transcript.end {
            Some((effective_start - transcript.end) as i64)
        } else {
            Some(0)
        }
    }

    fn is_in_coding_region(
        &self,
        var_start: u64,
        var_end: u64,
        coding_start: u64,
        coding_end: u64,
    ) -> bool {
        var_start <= coding_end && var_end >= coding_start
    }

    /// Does the variant INTERVAL overlap the 5'/3' UTR region (for the co-occurrence
    /// union — a variant spanning a UTR into the CDS is BOTH the UTR term and the
    /// coding term in Ensembl). The exon=1 include gate is applied separately.
    fn overlaps_5utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        let (s, e) = (var_start.min(var_end), var_start.max(var_end));
        let (cs, ce) = match (transcript.coding_region_start, transcript.coding_region_end) {
            (Some(a), Some(b)) => (a, b),
            _ => return false,
        };
        match transcript.strand {
            Strand::Forward => s <= cs.saturating_sub(1) && e >= transcript.start,
            Strand::Reverse => e >= ce + 1 && s <= transcript.end,
        }
    }
    fn overlaps_3utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        let (s, e) = (var_start.min(var_end), var_start.max(var_end));
        let (cs, ce) = match (transcript.coding_region_start, transcript.coding_region_end) {
            (Some(a), Some(b)) => (a, b),
            _ => return false,
        };
        match transcript.strand {
            Strand::Forward => e >= ce + 1 && s <= transcript.end,
            Strand::Reverse => s <= cs.saturating_sub(1) && e >= transcript.start,
        }
    }


}

fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
    }
}

impl Default for ConsequencePredictor {
    fn default() -> Self {
        Self::new(5000, 5000)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_genome::{Exon, Gene, Translation};

    fn make_coding_transcript() -> Transcript {
        // A simple protein-coding transcript on forward strand:
        // Exon1: 1000-1200 (UTR: 1000-1049, CDS: 1050-1200)
        // Intron: 1201-1999
        // Exon2: 2000-2300 (all CDS)
        // Intron: 2301-3999
        // Exon3: 4000-5000 (CDS: 4000-4500, UTR: 4501-5000)
        //
        // CDS length: 151 + 301 + 501 = 953 bases
        // translateable_seq: from cDNA pos 51 to 953+50=1003
        let translateable = "ATGGCTTCAAAGCCC".to_string() + &"A".repeat(938); // starts with ATG

        Transcript {
            stable_id: "ENST00000001".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG00000001".into(),
                symbol: Some("TESTGENE".into()),
                symbol_source: Some("HGNC".into()),
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 5000,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 5000,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 2000,
                    end: 2300,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 1,
                    rank: 2,
                },
                Exon {
                    stable_id: "E3".into(),
                    start: 4000,
                    end: 5000,
                    strand: Strand::Forward,
                    phase: 1,
                    end_phase: -1,
                    rank: 3,
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP00000001".into(),
                genomic_start: 1050,
                genomic_end: 4500,
                start_exon_rank: 1,
                start_exon_offset: 50,
                end_exon_rank: 3,
                end_exon_offset: 500,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(1003),
            coding_region_start: Some(1050),
            coding_region_end: Some(4500),
            spliced_seq: None,
            translateable_seq: Some(translateable),
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: None,
            ccds: None,
            protein_id: Some("ENSP00000001".into()),
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        }
    }

    fn make_noncoding_transcript() -> Transcript {
        Transcript {
            stable_id: "ENST_NC".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_NC".into(),
                symbol: Some("NCRNA1".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "lncRNA".into(),
                chromosome: "chr1".into(),
                start: 10000,
                end: 12000,
                strand: Strand::Forward,
            },
            biotype: "lncRNA".into(),
            chromosome: "chr1".into(),
            start: 10000,
            end: 12000,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 10000,
                    end: 10500,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: -1,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 11500,
                    end: 12000,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: -1,
                    rank: 2,
                },
            ],
            translation: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        }
    }

    #[test]
    fn test_upstream_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let pos = GenomicPosition::new("chr1", 500, 500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        assert_eq!(result.transcript_consequences.len(), 1);
        let tc = &result.transcript_consequences[0];
        assert_eq!(tc.allele_consequences.len(), 1);
        assert!(tc.allele_consequences[0]
            .consequences
            .contains(&Consequence::UpstreamGeneVariant));
        assert_eq!(tc.allele_consequences[0].distance, Some(500));
    }

    #[test]
    fn test_downstream_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let pos = GenomicPosition::new("chr1", 5500, 5500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac
            .consequences
            .contains(&Consequence::DownstreamGeneVariant));
        assert_eq!(ac.distance, Some(500));
    }

    #[test]
    fn test_intergenic_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Very far away
        let pos = GenomicPosition::new("chr1", 100000, 100000, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::IntergenicVariant));
    }

    #[test]
    fn test_5_prime_utr_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1020 is in exon1 (1000-1200), before CDS start (1050)
        let pos = GenomicPosition::new("chr1", 1020, 1020, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::FivePrimeUtrVariant));
    }

    #[test]
    fn test_3_prime_utr_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 4600 is in exon3 (4000-5000), after CDS end (4500)
        let pos = GenomicPosition::new("chr1", 4600, 4600, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::ThreePrimeUtrVariant));
    }

    #[test]
    fn test_intron_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1500 is in intron1 (1201-1999), away from splice sites
        let pos = GenomicPosition::new("chr1", 1500, 1500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_splice_donor_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1201 is first base of intron1 → splice donor
        let pos = GenomicPosition::new("chr1", 1201, 1201, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::SpliceDonorVariant));
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_splice_acceptor_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1999 is last base of intron1 → splice acceptor
        let pos = GenomicPosition::new("chr1", 1999, 1999, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac
            .consequences
            .contains(&Consequence::SpliceAcceptorVariant));
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_synonymous_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // First CDS position is at genomic 1050, which is cDNA pos 51, CDS pos 1
        // translateable_seq starts with "ATG" (Met)
        // CDS pos 3 (third base of first codon) - change G to A: ATA still codes for... wait
        // ATG -> Met. Let's change position 3 of ATG from G to something that's still Met: not possible
        // Let's use a different codon. CDS pos 4-6 is "GCT" (Ala). GCC also codes for Ala.
        // Genomic pos for CDS pos 4 = 1050 + 3 = 1053
        // Change T at CDS pos 6 to C: GCT -> GCC both = Ala
        let pos = GenomicPosition::new("chr1", 1055, 1055, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("T"),
            &[Allele::from_str("C")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::SynonymousVariant),
            "Expected synonymous, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::Low);
    }

    #[test]
    fn test_missense_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 4-6 is "GCT" (Ala). Change first base G to T: TCT = Ser (different!)
        // Genomic pos for CDS pos 4 = 1050 + 3 = 1053
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("T")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::MissenseVariant),
            "Expected missense, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::Moderate);
    }

    #[test]
    fn test_stop_gained() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 4-6 is "GCT". Change to "TAA" (stop) → need to change pos 4,5,6
        // For simplicity, change CDS pos 4: G->T, pos 5: C->A, pos 6: T->A
        // But our predictor works one SNV at a time. Let's pick a codon that's one base
        // away from a stop. "TCA" (Ser) → change C→A: TAA (stop). But we'd need that codon.
        // Actually, translateable_seq[6..9] = "TCA" (positions 7-9 in 1-based)
        // CDS pos 7 is at genomic 1050+6 = 1056
        // Change T to T (no), we need C at pos 8 to become something.
        // Let's just use translateable[3..6] = "GCT" and change pos 4 (G) to T: "TCT" = Ser
        // That's missense, not stop. Let's try another approach.
        // translateable[9..12] = "AAG" (Lys). Change A at pos 10 to T: TAG = stop!
        // CDS pos 10 is at genomic 1050+9 = 1059
        let pos = GenomicPosition::new("chr1", 1059, 1059, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("T")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::StopGained),
            "Expected stop_gained, got: {:?}. translateable[9..12]={:?}",
            ac.consequences,
            &tr.translateable_seq.as_ref().unwrap()[9..12]
        );
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_frameshift_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Deletion of 1 base in CDS → frameshift
        // CDS pos 4 at genomic 1053
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::Deletion],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::FrameshiftVariant),
            "Expected frameshift, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_inframe_deletion() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Deletion of 3 bases in CDS → inframe deletion
        let pos = GenomicPosition::new("chr1", 1053, 1055, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("GCT"),
            &[Allele::Deletion],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::InframeDeletion),
            "Expected inframe_deletion, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::Moderate);
    }

    #[test]
    fn test_noncoding_exon_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_noncoding_transcript();
        let pos = GenomicPosition::new("chr1", 10100, 10100, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant));
    }

    #[test]
    fn test_start_lost() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 1 (first base of ATG) is at genomic 1050
        // Change A to G: GTG is not a standard start codon
        let pos = GenomicPosition::new("chr1", 1050, 1050, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::StartLost),
            "Expected start_lost, got: {:?}",
            ac.consequences
        );
    }

    // A vertebrate-mito transcript whose annotated initiator is ATT (a valid mito
    // start that *translates* to Ile). The initiator always encodes Met, so a change
    // that no longer yields Met is start_lost even though ATT->ATC is a start->start
    // change — the old `is_start(ref) && !is_start(alt)` gate (standard-table only)
    // missed this. ATT->ATA stays Met (mito) -> start_retained, not lost.
    fn make_mito_start_transcript() -> Transcript {
        let mut tr = make_coding_transcript();
        tr.chromosome = "MT".into();
        tr.gene.chromosome = "MT".into();
        // CDS begins at genomic 1050 with ATT (vertebrate-mito initiator).
        tr.translateable_seq = Some("ATTGCTTCAAAGCCC".to_string() + &"A".repeat(938));
        tr
    }

    #[test]
    fn test_mito_start_lost_non_atg() {
        let predictor = ConsequencePredictor::default();
        let tr = make_mito_start_transcript();
        // 3rd base of the ATT initiator is genomic 1052; T->C => ATC (Ile, not Met).
        let pos = GenomicPosition::new("MT", 1052, 1052, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("T"),
            &[Allele::from_str("C")],
            &[&tr],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::StartLost),
            "ATT->ATC at the mito initiator should be start_lost, got: {:?}",
            ac.consequences
        );
    }

    #[test]
    fn test_mito_start_retained_to_met() {
        let predictor = ConsequencePredictor::default();
        let tr = make_mito_start_transcript();
        // ATT->ATA: ATA is Met in the vertebrate-mito table, so the initiator is kept.
        let pos = GenomicPosition::new("MT", 1052, 1052, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("T"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences
                .contains(&Consequence::StartRetainedVariant),
            "ATT->ATA (still Met in mito) should be start_retained, got: {:?}",
            ac.consequences
        );
    }

    // The haplotype-ready abstraction: two phased edits in the SAME codon are applied
    // together to the reference CDS before translation (haplosaurus / bcftools csq),
    // yielding the COMBINED amino acid — not two independent per-base calls. codon 2 of
    // make_coding_transcript is GCT (Ala); editing base 1 (G->T) and base 3 (T->A)
    // together gives TCA (Ser), which neither edit produces alone.
    #[test]
    fn test_haplotype_two_edits_one_codon() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let edits = vec![
            CdsEdit { cds_idx: 3, ref_bases: vec![b'G'], alt_bases: vec![b'T'] },
            CdsEdit { cds_idx: 5, ref_bases: vec![b'T'], alt_bases: vec![b'A'] },
        ];
        let ctx = predictor.build_coding_context(&tr, &edits).unwrap();
        assert_eq!(ctx.ref_pep, b"A", "reference codon GCT = Ala");
        assert_eq!(
            ctx.alt_pep, b"S",
            "combined edits give TCA = Ser (proves edits applied together, not independently)"
        );
        assert!(predictor
            .coding_consequence_terms(&ctx, &tr)
            .contains(&Consequence::MissenseVariant));
    }

    // An inframe delins that rearranges the codon sequence (alt is NOT a clean
    // prefix/suffix/internal match of ref) is protein_altering_variant, not
    // inframe_deletion — matching Ensembl. codon 2-3 = GCTTCA (Ala-Ser); replacing it
    // with TTT (Phe) deletes 3 bases inframe but keeps none of the ref at either end.
    #[test]
    fn test_inframe_delins_is_protein_altering() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let delins = vec![CdsEdit {
            cds_idx: 3,
            ref_bases: b"GCTTCA".to_vec(),
            alt_bases: b"TTT".to_vec(),
        }];
        let ctx = predictor.build_coding_context(&tr, &delins).unwrap();
        assert!(predictor
            .coding_consequence_terms(&ctx, &tr)
            .contains(&Consequence::ProteinAlteringVariant));
        // A clean deletion of TCA (ref keeps the GCT prefix) IS inframe_deletion.
        let clean = vec![CdsEdit {
            cds_idx: 3,
            ref_bases: b"GCTTCA".to_vec(),
            alt_bases: b"GCT".to_vec(),
        }];
        let ctx2 = predictor.build_coding_context(&tr, &clean).unwrap();
        assert!(predictor
            .coding_consequence_terms(&ctx2, &tr)
            .contains(&Consequence::InframeDeletion));
    }

    // Insertion stop_gained uses Ensembl's exact window: a boundary insertion's alt
    // codon is ONLY the inserted bases (no downstream base pulled in). Inserting a
    // whole stop codon -> stop_gained; inserting bases that don't form a stop codon ->
    // bare frameshift (no spurious stop from a trailing reference base).
    #[test]
    fn test_insertion_stop_gained_exact_window() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Boundary insertion (cds_idx 6) of TAATT: first whole codon TAA = stop.
        let with_stop = vec![CdsEdit { cds_idx: 6, ref_bases: vec![], alt_bases: b"TAATT".to_vec() }];
        let ctx = predictor.build_coding_context(&tr, &with_stop).unwrap();
        let terms = predictor.coding_consequence_terms(&ctx, &tr);
        assert!(terms.contains(&Consequence::StopGained), "TAATT inserts a stop codon");
        assert!(terms.contains(&Consequence::FrameshiftVariant));
        // Inserting TT (no whole codon) must NOT pull in a downstream base to form a stop.
        let no_stop = vec![CdsEdit { cds_idx: 6, ref_bases: vec![], alt_bases: b"TT".to_vec() }];
        let ctx2 = predictor.build_coding_context(&tr, &no_stop).unwrap();
        let terms2 = predictor.coding_consequence_terms(&ctx2, &tr);
        assert!(!terms2.contains(&Consequence::StopGained), "TT alone is not a stop");
        assert!(terms2.contains(&Consequence::FrameshiftVariant));
    }

    #[test]
    fn test_multiple_transcripts() {
        let predictor = ConsequencePredictor::default();
        let tr1 = make_coding_transcript();
        let tr2 = make_noncoding_transcript();

        // Position in tr1's intron, not overlapping tr2
        let pos = GenomicPosition::new("chr1", 1500, 1500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr1, &tr2],
            None,
        );

        assert_eq!(result.transcript_consequences.len(), 2);
        // tr1: intron variant
        assert!(result.transcript_consequences[0].allele_consequences[0]
            .consequences
            .contains(&Consequence::IntronVariant));
        // tr2: 8500bp away (>5000), so intergenic
        assert!(result.transcript_consequences[1].allele_consequences[0]
            .consequences
            .contains(&Consequence::IntergenicVariant));
    }

    #[test]
    fn test_most_severe_across_transcripts() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();

        // Splice donor is more severe than intron variant
        let pos = GenomicPosition::new("chr1", 1201, 1201, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );

        assert_eq!(result.most_severe, Some(Consequence::SpliceDonorVariant));
    }

    #[test]
    fn test_multi_allelic() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();

        // Two alt alleles at a coding position
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("T"), Allele::from_str("C")],
            &[&tr],
            None,
        );

        let tc = &result.transcript_consequences[0];
        assert_eq!(tc.allele_consequences.len(), 2);
    }
}
