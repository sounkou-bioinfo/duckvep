//! Structural variant consequence prediction.
//!
//! Predicts consequences for SVs (deletions, duplications, inversions,
//! CNVs, breakends) based on their overlap with transcripts.

use fastvep_core::{Allele, Consequence, Impact, Strand, VariantType};
use fastvep_genome::Transcript;
use std::sync::Arc;

use crate::predictor::{AlleleConsequenceResult, TranscriptConsequence};

/// Predict consequences of a structural variant against a set of transcripts.
///
/// The SV is defined by its genomic coordinates [start, end] and variant type.
/// For each transcript that overlaps, determines whether the SV causes
/// ablation, amplification, or partial disruption.
pub fn predict_sv_consequences(
    chrom: &str,
    start: u64,
    end: u64,
    variant_type: VariantType,
    alt_alleles: &[Allele],
    transcripts: &[&Transcript],
    upstream_distance: u64,
    downstream_distance: u64,
) -> Vec<TranscriptConsequence> {
    let mut results = Vec::new();

    for &transcript in transcripts {
        if *transcript.chromosome != *chrom {
            continue;
        }

        let allele_consequences: Vec<AlleleConsequenceResult> = alt_alleles
            .iter()
            .map(|allele| {
                predict_sv_for_transcript(
                    start,
                    end,
                    variant_type,
                    allele,
                    transcript,
                    upstream_distance,
                    downstream_distance,
                )
            })
            .collect();

        if !allele_consequences.is_empty() {
            results.push(TranscriptConsequence {
                transcript_id: Arc::clone(&transcript.stable_id),
                gene_id: transcript.gene.stable_id.clone(),
                gene_symbol: transcript.gene.symbol.clone(),
                biotype: Arc::clone(&transcript.biotype),
                allele_consequences,
                canonical: transcript.canonical,
                strand: transcript.strand,
            });
        }
    }

    results
}

fn predict_sv_for_transcript(
    sv_start: u64,
    sv_end: u64,
    variant_type: VariantType,
    allele: &Allele,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
) -> AlleleConsequenceResult {
    let tx_start = transcript.start;
    let tx_end = transcript.end;

    // Check if SV overlaps the transcript at all
    let overlaps = sv_start <= tx_end && sv_end >= tx_start;

    if !overlaps {
        // Check upstream/downstream
        let (up_dist, down_dist) = if transcript.strand == Strand::Forward {
            (upstream_distance, downstream_distance)
        } else {
            (downstream_distance, upstream_distance)
        };

        let distance = if sv_end < tx_start {
            // SV is before transcript
            if transcript.strand == Strand::Forward {
                Some(tx_start as i64 - sv_end as i64) // upstream
            } else {
                Some(tx_start as i64 - sv_end as i64) // downstream on reverse
            }
        } else if sv_start > tx_end {
            if transcript.strand == Strand::Forward {
                Some(sv_start as i64 - tx_end as i64) // downstream
            } else {
                Some(sv_start as i64 - tx_end as i64) // upstream on reverse
            }
        } else {
            None
        };

        let consequence = if let Some(d) = distance {
            let abs_d = d.unsigned_abs();
            if sv_end < tx_start {
                if (transcript.strand == Strand::Forward && abs_d <= up_dist)
                    || (transcript.strand == Strand::Reverse && abs_d <= down_dist)
                {
                    if transcript.strand == Strand::Forward {
                        Consequence::UpstreamGeneVariant
                    } else {
                        Consequence::DownstreamGeneVariant
                    }
                } else {
                    Consequence::IntergenicVariant
                }
            } else if (transcript.strand == Strand::Forward && abs_d <= down_dist)
                || (transcript.strand == Strand::Reverse && abs_d <= up_dist)
            {
                if transcript.strand == Strand::Forward {
                    Consequence::DownstreamGeneVariant
                } else {
                    Consequence::UpstreamGeneVariant
                }
            } else {
                Consequence::IntergenicVariant
            }
        } else {
            Consequence::IntergenicVariant
        };

        return AlleleConsequenceResult {
            allele: allele.clone(),
            consequences: vec![consequence],
            impact: consequence.impact(),
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

    // SV overlaps transcript — determine consequences
    let completely_contains = sv_start <= tx_start && sv_end >= tx_end;
    let consequences = determine_sv_overlap_consequences(
        sv_start,
        sv_end,
        variant_type,
        transcript,
        completely_contains,
    );

    let impact = consequences
        .iter()
        .map(|c| c.impact())
        .min()
        .unwrap_or(Impact::Modifier);

    AlleleConsequenceResult {
        allele: allele.clone(),
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
        distance: None,
    }
}

fn determine_sv_overlap_consequences(
    sv_start: u64,
    sv_end: u64,
    variant_type: VariantType,
    transcript: &Transcript,
    completely_contains: bool,
) -> Vec<Consequence> {
    let mut consequences = Vec::new();

    match variant_type {
        VariantType::CopyNumberLoss | VariantType::Deletion => {
            if completely_contains {
                consequences.push(Consequence::TranscriptAblation);
            } else {
                // Partial overlap: Ensembl emits feature_truncation PLUS the small-variant
                // interval terms over the SV span (coding_sequence / intron) — NOT a splice
                // term (a multi-kb deletion boundary lands in an intron, not the 2 bp splice
                // site). Reuse the interval predicates rather than the splice heuristic.
                consequences.push(Consequence::FeatureTruncation);
                if transcript.is_coding() && hits_coding_region(sv_start, sv_end, transcript) {
                    consequences.push(Consequence::CodingSequenceVariant);
                }
                if hits_intron(sv_start, sv_end, transcript) {
                    consequences.push(Consequence::IntronVariant);
                }
                // `copy_number_decrease` is a COPY-NUMBER allele term (<CN0> / copy_number_loss),
                // NOT a plain <DEL> — Ensembl does not attach it to a sequence deletion
                // (the `<DEL> != <CN0>` distinction).
                if variant_type == VariantType::CopyNumberLoss {
                    consequences.push(Consequence::CopyNumberDecrease);
                }
            }
        }
        VariantType::TandemDuplication | VariantType::CopyNumberGain => {
            if completely_contains {
                consequences.push(Consequence::TranscriptAmplification);
            } else {
                consequences.push(Consequence::FeatureElongation);
                consequences.push(Consequence::CopyNumberIncrease);
            }
        }
        VariantType::CopyNumberVariation => {
            if completely_contains {
                // Could be either ablation or amplification
                consequences.push(Consequence::TranscriptVariant);
            }
            consequences.push(Consequence::CopyNumberChange);
        }
        VariantType::Inversion => {
            if completely_contains {
                consequences.push(Consequence::TranscriptVariant);
            } else {
                // An inversion does NOT delete/truncate the feature — Ensembl gives it the
                // small-variant interval terms over its span (coding_sequence / intron / UTR),
                // not feature_truncation. Reuse the interval predicates.
                if transcript.is_coding() && hits_coding_region(sv_start, sv_end, transcript) {
                    consequences.push(Consequence::CodingSequenceVariant);
                }
                if hits_intron(sv_start, sv_end, transcript) {
                    consequences.push(Consequence::IntronVariant);
                }
                if hits_utr(sv_start, sv_end, transcript, true) {
                    consequences.push(Consequence::FivePrimeUtrVariant);
                }
                if hits_utr(sv_start, sv_end, transcript, false) {
                    consequences.push(Consequence::ThreePrimeUtrVariant);
                }
            }
        }
        VariantType::TranslocationBreakend => {
            // Breakend landing in a transcript
            consequences.push(Consequence::TranscriptVariant);
            if transcript.is_coding() && hits_coding_region(sv_start, sv_end, transcript) {
                consequences.push(Consequence::CodingSequenceVariant);
            }
        }
        VariantType::ShortTandemRepeatVariation => {
            consequences.push(Consequence::ShortTandemRepeatChange);
        }
        _ => {
            consequences.push(Consequence::TranscriptVariant);
        }
    }

    if consequences.is_empty() {
        consequences.push(Consequence::TranscriptVariant);
    }

    consequences
}

/// Check if the SV interval overlaps the CODING portion of any exon (the exon ∩ CDS span) —
/// VEP uses mapped CDS coordinates, so an intron-only SV that merely falls WITHIN the genomic
/// CDS bounds (between two coding exons) must NOT count as coding_sequence_variant.
fn hits_coding_region(sv_start: u64, sv_end: u64, transcript: &Transcript) -> bool {
    let (Some(cds_start), Some(cds_end)) =
        (transcript.coding_region_start, transcript.coding_region_end)
    else {
        return false;
    };
    let coding_lo = cds_start.min(cds_end);
    let coding_hi = cds_start.max(cds_end);
    transcript.exons.iter().any(|e| {
        let ex_cds_lo = e.start.max(coding_lo); // coding portion of this exon
        let ex_cds_hi = e.end.min(coding_hi);
        ex_cds_lo <= ex_cds_hi && sv_start <= ex_cds_hi && sv_end >= ex_cds_lo
    })
}

/// Check if the SV interval overlaps any intron (the gap between two consecutive exons).
fn hits_intron(sv_start: u64, sv_end: u64, transcript: &Transcript) -> bool {
    let mut exons: Vec<(u64, u64)> = transcript.exons.iter().map(|e| (e.start, e.end)).collect();
    exons.sort_unstable();
    for w in exons.windows(2) {
        let intron_lo = w[0].1 + 1; // one base past the previous exon end
        let intron_hi = w[1].0.saturating_sub(1); // one base before the next exon start
        if intron_lo <= intron_hi && sv_start <= intron_hi && sv_end >= intron_lo {
            return true;
        }
    }
    false
}

/// Check if the SV interval overlaps the 5' (`five_prime=true`) or 3' UTR genomic region —
/// the transcript span outside the coding region, strand-aware (5'UTR is upstream of the CDS
/// in transcript orientation: low-coordinate side on +, high-coordinate side on -).
fn hits_utr(sv_start: u64, sv_end: u64, transcript: &Transcript, five_prime: bool) -> bool {
    let (Some(crs), Some(cre)) =
        (transcript.coding_region_start, transcript.coding_region_end)
    else {
        return false;
    };
    let cds_lo = crs.min(cre);
    let cds_hi = crs.max(cre);
    // The UTR on the low-coordinate side is 5' on +strand, 3' on -strand (and vice-versa).
    let low_is_five_prime = transcript.strand == Strand::Forward;
    let (utr_lo, utr_hi) = if five_prime == low_is_five_prime {
        (transcript.start, cds_lo.saturating_sub(1)) // low-coordinate UTR
    } else {
        (cds_hi + 1, transcript.end) // high-coordinate UTR
    };
    if utr_lo > utr_hi {
        return false;
    }
    // The UTR is EXONIC — require overlap with the UTR portion of an exon, not just the genomic
    // UTR-side span, so an SV wholly in a UTR-side intron does not over-call the UTR term
    // (VEP's UTR predicates require cDNA / exon overlap).
    transcript.exons.iter().any(|e| {
        let ex_utr_lo = e.start.max(utr_lo);
        let ex_utr_hi = e.end.min(utr_hi);
        ex_utr_lo <= ex_utr_hi && sv_start <= ex_utr_hi && sv_end >= ex_utr_lo
    })
}

/// Check if the SV overlaps a splice site (2bp at exon boundaries).
#[allow(dead_code)] // retained for precise-breakpoint SV splice handling (interval-predicate WIP)
fn hits_splice_site(sv_start: u64, sv_end: u64, transcript: &Transcript) -> bool {
    for exon in &transcript.exons {
        // Donor site: 2bp after exon end (exon.end+1, exon.end+2)
        if sv_start <= exon.end + 2 && sv_end > exon.end {
            return true;
        }
        // Acceptor site: 2bp before exon start (exon.start-2, exon.start-1)
        if exon.start >= 3 && sv_start < exon.start && sv_end >= exon.start - 2 {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_genome::{Exon, Gene, Translation};

    fn make_coding_transcript(start: u64, end: u64) -> Transcript {
        Transcript {
            stable_id: Arc::from("ENST_TEST"),
            version: None,
            gene: Gene {
                stable_id: "ENSG_TEST".into(),
                symbol: Some("TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start,
                end,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start,
            end,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start,
                    end: start + 200,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: start + 500,
                    end: start + 700,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "P_TEST".into(),
                genomic_start: start + 50,
                genomic_end: start + 650,
                start_exon_rank: 1,
                start_exon_offset: 50,
                end_exon_rank: 2,
                end_exon_offset: 150,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(401),
            coding_region_start: Some(start + 50),
            coding_region_end: Some(start + 650),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
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
    fn test_deletion_completely_contains_transcript() {
        let tx = make_coding_transcript(5000, 6000);
        let results = predict_sv_consequences(
            "chr1",
            4000,
            7000,
            VariantType::CopyNumberLoss,
            &[Allele::Symbolic("<DEL>".into())],
            &[&tx],
            5000,
            5000,
        );

        assert_eq!(results.len(), 1);
        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::TranscriptAblation));
    }

    #[test]
    fn test_deletion_partial_overlap() {
        let tx = make_coding_transcript(5000, 6000);
        let results = predict_sv_consequences(
            "chr1",
            5300,
            5800,
            VariantType::CopyNumberLoss,
            &[Allele::Symbolic("<DEL>".into())],
            &[&tx],
            5000,
            5000,
        );

        assert_eq!(results.len(), 1);
        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::FeatureTruncation));
        assert!(cons.contains(&Consequence::CodingSequenceVariant));
    }

    // A partial SYMBOLIC <DEL> (VariantType::Deletion) matches Ensembl: feature_truncation +
    // the small-variant interval terms (coding_sequence + intron over the SV span), and NOT a
    // copy_number_decrease (that is a <CN0> / copy_number_loss term — the <DEL> != <CN0>
    // distinction) nor a splice_acceptor (the kb-scale boundary lands in the intron).
    #[test]
    fn test_partial_deletion_is_not_copy_number() {
        let tx = make_coding_transcript(5000, 6000); // exons 5000-5200, 5500-5700 (intron 5201-5499)
        let results = predict_sv_consequences(
            "chr1",
            5300,
            5800,
            VariantType::Deletion,
            &[Allele::Symbolic("<DEL>".into())],
            &[&tx],
            5000,
            5000,
        );
        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::FeatureTruncation));
        assert!(cons.contains(&Consequence::CodingSequenceVariant));
        assert!(cons.contains(&Consequence::IntronVariant), "SV spans the intron");
        assert!(!cons.contains(&Consequence::CopyNumberDecrease), "<DEL> is not <CN0>");
        assert!(!cons.contains(&Consequence::SpliceAcceptorVariant), "kb boundary is intronic");
    }

    // A partial <INV> reuses the interval predicates over its span (coding_sequence + intron +
    // 5'UTR here) and does NOT emit feature_truncation (an inversion does not truncate the
    // feature) — matching Ensembl. tx 5000-6000: 5'UTR 5000-5049, CDS 5050-5650, intron
    // 5201-5499, 3'UTR 5651-6000; the SV 5000-5550 hits the first three only.
    #[test]
    fn test_partial_inversion_interval_predicates() {
        let tx = make_coding_transcript(5000, 6000);
        let results = predict_sv_consequences(
            "chr1",
            5000,
            5550,
            VariantType::Inversion,
            &[Allele::Symbolic("<INV>".into())],
            &[&tx],
            5000,
            5000,
        );
        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::CodingSequenceVariant));
        assert!(cons.contains(&Consequence::IntronVariant));
        assert!(cons.contains(&Consequence::FivePrimeUtrVariant));
        assert!(!cons.contains(&Consequence::ThreePrimeUtrVariant), "SV does not reach the 3'UTR");
        assert!(!cons.contains(&Consequence::FeatureTruncation), "an inversion does not truncate");
    }

    // An SV entirely within an INTRON but inside the genomic CDS bounds must NOT be
    // coding_sequence_variant (VEP uses mapped/exonic CDS coords). tx 5000-6000: intron
    // 5201-5499 sits between coding exons; an SV 5250-5450 is intron-only.
    #[test]
    fn test_intron_only_sv_inside_cds_is_not_coding() {
        let tx = make_coding_transcript(5000, 6000);
        let results = predict_sv_consequences(
            "chr1",
            5250,
            5450,
            VariantType::Deletion,
            &[Allele::Symbolic("<DEL>".into())],
            &[&tx],
            5000,
            5000,
        );
        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::IntronVariant));
        assert!(
            !cons.contains(&Consequence::CodingSequenceVariant),
            "intron-only SV inside CDS bounds is not coding"
        );
    }

    // hits_utr is exon-gated: an SV on the UTR-side of the CDS but NOT overlapping a UTR exon
    // must not over-call the UTR term. tx 5000-6000: 3'UTR span 5651-6000, but exon2 ends at
    // 5700, so 5750-5900 is a non-exonic UTR-side region.
    #[test]
    fn test_non_exonic_utr_side_sv_is_not_utr() {
        let tx = make_coding_transcript(5000, 6000);
        let results = predict_sv_consequences(
            "chr1",
            5750,
            5900,
            VariantType::Inversion,
            &[Allele::Symbolic("<INV>".into())],
            &[&tx],
            5000,
            5000,
        );
        let cons = &results[0].allele_consequences[0].consequences;
        assert!(
            !cons.contains(&Consequence::ThreePrimeUtrVariant),
            "SV in a non-exonic UTR-side region is not a UTR variant"
        );
    }

    #[test]
    fn test_duplication_completely_contains() {
        let tx = make_coding_transcript(5000, 6000);
        let results = predict_sv_consequences(
            "chr1",
            4000,
            7000,
            VariantType::TandemDuplication,
            &[Allele::Symbolic("<DUP>".into())],
            &[&tx],
            5000,
            5000,
        );

        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::TranscriptAmplification));
    }

    #[test]
    fn test_sv_no_overlap_upstream() {
        let tx = make_coding_transcript(10000, 11000);
        let results = predict_sv_consequences(
            "chr1",
            5000,
            6000,
            VariantType::CopyNumberLoss,
            &[Allele::Symbolic("<DEL>".into())],
            &[&tx],
            5000,
            5000,
        );

        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::UpstreamGeneVariant));
    }
}
