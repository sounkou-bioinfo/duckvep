use fastvep_core::Strand;
use fastvep_genome::Transcript;

/// Splice site boundaries relative to the exon-intron junction.
///
/// Donor site (5' end of intron on forward strand):
///   exon ...XY | GT...  intron
///   splice_donor: positions 1-2 into intron (GT)
///   splice_donor_5th_base: position 5 into intron
///   splice_donor_region: positions 3-6 into intron
///   splice_region (exonic side): 3 bases into exon from boundary
///   splice_region (intronic side): 3-8 bases into intron
///
/// Acceptor site (3' end of intron on forward strand):
///   intron  ...AG | XY...  exon
///   splice_acceptor: last 2 bases of intron (AG)
///   splice_polypyrimidine: 3-17 bases from end of intron
///   splice_region (exonic side): 1-3 bases into exon
///   splice_region (intronic side): 3-8 bases into intron

/// Interval overlap: does `[a, b]` intersect `[lo, hi]`? All splice predicates
/// take the variant INTERVAL `[start, end]` (not a single position), so that
/// indels / MNVs / haplotypes spanning a splice boundary are classified correctly,
/// matching Ensembl's `overlap($r_start, $r_end, ...)`. For an SNV `start == end`.
fn ov(a: u64, b: u64, lo: u64, hi: u64) -> bool {
    a <= hi && b >= lo
}

/// Does the variant interval overlap a splice donor site (first 2 intronic bases at 5' of intron)?
pub fn is_splice_donor(transcript: &Transcript, start: u64, end: u64) -> bool {
    for_each_intron_boundary(
        transcript,
        |donor_start, donor_end, _acc_start, _acc_end| ov(start, end, donor_start, donor_end),
    )
}

/// Does the variant interval overlap a splice acceptor site (last 2 intronic bases at 3' of intron)?
pub fn is_splice_acceptor(transcript: &Transcript, start: u64, end: u64) -> bool {
    for_each_intron_boundary(
        transcript,
        |_donor_start, _donor_end, acc_start, acc_end| ov(start, end, acc_start, acc_end),
    )
}

/// Does the variant interval overlap the 5th base of the donor site?
pub fn is_splice_donor_5th_base(transcript: &Transcript, start: u64, end: u64) -> bool {
    for_each_intron_boundary_extended(transcript, |intron_start, intron_end, is_donor_at_start| {
        let p = if is_donor_at_start {
            intron_start + 4
        } else {
            intron_end - 4
        };
        ov(start, end, p, p)
    })
}

/// Does the variant interval overlap the splice donor region (positions 3-6 of intron)?
pub fn is_splice_donor_region(transcript: &Transcript, start: u64, end: u64) -> bool {
    for_each_intron_boundary_extended(transcript, |intron_start, intron_end, is_donor_at_start| {
        if is_donor_at_start {
            ov(start, end, intron_start + 2, intron_start + 5)
        } else {
            ov(start, end, intron_end - 5, intron_end - 2)
        }
    })
}

/// Check if position is in the splice polypyrimidine tract (3-17 bases from acceptor).
///
/// VEP defines this as positions `intron_end-16` to `intron_end-2` for the forward strand
/// (i.e., 3 to 17 bases from the 3' end of the intron), matching the Ensembl definition
/// of acceptor -3 to acceptor -17.
pub fn is_splice_polypyrimidine_tract(transcript: &Transcript, start: u64, end: u64) -> bool {
    // VEP swaps the affected interval to (min,max) before the polypyrimidine overlap
    // (BaseTranscriptVariationAllele::_intron_effects, first loop). For insertions
    // `normalized_interval` yields start=end+1, so without this swap an insertion at
    // the PPT edge (one base inside the window) is missed. For SNV/del start<=end so
    // the swap is a no-op.
    let (start, end) = (start.min(end), start.max(end));
    for_each_intron_boundary_extended(transcript, |intron_start, intron_end, is_donor_at_start| {
        // Polypyrimidine tract is near the acceptor end (Ensembl: intron_end-16 ..
        // intron_end-2 on the forward strand; mirrored on the reverse strand).
        if is_donor_at_start {
            let acc_region_start = if intron_end >= 16 {
                intron_end - 16
            } else {
                intron_start
            };
            ov(start, end, acc_region_start, intron_end.saturating_sub(2))
        } else {
            let acc_region_end = (intron_start + 16).min(intron_end);
            ov(start, end, intron_start + 2, acc_region_end)
        }
    })
}

/// Does the variant interval overlap a splice region (SO:0001630): within 1-3 bases
/// of the exon or 3-8 bases of the intron, at either junction. The bands are plain
/// genomic intervals (strand-independent for overlap), so this is correct for
/// indels/MNVs/haplotypes spanning the boundary.
pub fn is_splice_region(transcript: &Transcript, start: u64, end: u64) -> bool {
    // Introns lie in the gaps between position-sorted exons (strand-agnostic).
    // For insertions `normalized_interval` yields start=end+1; VEP's _intron_overlap
    // adds explicit boundary-touch cases for that zero-width interval.
    let insertion = start > end;
    let mut exons: Vec<_> = transcript.exons.iter().collect();
    exons.sort_by_key(|e| e.start);
    for w in exons.windows(2) {
        let intron_start = w[0].end + 1; // first intronic base (5' boundary side)
        let intron_end = w[1].start - 1; // last intronic base (3' boundary side)
        if intron_start > intron_end {
            continue;
        }
        // 3-8 bases into the intron from each boundary + 1-3 bases into each
        // flanking exon (the left exon ends at intron_start-1, the right exon
        // starts at intron_end+1).
        if ov(start, end, intron_start + 2, intron_start + 7)   // intronic, 5' band
            || ov(start, end, intron_end - 7, intron_end - 2)   // intronic, 3' band
            || ov(start, end, intron_start.saturating_sub(3), intron_start - 1) // left exon 1-3
            || ov(start, end, intron_end + 1, intron_end + 3)   // right exon 1-3
            // Insertion boundary-touch cases (VEP _intron_overlap insertion clause):
            // an insertion exactly at the exon/intron junction or the +2/-2 boundary.
            || (insertion
                && (start == intron_start
                    || end == intron_end
                    || start == intron_start + 2
                    || end == intron_end - 2))
        {
            return true;
        }
    }
    false
}

/// Helper: iterate intron boundaries and check a condition.
/// Calls `check(donor_start, donor_end, acceptor_start, acceptor_end)`.
fn for_each_intron_boundary<F>(transcript: &Transcript, check: F) -> bool
where
    F: Fn(u64, u64, u64, u64) -> bool,
{
    let sorted = sorted_exons(transcript);
    let n = sorted.len();
    if n < 2 {
        return false;
    }

    for i in 0..n - 1 {
        // Compute intron genomic coordinates based on strand
        // For forward: intron is between sorted[i].end and sorted[i+1].start
        // For reverse: sorted exons are in descending genomic order,
        //   so intron is between sorted[i+1].end and sorted[i].start
        let (intron_start, intron_end) = match transcript.strand {
            Strand::Forward => (sorted[i].end + 1, sorted[i + 1].start - 1),
            Strand::Reverse => (sorted[i + 1].end + 1, sorted[i].start - 1),
        };

        if intron_start > intron_end {
            continue;
        }

        // On forward strand: donor at intron_start, acceptor at intron_end
        // On reverse strand: donor at intron_end, acceptor at intron_start
        // Since we sort exons in transcript order, the first exon boundary
        // is always the donor side in transcript terms
        let (donor_start, donor_end, acc_start, acc_end) = match transcript.strand {
            Strand::Forward => (
                intron_start,
                (intron_start + 1).min(intron_end),
                if intron_end >= 1 {
                    intron_end - 1
                } else {
                    intron_start
                },
                intron_end,
            ),
            Strand::Reverse => (
                if intron_end >= 1 {
                    intron_end - 1
                } else {
                    intron_start
                },
                intron_end,
                intron_start,
                (intron_start + 1).min(intron_end),
            ),
        };

        if check(donor_start, donor_end, acc_start, acc_end) {
            return true;
        }
    }

    false
}

/// Extended intron boundary helper that provides full intron coords.
fn for_each_intron_boundary_extended<F>(transcript: &Transcript, check: F) -> bool
where
    F: Fn(u64, u64, bool) -> bool,
{
    let sorted = sorted_exons(transcript);
    let n = sorted.len();
    if n < 2 {
        return false;
    }

    for i in 0..n - 1 {
        let (intron_start, intron_end) = match transcript.strand {
            Strand::Forward => (sorted[i].end + 1, sorted[i + 1].start - 1),
            Strand::Reverse => (sorted[i + 1].end + 1, sorted[i].start - 1),
        };

        if intron_start > intron_end {
            continue;
        }

        // is_donor_at_start: true for forward strand (donor=5' end=start of intron)
        let is_donor_at_start = transcript.strand == Strand::Forward;

        if check(intron_start, intron_end, is_donor_at_start) {
            return true;
        }
    }

    false
}

fn sorted_exons(transcript: &Transcript) -> Vec<&fastvep_genome::Exon> {
    let mut exons: Vec<&fastvep_genome::Exon> = transcript.exons.iter().collect();
    match transcript.strand {
        Strand::Forward => exons.sort_by_key(|e| e.start),
        Strand::Reverse => exons.sort_by(|a, b| b.start.cmp(&a.start)),
    }
    exons
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_genome::{Exon, Gene, Transcript, Translation};

    fn make_forward_transcript() -> Transcript {
        // Exon1: 1000-1200, Intron1: 1201-1999, Exon2: 2000-2300
        Transcript {
            stable_id: "ENST_TEST".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_TEST".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 2300,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 2300,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 2000,
                    end: 2300,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "P1".into(),
                genomic_start: 1000,
                genomic_end: 2300,
                start_exon_rank: 1,
                start_exon_offset: 0,
                end_exon_rank: 2,
                end_exon_offset: 300,
            }),
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(502),
            coding_region_start: Some(1000),
            coding_region_end: Some(2300),
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
    fn test_splice_donor() {
        let tr = make_forward_transcript();
        // Intron: 1201-1999. Donor = first 2 bases: 1201, 1202
        assert!(is_splice_donor(&tr, 1201, 1201));
        assert!(is_splice_donor(&tr, 1202, 1202));
        assert!(!is_splice_donor(&tr, 1203, 1203));
        assert!(!is_splice_donor(&tr, 1200, 1200)); // exonic
    }

    #[test]
    fn test_splice_acceptor() {
        let tr = make_forward_transcript();
        // Intron: 1201-1999. Acceptor = last 2 bases: 1998, 1999
        assert!(is_splice_acceptor(&tr, 1998, 1998));
        assert!(is_splice_acceptor(&tr, 1999, 1999));
        assert!(!is_splice_acceptor(&tr, 1997, 1997));
        assert!(!is_splice_acceptor(&tr, 2000, 2000)); // exonic
    }

    #[test]
    fn test_splice_region() {
        let tr = make_forward_transcript();
        // Exonic splice region: last 3 bases of exon1 (1198, 1199, 1200)
        assert!(is_splice_region(&tr, 1198, 1198));
        assert!(is_splice_region(&tr, 1200, 1200));
        // Exonic splice region: first 3 bases of exon2 (2000, 2001, 2002)
        assert!(is_splice_region(&tr, 2000, 2000));
        assert!(is_splice_region(&tr, 2002, 2002));
        // Intronic splice region: 3-8 bases from donor (1203-1208)
        assert!(is_splice_region(&tr, 1203, 1203));
        assert!(is_splice_region(&tr, 1208, 1208));
        assert!(!is_splice_region(&tr, 1209, 1209));
        // Mid-intron: not splice region
        assert!(!is_splice_region(&tr, 1500, 1500));
    }

    #[test]
    fn test_polypyrimidine_forward() {
        let tr = make_forward_transcript();
        // Intron: 1201-1999. Acceptor at end (1999).
        // Polypyrimidine tract: 3-17 bases from acceptor = positions 1982-1996
        // (intron_end - 17 = 1982, intron_end - 3 = 1996)
        // VEP definition: intron_end-16 to intron_end-2 (1983-1997)

        // Check boundaries
        for pos in 1980..=2000 {
            let in_ppt = is_splice_polypyrimidine_tract(&tr, pos, pos);
            eprintln!(
                "  pos {} (dist from end = {}): ppt={}",
                pos,
                1999u64.saturating_sub(pos),
                in_ppt
            );
        }

        // Distance 17 from intron_end(1999) = 1982 = intron_end - 17
        // Distance 3 from intron_end(1999) = 1996 = intron_end - 3
        // But VEP measures from exon boundary (2000):
        //   dist 17 from exon = 2000-17 = 1983 = intron_end - 16
        //   dist 3 from exon = 2000-3 = 1997 = intron_end - 2
        assert!(
            is_splice_polypyrimidine_tract(&tr, 1983, 1983),
            "pos 1983 (dist 17 from exon) should be PPT"
        );
        assert!(
            is_splice_polypyrimidine_tract(&tr, 1997, 1997),
            "pos 1997 (dist 3 from exon) should be PPT"
        );
        assert!(
            !is_splice_polypyrimidine_tract(&tr, 1982, 1982),
            "pos 1982 (dist 18 from exon) should NOT be PPT"
        );
        assert!(
            !is_splice_polypyrimidine_tract(&tr, 1998, 1998),
            "pos 1998 (dist 2 from exon) should NOT be PPT - it's the acceptor site"
        );
    }

    fn make_reverse_transcript() -> Transcript {
        // Reverse strand: exons sorted in descending genomic order
        // Exon1 (rank 1, 5'): 2000-2300 (higher coords)
        // Exon2 (rank 2, 3'): 1000-1200 (lower coords)
        // Intron: 1201-1999
        // For reverse: donor at intron_end (1999), acceptor at intron_start (1201)
        Transcript {
            stable_id: "ENST_REV".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_REV".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 2300,
                strand: Strand::Reverse,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 2300,
            strand: Strand::Reverse,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 2000,
                    end: 2300,
                    strand: Strand::Reverse,
                    phase: 0,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Reverse,
                    phase: 0,
                    end_phase: 0,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "P1".into(),
                genomic_start: 1000,
                genomic_end: 2300,
                start_exon_rank: 1,
                start_exon_offset: 0,
                end_exon_rank: 2,
                end_exon_offset: 200,
            }),
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(502),
            coding_region_start: Some(1000),
            coding_region_end: Some(2300),
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
    fn test_polypyrimidine_reverse() {
        let tr = make_reverse_transcript();
        // Intron: 1201-1999. On reverse strand, acceptor at intron_start (1201).
        // Polypyrimidine tract: 3-17 bases from acceptor
        // VEP: distance measured from exon boundary (1200)
        //   dist 3: 1200+3 = 1203 = intron_start + 2
        //   dist 17: 1200+17 = 1217 = intron_start + 16

        for pos in 1199..=1220 {
            let in_ppt = is_splice_polypyrimidine_tract(&tr, pos, pos);
            eprintln!(
                "  REV pos {} (dist from intron_start=1201: {}): ppt={}",
                pos,
                pos as i64 - 1201,
                in_ppt
            );
        }

        // c.X-17 = 17 bases from exon boundary = 1200+17 = 1217 = intron_start + 16
        assert!(
            is_splice_polypyrimidine_tract(&tr, 1217, 1217),
            "pos 1217 (c.X-17, dist 16 from intron_start) should be PPT"
        );
        // c.X-3 = 3 bases from exon boundary = 1200+3 = 1203 = intron_start + 2
        assert!(
            is_splice_polypyrimidine_tract(&tr, 1203, 1203),
            "pos 1203 (c.X-3, dist 2 from intron_start) should be PPT"
        );
        // c.X-18 = 18 bases = 1218 = intron_start + 17
        assert!(
            !is_splice_polypyrimidine_tract(&tr, 1218, 1218),
            "pos 1218 (c.X-18) should NOT be PPT"
        );
        // c.X-2 = 2 bases = 1202 = intron_start + 1 (acceptor site)
        assert!(
            !is_splice_polypyrimidine_tract(&tr, 1202, 1202),
            "pos 1202 (c.X-2, acceptor site) should NOT be PPT"
        );
    }
}
