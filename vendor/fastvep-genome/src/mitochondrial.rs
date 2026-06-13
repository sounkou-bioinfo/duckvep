//! Mitochondrial genome-specific handling.
//!
//! Handles circular coordinate wrapping and the mitochondrial codon table.

use crate::CodonTable;

/// Length of the human mitochondrial genome (rCRS reference).
pub const MT_LENGTH: u64 = 16569;

/// Returns true if the chromosome name indicates mitochondrial DNA.
pub fn is_mitochondrial(chrom: &str) -> bool {
    let c = chrom.to_lowercase();
    c == "mt" || c == "chrm" || c == "chrmt" || c == "m"
}

/// Wrap a position around the circular mitochondrial genome.
/// Positions > MT_LENGTH wrap to the beginning.
pub fn wrap_position(pos: u64) -> u64 {
    if pos == 0 {
        return 0;
    }
    ((pos - 1) % MT_LENGTH) + 1
}

/// The vertebrate mitochondrial codon table (NCBI translation table 2).
///
/// Differences from standard table:
/// - AGA, AGG = Stop (not Arg)
/// - ATA = Met (not Ile)
/// - TGA = Trp (not Stop)
pub fn mitochondrial_codon_table() -> CodonTable {
    CodonTable::from_ncbi_table(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_mitochondrial() {
        assert!(is_mitochondrial("MT"));
        assert!(is_mitochondrial("chrM"));
        assert!(is_mitochondrial("chrMT"));
        assert!(is_mitochondrial("M"));
        assert!(!is_mitochondrial("chr1"));
        assert!(!is_mitochondrial("chrX"));
    }

    #[test]
    fn test_wrap_position() {
        assert_eq!(wrap_position(1), 1);
        assert_eq!(wrap_position(16569), 16569);
        assert_eq!(wrap_position(16570), 1);
        assert_eq!(wrap_position(16571), 2);
        assert_eq!(wrap_position(0), 0);
    }

    #[test]
    fn test_mt_codon_table() {
        let table = mitochondrial_codon_table();
        // ATA = Met in MT (Ile in standard)
        assert_eq!(table.translate(b"ATA"), b'M');
        // TGA = Trp in MT (Stop in standard)
        assert_eq!(table.translate(b"TGA"), b'W');
        // AGA = Stop in MT (Arg in standard)
        assert_eq!(table.translate(b"AGA"), b'*');
        // Normal codons unchanged
        assert_eq!(table.translate(b"ATG"), b'M');
        assert_eq!(table.translate(b"TAA"), b'*');
    }
}
