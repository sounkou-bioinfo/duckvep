pub mod codon;
pub mod mitochondrial;
mod transcript;

pub use codon::CodonTable;
pub use transcript::{Exon, Gene, Transcript, Translation};
