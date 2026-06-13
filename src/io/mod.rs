//! Genomics IO table functions, backed by noodles.
//!
//! Each reader is a DuckDB `VTab` that materializes records into result
//! columns. Ported in spirit from the duckhts reader surface (read_vcf,
//! read_gff, read_fasta, …), but pure Rust via noodles.

pub mod gff;
pub mod vcf;
