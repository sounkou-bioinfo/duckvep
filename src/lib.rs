//! duckvep — DuckDB-native variant effect prediction.
//!
//! See docs/DESIGN.md. This crate registers genomics IO table functions (noodles)
//! and, in later phases, the VEP consequence / HGVS / ACMG UDFs.

mod io;
mod vec_util;
mod vep;

use duckdb::{duckdb_entrypoint_c_api, Connection, Result};
use std::error::Error;

use crate::io::gff::ReadGffTranscripts;
use crate::io::vcf::{ReadVcf, VcfSamples};
use crate::vep::annotate::VepAnnotate;
use crate::vep::consequence::{VepConsequence, VepLoadCache};
use crate::vep::normalize::NormalizeVariant;

#[duckdb_entrypoint_c_api()]
pub unsafe fn extension_entrypoint(con: Connection) -> Result<(), Box<dyn Error>> {
    con.register_table_function::<ReadVcf>("read_vcf")
        .expect("failed to register read_vcf");
    con.register_table_function::<VcfSamples>("vcf_samples")
        .expect("failed to register vcf_samples");
    con.register_table_function::<VepAnnotate>("vep_annotate")
        .expect("failed to register vep_annotate");
    con.register_table_function::<ReadGffTranscripts>("read_gff_transcripts")
        .expect("failed to register read_gff_transcripts");
    con.register_scalar_function::<VepLoadCache>("vep_load_cache")
        .expect("failed to register vep_load_cache");
    con.register_scalar_function::<VepConsequence>("vep_consequence")
        .expect("failed to register vep_consequence");
    con.register_scalar_function::<NormalizeVariant>("normalize_variant")
        .expect("failed to register normalize_variant");
    Ok(())
}
