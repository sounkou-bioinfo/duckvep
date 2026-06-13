//! duckvep — DuckDB-native variant effect prediction.
//!
//! See DESIGN.md. This crate registers genomics IO table functions (noodles)
//! and, in later phases, the VEP consequence / HGVS / ACMG UDFs.

mod io;

use duckdb::{duckdb_entrypoint_c_api, Connection, Result};
use std::error::Error;

use crate::io::vcf::ReadVcf;

#[duckdb_entrypoint_c_api()]
pub unsafe fn extension_entrypoint(con: Connection) -> Result<(), Box<dyn Error>> {
    con.register_table_function::<ReadVcf>("read_vcf")
        .expect("failed to register read_vcf");
    Ok(())
}
