//! `normalize_variant(pos, ref, alt) -> STRUCT(pos BIGINT, ref VARCHAR, alt VARCHAR)`
//!
//! Canonical minimal variant representation (right-trim shared suffix, then
//! left-trim shared prefix advancing `pos`) — the same form Ensembl VEP and
//! fastVEP emit (deletion `CT>C @P` -> `(P+1, "T", "")`, insertion `C>CA @P` ->
//! `(P+1, "", "A")`, SNV unchanged). Empty ref/alt is the minimal indel form
//! (map to `-` for VEP-style display). This is the load-bearing join key: any two
//! variant sources (VCF-anchored vs VEP-trimmed) compare correctly only once
//! normalized to this — needed for the concordance harness AND every future
//! supplementary-annotation join (ClinVar/gnomAD/dbSNP).

use duckdb::arrow::array::{Array, AsArray};
use duckdb::arrow::datatypes::Int64Type;
use duckdb::core::{DataChunkHandle, Inserter, LogicalTypeHandle, LogicalTypeId};
use duckdb::vscalar::{ScalarFunctionSignature, VScalar};
use duckdb::vtab::arrow::{data_chunk_to_arrow, WritableVector};
use std::error::Error;

fn varchar() -> LogicalTypeHandle {
    LogicalTypeHandle::from(LogicalTypeId::Varchar)
}
fn bigint() -> LogicalTypeHandle {
    LogicalTypeHandle::from(LogicalTypeId::Bigint)
}

/// Right-trim shared suffix, then left-trim shared prefix (advancing pos). Keeps
/// at least one base across the two alleles (so an SNV is untouched).
pub(crate) fn normalize(pos: i64, r: &str, a: &str) -> (i64, String, String) {
    let mut rb = r.as_bytes().to_vec();
    let mut ab = a.as_bytes().to_vec();
    let mut p = pos;
    while !rb.is_empty() && !ab.is_empty() && rb.last() == ab.last() && rb.len() + ab.len() > 2 {
        rb.pop();
        ab.pop();
    }
    let mut i = 0;
    while i < rb.len() && i < ab.len() && rb[i] == ab[i] && (rb.len() - i) + (ab.len() - i) > 2 {
        i += 1;
        p += 1;
    }
    (
        p,
        String::from_utf8_lossy(&rb[i..]).into_owned(),
        String::from_utf8_lossy(&ab[i..]).into_owned(),
    )
}

pub struct NormalizeVariant;

impl VScalar for NormalizeVariant {
    type State = ();

    unsafe fn invoke(
        _: &Self::State,
        input: &mut DataChunkHandle,
        output: &mut dyn WritableVector,
    ) -> Result<(), Box<dyn Error>> {
        let batch = data_chunk_to_arrow(input)?;
        let pos = batch.column(0).as_primitive::<Int64Type>();
        let refa = batch.column(1).as_string::<i32>();
        let alta = batch.column(2).as_string::<i32>();
        let n = input.len();

        let sv = output.struct_vector();
        let mut pv = sv.child(0, n);
        let pslice = pv.as_mut_slice::<i64>();
        let rv = sv.child(1, n);
        let av = sv.child(2, n);
        for i in 0..n {
            let (np, nr, na) = normalize(pos.value(i), refa.value(i), alta.value(i));
            pslice[i] = np;
            rv.insert(i, nr.as_str());
            av.insert(i, na.as_str());
        }
        Ok(())
    }

    fn signatures() -> Vec<ScalarFunctionSignature> {
        vec![ScalarFunctionSignature::exact(
            vec![bigint(), varchar(), varchar()],
            LogicalTypeHandle::struct_type(&[
                ("pos", bigint()),
                ("ref", varchar()),
                ("alt", varchar()),
            ]),
        )]
    }
}

#[cfg(test)]
mod tests {
    use super::normalize;
    #[test]
    fn vep_minimal_forms() {
        assert_eq!(normalize(100, "CT", "C"), (101, "T".into(), "".into())); // del
        assert_eq!(normalize(100, "C", "CA"), (101, "".into(), "A".into())); // ins
        assert_eq!(normalize(100, "A", "G"), (100, "A".into(), "G".into())); // snv
        assert_eq!(normalize(100, "AT", "GC"), (100, "AT".into(), "GC".into())); // mnv
        assert_eq!(normalize(100, "CTG", "CG"), (101, "T".into(), "".into())); // non-minimal del
    }
}
