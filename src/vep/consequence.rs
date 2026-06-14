//! Scan-driven scalar UDFs — the idiomatic, performant path (docs/DESIGN.md §3.2).
//!
//! `vep_load_cache(gff3, fasta)` builds the engine once into a global (set once,
//! read by all worker threads). `vep_consequence(chrom,pos,ref,alt)` is a pure
//! scalar over rows DuckDB scans (`read_vcf`, `read_parquet`, any source) — DuckDB
//! owns the scan; the kernel only computes per variant.
//!
//! Output is a native `LIST<STRUCT>` written **directly to DuckDB vectors** (the
//! same core API the VTabs use) — no Arrow on the output path. Arrow is used only
//! to *read* the input batch (`data_chunk_to_arrow`), the safe non-nested
//! direction; the Arrow *output* conversion mangles nested structs for multi-row
//! batches, which is why we avoid it.

use crate::vep::engine::{build_context, AnnotatedRow, EngineContext, DEFAULT_DISTANCE};
use arc_swap::ArcSwapOption;
use duckdb::arrow::array::{Array, AsArray};
use duckdb::arrow::datatypes::Int64Type;
use duckdb::core::{DataChunkHandle, Inserter, LogicalTypeHandle, LogicalTypeId};
use duckdb::vscalar::{ScalarFunctionSignature, VScalar};
use duckdb::vtab::arrow::{data_chunk_to_arrow, WritableVector};
use std::error::Error;
use std::sync::{Arc, OnceLock};

/// Load-once engine, shared across DuckDB worker threads. `ArcSwapOption` gives
/// **lock-free** reads (an atomic pointer load + Arc clone) on the hot scalar
/// path — no reader-counter contention, no guard held across the compute — while
/// `vep_load_cache` can still atomically swap in a new engine.
static CACHE: OnceLock<ArcSwapOption<EngineContext>> = OnceLock::new();
fn cache() -> &'static ArcSwapOption<EngineContext> {
    CACHE.get_or_init(ArcSwapOption::empty)
}

fn varchar() -> LogicalTypeHandle {
    LogicalTypeHandle::from(LogicalTypeId::Varchar)
}

/// `vep_load_cache(gff3, fasta) -> VARCHAR` — build + store the engine.
pub struct VepLoadCache;

impl VScalar for VepLoadCache {
    type State = ();

    unsafe fn invoke(
        _: &Self::State,
        input: &mut DataChunkHandle,
        output: &mut dyn WritableVector,
    ) -> Result<(), Box<dyn Error>> {
        let batch = data_chunk_to_arrow(input)?;
        let gff3 = batch.column(0).as_string::<i32>();
        let fasta = batch.column(1).as_string::<i32>();
        let fastav = if fasta.is_null(0) || fasta.value(0).is_empty() {
            None
        } else {
            Some(fasta.value(0))
        };
        let ctx = build_context(gff3.value(0), fastav, DEFAULT_DISTANCE)?;
        cache().store(Some(Arc::new(ctx)));

        let v = output.flat_vector();
        for i in 0..input.len() {
            v.insert(i, "loaded");
        }
        Ok(())
    }

    fn signatures() -> Vec<ScalarFunctionSignature> {
        vec![ScalarFunctionSignature::exact(
            vec![varchar(), varchar()],
            varchar(),
        )]
    }
}

/// The per-transcript struct type returned inside the list.
fn consequence_struct() -> LogicalTypeHandle {
    LogicalTypeHandle::struct_type(&[
        ("transcript_id", varchar()),
        ("gene_id", varchar()),
        ("gene_symbol", varchar()),
        ("biotype", varchar()),
        ("canonical", LogicalTypeHandle::from(LogicalTypeId::Boolean)),
        // SO terms as a proper LIST<VARCHAR> (one variant can hit several, e.g.
        // splice_region_variant & intron_variant). This is a 3-level nesting
        // (LIST<STRUCT<…LIST<VARCHAR>…>>); see `invoke` for how the inner list's
        // entry array is written past duckdb-rs's per-vector wrapper bound.
        ("consequence", LogicalTypeHandle::list(&varchar())),
        ("impact", varchar()),
        ("amino_acids", varchar()),
        ("codons", varchar()),
        (
            "protein_pos",
            LogicalTypeHandle::from(LogicalTypeId::Bigint),
        ),
        ("alt", varchar()),
        ("hgvsc", varchar()),
        ("hgvsp", varchar()),
        ("hgvsg", varchar()),
    ])
}

/// `vep_consequence(chrom, pos, ref, alt) -> LIST<STRUCT>` — per-variant kernel.
pub struct VepConsequence;

impl VScalar for VepConsequence {
    type State = ();

    unsafe fn invoke(
        _: &Self::State,
        input: &mut DataChunkHandle,
        output: &mut dyn WritableVector,
    ) -> Result<(), Box<dyn Error>> {
        let batch = data_chunk_to_arrow(input)?;
        // Two arities: (chrom,pos,ref,alt) and the END-aware (chrom,pos,end,ref,alt)
        // used to span structural variants (symbolic <DEL>/<CNV>/… and breakends).
        let spanned = batch.num_columns() == 5;
        let chrom = batch.column(0).as_string::<i32>();
        let pos = batch.column(1).as_primitive::<Int64Type>();
        let end = spanned.then(|| batch.column(2).as_primitive::<Int64Type>());
        let refa = batch.column(if spanned { 3 } else { 2 }).as_string::<i32>();
        let alta = batch.column(if spanned { 4 } else { 3 }).as_string::<i32>();
        let nrows = input.len();

        // Compute all rows, holding the cache only for the compute; results are
        // owned, so we release the lock before writing vectors.
        let mut flat: Vec<AnnotatedRow> = Vec::new();
        let mut list_offsets: Vec<usize> = Vec::with_capacity(nrows + 1);
        list_offsets.push(0);
        {
            // Lock-free snapshot: atomic load + Arc clone, released for reload.
            let ctx = cache()
                .load_full()
                .ok_or("vep_consequence: call vep_load_cache(gff3, fasta) first")?;
            for i in 0..nrows {
                let p = pos.value(i) as u64;
                let e = end.map(|c| c.value(i) as u64).unwrap_or(p);
                flat.extend(ctx.annotate_variant_spanned(
                    chrom.value(i),
                    p,
                    e,
                    refa.value(i),
                    alta.value(i),
                ));
                list_offsets.push(flat.len());
            }
        }
        let total = flat.len();

        // Flatten the nested consequence (LIST<VARCHAR>) child: one inner list per
        // struct row, `so_offsets[j]..so_offsets[j+1]` slicing the shared values.
        let mut so_flat: Vec<&str> = Vec::new();
        let mut so_offsets: Vec<usize> = Vec::with_capacity(total + 1);
        so_offsets.push(0);
        for r in &flat {
            so_flat.extend(r.consequence.iter().map(String::as_str));
            so_offsets.push(so_flat.len());
        }

        let mut list = output.list_vector();
        {
            let sv = list.struct_child(total);
            macro_rules! vcol {
                ($i:expr, $f:ident) => {{
                    let v = sv.child($i, total);
                    for (j, r) in flat.iter().enumerate() {
                        v.insert(j, &*r.$f);
                    }
                }};
            }
            vcol!(0, transcript_id);
            vcol!(1, gene_id);
            vcol!(2, gene_symbol);
            vcol!(3, biotype);
            {
                let mut v = sv.child(4, total);
                let s = v.as_mut_slice::<bool>();
                for (j, r) in flat.iter().enumerate() {
                    s[j] = r.canonical;
                }
            }
            {
                // consequence: a LIST<VARCHAR> child of the struct. Write the
                // inner VARCHAR values (reserves the inner child to so_flat.len),
                // set its size, then fill the per-row entry array. We must NOT use
                // `ListVector::set_entry` here: its slice is bounded by the wrapper's
                // hardcoded vector size (2048), but `total` can exceed that. Instead
                // we view the inner list vector's own data — the contiguous
                // `duckdb_list_entry` (offset,length) array — as a FlatVector sized
                // to `total` (the struct's reserved capacity) and write it directly.
                {
                    let inner = sv.list_vector_child(5);
                    {
                        let cc = inner.child(so_flat.len().max(1));
                        for (k, s) in so_flat.iter().enumerate() {
                            cc.insert(k, *s);
                        }
                    }
                    inner.set_len(so_flat.len());
                }
                let mut entry_vec = sv.child(5, total);
                let entries = entry_vec.as_mut_slice::<duckdb::ffi::duckdb_list_entry>();
                for j in 0..total {
                    entries[j].offset = so_offsets[j] as u64;
                    entries[j].length = (so_offsets[j + 1] - so_offsets[j]) as u64;
                }
            }
            vcol!(6, impact);
            vcol!(7, amino_acids);
            vcol!(8, codons);
            {
                let mut v = sv.child(9, total);
                {
                    let s = v.as_mut_slice::<i64>();
                    for (j, r) in flat.iter().enumerate() {
                        s[j] = r.protein_pos.unwrap_or(0);
                    }
                }
                for (j, r) in flat.iter().enumerate() {
                    if r.protein_pos.is_none() {
                        v.set_null(j);
                    }
                }
            }
            vcol!(10, alt);
            vcol!(11, hgvsc);
            vcol!(12, hgvsp);
            vcol!(13, hgvsg);
        }
        for i in 0..nrows {
            list.set_entry(i, list_offsets[i], list_offsets[i + 1] - list_offsets[i]);
        }
        list.set_len(total);
        Ok(())
    }

    fn signatures() -> Vec<ScalarFunctionSignature> {
        let bigint = || LogicalTypeHandle::from(LogicalTypeId::Bigint);
        let ret = || LogicalTypeHandle::list(&consequence_struct());
        vec![
            // (chrom, pos, ref, alt)
            ScalarFunctionSignature::exact(vec![varchar(), bigint(), varchar(), varchar()], ret()),
            // (chrom, pos, end, ref, alt) — END spans structural variants.
            ScalarFunctionSignature::exact(
                vec![varchar(), bigint(), bigint(), varchar(), varchar()],
                ret(),
            ),
        ]
    }
}
