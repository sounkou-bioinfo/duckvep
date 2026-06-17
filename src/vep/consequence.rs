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
use crate::vep::tcache;
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

/// Lock-free snapshot of the loaded engine, for other modules (e.g. the
/// `vep_transcripts()` table function) that need to read the shared gene model.
pub(crate) fn loaded_engine() -> Option<std::sync::Arc<EngineContext>> {
    cache().load_full()
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
        // up/downstream window is an EXPLICIT parameter (VEP's --distance), defaulting
        // to the VEP-compatible 5000 — not a buried constant. The 3-arg overload
        // `vep_load_cache(gff3, fasta, distance)` sets it.
        let distance = if input.num_columns() >= 3 {
            let d = batch.column(2).as_primitive::<Int64Type>();
            if d.is_null(0) {
                DEFAULT_DISTANCE
            } else {
                d.value(0).max(0) as u64
            }
        } else {
            DEFAULT_DISTANCE
        };
        // zstd level for the columnar transcript cache — an EXPLICIT space/build-time
        // tradeoff (1 = fastest build/largest file … 22 = smallest/slowest), not a
        // buried constant. Only consulted on a COLD build; ignored when an existing
        // cache is loaded. The 4-arg overload `vep_load_cache(gff3, fasta, distance,
        // zstd_level)` sets it; otherwise the default applies.
        let zstd = if input.num_columns() >= 4 {
            let z = batch.column(3).as_primitive::<Int64Type>();
            if z.is_null(0) {
                tcache::DEFAULT_ZSTD
            } else {
                z.value(0) as i32
            }
        } else {
            tcache::DEFAULT_ZSTD
        };
        let ctx = build_context(gff3.value(0), fastav, distance, zstd)?;
        cache().store(Some(Arc::new(ctx)));

        let v = output.flat_vector();
        for i in 0..input.len() {
            v.insert(i, "loaded");
        }
        Ok(())
    }

    fn signatures() -> Vec<ScalarFunctionSignature> {
        let bigint = || LogicalTypeHandle::from(LogicalTypeId::Bigint);
        vec![
            // vep_load_cache(gff3, fasta)            -> distance defaults to 5000 (VEP)
            ScalarFunctionSignature::exact(vec![varchar(), varchar()], varchar()),
            // vep_load_cache(gff3, fasta, distance)  -> explicit up/downstream window
            ScalarFunctionSignature::exact(vec![varchar(), varchar(), bigint()], varchar()),
            // vep_load_cache(gff3, fasta, distance, zstd_level) -> + cold-build zstd level
            ScalarFunctionSignature::exact(
                vec![varchar(), varchar(), bigint(), bigint()],
                varchar(),
            ),
        ]
    }
}

/// `vep_consequence_pair_idx(chrom, pos, [end_pos,] ref, alt, tx_idx) -> STRUCT`
/// — the ordinal-keyed per-pair kernel. Identical to `vep_consequence_pair` but
/// the transcript is named by its compact `tx_idx` (UINTEGER) ordinal, so the
/// candidate join carries an integer (not a `transcript_id` string) and lookup is
/// a direct array index. The no-regret hot-path win (docs/kernel-algorithm.md §9).
pub struct VepConsequencePairIdx;

impl VScalar for VepConsequencePairIdx {
    type State = ();

    unsafe fn invoke(
        _: &Self::State,
        input: &mut DataChunkHandle,
        output: &mut dyn WritableVector,
    ) -> Result<(), Box<dyn Error>> {
        let batch = data_chunk_to_arrow(input)?;
        // Arities: (chrom,pos,ref,alt,tx_idx) and END-aware
        // (chrom,pos,end_pos,ref,alt,tx_idx).
        let has_end = batch.num_columns() == 6;
        let chrom = batch.column(0).as_string::<i32>();
        let pos = batch.column(1).as_primitive::<Int64Type>();
        let end = has_end.then(|| batch.column(2).as_primitive::<Int64Type>());
        let refa = batch.column(if has_end { 3 } else { 2 }).as_string::<i32>();
        let alta = batch.column(if has_end { 4 } else { 3 }).as_string::<i32>();
        let tidx = batch.column(if has_end { 5 } else { 4 }).as_primitive::<Int64Type>();
        let nrows = input.len();

        let mut rows: Vec<Option<AnnotatedRow>> = Vec::with_capacity(nrows);
        {
            let ctx = cache()
                .load_full()
                .ok_or("vep_consequence_pair_idx: call vep_load_cache(gff3, fasta) first")?;
            for i in 0..nrows {
                let p = pos.value(i) as u64;
                let e = end.map(|c| c.value(i) as u64).unwrap_or(p);
                let mut v = ctx.annotate_pair_idx(
                    chrom.value(i),
                    p,
                    e,
                    refa.value(i),
                    alta.value(i),
                    tidx.value(i).max(0) as usize,
                );
                rows.push(v.pop());
            }
        }
        write_consequence_struct(output, &rows, nrows);
        Ok(())
    }

    fn signatures() -> Vec<ScalarFunctionSignature> {
        let bigint = || LogicalTypeHandle::from(LogicalTypeId::Bigint);
        let ret = consequence_struct;
        vec![
            // (chrom, pos, ref, alt, tx_idx)
            ScalarFunctionSignature::exact(
                vec![varchar(), bigint(), varchar(), varchar(), bigint()],
                ret(),
            ),
            // (chrom, pos, end_pos, ref, alt, tx_idx) — END-aware (SV span)
            ScalarFunctionSignature::exact(
                vec![varchar(), bigint(), bigint(), varchar(), varchar(), bigint()],
                ret(),
            ),
        ]
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
        write_consequence_list(output, &flat, &list_offsets, nrows);
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

/// Write a flat `[AnnotatedRow]` (grouped per input row by `list_offsets`) into
/// `output` as `LIST<STRUCT>` matching `consequence_struct()`. Shared by the
/// per-variant (`vep_consequence`) and per-pair (`vep_consequence_pair`) kernels
/// so the intricate nested-vector marshalling lives in exactly one place.
unsafe fn write_consequence_list(
    output: &mut dyn WritableVector,
    flat: &[AnnotatedRow],
    list_offsets: &[usize],
    nrows: usize,
) {
    let total = flat.len();

    // Flatten the nested consequence (LIST<VARCHAR>) child: one inner list per
    // struct row, `so_offsets[j]..so_offsets[j+1]` slicing the shared values.
    let mut so_flat: Vec<&str> = Vec::new();
    let mut so_offsets: Vec<usize> = Vec::with_capacity(total + 1);
    so_offsets.push(0);
    for r in flat {
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
}

/// Write one `consequence_struct()` per input row into `output` as a nullable
/// top-level `STRUCT` (NULL where the pair yields no consequence). Used by the
/// per-pair kernel: the range join already selected exactly one transcript, so
/// the result is 0-or-1 — a `STRUCT` is the right shape, not a 1-element list.
/// The inner `consequence` field stays `LIST<VARCHAR>` (one pair can hit several
/// SO terms, e.g. splice_region & intron).
unsafe fn write_consequence_struct(
    output: &mut dyn WritableVector,
    rows: &[Option<AnnotatedRow>],
    nrows: usize,
) {
    // Per struct row, the slice of SO terms; None rows contribute an empty list.
    let mut so_flat: Vec<&str> = Vec::new();
    let mut so_offsets: Vec<usize> = Vec::with_capacity(nrows + 1);
    so_offsets.push(0);
    for r in rows {
        if let Some(r) = r {
            so_flat.extend(r.consequence.iter().map(String::as_str));
        }
        so_offsets.push(so_flat.len());
    }

    let mut sv = output.struct_vector();
    // varchar/scalar children are written for every row (empty for None rows);
    // the struct validity is nulled at the end so DuckDB ignores those slots.
    macro_rules! vcol {
        ($i:expr, $f:ident) => {{
            let v = sv.child($i, nrows);
            for (j, r) in rows.iter().enumerate() {
                match r {
                    Some(r) => v.insert(j, &*r.$f),
                    None => v.insert(j, ""),
                }
            }
        }};
    }
    vcol!(0, transcript_id);
    vcol!(1, gene_id);
    vcol!(2, gene_symbol);
    vcol!(3, biotype);
    {
        let mut v = sv.child(4, nrows);
        let s = v.as_mut_slice::<bool>();
        for (j, r) in rows.iter().enumerate() {
            s[j] = r.as_ref().map(|r| r.canonical).unwrap_or(false);
        }
    }
    {
        // consequence: LIST<VARCHAR> child of the struct (see write_consequence_list
        // for why we write the inner entry array directly rather than set_entry).
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
        let mut entry_vec = sv.child(5, nrows);
        let entries = entry_vec.as_mut_slice::<duckdb::ffi::duckdb_list_entry>();
        for j in 0..nrows {
            entries[j].offset = so_offsets[j] as u64;
            entries[j].length = (so_offsets[j + 1] - so_offsets[j]) as u64;
        }
    }
    vcol!(6, impact);
    vcol!(7, amino_acids);
    vcol!(8, codons);
    {
        let mut v = sv.child(9, nrows);
        {
            let s = v.as_mut_slice::<i64>();
            for (j, r) in rows.iter().enumerate() {
                s[j] = r.as_ref().and_then(|r| r.protein_pos).unwrap_or(0);
            }
        }
        for (j, r) in rows.iter().enumerate() {
            if r.as_ref().and_then(|r| r.protein_pos).is_none() {
                v.set_null(j);
            }
        }
    }
    vcol!(10, alt);
    vcol!(11, hgvsc);
    vcol!(12, hgvsp);
    vcol!(13, hgvsg);

    // Null the whole struct for rows that produced no consequence.
    for (j, r) in rows.iter().enumerate() {
        if r.is_none() {
            sv.set_null(j);
        }
    }
}

/// `vep_consequence_pair(chrom, pos, [end_pos,] ref, alt, transcript_id) -> STRUCT`
/// — the **per-pair** kernel. Annotates a variant against ONE named transcript
/// (looked up by id), so DuckDB can drive annotation from a parallel range join
/// instead of the serial per-variant spatial lookup. Returns a nullable `STRUCT`
/// (NULL when the pair yields no consequence — filter with `WHERE c IS NOT NULL`),
/// not a 1-element `LIST`, since the join already picked exactly one transcript.
/// The 6-arg form carries `end_pos` (the full ref span / `INFO/END`) so symbolic
/// structural alleles span their interval. See docs/kernel-algorithm.md §6/§8.
pub struct VepConsequencePair;

impl VScalar for VepConsequencePair {
    type State = ();

    unsafe fn invoke(
        _: &Self::State,
        input: &mut DataChunkHandle,
        output: &mut dyn WritableVector,
    ) -> Result<(), Box<dyn Error>> {
        let batch = data_chunk_to_arrow(input)?;
        // Arities: (chrom,pos,ref,alt,tid) and the END-aware
        // (chrom,pos,end_pos,ref,alt,tid) that spans structural variants.
        let has_end = batch.num_columns() == 6;
        let chrom = batch.column(0).as_string::<i32>();
        let pos = batch.column(1).as_primitive::<Int64Type>();
        let end = has_end.then(|| batch.column(2).as_primitive::<Int64Type>());
        let refa = batch.column(if has_end { 3 } else { 2 }).as_string::<i32>();
        let alta = batch.column(if has_end { 4 } else { 3 }).as_string::<i32>();
        let tid = batch.column(if has_end { 5 } else { 4 }).as_string::<i32>();
        let nrows = input.len();

        let mut rows: Vec<Option<AnnotatedRow>> = Vec::with_capacity(nrows);
        {
            let ctx = cache()
                .load_full()
                .ok_or("vep_consequence_pair: call vep_load_cache(gff3, fasta) first")?;
            for i in 0..nrows {
                let p = pos.value(i) as u64;
                let e = end.map(|c| c.value(i) as u64).unwrap_or(p);
                // One (variant, transcript, allele) yields 0 or 1 consequence row.
                let mut v =
                    ctx.annotate_pair(chrom.value(i), p, e, refa.value(i), alta.value(i), tid.value(i));
                rows.push(v.pop());
            }
        }
        write_consequence_struct(output, &rows, nrows);
        Ok(())
    }

    fn signatures() -> Vec<ScalarFunctionSignature> {
        let bigint = || LogicalTypeHandle::from(LogicalTypeId::Bigint);
        let ret = consequence_struct;
        vec![
            // (chrom, pos, ref, alt, transcript_id) — end defaults to pos
            ScalarFunctionSignature::exact(
                vec![varchar(), bigint(), varchar(), varchar(), varchar()],
                ret(),
            ),
            // (chrom, pos, end_pos, ref, alt, transcript_id) — END-aware (SV span)
            ScalarFunctionSignature::exact(
                vec![varchar(), bigint(), bigint(), varchar(), varchar(), varchar()],
                ret(),
            ),
        ]
    }
}

/// `vep_haplotype_consequence(chrom, transcript_id, variants) -> VARCHAR`
///
/// Combine a set of PHASED variants on ONE transcript into a single haplotype coding
/// consequence (the bcftools-csq / haplosaurus model: co-located variants are merged in
/// transcript/CDS coordinates and translated together, so an in-codon SNV pair becomes
/// one MNV and a frameshift + restoring indel cancel — the capability fastVEP lacks).
///
/// `variants` is a `;`-delimited list of `pos:ref:alt` (e.g. `'1053:G:T;1055:T:A'`),
/// built in SQL with `string_agg(pos||':'||ref||':'||alt, ';')` after
/// `GROUP BY sample, haplotype, transcript_id` — the grouping layer stays in SQL.
/// Returns the `&`-joined sorted SO terms (empty string if the transcript is out of
/// range or no variant is coding).
pub struct VepHaplotypeConsequence;

impl VScalar for VepHaplotypeConsequence {
    type State = ();

    unsafe fn invoke(
        _: &Self::State,
        input: &mut DataChunkHandle,
        output: &mut dyn WritableVector,
    ) -> Result<(), Box<dyn Error>> {
        let batch = data_chunk_to_arrow(input)?;
        let chrom = batch.column(0).as_string::<i32>();
        let tid = batch.column(1).as_string::<i32>();
        let vars = batch.column(2).as_string::<i32>();
        let nrows = input.len();

        let mut out_strings: Vec<String> = Vec::with_capacity(nrows);
        {
            let ctx = cache()
                .load_full()
                .ok_or("vep_haplotype_consequence: call vep_load_cache(gff3, fasta) first")?;
            for i in 0..nrows {
                if chrom.is_null(i) || tid.is_null(i) || vars.is_null(i) {
                    out_strings.push(String::new());
                    continue;
                }
                let variants = parse_haplotype_variants(vars.value(i));
                let terms = ctx.haplotype_consequence(chrom.value(i), tid.value(i), &variants);
                out_strings.push(terms.join("&"));
            }
        }

        let v = output.flat_vector();
        for (i, s) in out_strings.iter().enumerate() {
            v.insert(i, s.as_str());
        }
        Ok(())
    }

    fn signatures() -> Vec<ScalarFunctionSignature> {
        vec![ScalarFunctionSignature::exact(
            vec![varchar(), varchar(), varchar()],
            varchar(),
        )]
    }
}

/// Parse a `;`-delimited `pos:ref:alt` list into `(pos, ref, alt)` triples. Tokens that
/// don't have all three colon-separated fields (or a non-numeric pos) are skipped.
fn parse_haplotype_variants(s: &str) -> Vec<(u64, String, String)> {
    s.split(';')
        .filter_map(|tok| {
            let tok = tok.trim();
            if tok.is_empty() {
                return None;
            }
            let mut it = tok.split(':');
            let pos = it.next()?.trim().parse::<u64>().ok()?;
            let r = it.next()?.trim().to_string();
            let a = it.next()?.trim().to_string();
            Some((pos, r, a))
        })
        .collect()
}
