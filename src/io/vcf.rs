//! `read_vcf(path [, region := 'chr:start-end'])` — a DuckDB table function that
//! reads VCF/BCF(.gz) via noodles and emits one row per record.
//!
//! v1 eagerly materializes rows in `bind` and serves them from `func`; region
//! filtering is done in memory. Indexed (tabix/csi) streaming and `VARIANT`
//! INFO are tracked follow-ups (see DESIGN.md §3.1, §8).
//!
//! `alt` and `filter` are `LIST<VARCHAR>` so multiallelic (`A,AT`), symbolic
//! (`<DEL>`, `<CNV>`), and breakend alleles are first-class.

use duckdb::core::{DataChunkHandle, Inserter, LogicalTypeHandle, LogicalTypeId};
use duckdb::vtab::{BindInfo, InitInfo, TableFunctionInfo, VTab};
use noodles::vcf;
use std::error::Error;
use std::sync::atomic::{AtomicUsize, Ordering};

/// One materialized VCF record, owned (so `BindData` is `Send + Sync`).
pub struct VcfRow {
    chrom: String,
    pos: i64,
    /// End coordinate of the variant interval. `INFO/END` for SV/CNV with
    /// symbolic ALTs, otherwise `pos + len(ref) - 1` for precise variants.
    end: i64,
    id: String,
    reference: String,
    alt: Vec<String>,
    qual: Option<f64>,
    filter: Vec<String>,
    info: String,
    /// Per-sample GT strings in header sample order (phasing kept as `0|1`).
    /// Empty when the VCF has no genotype columns. Sample names come from
    /// `vcf_samples(path)` (annotation is site-wise, so genotypes ride along).
    gt: Vec<String>,
}

pub struct ReadVcf;

pub struct VcfBind {
    rows: Vec<VcfRow>,
}

pub struct VcfInit {
    cursor: AtomicUsize,
}

/// DuckDB's standard vector size.
const VECTOR_SIZE: usize = 2048;

impl VTab for ReadVcf {
    type BindData = VcfBind;
    type InitData = VcfInit;

    fn bind(bind: &BindInfo) -> Result<Self::BindData, Box<dyn Error>> {
        let varchar = || LogicalTypeHandle::from(LogicalTypeId::Varchar);
        bind.add_result_column("chrom", varchar());
        bind.add_result_column("pos", LogicalTypeHandle::from(LogicalTypeId::Bigint));
        // `end` is a SQL reserved word; expose as `end_pos` to stay unquoted.
        bind.add_result_column("end_pos", LogicalTypeHandle::from(LogicalTypeId::Bigint));
        bind.add_result_column("id", varchar());
        bind.add_result_column("ref", varchar());
        bind.add_result_column("alt", LogicalTypeHandle::list(&varchar()));
        bind.add_result_column("qual", LogicalTypeHandle::from(LogicalTypeId::Double));
        bind.add_result_column("filter", LogicalTypeHandle::list(&varchar()));
        bind.add_result_column("info", varchar());
        bind.add_result_column("gt", LogicalTypeHandle::list(&varchar()));

        let path = bind.get_parameter(0).to_string();
        let region = bind
            .get_named_parameter("region")
            .map(|v| v.to_string())
            .filter(|s| !s.is_empty());

        let rows = read_all(&path, region.as_deref())?;
        Ok(VcfBind { rows })
    }

    fn init(_: &InitInfo) -> Result<Self::InitData, Box<dyn Error>> {
        Ok(VcfInit {
            cursor: AtomicUsize::new(0),
        })
    }

    fn func(
        func: &TableFunctionInfo<Self>,
        output: &mut DataChunkHandle,
    ) -> Result<(), Box<dyn Error>> {
        let bind = func.get_bind_data();
        let init = func.get_init_data();

        let start = init.cursor.load(Ordering::Relaxed);
        let total = bind.rows.len();
        let n = (total - start).min(VECTOR_SIZE);
        let rows = &bind.rows[start..start + n];

        // Scalar columns, one vector borrow live at a time.
        {
            let v = output.flat_vector(0);
            for (i, r) in rows.iter().enumerate() {
                v.insert(i, r.chrom.as_str());
            }
        }
        {
            let mut v = output.flat_vector(1);
            let s = unsafe { v.as_mut_slice::<i64>() };
            for (i, r) in rows.iter().enumerate() {
                s[i] = r.pos;
            }
        }
        {
            let mut v = output.flat_vector(2);
            let s = unsafe { v.as_mut_slice::<i64>() };
            for (i, r) in rows.iter().enumerate() {
                s[i] = r.end;
            }
        }
        {
            let v = output.flat_vector(3);
            for (i, r) in rows.iter().enumerate() {
                v.insert(i, r.id.as_str());
            }
        }
        {
            let v = output.flat_vector(4);
            for (i, r) in rows.iter().enumerate() {
                v.insert(i, r.reference.as_str());
            }
        }
        fill_string_list(output, 5, rows, |r| &r.alt);
        {
            let mut v = output.flat_vector(6);
            {
                let s = unsafe { v.as_mut_slice::<f64>() };
                for (i, r) in rows.iter().enumerate() {
                    s[i] = r.qual.unwrap_or(0.0);
                }
            }
            for (i, r) in rows.iter().enumerate() {
                if r.qual.is_none() {
                    v.set_null(i);
                }
            }
        }
        fill_string_list(output, 7, rows, |r| &r.filter);
        {
            let v = output.flat_vector(8);
            for (i, r) in rows.iter().enumerate() {
                v.insert(i, r.info.as_str());
            }
        }
        fill_string_list(output, 9, rows, |r| &r.gt);

        init.cursor.store(start + n, Ordering::Relaxed);
        output.set_len(n);
        Ok(())
    }

    fn parameters() -> Option<Vec<LogicalTypeHandle>> {
        Some(vec![LogicalTypeHandle::from(LogicalTypeId::Varchar)])
    }

    fn named_parameters() -> Option<Vec<(String, LogicalTypeHandle)>> {
        Some(vec![(
            "region".to_string(),
            LogicalTypeHandle::from(LogicalTypeId::Varchar),
        )])
    }
}

/// Fill a `LIST<VARCHAR>` output column from a per-row slice of strings.
fn fill_string_list(
    output: &mut DataChunkHandle,
    col: usize,
    rows: &[VcfRow],
    get: impl Fn(&VcfRow) -> &[String],
) {
    let total: usize = rows.iter().map(|r| get(r).len()).sum();
    let mut list = output.list_vector(col);
    {
        let child = list.child(total.max(1));
        let mut off = 0usize;
        for r in rows {
            for (j, s) in get(r).iter().enumerate() {
                child.insert(off + j, s.as_str());
            }
            off += get(r).len();
        }
    }
    let mut off = 0usize;
    for (i, r) in rows.iter().enumerate() {
        let len = get(r).len();
        list.set_entry(i, off, len);
        off += len;
    }
    list.set_len(total);
}

/// Parsed `chrom[:start-end]` region filter.
struct Region {
    chrom: String,
    start: i64,
    end: i64,
}

fn parse_region(s: &str) -> Region {
    match s.split_once(':') {
        None => Region {
            chrom: s.to_string(),
            start: i64::MIN,
            end: i64::MAX,
        },
        Some((chrom, range)) => {
            let (start, end) = range
                .split_once('-')
                .map(|(a, b)| {
                    (
                        a.replace(',', "").parse().unwrap_or(i64::MIN),
                        b.replace(',', "").parse().unwrap_or(i64::MAX),
                    )
                })
                .unwrap_or((i64::MIN, i64::MAX));
            Region {
                chrom: chrom.to_string(),
                start,
                end,
            }
        }
    }
}

/// Variant END coordinate. `INFO/END` (SV/CNV with symbolic ALTs) takes
/// precedence; otherwise the precise interval `pos + len(ref) - 1`.
fn compute_end(pos: i64, reference: &str, info_raw: &str) -> i64 {
    for field in info_raw.split(';') {
        if let Some(v) = field.strip_prefix("END=") {
            if let Ok(e) = v.parse::<i64>() {
                return e;
            }
        }
    }
    pos + (reference.len() as i64 - 1).max(0)
}

/// Per-sample GT strings (header sample order), phasing preserved as `0|1`.
/// `raw` is the noodles samples region: `"<FORMAT keys>\t<sample1>\t<sample2>"`.
fn parse_gts(raw: &str) -> Vec<String> {
    if raw.is_empty() {
        return Vec::new();
    }
    let mut fields = raw.split('\t');
    let keys = fields.next().unwrap_or("");
    let gt_idx = keys.split(':').position(|k| k == "GT");
    fields
        .map(|sample| match gt_idx {
            Some(gi) => sample.split(':').nth(gi).unwrap_or(".").to_string(),
            None => ".".to_string(),
        })
        .collect()
}

/// Split a VCF list field, treating `.`/empty as an empty list.
fn split_field(raw: &str, sep: char) -> Vec<String> {
    if raw.is_empty() || raw == "." {
        Vec::new()
    } else {
        raw.split(sep).map(|s| s.to_string()).collect()
    }
}

fn read_all(path: &str, region: Option<&str>) -> Result<Vec<VcfRow>, Box<dyn Error>> {
    let filt = region.map(parse_region);
    let mut reader = vcf::io::reader::Builder::default().build_from_path(path)?;
    let _header = reader.read_header()?;

    let mut record = vcf::Record::default();
    let mut rows = Vec::new();
    while reader.read_record(&mut record)? != 0 {
        let chrom = record.reference_sequence_name().to_string();
        let pos = record
            .variant_start()
            .transpose()?
            .map(usize::from)
            .unwrap_or(0) as i64;

        if let Some(r) = &filt {
            if chrom != r.chrom || pos < r.start || pos > r.end {
                continue;
            }
        }

        let reference = record.reference_bases().to_string();
        let info = record.info().as_ref().to_string();
        let end = compute_end(pos, &reference, &info);
        let gt = parse_gts(record.samples().as_ref());
        rows.push(VcfRow {
            chrom,
            pos,
            end,
            id: record.ids().as_ref().to_string(),
            reference,
            alt: split_field(record.alternate_bases().as_ref(), ','),
            qual: record.quality_score().transpose()?.map(|q| q as f64),
            filter: split_field(record.filters().as_ref(), ';'),
            info,
            gt,
        });
    }
    Ok(rows)
}
