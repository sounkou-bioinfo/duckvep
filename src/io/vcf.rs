//! `read_vcf(path [, region := 'chr:start-end'])` — a DuckDB table function that
//! reads VCF/BCF(.gz) via noodles and emits one row per record.
//!
//! **Streams**: the reader is opened lazily and `func` reads one ~2048-row chunk
//! per call, so memory is bounded to a chunk regardless of file size (full GIAB
//! 4M variants: ~2 GB vs ~7 GB eager). Region filtering is applied per record.
//! Indexed (tabix/csi) seeking and `VARIANT` INFO are follow-ups (docs/DESIGN.md §3.1).
//!
//! `alt` and `filter` are `LIST<VARCHAR>` so multiallelic (`A,AT`), symbolic
//! (`<DEL>`, `<CNV>`), and breakend alleles are first-class.

use crate::vec_util::fill_string_list;
use duckdb::core::{DataChunkHandle, Inserter, LogicalTypeHandle, LogicalTypeId};
use duckdb::vtab::{BindInfo, InitInfo, TableFunctionInfo, VTab};
use noodles::{bgzf, vcf};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

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
    path: String,
    region: Option<Region>,
}

/// A `Send` VCF reader (plain or bgzf), opened lazily in `func`. Held behind a
/// Mutex so `InitData` is `Send + Sync`; `read_vcf` scans single-threaded.
enum VcfRdr {
    Plain(vcf::io::Reader<BufReader<File>>),
    Bgzf(vcf::io::Reader<bgzf::Reader<File>>),
}

impl VcfRdr {
    fn open(path: &str) -> io::Result<Self> {
        let mut file = File::open(path)?;
        let mut magic = [0u8; 2];
        let n = file.read(&mut magic)?;
        file.seek(SeekFrom::Start(0))?;
        if n == 2 && magic == [0x1f, 0x8b] {
            let mut r = vcf::io::Reader::new(bgzf::Reader::new(file));
            r.read_header()?;
            Ok(VcfRdr::Bgzf(r))
        } else {
            let mut r = vcf::io::Reader::new(BufReader::new(file));
            r.read_header()?;
            Ok(VcfRdr::Plain(r))
        }
    }

    fn read_record(&mut self, rec: &mut vcf::Record) -> io::Result<usize> {
        match self {
            VcfRdr::Plain(r) => r.read_record(rec),
            VcfRdr::Bgzf(r) => r.read_record(rec),
        }
    }
}

pub struct VcfInit {
    reader: Mutex<Option<VcfRdr>>,
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
            .filter(|s| !s.is_empty())
            .map(|s| parse_region(&s))
            .transpose()?;

        Ok(VcfBind { path, region })
    }

    fn init(_: &InitInfo) -> Result<Self::InitData, Box<dyn Error>> {
        Ok(VcfInit {
            reader: Mutex::new(None),
        })
    }

    fn func(
        func: &TableFunctionInfo<Self>,
        output: &mut DataChunkHandle,
    ) -> Result<(), Box<dyn Error>> {
        let bind = func.get_bind_data();
        let init = func.get_init_data();

        // Stream one chunk: read up to VECTOR_SIZE region-matching records.
        let mut guard = init.reader.lock().map_err(|_| "reader lock poisoned")?;
        if guard.is_none() {
            *guard = Some(VcfRdr::open(&bind.path)?);
        }
        let rdr = guard.as_mut().unwrap();
        let mut record = vcf::Record::default();
        let mut chunk: Vec<VcfRow> = Vec::with_capacity(VECTOR_SIZE);
        while chunk.len() < VECTOR_SIZE {
            if rdr.read_record(&mut record)? == 0 {
                break; // EOF
            }
            if let Some(row) = record_to_row(&record, bind.region.as_ref())? {
                chunk.push(row);
            }
        }
        let rows = &chunk[..];

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

        output.set_len(chunk.len());
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

/// `vcf_samples(path)` — one row per sample in header order: `(idx, sample)`.
///
/// `idx` is 1-based so it lines up with DuckDB's `UNNEST(... WITH ORDINALITY)`,
/// letting the positional `gt` list from [`ReadVcf`] be exploded into tidy
/// per-sample genotypes (annotation stays site-wise; see docs/DESIGN.md §3.0).
pub struct VcfSamples;

pub struct VcfSamplesBind {
    names: Vec<String>,
}

pub struct VcfSamplesInit {
    cursor: AtomicUsize,
}

impl VTab for VcfSamples {
    type BindData = VcfSamplesBind;
    type InitData = VcfSamplesInit;

    fn bind(bind: &BindInfo) -> Result<Self::BindData, Box<dyn Error>> {
        bind.add_result_column("idx", LogicalTypeHandle::from(LogicalTypeId::Bigint));
        bind.add_result_column("sample", LogicalTypeHandle::from(LogicalTypeId::Varchar));

        let path = bind.get_parameter(0).to_string();
        let mut reader = vcf::io::reader::Builder::default().build_from_path(&path)?;
        let header = reader.read_header()?;
        let names = header.sample_names().iter().cloned().collect();
        Ok(VcfSamplesBind { names })
    }

    fn init(_: &InitInfo) -> Result<Self::InitData, Box<dyn Error>> {
        Ok(VcfSamplesInit {
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
        let total = bind.names.len();
        let n = (total - start).min(VECTOR_SIZE);
        let names = &bind.names[start..start + n];

        {
            let mut v = output.flat_vector(0);
            let s = unsafe { v.as_mut_slice::<i64>() };
            for (i, _) in names.iter().enumerate() {
                s[i] = (start + i + 1) as i64;
            }
        }
        {
            let v = output.flat_vector(1);
            for (i, name) in names.iter().enumerate() {
                v.insert(i, name.as_str());
            }
        }

        init.cursor.store(start + n, Ordering::Relaxed);
        output.set_len(n);
        Ok(())
    }

    fn parameters() -> Option<Vec<LogicalTypeHandle>> {
        Some(vec![LogicalTypeHandle::from(LogicalTypeId::Varchar)])
    }
}

/// Parsed `chrom[:start-end]` region filter.
struct Region {
    chrom: String,
    start: i64,
    end: i64,
}

/// Normalize a contig name for comparison so `chr1` and `1` match (a common
/// VCF/region naming mismatch). Only the leading `chr` is stripped.
fn normalize_chrom(c: &str) -> &str {
    c.strip_prefix("chr").unwrap_or(c)
}

fn parse_region(s: &str) -> Result<Region, String> {
    let bad = |what: &str| format!("invalid region '{s}': bad {what}");
    match s.split_once(':') {
        None => Ok(Region {
            chrom: s.to_string(),
            start: i64::MIN,
            end: i64::MAX,
        }),
        Some((chrom, range)) => {
            let (a, b) = range.split_once('-').ok_or_else(|| bad("range"))?;
            let start = a.replace(',', "").parse().map_err(|_| bad("start"))?;
            let end = b.replace(',', "").parse().map_err(|_| bad("end"))?;
            Ok(Region {
                chrom: chrom.to_string(),
                start,
                end,
            })
        }
    }
}

/// Field value lookup in a raw VCF INFO string (`KEY=value;KEY2=...`).
fn info_value<'a>(info_raw: &'a str, key: &str) -> Option<&'a str> {
    info_raw
        .split(';')
        .find_map(|f| f.strip_prefix(key).filter(|r| r.starts_with('=')))
        .map(|r| &r[1..])
}

/// Variant END coordinate. `INFO/END` wins (SV/CNV with symbolic ALTs); else
/// `INFO/SVLEN` (insertions stay a reference point); else the precise interval
/// `pos + len(ref) - 1`.
fn compute_end(pos: i64, reference: &str, info_raw: &str) -> i64 {
    if let Some(e) = info_value(info_raw, "END").and_then(|v| v.parse::<i64>().ok()) {
        return e;
    }
    if let Some(svlen) = info_value(info_raw, "SVLEN")
        .and_then(|v| v.split(',').next())
        .and_then(|v| v.parse::<i64>().ok())
    {
        // Insertions add sequence between pos and pos+1; they don't span the
        // reference, so their interval end stays at pos.
        if info_value(info_raw, "SVTYPE") == Some("INS") {
            return pos;
        }
        return pos + svlen.abs();
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

/// Build a `VcfRow` from one record, or `None` if it is filtered out by `region`.
fn record_to_row(
    record: &vcf::Record,
    region: Option<&Region>,
) -> Result<Option<VcfRow>, Box<dyn Error>> {
    let chrom = record.reference_sequence_name().to_string();
    // A record with no POS is malformed; skip it (Ok(None)) rather than surfacing it at the
    // invalid coordinate 0 (`unwrap_or(0)` would have masked it). `?` propagates a parse Err.
    let pos = match record.variant_start().transpose()? {
        Some(p) => usize::from(p) as i64,
        None => return Ok(None),
    };
    let reference = record.reference_bases().to_string();
    let info = record.info().as_ref().to_string();
    let end = compute_end(pos, &reference, &info);

    if let Some(r) = region {
        // Keep records whose interval [pos, end] overlaps the region, so
        // SVs/indels starting before the region but spanning into it match.
        if normalize_chrom(&chrom) != normalize_chrom(&r.chrom) || pos > r.end || end < r.start {
            return Ok(None);
        }
    }

    let gt = parse_gts(record.samples().as_ref());
    Ok(Some(VcfRow {
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
    }))
}
