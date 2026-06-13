//! The annotation engine core — **DuckDB-free**.
//!
//! This module wraps fastVEP's consequence engine and produces plain owned
//! rows. It has no `duckdb` dependency, which is the point: the DuckDB
//! extension (`super::annotate`), a future CLI binary, and a future C API are
//! all thin frontends over this one function. Only the extension links the
//! DuckDB C API; the CLI/C-API do not link libduckdb at all.

use fastvep_cache::fasta::{FastaReader, MmapFastaReader};
use fastvep_cache::gff::parse_gff3;
use fastvep_cache::providers::{
    FastaSequenceProvider, IndexedTranscriptProvider, MmapFastaSequenceProvider, SequenceProvider,
    TranscriptProvider,
};
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::{Allele, GenomicPosition, Strand};
use noodles::vcf;
use std::error::Error;
use std::fs::File;
use std::path::Path;

pub const DEFAULT_DISTANCE: u64 = 5000;

type SeqProvider = Box<dyn SequenceProvider + Send + Sync>;

/// The lean consequence context: transcript index + reference + predictor.
/// Built directly from the cache/consequence crates, deliberately bypassing
/// fastVEP's `AnnotationContext` (which also drags in SA/ACMG/HGVS/IO).
///
/// `Send + Sync` so it can live in the extension's load-once cache and be shared
/// across DuckDB worker threads (`super::consequence`).
pub(crate) struct EngineContext {
    transcripts: IndexedTranscriptProvider,
    seq: Option<SeqProvider>,
    predictor: ConsequencePredictor,
    distance: u64,
}

pub(crate) fn build_context(
    gff3: &str,
    fasta: Option<&str>,
    distance: u64,
) -> Result<EngineContext, Box<dyn Error>> {
    let gff_file = File::open(gff3)?;
    let mut transcripts = if gff3.ends_with(".gz") || gff3.ends_with(".bgz") {
        parse_gff3(flate2::read::MultiGzDecoder::new(gff_file))?
    } else {
        parse_gff3(gff_file)?
    };

    let seq: Option<SeqProvider> = match fasta {
        Some(path) if Path::new(&format!("{path}.fai")).exists() => Some(Box::new(
            MmapFastaSequenceProvider::new(MmapFastaReader::open(Path::new(path))?),
        )),
        Some(path) => Some(Box::new(FastaSequenceProvider::new(
            FastaReader::from_reader(File::open(path)?)?,
        ))),
        None => None,
    };

    // Build spliced cDNA / protein sequences so coding consequences (missense,
    // synonymous, amino acids) are exact. Requires a reference.
    if let Some(sp) = &seq {
        for t in &mut transcripts {
            let _ = t.build_sequences(|chrom, start, end| {
                sp.fetch_sequence(chrom, start, end)
                    .map_err(|e| e.to_string())
            });
        }
    }

    Ok(EngineContext {
        transcripts: IndexedTranscriptProvider::new(transcripts),
        seq,
        predictor: ConsequencePredictor::new(distance, distance),
        distance,
    })
}

impl EngineContext {
    /// Annotate a single variant against overlapping transcripts — the pure
    /// per-variant kernel. No IO of variants; DuckDB feeds the rows. Used by
    /// both the streaming `annotate_all` and the scalar `vep_consequence`.
    pub(crate) fn annotate_variant(
        &self,
        chrom: &str,
        pos: u64,
        ref_str: &str,
        alt_raw: &str,
    ) -> Vec<AnnotatedRow> {
        let mut rows = Vec::new();
        if alt_raw.is_empty() || alt_raw == "." {
            return rows;
        }
        let end = pos + (ref_str.len() as u64).saturating_sub(1);
        let ref_allele = Allele::from_str(ref_str);
        let alt_alleles: Vec<Allele> = alt_raw.split(',').map(Allele::from_str).collect();

        let position = GenomicPosition::new(chrom.to_string(), pos, end, Strand::Forward);
        let query_start = pos.saturating_sub(self.distance).max(1);
        let query_end = end + self.distance;
        let overlapping = self
            .transcripts
            .get_transcripts(chrom, query_start, query_end)
            .unwrap_or_default();
        if overlapping.is_empty() {
            return rows;
        }
        let ref_seq = self
            .seq
            .as_ref()
            .and_then(|sp| sp.fetch_sequence(chrom, query_start, query_end).ok());

        let result = self.predictor.predict(
            &position,
            &ref_allele,
            &alt_alleles,
            &overlapping,
            ref_seq.as_deref(),
        );

        for tc in &result.transcript_consequences {
            for ac in &tc.allele_consequences {
                rows.push(AnnotatedRow {
                    chrom: chrom.to_string(),
                    pos: pos as i64,
                    reference: ref_str.to_string(),
                    alt: ac.allele.to_string(),
                    gene_id: tc.gene_id.to_string(),
                    gene_symbol: tc.gene_symbol.as_deref().unwrap_or("").to_string(),
                    transcript_id: tc.transcript_id.to_string(),
                    biotype: tc.biotype.to_string(),
                    canonical: tc.canonical,
                    consequence: ac
                        .consequences
                        .iter()
                        .map(|c| c.so_term().to_string())
                        .collect(),
                    impact: ac.impact.as_str().to_string(),
                    amino_acids: pair_to_string(&ac.amino_acids),
                    codons: pair_to_string(&ac.codons),
                    protein_pos: ac.protein_start.map(|p| p as i64),
                });
            }
        }
        rows
    }
}

/// One annotated (variant, transcript, allele) row. Plain owned data so it can
/// cross any frontend boundary (DuckDB vectors, CLI output, C API).
pub(crate) struct AnnotatedRow {
    pub chrom: String,
    pub pos: i64,
    pub reference: String,
    pub alt: String,
    pub gene_id: String,
    pub gene_symbol: String,
    pub transcript_id: String,
    pub biotype: String,
    pub canonical: bool,
    pub consequence: Vec<String>,
    pub impact: String,
    pub amino_acids: String,
    pub codons: String,
    pub protein_pos: Option<i64>,
}

/// Annotate every variant in `vcf_path` against the `gff3` gene model (optional
/// `fasta` for sequence-dependent consequences). The single library entry point.
pub(crate) fn annotate(
    vcf_path: &str,
    gff3: &str,
    fasta: Option<&str>,
    distance: u64,
) -> Result<Vec<AnnotatedRow>, Box<dyn Error>> {
    let ctx = build_context(gff3, fasta, distance)?;
    annotate_all(&ctx, vcf_path, distance)
}

fn pair_to_string(p: &Option<(String, String)>) -> String {
    match p {
        Some((a, b)) => format!("{a}/{b}"),
        None => String::new(),
    }
}

fn annotate_all(
    ctx: &EngineContext,
    vcf_path: &str,
    _distance: u64,
) -> Result<Vec<AnnotatedRow>, Box<dyn Error>> {
    let mut reader = vcf::io::reader::Builder::default().build_from_path(vcf_path)?;
    let _header = reader.read_header()?;
    let mut record = vcf::Record::default();
    let mut rows = Vec::new();

    while reader.read_record(&mut record)? != 0 {
        let chrom = record.reference_sequence_name().to_string();
        let pos = record
            .variant_start()
            .transpose()?
            .map(usize::from)
            .unwrap_or(0) as u64;
        let ref_str = record.reference_bases().to_string();
        let alt_raw = record.alternate_bases().as_ref().to_string();
        rows.extend(ctx.annotate_variant(&chrom, pos, &ref_str, &alt_raw));
    }
    Ok(rows)
}
