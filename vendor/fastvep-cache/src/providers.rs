use anyhow::Result;
use fastvep_genome::Transcript;
use std::collections::HashMap;
use std::sync::Arc;

use crate::info::CacheInfo;
use crate::variation::{self, VariationTabixReader};

/// Trait for providing transcript annotations for a genomic region.
pub trait TranscriptProvider: Send + Sync {
    /// Return all transcripts that overlap the given genomic region.
    fn get_transcripts(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<&Transcript>>;

    /// Return all transcripts on a chromosome.
    fn get_transcripts_by_chrom(&self, chrom: &str) -> Result<Vec<&Transcript>>;
}

/// Trait for providing reference sequences.
pub trait SequenceProvider: Send + Sync {
    /// Fetch reference sequence for a region (1-based, inclusive).
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>>;

    /// Fetch a reference sequence slice. Default delegates to fetch_sequence.
    fn fetch_sequence_slice(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        self.fetch_sequence(chrom, start, end)
    }
}

/// In-memory transcript provider backed by a Vec<Transcript>.
pub struct MemoryTranscriptProvider {
    transcripts: Vec<Transcript>,
}

impl MemoryTranscriptProvider {
    pub fn new(transcripts: Vec<Transcript>) -> Self {
        Self { transcripts }
    }

    pub fn transcript_count(&self) -> usize {
        self.transcripts.len()
    }
}

impl TranscriptProvider for MemoryTranscriptProvider {
    fn get_transcripts(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<&Transcript>> {
        Ok(self
            .transcripts
            .iter()
            .filter(|t| &*t.chromosome == chrom && t.start <= end && t.end >= start)
            .collect())
    }

    fn get_transcripts_by_chrom(&self, chrom: &str) -> Result<Vec<&Transcript>> {
        Ok(self
            .transcripts
            .iter()
            .filter(|t| &*t.chromosome == chrom)
            .collect())
    }
}

/// High-performance transcript provider using per-chromosome sorted arrays,
/// binary search, and suffix-max-end for O(log n + k) lookups with early termination.
///
/// Each chromosome's transcripts are sorted by start position. A parallel
/// `suffix_max_end` array stores the maximum `end` value from index `i` to the
/// end of the array, enabling early termination when scanning backwards:
/// if `suffix_max_end[i] < query_start`, no transcript at index <= i can overlap.
pub struct IndexedTranscriptProvider {
    /// ALL transcripts in one flat array, globally ordered by `(chrom, start,
    /// stable_id)`. The array index **is** the transcript ordinal (`tx_idx`) — the
    /// compact integer key DuckDB carries through the candidate join instead of a
    /// `transcript_id` string (see `get_by_idx`). This is also the first brick of
    /// the SoA execution model. Because the order is `(chrom, start, …)`, every
    /// chromosome's transcripts occupy a CONTIGUOUS ordinal range.
    all: Vec<Transcript>,
    /// Suffix-max-end, aligned 1:1 with `all`; within each chromosome's `[lo,hi)`
    /// segment, `sme[i] = max(all[i..hi].end)` — enables backward-scan early-out.
    sme: Vec<u64>,
    /// Chromosome → its contiguous `[lo, hi)` ordinal range in `all`.
    by_chrom: HashMap<Arc<str>, (usize, usize)>,
    /// stable_id → `tx_idx`.
    by_id: HashMap<Arc<str>, usize>,
}

impl IndexedTranscriptProvider {
    pub fn new(mut transcripts: Vec<Transcript>) -> Self {
        // Global total order — IDENTICAL to the cache writer's (vep::tcache::save)
        // sort, so the on-disk `tx_idx` (parquet row order) equals this ordinal.
        // stable_id breaks ties so the order is total and reproducible.
        transcripts.sort_by(|a, b| {
            a.chromosome
                .cmp(&b.chromosome)
                .then(a.start.cmp(&b.start))
                .then_with(|| a.stable_id.cmp(&b.stable_id))
        });
        let all = transcripts;
        let n = all.len();

        // Contiguous per-chromosome ranges + per-segment suffix-max-end.
        let mut by_chrom: HashMap<Arc<str>, (usize, usize)> = HashMap::new();
        let mut sme = vec![0u64; n];
        let mut lo = 0;
        while lo < n {
            let chrom = Arc::clone(&all[lo].chromosome);
            let mut hi = lo + 1;
            while hi < n && all[hi].chromosome == chrom {
                hi += 1;
            }
            // suffix-max-end within [lo, hi)
            sme[hi - 1] = all[hi - 1].end;
            for i in (lo..hi - 1).rev() {
                sme[i] = all[i].end.max(sme[i + 1]);
            }
            by_chrom.insert(chrom, (lo, hi));
            lo = hi;
        }

        let mut by_id = HashMap::with_capacity(n);
        for (i, t) in all.iter().enumerate() {
            by_id.insert(Arc::clone(&t.stable_id), i);
        }
        Self {
            all,
            sme,
            by_chrom,
            by_id,
        }
    }

    pub fn transcript_count(&self) -> usize {
        self.all.len()
    }

    /// Resolve a transcript by its compact ordinal (`tx_idx`) in O(1) — the
    /// per-pair kernel's lookup once DuckDB's range join has named the pair by
    /// integer ordinal rather than `transcript_id` string. Index = position in the
    /// `(chrom,start,stable_id)` order, matching the cache's `tx_idx` column.
    pub fn get_by_idx(&self, tx_idx: usize) -> Option<&Transcript> {
        self.all.get(tx_idx)
    }

    /// Resolve a transcript by stable id in O(1).
    pub fn get_by_id(&self, stable_id: &str) -> Option<&Transcript> {
        self.all.get(*self.by_id.get(stable_id)?)
    }

    /// Iterate every transcript in `tx_idx` order — backs `vep_transcripts()`.
    pub fn iter(&self) -> impl Iterator<Item = &Transcript> {
        self.all.iter()
    }

    /// Resolve a query chromosome to its stored `(transcripts, suffix_max_end)`
    /// slices, tolerating a `chr` prefix mismatch between the VCF and the gene
    /// model — VEP normalizes `chr1`↔`1` and `chrM`↔`MT`. Without this, a
    /// `chr`-prefixed VCF against an Ensembl cache (`1`,`2`,…,`MT`) annotates
    /// **nothing**. The exact match is tried first and is allocation-free; the
    /// fallbacks only run on a miss.
    fn buckets(&self, chrom: &str) -> Option<(&[Transcript], &[u64])> {
        let get = |k: &str| -> Option<(&[Transcript], &[u64])> {
            let &(lo, hi) = self.by_chrom.get(k)?;
            Some((&self.all[lo..hi], &self.sme[lo..hi]))
        };
        if let Some(x) = get(chrom) {
            return Some(x); // fast path: exact match, no allocation
        }
        match chrom.strip_prefix("chr") {
            // `chr1` → `1`, `chrM`/`chrMT` → `MT` (both alloc-free `&str` lookups)
            Some(core) => get(core).or_else(|| match core {
                "M" | "MT" => get("MT"),
                _ => None,
            }),
            // `1` → `chr1`, `MT`/`M` → `chrM` (allocates only on this miss path)
            None => get(&format!("chr{chrom}")).or_else(|| match chrom {
                "MT" | "M" => get("chrM").or_else(|| get("MT")),
                _ => None,
            }),
        }
    }
}

impl TranscriptProvider for IndexedTranscriptProvider {
    fn get_transcripts(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<&Transcript>> {
        let (trs, sme) = match self.buckets(chrom) {
            Some(x) => x,
            None => return Ok(Vec::new()),
        };

        // Binary search: find the first transcript whose start > end (query end).
        // All transcripts that could overlap must have start <= end, so they're in [0..upper).
        let upper = trs.partition_point(|t| t.start <= end);

        // From [0..upper), filter those whose end >= start (query start).
        // Use suffix_max_end for early termination: if the max end from index i
        // onwards is less than query start, no transcript at i or earlier can overlap.
        let mut results = Vec::new();
        for i in (0..upper).rev() {
            if sme[i] < start {
                break; // No transcript from [0..=i] can reach query start
            }
            if trs[i].end >= start {
                results.push(&trs[i]);
            }
        }
        results.reverse(); // Restore start-position order
        Ok(results)
    }

    fn get_transcripts_by_chrom(&self, chrom: &str) -> Result<Vec<&Transcript>> {
        match self.buckets(chrom) {
            Some((trs, _)) => Ok(trs.iter().collect()),
            None => Ok(Vec::new()),
        }
    }
}

/// Sequence provider backed by a FASTA reader.
pub struct FastaSequenceProvider {
    reader: crate::fasta::FastaReader,
}

impl FastaSequenceProvider {
    pub fn new(reader: crate::fasta::FastaReader) -> Self {
        Self { reader }
    }
}

impl SequenceProvider for FastaSequenceProvider {
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        self.reader.fetch(chrom, start, end)
    }
}

/// Sequence provider backed by a memory-mapped FASTA reader.
/// Uses .fai index for random access without loading the full file into RAM.
pub struct MmapFastaSequenceProvider {
    reader: crate::fasta::MmapFastaReader,
}

impl MmapFastaSequenceProvider {
    pub fn new(reader: crate::fasta::MmapFastaReader) -> Self {
        Self { reader }
    }
}

impl SequenceProvider for MmapFastaSequenceProvider {
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        self.reader.fetch(chrom, start, end)
    }
}

/// A matched known variant with its allele-specific frequency data.
#[derive(Debug, Clone)]
pub struct MatchedVariant {
    pub name: String,
    pub matched_allele: String,
    pub minor_allele: Option<String>,
    pub minor_allele_freq: Option<f64>,
    pub clin_sig: Option<String>,
    pub somatic: bool,
    pub phenotype_or_disease: bool,
    pub pubmed: Vec<String>,
    /// Population → allele-specific frequency for the matched allele.
    pub frequencies: HashMap<String, f64>,
}

/// Trait for providing co-located known variant annotations.
pub trait VariationProvider {
    /// Look up known variants overlapping a position that match the given alleles.
    fn get_matched_variants(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Vec<MatchedVariant>>;
}

/// Variation provider backed by VEP's tabix-indexed cache files.
pub struct TabixVariationProvider {
    reader: VariationTabixReader,
}

impl TabixVariationProvider {
    /// Create a provider from a VEP cache directory.
    ///
    /// Reads `info.txt` for column definitions and valid chromosomes.
    pub fn new(cache_dir: &std::path::Path, cache_info: &CacheInfo) -> Result<Self> {
        let reader = VariationTabixReader::new(
            cache_dir,
            &cache_info.variation_cols,
            &cache_info.valid_chromosomes,
        )?;
        Ok(Self { reader })
    }
}

impl VariationProvider for TabixVariationProvider {
    fn get_matched_variants(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Vec<MatchedVariant>> {
        let records = self.reader.query(chrom, start, end)?;
        let mut matched = Vec::new();

        for record in &records {
            // Skip failed variants
            if record.failed {
                continue;
            }

            // Check allele match
            if let Some(matched_alt) =
                variation::match_alleles(ref_allele, alt_allele, start, end, record)
            {
                // Extract per-allele frequencies for the matched allele
                let mut freqs = HashMap::new();
                for (pop, freq_str) in &record.frequencies {
                    if let Some(f) = variation::get_allele_freq(freq_str, &matched_alt) {
                        freqs.insert(pop.clone(), f);
                    }
                }

                // Also include MAF if minor_allele matches
                if let (Some(ref ma), Some(maf)) = (&record.minor_allele, record.minor_allele_freq)
                {
                    if ma.eq_ignore_ascii_case(&matched_alt) {
                        freqs.entry("minor_allele_freq".into()).or_insert(maf);
                    }
                }

                matched.push(MatchedVariant {
                    name: record.variation_name.clone(),
                    matched_allele: matched_alt,
                    minor_allele: record.minor_allele.clone(),
                    minor_allele_freq: record.minor_allele_freq,
                    clin_sig: record.clin_sig.clone(),
                    somatic: record.somatic,
                    phenotype_or_disease: record.phenotype_or_disease,
                    pubmed: record.pubmed.clone(),
                    frequencies: freqs,
                });
            }
        }

        Ok(matched)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_core::Strand;
    use fastvep_genome::{Exon, Gene};

    fn make_transcript(chrom: &str, start: u64, end: u64) -> Transcript {
        Transcript {
            stable_id: Arc::from(format!("ENST_{}", start).as_str()),
            version: None,
            gene: Gene {
                stable_id: "ENSG_1".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: chrom.into(),
                start,
                end,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: chrom.into(),
            start,
            end,
            strand: Strand::Forward,
            exons: vec![Exon {
                stable_id: "ENSE_1".into(),
                start,
                end,
                strand: Strand::Forward,
                phase: 0,
                end_phase: 0,
                rank: 1,
            }],
            translation: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        }
    }

    #[test]
    fn test_memory_transcript_provider() {
        let provider = MemoryTranscriptProvider::new(vec![
            make_transcript("chr1", 1000, 2000),
            make_transcript("chr1", 3000, 4000),
            make_transcript("chr2", 1000, 2000),
        ]);

        // Overlapping query
        let results = provider.get_transcripts("chr1", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 1000);

        // Non-overlapping
        let results = provider.get_transcripts("chr1", 2500, 2600).unwrap();
        assert_eq!(results.len(), 0);

        // Different chromosome
        let results = provider.get_transcripts("chr2", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);

        // By chromosome
        let results = provider.get_transcripts_by_chrom("chr1").unwrap();
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_indexed_transcript_provider() {
        let provider = IndexedTranscriptProvider::new(vec![
            make_transcript("chr1", 1000, 2000),
            make_transcript("chr1", 3000, 4000),
            make_transcript("chr2", 1000, 2000),
        ]);

        // Overlapping query
        let results = provider.get_transcripts("chr1", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 1000);

        // Non-overlapping
        let results = provider.get_transcripts("chr1", 2500, 2600).unwrap();
        assert_eq!(results.len(), 0);

        // Different chromosome
        let results = provider.get_transcripts("chr2", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);

        // By chromosome
        let results = provider.get_transcripts_by_chrom("chr1").unwrap();
        assert_eq!(results.len(), 2);

        // Missing chromosome
        let results = provider.get_transcripts("chr99", 1, 100).unwrap();
        assert_eq!(results.len(), 0);

        // Query spanning two transcripts
        let results = provider.get_transcripts("chr1", 1500, 3500).unwrap();
        assert_eq!(results.len(), 2);

        // Query at exact boundaries
        let results = provider.get_transcripts("chr1", 2000, 2000).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 1000);

        // Transcript count
        assert_eq!(provider.transcript_count(), 3);
    }

    #[test]
    fn test_tx_idx_ordinal_matches_chrom_start_order() {
        // The tx_idx ordinal must equal position in (chrom, start, stable_id)
        // order regardless of input order — this is what makes the on-disk
        // tx_idx column (same sort in tcache::save) agree with the engine, so
        // get_by_idx and get_by_id must be consistent.
        let provider = IndexedTranscriptProvider::new(vec![
            make_transcript("chr2", 5000, 6000),
            make_transcript("chr1", 3000, 4000),
            make_transcript("chr1", 1000, 2000),
        ]);
        // ordinal 0,1,2 = chr1@1000, chr1@3000, chr2@5000
        assert_eq!(provider.get_by_idx(0).unwrap().start, 1000);
        assert_eq!(&*provider.get_by_idx(0).unwrap().chromosome, "chr1");
        assert_eq!(provider.get_by_idx(1).unwrap().start, 3000);
        assert_eq!(provider.get_by_idx(2).unwrap().start, 5000);
        assert_eq!(&*provider.get_by_idx(2).unwrap().chromosome, "chr2");
        assert!(provider.get_by_idx(3).is_none());
        // get_by_id and get_by_idx resolve to the SAME transcript
        let by_id = provider.get_by_id("ENST_3000").unwrap();
        let by_idx = provider.get_by_idx(1).unwrap();
        assert_eq!(by_id.start, by_idx.start);
        assert_eq!(by_id.stable_id, by_idx.stable_id);
    }

    #[test]
    fn test_chr_prefix_and_mt_normalization() {
        // Regression: a `chr`-prefixed VCF against an Ensembl-style gene model
        // (`1`,`2`,…,`MT`) used to annotate NOTHING because lookups were exact.
        // VEP normalizes `chr1`↔`1` and `chrM`↔`MT`; so must we.
        let provider = IndexedTranscriptProvider::new(vec![
            make_transcript("1", 1000, 2000),  // Ensembl-style (no prefix)
            make_transcript("MT", 100, 500),   // mitochondrial, Ensembl spelling
        ]);

        // chr-prefixed query must find the non-prefixed model
        assert_eq!(provider.get_transcripts("chr1", 1500, 1600).unwrap().len(), 1);
        assert_eq!(provider.get_transcripts_by_chrom("chr1").unwrap().len(), 1);
        // exact (no-prefix) query still works
        assert_eq!(provider.get_transcripts("1", 1500, 1600).unwrap().len(), 1);
        // mitochondrial: chrM / M both resolve to the stored `MT`
        assert_eq!(provider.get_transcripts("chrM", 200, 300).unwrap().len(), 1);
        assert_eq!(provider.get_transcripts("M", 200, 300).unwrap().len(), 1);
        assert_eq!(provider.get_transcripts("MT", 200, 300).unwrap().len(), 1);
        // a genuinely absent contig still yields nothing
        assert_eq!(provider.get_transcripts("chr99", 1, 100).unwrap().len(), 0);

        // The reverse direction: model stored WITH `chr`, query without it.
        let pfx = IndexedTranscriptProvider::new(vec![
            make_transcript("chr7", 1000, 2000),
            make_transcript("chrM", 100, 500),
        ]);
        assert_eq!(pfx.get_transcripts("7", 1500, 1600).unwrap().len(), 1);
        assert_eq!(pfx.get_transcripts("MT", 200, 300).unwrap().len(), 1);
    }

    #[test]
    fn test_indexed_provider_overlapping_transcripts() {
        // Test with overlapping transcripts (common in real genomes)
        let provider = IndexedTranscriptProvider::new(vec![
            make_transcript("chr1", 1000, 5000),
            make_transcript("chr1", 2000, 3000),
            make_transcript("chr1", 4000, 6000),
        ]);

        // Query in the overlap region
        let results = provider.get_transcripts("chr1", 2500, 2600).unwrap();
        assert_eq!(results.len(), 2); // Both 1000-5000 and 2000-3000

        // Query spanning all three
        let results = provider.get_transcripts("chr1", 1000, 6000).unwrap();
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_fasta_sequence_provider() {
        let fasta = ">chr1\nACGTACGTAAAACCCC\n";
        let reader = crate::fasta::FastaReader::from_reader(fasta.as_bytes()).unwrap();
        let provider = FastaSequenceProvider::new(reader);
        let seq = provider.fetch_sequence("chr1", 1, 4).unwrap();
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn test_fasta_fetch_slice_matches_fetch() {
        let fasta = ">chr1\nacgtACGTaaaa\n>chr2\nTTTTgggg\n";
        let reader = crate::fasta::FastaReader::from_reader(fasta.as_bytes()).unwrap();
        let provider = FastaSequenceProvider::new(reader);

        // fetch_slice returns same data as fetch (both uppercase)
        let slice = provider.fetch_sequence_slice("chr1", 1, 4).unwrap();
        let vec = provider.fetch_sequence("chr1", 1, 4).unwrap();
        assert_eq!(slice, vec.as_slice());
        assert_eq!(slice, b"ACGT");

        // Lowercase input is uppercased at load time
        let slice = provider.fetch_sequence_slice("chr2", 5, 8).unwrap();
        assert_eq!(slice, b"GGGG");
    }
}
