// src/utils/mod.rs

use noodles_bam as bam;
use std::fs::File;
use std::io::BufReader;

pub struct ReadStats {
    pub n50: u32,
    pub mean_len: f64,
    pub median_len: f64,
    pub mean_qual: f64,
    pub median_qual: f64,
    pub num_reads: u64,
    pub num_bases: u64,
    pub lengths: Option<Vec<u32>>,
}

pub fn extract_read_stats(bam_path: &std::path::Path) -> Result<ReadStats, Box<dyn std::error::Error>> {
    let mut reader = bam::io::Reader::new(BufReader::new(File::open(bam_path)?));
    let _header = reader.read_header()?;
    let mut lengths = Vec::new();
    let mut quals = Vec::new();
    let mut num_reads: u64 = 0;
    let mut num_bases: u64 = 0;

    for result in reader.records() {
        let record = result?;
        let len = record.sequence().len() as u32;
        lengths.push(len);
        num_reads += 1;
        num_bases += len as u64;

        // Collect mean quality per read (if available)
        let qual = record.quality_scores();
        let qual_slice = qual.as_ref();
        if !qual_slice.is_empty() {
            let q: f64 = qual_slice.iter().map(|&q| q as u32).sum::<u32>() as f64 / qual_slice.len() as f64;
            quals.push(q);
        }
    }

    // N50 calculation
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let total: u64 = lengths.iter().map(|&l| l as u64).sum();
    let mut acc = 0;
    let mut n50 = 0;
    for &l in &lengths {
        acc += l as u64;
        if acc >= total / 2 {
            n50 = l;
            break;
        }
    }

    // Mean/median length
    let mean_len = if lengths.is_empty() { 0.0 } else { lengths.iter().sum::<u32>() as f64 / lengths.len() as f64 };
    let median_len = if lengths.is_empty() {
        0.0
    } else {
        let mid = lengths.len() / 2;
        if lengths.len() % 2 == 0 {
            (lengths[mid - 1] + lengths[mid]) as f64 / 2.0
        } else {
            lengths[mid] as f64
        }
    };

    // Mean/median quality
    let mean_qual = if quals.is_empty() { 0.0 } else { quals.iter().sum::<f64>() / quals.len() as f64 };
    let median_qual = if quals.is_empty() {
        0.0
    } else {
        let mut sorted = quals.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mid = sorted.len() / 2;
        if sorted.len() % 2 == 0 {
            (sorted[mid - 1] + sorted[mid]) / 2.0
        } else {
            sorted[mid]
        }
    };

    Ok(ReadStats {
        n50,
        mean_len,
        median_len,
        mean_qual,
        median_qual,
        num_reads,
        num_bases,
        lengths: Some(lengths),
    })
}
// Utility functions for nanocov

// Place for helpers, region job creation, merging, etc.
// (Implementations will be moved in the next step)
