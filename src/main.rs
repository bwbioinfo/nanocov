use std::{collections::HashMap, path::PathBuf};
use std::fs::File;
use std::io::Write;
use std::io::BufReader;

use clap::Parser;
use noodles_bam as bam;
use noodles_core::Region;
use noodles_core::Position;
use nanocov::parse_bed;
use rayon::prelude::*;

/// Simple BAM coverage calculator
#[derive(Parser, Debug)]
#[command(name = "bam-coverage")]
#[command(about = "Calculates per-base coverage from a BAM file", long_about = None)]
struct Cli {
    /// Input BAM file
    #[arg(short, long)]
    input: PathBuf,

    /// BED file with regions to include (chrom, start, end)
    #[arg(short = 'b', long = "bed")]
    bed: Option<PathBuf>,

    /// Number of threads to use (default: half of available cores)
    #[arg(short = 't', long = "threads")]
    threads: Option<usize>,

    /// Output file path (default: coverage.tsv)
    #[arg(short = 'o', long = "output", default_value = "coverage.tsv")]
    output: PathBuf,

    /// Chunk size for parallel processing (default: 10,000)
    #[arg(short = 'c', long = "chunk-size", default_value_t = 10_000)]
    chunk_size: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments using clap's derive API
    let cli = Cli::parse();

    // Prepare BAM and BAI file paths
    let bam_path = &cli.input;
    let bai_path = bam_path.with_extension("bam.bai");
    // Check for BAM index file (.bai); if missing, instruct user to create it
    if !bai_path.exists() {
        eprintln!("BAM index not found at {:?}. Please run 'samtools index {:?}' to create it.", bai_path, bam_path);
        std::process::exit(1);
    }
    // Open BAM file using noodles_bam::io::reader::Builder
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_path)?;
    let header = reader.read_header()?;
    let reference_sequences = header.reference_sequences();

    // Parse BED file if provided; returns Option<HashMap<String, Vec<(u32, u32)>>> mapping chrom to regions
    let bed_regions = if let Some(bed_path) = &cli.bed {
        Some(parse_bed(bed_path)?)
    } else {
        None
    };

    let coverage_results: Vec<(HashMap<String, HashMap<u32, u32>>, HashMap<String, f64>)> = if let Some(ref bed) = bed_regions {
        // Collect all (chrom, region) pairs
        let mut region_jobs = Vec::new();
        for (chrom, regions) in bed.iter() {
            for &(start, end) in regions {
                region_jobs.push((chrom.clone(), start, end));
            }
        }
        // Parallelize over regions
        region_jobs.par_iter().map(|(chrom, start, end)| {
            let mut reader = bam::io::Reader::new(BufReader::new(File::open(bam_path).unwrap()));
            let mut bai_reader = noodles_bam::bai::Reader::new(File::open(&bai_path).unwrap());
            let index = bai_reader.read_index().unwrap();
            let header = reader.read_header().unwrap();
            let reference_sequences = header.reference_sequences();
            let ref_names: Vec<String> = reference_sequences.keys().map(|b| b.to_string()).collect();
            let region = Region::new(
                chrom.clone(),
                Position::try_from((*start + 1) as usize).unwrap()..=Position::try_from(*end as usize).unwrap()
            );
            let mut local_coverage: HashMap<String, HashMap<u32, u32>> = HashMap::new();
            let mut local_averages: HashMap<String, f64> = HashMap::new();
            let query = reader.query(&header, &index, &region).unwrap();
            for result in query {
                let record = result.unwrap();
                if record.flags().is_unmapped() {
                    continue;
                }
                let ref_id = match record.reference_sequence_id() {
                    Some(Ok(id)) => id,
                    _ => continue,
                };
                let ref_name = ref_names.get(ref_id)
                    .map(|name| name.clone())
                    .unwrap_or_else(|| "unknown".to_string());
                let start = match record.alignment_start() {
                    Some(Ok(pos)) => pos.get() as u32,
                    _ => continue,
                };
                let len = record.cigar().len() as u32;
                let region_coverage = local_coverage.entry(ref_name.clone()).or_default();
                for pos in start..start + len {
                    *region_coverage.entry(pos).or_insert(0) += 1;
                }
            }
            for (ref_name, region_coverage) in &local_coverage {
                let (total, count) = region_coverage.values().fold((0u64, 0u64), |(t, c), v| (t + *v as u64, c + 1));
                if count > 0 {
                    let avg = total as f64 / count as f64;
                    local_averages.insert(ref_name.clone(), avg);
                }
            }
            (local_coverage, local_averages)
        }).collect()
    } else {
        // No BED: parallelize by chromosome
        let ref_names: Vec<String> = reference_sequences.keys().map(|b| b.to_string()).collect();
        let chrom_jobs: Vec<String> = ref_names.clone();
        chrom_jobs.par_iter().map(|chrom| {
            let mut reader = bam::io::Reader::new(BufReader::new(File::open(bam_path).unwrap()));
            let mut bai_reader = noodles_bam::bai::Reader::new(File::open(&bai_path).unwrap());
            let index = bai_reader.read_index().unwrap();
            let header = reader.read_header().unwrap();
            let reference_sequences = header.reference_sequences();
            // Find the reference sequence for this chromosome
            let ref_seq = reference_sequences.iter().find(|(k, _)| k.to_string() == *chrom).map(|(_, v)| v).unwrap();
            let len = ref_seq.length();
            let region = Region::new(
                chrom.clone(),
                Position::try_from(1usize).unwrap()..=Position::try_from(len.get()).unwrap()
            );
            let mut local_coverage: HashMap<String, HashMap<u32, u32>> = HashMap::new();
            let mut local_averages: HashMap<String, f64> = HashMap::new();
            let query = reader.query(&header, &index, &region).unwrap();
            for result in query {
                let record = result.unwrap();
                if record.flags().is_unmapped() {
                    continue;
                }
                let ref_id = match record.reference_sequence_id() {
                    Some(Ok(id)) => id,
                    _ => continue,
                };
                let ref_name = ref_names.get(ref_id)
                    .map(|name| name.clone())
                    .unwrap_or_else(|| "unknown".to_string());
                let start = match record.alignment_start() {
                    Some(Ok(pos)) => pos.get() as u32,
                    _ => continue,
                };
                let len = record.cigar().len() as u32;
                let region_coverage = local_coverage.entry(ref_name.clone()).or_default();
                for pos in start..start + len {
                    *region_coverage.entry(pos).or_insert(0) += 1;
                }
            }
            for (ref_name, region_coverage) in &local_coverage {
                let (total, count) = region_coverage.values().fold((0u64, 0u64), |(t, c), v| (t + *v as u64, c + 1));
                if count > 0 {
                    let avg = total as f64 / count as f64;
                    local_averages.insert(ref_name.clone(), avg);
                }
            }
            (local_coverage, local_averages)
        }).collect()
    };

    // Merge results from all jobs
    let mut coverage: HashMap<String, HashMap<u32, u32>> = HashMap::new();
    let mut per_chrom_averages: Vec<f64> = Vec::new();
    for (local_coverage, local_averages) in coverage_results {
        for (ref_name, region_coverage) in local_coverage {
            let entry = coverage.entry(ref_name).or_default();
            for (pos, count) in region_coverage {
                *entry.entry(pos).or_insert(0) += count;
            }
        }
        for (_ref_name, avg) in local_averages {
            per_chrom_averages.push(avg);
        }
    }

    // Write per-position coverage to a file (TSV format)
    let mut out = File::create(&cli.output)?;
    writeln!(out, "#chromosome\tposition\tcount")?;
    for (ref_name, region_coverage) in &coverage {
        let mut positions: Vec<_> = region_coverage.iter().collect();
        positions.sort_by_key(|&(pos, _)| *pos);
        for (pos, count) in positions {
            writeln!(out, "{}\t{}\t{}", ref_name, pos, count)?;
        }
    }

    // Print per-chromosome averages and global average to stdout
    if !per_chrom_averages.is_empty() {
        for (ref_name, region_coverage) in &coverage {
            let (total, count) = region_coverage.values().fold((0u64, 0u64), |(t, c), v| (t + *v as u64, c + 1));
            if count > 0 {
                let avg = total as f64 / count as f64;
                println!("{} average coverage: {:.2}", ref_name, avg);
            }
        }
        let global_avg = per_chrom_averages.iter().sum::<f64>() / per_chrom_averages.len() as f64;
        println!("Global average coverage: {:.2}", global_avg);
    } else {
        println!("No coverage data found.");
    }
    Ok(())
}
