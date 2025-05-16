use std::{collections::HashMap, path::PathBuf};
use std::fs::File;
use std::io::Write;
use std::io::BufReader;

use clap::Parser;
use noodles_bam as bam;
use noodles_core::Region;
use noodles_core::Position;
use nanocov::parse_bed;

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
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
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
    // This returns a noodles_bam::io::Reader, which wraps a BGZF-compressed BAM file
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_path)?;
    // Read the BAM header (noodles_sam::Header)
    let header = reader.read_header()?;
    // Get reference sequence dictionary (HashMap<String, ReferenceSequence>)
    let reference_sequences = header.reference_sequences();
    // Collect reference names as Vec<String>
    let ref_names: Vec<_> = reference_sequences.keys().cloned().collect();

    // Parse BED file if provided; returns Option<HashMap<String, Vec<(u32, u32)>>> mapping chrom to regions
    let bed_regions = if let Some(bed_path) = &cli.bed {
        Some(parse_bed(bed_path)?)
    } else {
        None
    };

    // Determine number of threads to use for chunked processing
    let n_threads = cli.threads.unwrap_or_else(|| std::cmp::max(1, num_cpus::get() / 2));

    // Shared state for coverage and per-chromosome averages, protected by Mutex for thread safety
    use std::sync::{Arc, Mutex};
    let coverage = Arc::new(Mutex::new(HashMap::new()));
    let per_chrom_averages = Arc::new(Mutex::new(Vec::new()));

    // Use index-backed region queries if BED is provided
    let mut records = Vec::new();
    if let Some(ref bed) = bed_regions {
        // Open BAM file as std::fs::File and wrap in BufReader for noodles_bam::io::Reader
        let bam_file = File::open(bam_path)?;
        let mut reader = bam::io::Reader::new(BufReader::new(bam_file));
        // Open BAI index file using noodles_bam::bai::Reader
        let bai_file = File::open(&bai_path)?;
        let mut bai_reader = noodles_bam::bai::Reader::new(bai_file);
        // Load the BAI index using noodles_bam::bai::Reader::read_index()
        let index = bai_reader.read_index()?;
        // Read the BAM header again for this reader
        let header = reader.read_header()?;
        let reference_sequences = header.reference_sequences();
        let ref_names: Vec<_> = reference_sequences.keys().cloned().collect();
        // Iterate over each chromosome and its regions from the BED file
        for (chrom, regions) in bed.iter() {
            // Find the reference sequence index for this chromosome
            let ref_id = ref_names.iter().position(|r| r == chrom);
            if ref_id.is_none() {
                continue;
            }
            for &(start, end) in regions {
                // BED is 0-based, half-open; BAM is 1-based, closed
                // noodles_core::Region is constructed with a chromosome name and a RangeInclusive<Position>
                // Position is a 1-based coordinate (usize)
                let region = Region::new(
                    chrom.clone(),
                    Position::try_from((start + 1) as usize)?..=Position::try_from(end as usize)?
                );
                // Query the BAM file for records overlapping this region using the index
                // reader.query returns an iterator over Result<Record, _>
                let query = reader.query(&header, &index, &region)?;
                for result in query {
                    let record = result?; // noodles_bam::Record
                    records.push(record);
                }
            }
        }
    } else {
        // No BED: read all records sequentially from the BAM file
        // reader.records() returns an iterator over Result<Record, _>
        records = reader.records().collect::<Result<_, _>>()?;
    }

    // Calculate per-chromosome coverage and average, then global average
    // Split records into chunks for parallel processing using tokio tasks
    let mut handles = Vec::new();
    for chunk in records.chunks(n_threads) {
        let chunk = chunk.to_vec();
        let ref_names = ref_names.clone();
        let bed_regions = bed_regions.clone();
        // Spawn a tokio task for each chunk
        let handle = tokio::spawn(async move {
            // local_coverage: HashMap<String, HashMap<u32, u32>>
            //   maps chromosome -> position -> count
            let mut local_coverage: HashMap<String, HashMap<u32, u32>> = HashMap::new();
            // local_averages: HashMap<String, f64> for per-chromosome averages
            let mut local_averages: HashMap<String, f64> = HashMap::new();
            for record in chunk {
                // Skip unmapped records
                if record.flags().is_unmapped() {
                    continue;
                }
                // Get reference sequence index for this record
                let ref_id = match record.reference_sequence_id() {
                    Some(Ok(id)) => id,
                    _ => continue,
                };
                // Convert reference index to chromosome name (String)
                let ref_name = ref_names.get(ref_id)
                    .map(|bstr| String::from_utf8_lossy(&bstr.to_vec()).into_owned())
                    .unwrap_or_else(|| "unknown".to_string());
                // If BED regions are present, check if this record overlaps any region
                if let Some(ref regions) = bed_regions {
                    if let Some(region_list) = regions.get(&ref_name) {
                        let start = match record.alignment_start() {
                            Some(Ok(pos)) => pos.get() as u32,
                            _ => continue,
                        };
                        let len = record.cigar().len() as u32;
                        let end = start + len;
                        // Only count records overlapping any region in the BED
                        if !region_list.iter().any(|&(r_start, r_end)| start < r_end && end > r_start) {
                            continue;
                        }
                    } else {
                        continue;
                    }
                }
                // Calculate per-base coverage for this record
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
            // Calculate per-chromosome averages for this chunk
            for (ref_name, region_coverage) in &local_coverage {
                let (total, count) = region_coverage.values().fold((0u64, 0u64), |(t, c), v| (t + *v as u64, c + 1));
                if count > 0 {
                    let avg = total as f64 / count as f64;
                    local_averages.insert(ref_name.clone(), avg);
                }
            }
            (local_coverage, local_averages)
        });
        handles.push(handle);
    }
    // Await all jobs and merge results
    for handle in handles {
        let (local_coverage, local_averages) = handle.await?;
        {
            let mut cov = coverage.lock().unwrap();
            let cov: &mut HashMap<String, HashMap<u32, u32>> = &mut *cov;
            for (ref_name, region_coverage) in local_coverage {
                let entry: &mut HashMap<u32, u32> = cov.entry(ref_name).or_default();
                for (pos, count) in region_coverage {
                    *entry.entry(pos).or_insert(0) += count;
                }
            }
        }
        {
            let mut avgs = per_chrom_averages.lock().unwrap();
            for (_ref_name, avg) in local_averages {
                avgs.push(avg);
            }
        }
    }
    let coverage = Arc::try_unwrap(coverage).unwrap().into_inner().unwrap();
    let per_chrom_averages = Arc::try_unwrap(per_chrom_averages).unwrap().into_inner().unwrap();

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
