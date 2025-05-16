use std::{collections::HashMap, path::PathBuf, io::BufRead};
use std::fs::File;
use std::io::Write;

use clap::Parser;
use noodles_bam as bam;
use tokio::task;
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
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let mut reader = bam::io::reader::Builder::default().build_from_path(&cli.input)?;
    let header = reader.read_header()?;
    let reference_sequences = header.reference_sequences();
    let ref_names: Vec<_> = reference_sequences.keys().cloned().collect();

    let bed_regions = if let Some(bed_path) = &cli.bed {
        Some(parse_bed(bed_path)?)
    } else {
        None
    };

    let records: Vec<_> = reader.records().collect::<Result<_, _>>()?;
    let n_threads = num_cpus::get();
    let chunk_size = (records.len() + n_threads - 1) / n_threads;
    let coverage_results = futures::future::join_all(
        records
            .chunks(chunk_size)
            .map(|chunk| {
                let chunk = chunk.to_vec();
                let ref_names = ref_names.clone();
                let bed_regions = bed_regions.clone();
                task::spawn_blocking(move || {
                    let mut coverage: HashMap<String, HashMap<u32, u32>> = HashMap::new();
                    for record in chunk {
                        if record.flags().is_unmapped() {
                            continue;
                        }
                        let ref_id = match record.reference_sequence_id() {
                            Some(Ok(id)) => id,
                            _ => continue,
                        };
                        let ref_name = ref_names.get(ref_id)
                            .map(|bstr| String::from_utf8_lossy(&bstr.to_vec()).into_owned())
                            .unwrap_or_else(|| "unknown".to_string());
                        if let Some(ref regions) = bed_regions {
                            if let Some(region_list) = regions.get(&ref_name) {
                                let start = match record.alignment_start() {
                                    Some(Ok(pos)) => pos.get() as u32,
                                    _ => continue,
                                };
                                let len = record.cigar().len() as u32;
                                let end = start + len;
                                if !region_list.iter().any(|&(r_start, r_end)| start < r_end && end > r_start) {
                                    continue;
                                }
                            } else {
                                continue;
                            }
                        }
                        let start = match record.alignment_start() {
                            Some(Ok(pos)) => pos.get() as u32,
                            _ => continue,
                        };
                        let len = record.cigar().len() as u32;
                        let region_coverage = coverage.entry(ref_name).or_default();
                        for pos in start..start + len {
                            *region_coverage.entry(pos).or_insert(0) += 1;
                        }
                    }
                    coverage
                })
            })
    ).await;

    let mut coverage: HashMap<String, HashMap<u32, u32>> = HashMap::new();
    for result in coverage_results {
        let partial = result?;
        for (ref_name, region_coverage) in partial {
            let entry = coverage.entry(ref_name).or_default();
            for (pos, count) in region_coverage {
                *entry.entry(pos).or_insert(0) += count;
            }
        }
    }

    // Write per-position coverage to a file
    let mut out = File::create("coverage.tsv")?;
    for (ref_name, region_coverage) in &coverage {
        writeln!(out, "# {}", ref_name)?;
        let mut positions: Vec<_> = region_coverage.iter().collect();
        positions.sort_by_key(|&(pos, _)| *pos);
        for (pos, count) in positions {
            writeln!(out, "{}\t{}", pos, count)?;
        }
    }

    // Only print summary and average coverage to stdout
    let mut total_coverage = 0u64;
    let mut total_positions = 0u64;
    for (_ref_name, region_coverage) in &coverage {
        for (_pos, count) in region_coverage {
            total_coverage += *count as u64;
            total_positions += 1;
        }
    }
    if total_positions > 0 {
        let avg_coverage = total_coverage as f64 / total_positions as f64;
        println!("Average coverage: {:.2}", avg_coverage);
    } else {
        println!("No coverage data found.");
    }
    Ok(())
}
