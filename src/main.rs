use std::{collections::HashMap, path::PathBuf};
use std::fs::File;
use std::io::Write;

use clap::Parser;
use noodles_bam as bam;
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
    // Determine number of threads
    let n_threads = cli.threads.unwrap_or_else(|| std::cmp::max(1, num_cpus::get() / 2));

    // Calculate per-chromosome coverage and average, then global average
    use std::sync::{Arc, Mutex};
    let coverage = Arc::new(Mutex::new(HashMap::new()));
    let per_chrom_averages = Arc::new(Mutex::new(Vec::new()));

    let mut handles = Vec::new();
    for chunk in records.chunks(n_threads) {
        let chunk = chunk.to_vec();
        let ref_names = ref_names.clone();
        let bed_regions = bed_regions.clone();
        let handle = tokio::spawn(async move {
            let mut local_coverage: HashMap<String, HashMap<u32, u32>> = HashMap::new();
            let mut local_averages: HashMap<String, f64> = HashMap::new();
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

    // Write per-position coverage to a file
    let mut out = File::create(&cli.output)?;
    writeln!(out, "#chromosome\tposition\tcount")?;
    for (ref_name, region_coverage) in &coverage {
        let mut positions: Vec<_> = region_coverage.iter().collect();
        positions.sort_by_key(|&(pos, _)| *pos);
        for (pos, count) in positions {
            writeln!(out, "{}\t{}\t{}", ref_name, pos, count)?;
        }
    }

    // Print per-chromosome averages and global average
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
