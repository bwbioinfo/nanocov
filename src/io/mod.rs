// src/io/mod.rs
// IO module for nanocov: BAM/BED reading, coverage writing

pub mod cramino;

use crate::cli::Cli;

// Functions for BAM/BED reading and coverage writing will be moved here from main.rs
// (Implementations will be moved in the next step)

use crate::utils::ReadStats;

pub fn run_coverage(cli: &Cli, read_stats: Option<ReadStats>) -> Result<(), Box<dyn std::error::Error>> {
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::Write;
    use std::io::BufReader;
    use noodles_bam as bam;
    use noodles_core::Region;
    use noodles_core::Position;
    use nanocov::parse_bed;
    use rayon::prelude::*;

    // Open BAM file using noodles_bam::io::reader::Builder
    let mut reader = bam::io::reader::Builder::default().build_from_path(&cli.input)?;
    let header = reader.read_header()?;
    let reference_sequences = header.reference_sequences();

    // Parse BED file if provided; returns Option<HashMap<String, Vec<(u32, u32)>>> mapping chrom to regions
    let bed_regions = if let Some(bed_path) = &cli.bed {
        Some(parse_bed(bed_path)?)
    } else {
        None
    };

    // Parse chrom_bed file if provided; returns Option<HashMap<String, Vec<(u32, u32)>>> mapping chrom to full ranges
    let chrom_bed_regions = if let Some(chrom_bed_path) = &cli.chrom_bed {
        Some(parse_bed(chrom_bed_path)?)
    } else {
        None
    };

    let coverage_results: Vec<(HashMap<String, HashMap<u32, u32>>, HashMap<String, f64>)> = 
        if let Some(ref bed) = bed_regions {
            // Collect all (chrom, region) pairs
            let mut region_jobs = Vec::new();
            for (chrom, regions) in bed.iter() {
                for &(start, end) in regions {
                    region_jobs.push((chrom.clone(), start, end));
                }
            }
            // Parallelize over regions
            region_jobs.par_iter().map(|(chrom, start, end)| {
                let mut reader = bam::io::Reader::new(BufReader::new(File::open(&cli.input).unwrap()));
                let mut bai_reader = noodles_bam::bai::Reader::new(File::open(&cli.input.with_extension("bam.bai")).unwrap());
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
        } else if let Some(ref chrom_bed) = chrom_bed_regions {
            // Use full chromosome ranges from chrom_bed
            let mut chrom_jobs = Vec::new();
            for (chrom, regions) in chrom_bed.iter() {
                for &(start, end) in regions {
                    chrom_jobs.push((chrom.clone(), start, end));
                }
            }
            chrom_jobs.par_iter().map(|(chrom, start, end)| {
                let mut reader = bam::io::Reader::new(BufReader::new(File::open(&cli.input).unwrap()));
                let mut bai_reader = noodles_bam::bai::Reader::new(File::open(&cli.input.with_extension("bam.bai")).unwrap());
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
                let mut reader = bam::io::Reader::new(BufReader::new(File::open(&cli.input).unwrap()));
                let mut bai_reader = noodles_bam::bai::Reader::new(File::open(&cli.input.with_extension("bam.bai")).unwrap());
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

    // Write per-position coverage to a file (TSV format) using BufWriter and parallel formatting
    use std::io::BufWriter;
    let mut out = BufWriter::new(File::create(&cli.output)?);
    writeln!(out, "#chromosome\tposition\tcount")?;
    // Parallelize formatting of coverage lines per chromosome
    let chrom_blocks: Vec<(String, String)> = coverage.par_iter()
        .map(|(ref_name, region_coverage)| {
            let mut positions: Vec<_> = region_coverage.iter().collect();
            positions.sort_by_key(|&(pos, _)| *pos);
            let mut block = String::with_capacity(positions.len() * 24); // estimate
            for (pos, count) in positions {
                use std::fmt::Write as _;
                let _ = write!(block, "{}\t{}\t{}\n", ref_name, pos, count);
            }
            (ref_name.clone(), block)
        })
        .collect();
    let mut chrom_blocks = chrom_blocks;
    chrom_blocks.sort_by(|a, b| a.0.cmp(&b.0));
    for (_ref_name, block) in chrom_blocks {
        out.write_all(block.as_bytes())?;
    }
    out.flush()?;

    // Print per-chromosome averages and global average to stdout
    if !per_chrom_averages.is_empty() {
        // Create a HashMap for multi-chromosome plotting
        let mut chrom_coverages: HashMap<String, &HashMap<u32, u32>> = HashMap::new();
        
        for (ref_name, region_coverage) in &coverage {
            let (total, count) = region_coverage.values().fold((0u64, 0u64), |(t, c), v| (t + *v as u64, c + 1));
            if count > 0 {
                let avg = total as f64 / count as f64;
                println!("{} average coverage: {:.2}", ref_name, avg);
                chrom_coverages.insert(ref_name.clone(), region_coverage);
            }
            // Call plotting for each chromosome
            let output_stem = cli.output.file_stem().unwrap_or_default().to_string_lossy();
            let output_dir = cli.output.parent().unwrap_or_else(|| std::path::Path::new("."));
            
            // Determine output format - support both PNG and SVG
            let file_format = if cli.svg_output {
                "svg"
            } else {
                "png"
            };
            
            let plot_path = output_dir.join(format!("{}.{}.{}", output_stem, ref_name, file_format));
            
            // Determine plot range: chrom_bed > bed > coverage
            let (plot_start, plot_end) = if let Some(ref chrom_bed) = chrom_bed_regions {
                if let Some(regions) = chrom_bed.get(ref_name) {
                    let min_start = regions.iter().map(|(s, _)| *s).min().unwrap_or(0);
                    let max_end = regions.iter().map(|(_, e)| *e).max().unwrap_or(0);
                    (min_start, max_end)
                } else {
                    (0, 0)
                }
            } else if let Some(ref bed) = bed_regions {
                if let Some(regions) = bed.get(ref_name) {
                    let min_start = regions.iter().map(|(s, _)| *s).min().unwrap_or(0);
                    let max_end = regions.iter().map(|(_, e)| *e).max().unwrap_or(0);
                    (min_start, max_end)
                } else {
                    (0, 0)
                }
            } else {
                // fallback to coverage min/max if no BED
                let min_pos = region_coverage.keys().min().copied().unwrap_or(0);
                let max_pos = region_coverage.keys().max().copied().unwrap_or(0);
                (min_pos, max_pos)
            };
            
            // Apply theme if specified
            if let Some(theme) = &cli.theme {
                crate::plotting::set_theme(theme);
            }
            
            crate::plotting::plot_per_base_coverage_with_range(
                ref_name,
                region_coverage,
                plot_path.to_str().unwrap(),
                plot_start,
                plot_end,
                read_stats.as_ref(),
                cli.show_zero_regions,
            )?;
        }
        let global_avg = per_chrom_averages.iter().sum::<f64>() / per_chrom_averages.len() as f64;
        println!("Global average coverage: {:.2}", global_avg);
        
        // Generate multi-chromosome plot if we have data from multiple chromosomes
        if chrom_coverages.len() > 1 && !cli.skip_plotting && !cli.skip_multi_plot {
            let output_stem = cli.output.file_stem().unwrap_or_default().to_string_lossy();
            let output_dir = cli.output.parent().unwrap_or_else(|| std::path::Path::new("."));
            
            // Get current theme
            let theme = unsafe { crate::plotting::CURRENT_THEME };
            
            // Plot with appropriate scale based on CLI option
            let plot_path = output_dir.join(format!("{}.multi_chrom.png", output_stem));
            crate::plotting::plot_all_chromosomes(
                &chrom_coverages,
                plot_path.to_str().unwrap(),
                cli.log_scale, // Use log scale if requested
                read_stats.as_ref(),
                theme,
            )?;
            
            println!("Generated multi-chromosome plot ({}): {}", 
                if cli.log_scale { "log scale" } else { "linear scale" }, 
                plot_path.display());
        }
    } else {
        println!("No coverage data found.");
    }

    // DEBUG: Print chromosome names and coverage counts
    eprintln!("[DEBUG] Chromosomes in coverage: {:?}", coverage.keys().collect::<Vec<_>>());
    for chrom in coverage.keys() {
        eprintln!("[DEBUG] {}: {} positions", chrom, coverage[chrom].len());
    }

    Ok(())
}
