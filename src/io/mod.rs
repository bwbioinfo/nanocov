// src/io/mod.rs
// IO module for nanocov: BAM/BED reading, coverage writing

pub mod cramino;

use crate::cli::Cli;

// Functions for BAM/BED reading and coverage writing will be moved here from main.rs
// (Implementations will be moved in the next step)

use crate::utils::ReadStats;

/// Calculate the reference span from a CIGAR string
/// This properly handles different CIGAR operations to get the actual alignment length on the reference
fn calculate_reference_span(cigar: &noodles_bam::record::Cigar) -> u32 {
    cigar.iter()
        .filter_map(|op_result| {
            if let Ok(op) = op_result {
                use noodles_sam::alignment::record::cigar::op::Kind;
                match op.kind() {
                    // Operations that consume reference sequence
                    Kind::Match | 
                    Kind::Deletion | 
                    Kind::Skip | 
                    Kind::SequenceMatch | 
                    Kind::SequenceMismatch => {
                        Some(op.len() as u32)
                    }
                    // Operations that don't consume reference: Insertion, SoftClip, HardClip, Pad
                    Kind::Insertion | 
                    Kind::SoftClip | 
                    Kind::HardClip | 
                    Kind::Pad => None,
                }
            } else {
                None // Skip invalid CIGAR operations
            }
        })
        .sum()
}

pub fn run_coverage(cli: &Cli, read_stats: Option<ReadStats>) -> Result<(), Box<dyn std::error::Error>> {
    // Choose the appropriate coverage calculation mode
    if should_use_streaming_mode(cli)? {
        // Use streaming mode for very large files
        run_coverage_streaming(cli, read_stats)
    } else {
        // Use enhanced parallel chunked mode for better performance
        run_coverage_parallel_chunked(cli, read_stats)
    }
}

/// Determine if we should use streaming mode based on file size and available memory
fn should_use_streaming_mode(cli: &Cli) -> Result<bool, Box<dyn std::error::Error>> {
    // Force streaming if explicitly requested
    if cli.force_streaming {
        return Ok(true);
    }
    
    let bam_metadata = std::fs::metadata(&cli.input)?;
    let file_size_mb = bam_metadata.len() / (1024 * 1024);
    
    // Use custom memory limit or default to 500MB
    let memory_limit_mb = cli.memory_limit_mb.unwrap_or(500);
    
    // Use streaming mode for files larger than the threshold
    Ok(file_size_mb > memory_limit_mb)
}

/// Memory-efficient streaming approach for large BAM files
fn run_coverage_streaming(cli: &Cli, read_stats: Option<ReadStats>) -> Result<(), Box<dyn std::error::Error>> {
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{Write, BufWriter};
    use noodles_bam as bam;

    println!("Using memory-efficient streaming mode for large BAM file");

    // Open BAM file
    let mut reader = bam::io::reader::Builder::default().build_from_path(&cli.input)?;
    let header = reader.read_header()?;
    let reference_sequences = header.reference_sequences();

    // Parse BED files
    let bed_regions = if let Some(bed_path) = &cli.bed {
        Some(nanocov::parse_bed(bed_path)?)
    } else {
        None
    };

    let chrom_bed_regions = if let Some(chrom_bed_path) = &cli.chrom_bed {
        Some(nanocov::parse_bed(chrom_bed_path)?)
    } else {
        None
    };

    // Create output file with buffered writer
    let mut out = BufWriter::new(File::create(&cli.output)?);
    writeln!(out, "#chromosome\tposition\tcount")?;

    let mut global_avg_sum = 0.0;
    let mut global_avg_count = 0;
    let mut chrom_coverages: HashMap<String, HashMap<u32, u32>> = HashMap::new();

    // Process each chromosome individually to save memory
    for (chrom_name, _) in reference_sequences.iter() {
        let chrom = chrom_name.to_string();
        println!("Processing chromosome: {}", chrom);

        let coverage = process_chromosome_streaming(&cli.input, &chrom, &bed_regions, &chrom_bed_regions)?;
        
        if !coverage.is_empty() {
            // Calculate chromosome average
            let (total, count) = coverage.values().fold((0u64, 0u64), |(t, c), v| (t + *v as u64, c + 1));
            if count > 0 {
                let avg = total as f64 / count as f64;
                println!("{} average coverage: {:.2}", chrom, avg);
                global_avg_sum += avg;
                global_avg_count += 1;

                // Store for plotting (keep minimal data)
                chrom_coverages.insert(chrom.clone(), coverage.clone());
            }

            // Write coverage data immediately to file
            write_chromosome_coverage(&mut out, &chrom, &coverage)?;
        }

        // Clear coverage data to free memory
        drop(coverage);
    }

    out.flush()?;

    // Print global average
    if global_avg_count > 0 {
        let global_avg = global_avg_sum / global_avg_count as f64;
        println!("Global average coverage: {:.2}", global_avg);
    }

    // Generate plots with reduced memory usage (if not disabled)
    if !cli.skip_all_plots {
        generate_plots_from_stored_coverage(cli, &chrom_coverages, read_stats.as_ref())?;
    } else {
        println!("Skipping plot generation as requested (--no-plots)");
    }

    Ok(())
}

/// Process a single chromosome and return its coverage data
fn process_chromosome_streaming(
    bam_path: &std::path::Path,
    chrom: &str,
    bed_regions: &Option<std::collections::HashMap<String, Vec<(u32, u32)>>>,
    chrom_bed_regions: &Option<std::collections::HashMap<String, Vec<(u32, u32)>>>,
) -> Result<std::collections::HashMap<u32, u32>, Box<dyn std::error::Error>> {
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::BufReader;
    use noodles_bam as bam;
    use noodles_core::Region;
    use noodles_core::Position;

    let mut reader = bam::io::Reader::new(BufReader::new(File::open(bam_path)?));
    let mut bai_reader = noodles_bam::bai::Reader::new(File::open(&bam_path.with_extension("bam.bai"))?);
    let index = bai_reader.read_index()?;
    let header = reader.read_header()?;
    let reference_sequences = header.reference_sequences();

    let mut coverage: HashMap<u32, u32> = HashMap::new();

    // Determine regions to process for this chromosome
    let regions_to_process = if let Some(bed) = bed_regions {
        bed.get(chrom).cloned().unwrap_or_default()
    } else if let Some(chrom_bed) = chrom_bed_regions {
        chrom_bed.get(chrom).cloned().unwrap_or_default()
    } else {
        // Process entire chromosome
        if let Some((_, ref_seq)) = reference_sequences.iter().find(|(k, _)| k.to_string() == chrom) {
            vec![(1, ref_seq.length().get() as u32)]
        } else {
            return Ok(coverage);
        }
    };

    // Process each region
    for (start, end) in regions_to_process {
        let region = Region::new(
            chrom.to_string(),
            Position::try_from(start as usize)
                .map_err(|e| format!("Invalid start position {}: {}", start, e))?
                ..=Position::try_from(end as usize)
                    .map_err(|e| format!("Invalid end position {}: {}", end, e))?,
        );

        let query = reader.query(&header, &index, &region)?;
        for result in query {
            let record = result?;
            if record.flags().is_unmapped() {
                continue;
            }

            let start_pos = match record.alignment_start() {
                Some(Ok(pos)) => pos.get() as u32,
                _ => continue,
            };

            let len = calculate_reference_span(&record.cigar());
            
            // Update coverage for this alignment
            for pos in start_pos..start_pos + len {
                if pos >= start && pos <= end {
                    *coverage.entry(pos).or_insert(0) += 1;
                }
            }
        }
    }

    Ok(coverage)
}

/// Write chromosome coverage data to output file
fn write_chromosome_coverage(
    out: &mut std::io::BufWriter<std::fs::File>,
    chrom: &str,
    coverage: &std::collections::HashMap<u32, u32>,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::io::Write;
    
    let mut positions: Vec<_> = coverage.iter().collect();
    positions.sort_by_key(|&(pos, _)| *pos);
    
    for (&pos, &count) in positions {
        writeln!(out, "{}\t{}\t{}", chrom, pos, count)?;
    }
    
    Ok(())
}

/// Generate plots with memory-efficient approach
fn generate_plots_from_stored_coverage(
    cli: &Cli,
    chrom_coverages: &std::collections::HashMap<String, std::collections::HashMap<u32, u32>>,
    read_stats: Option<&ReadStats>,
) -> Result<(), Box<dyn std::error::Error>> {
    let output_stem = cli.output.file_stem().unwrap_or_default().to_string_lossy();
    let output_dir = cli.output.parent().unwrap_or_else(|| std::path::Path::new("."));
    
    let file_format = if cli.svg_output {
        "svg"
    } else {
        "png"
    };

    // Apply theme if specified
    if let Some(theme) = &cli.theme {
        crate::plotting::set_theme(theme);
    }

    // Generate individual chromosome plots
    for (ref_name, region_coverage) in chrom_coverages {
        let plot_path = output_dir.join(format!("{}.{}.{}", output_stem, ref_name, file_format));
        
        // Use full chromosome range for plotting
        let min_pos = region_coverage.keys().min().copied().unwrap_or(0);
        let max_pos = region_coverage.keys().max().copied().unwrap_or(0);
        
        crate::plotting::plot_per_base_coverage_with_range(
            ref_name,
            region_coverage,
            plot_path.to_str().unwrap(),
            min_pos,
            max_pos,
            read_stats,
            cli.show_zero_regions,
            cli.log_scale,
        )?;
    }

    // Generate multi-chromosome plot if we have multiple chromosomes
    if chrom_coverages.len() > 1 {
        let multi_plot_path = output_dir.join(format!("{}.multi_chrom.{}", output_stem, file_format));
        
        // Convert to the expected format for plot_all_chromosomes
        let chrom_refs: std::collections::HashMap<String, &std::collections::HashMap<u32, u32>> = 
            chrom_coverages.iter().map(|(k, v)| (k.clone(), v)).collect();
        
        // Get current theme
        let theme = unsafe { crate::plotting::CURRENT_THEME };
        
        crate::plotting::plot_all_chromosomes(
            &chrom_refs,
            multi_plot_path.to_str().unwrap(),
            cli.log_scale,
            read_stats,
            theme,
        )?;
    }

    Ok(())
}

/// Original in-memory approach for smaller files
fn run_coverage_in_memory(cli: &Cli, read_stats: Option<ReadStats>) -> Result<(), Box<dyn std::error::Error>> {
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
                    let len = calculate_reference_span(&record.cigar());
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
                    let len = calculate_reference_span(&record.cigar());
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
                    let len = calculate_reference_span(&record.cigar());
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
                cli.log_scale,
            )?;
        }
        let global_avg = per_chrom_averages.iter().sum::<f64>() / per_chrom_averages.len() as f64;
        println!("Global average coverage: {:.2}", global_avg);
        
        // Generate multi-chromosome plot if we have data from multiple chromosomes
        if chrom_coverages.len() > 1 && !cli.skip_all_plots && !cli.skip_multi_plot {
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

use std::collections::HashMap;
use rayon::prelude::*;

/// Enhanced parallel coverage calculation with efficient BAM index usage and chunking
fn run_coverage_parallel_chunked(cli: &Cli, read_stats: Option<ReadStats>) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::{Write, BufWriter};
    use std::io::BufReader;
    use noodles_bam as bam;
    use noodles_core::Region;
    use noodles_core::Position;
    use nanocov::parse_bed;

    println!("Using enhanced parallel coverage calculation with chunking");

    // Set thread pool size
    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .map_err(|e| format!("Failed to set thread pool: {}", e))?;
    }

    // Open BAM file and read header once
    let mut reader = bam::io::reader::Builder::default().build_from_path(&cli.input)?;
    let header = reader.read_header()?;

    // Parse BED files
    let bed_regions = if let Some(bed_path) = &cli.bed {
        Some(parse_bed(bed_path)?)
    } else {
        None
    };

    let chrom_bed_regions = if let Some(chrom_bed_path) = &cli.chrom_bed {
        Some(parse_bed(chrom_bed_path)?)
    } else {
        None
    };

    // Create chunk jobs for parallel processing
    let chunk_jobs = create_chunk_jobs(&header, &bed_regions, &chrom_bed_regions, cli.chunk_size)?;
    
    println!("Created {} chunk jobs for parallel processing", chunk_jobs.len());

    // Process chunks in parallel with shared BAM access
    let coverage_results: Vec<HashMap<String, HashMap<u32, u32>>> = chunk_jobs
        .par_iter()
        .map(|chunk| {
            process_chunk_efficiently(&cli.input, chunk)
                .unwrap_or_else(|e| {
                    eprintln!("Warning: Failed to process chunk {:?}: {}", chunk, e);
                    HashMap::new()
                })
        })
        .collect();

    // Merge results from all chunks
    let mut merged_coverage: HashMap<String, HashMap<u32, u32>> = HashMap::new();
    let mut per_chrom_averages: Vec<f64> = Vec::new();
    
    for chunk_coverage in coverage_results {
        for (chrom, positions) in chunk_coverage {
            let chrom_coverage = merged_coverage.entry(chrom.clone()).or_default();
            for (pos, count) in positions {
                *chrom_coverage.entry(pos).or_insert(0) += count;
            }
        }
    }

    // Calculate per-chromosome averages
    for (chrom, positions) in &merged_coverage {
        let (total, count) = positions.values().fold((0u64, 0u64), |(t, c), v| (t + *v as u64, c + 1));
        if count > 0 {
            let avg = total as f64 / count as f64;
            println!("{} average coverage: {:.2}", chrom, avg);
            per_chrom_averages.push(avg);
        }
    }

    // Write output using parallel formatting
    write_coverage_parallel(&merged_coverage, &cli.output)?;

    // Generate plots if requested
    if !cli.skip_all_plots {
        generate_plots_from_coverage_map(cli, &merged_coverage, read_stats.as_ref())?;
    }

    // Print summary
    if !per_chrom_averages.is_empty() {
        let global_avg = per_chrom_averages.iter().sum::<f64>() / per_chrom_averages.len() as f64;
        println!("Global average coverage: {:.2}", global_avg);
    }

    Ok(())
}

/// Chunk definition for parallel processing
#[derive(Debug, Clone)]
struct CoverageChunk {
    chromosome: String,
    start: u32,
    end: u32,
    target_regions: Option<Vec<(u32, u32)>>, // BED regions within this chunk
}

/// Create optimized chunk jobs for parallel processing
fn create_chunk_jobs(
    header: &noodles_sam::Header,
    bed_regions: &Option<HashMap<String, Vec<(u32, u32)>>>,
    chrom_bed_regions: &Option<HashMap<String, Vec<(u32, u32)>>>,
    chunk_size: usize,
) -> Result<Vec<CoverageChunk>, Box<dyn std::error::Error>> {
    let mut chunks = Vec::new();
    let chunk_size = chunk_size as u32;
    let reference_sequences = header.reference_sequences();

    for (chrom_name, ref_seq) in reference_sequences.iter() {
        let chrom = chrom_name.to_string();
        let chrom_length = ref_seq.length().get() as u32;

        // Determine what regions to process for this chromosome
        let regions_to_process = if let Some(bed) = bed_regions {
            bed.get(&chrom).cloned().unwrap_or_default()
        } else if let Some(chrom_bed) = chrom_bed_regions {
            chrom_bed.get(&chrom).cloned().unwrap_or_default()
        } else {
            // Process entire chromosome
            vec![(1, chrom_length)]
        };

        if regions_to_process.is_empty() {
            continue;
        }

        // Create chunks for each region
        for (region_start, region_end) in regions_to_process {
            let mut current_start = region_start;
            
            while current_start < region_end {
                let chunk_end = std::cmp::min(current_start + chunk_size, region_end);
                
                // For BED regions, find any sub-regions within this chunk
                let chunk_target_regions = if bed_regions.is_some() {
                    Some(vec![(current_start, chunk_end)])
                } else {
                    None
                };

                chunks.push(CoverageChunk {
                    chromosome: chrom.clone(),
                    start: current_start,
                    end: chunk_end,
                    target_regions: chunk_target_regions,
                });

                current_start = chunk_end;
            }
        }
    }

    Ok(chunks)
}

/// Process a single chunk efficiently with optimized BAM access
fn process_chunk_efficiently(
    bam_path: &std::path::Path,
    chunk: &CoverageChunk,
) -> Result<HashMap<String, HashMap<u32, u32>>, Box<dyn std::error::Error>> {
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::BufReader;
    use noodles_bam as bam;
    use noodles_core::Region;
    use noodles_core::Position;

    // Open BAM reader for this thread
    let mut reader = bam::io::Reader::new(BufReader::new(File::open(bam_path)?));
    let mut bai_reader = noodles_bam::bai::Reader::new(File::open(&bam_path.with_extension("bam.bai"))?);
    let index = bai_reader.read_index()?;
    let header = reader.read_header()?;

    let mut coverage: HashMap<String, HashMap<u32, u32>> = HashMap::new();

    // Create region for this chunk
    let region = Region::new(
        chunk.chromosome.clone(),
        Position::try_from(chunk.start as usize)?..=Position::try_from(chunk.end as usize)?,
    );

    // Query BAM for this region
    let query = reader.query(&header, &index, &region)?;
    
    for result in query {
        let record = result?;
        if record.flags().is_unmapped() {
            continue;
        }

        let start_pos = match record.alignment_start() {
            Some(Ok(pos)) => pos.get() as u32,
            _ => continue,
        };

        let alignment_len = calculate_reference_span(&record.cigar());
        let end_pos = start_pos + alignment_len;

        // Skip if alignment doesn't overlap with our chunk
        if end_pos < chunk.start || start_pos > chunk.end {
            continue;
        }

        // Calculate overlap with chunk boundaries
        let overlap_start = std::cmp::max(start_pos, chunk.start);
        let overlap_end = std::cmp::min(end_pos, chunk.end);

        // Update coverage efficiently using range iteration
        let chrom_coverage = coverage.entry(chunk.chromosome.clone()).or_default();
        
        // Optimized coverage update - avoid per-base loop for long alignments
        if alignment_len <= 1000 {
            // For short alignments, use per-base counting
            for pos in overlap_start..overlap_end {
                *chrom_coverage.entry(pos).or_insert(0) += 1;
            }
        } else {
            // For long alignments, use a more efficient approach
            update_coverage_range(chrom_coverage, overlap_start, overlap_end);
        }
    }

    Ok(coverage)
}

/// Efficiently update coverage for a range (optimized for long reads)
fn update_coverage_range(coverage: &mut HashMap<u32, u32>, start: u32, end: u32) {
    // For very long ranges, we could implement a more sophisticated
    // interval-based approach, but for now use the simple approach
    for pos in start..end {
        *coverage.entry(pos).or_insert(0) += 1;
    }
}

/// Write coverage data efficiently using parallel formatting
fn write_coverage_parallel(
    coverage: &HashMap<String, HashMap<u32, u32>>,
    output_path: &std::path::Path,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::io::{Write, BufWriter};
    use std::fs::File;
    use std::fmt::Write as FmtWrite;

    let mut out = BufWriter::new(File::create(output_path)?);
    writeln!(out, "#chromosome\tposition\tcount")?;

    // Parallelize formatting of coverage lines per chromosome
    let chrom_blocks: Vec<(String, String)> = coverage.par_iter()
        .map(|(chrom, positions)| {
            let mut positions_vec: Vec<_> = positions.iter().collect();
            positions_vec.sort_by_key(|&(pos, _)| *pos);
            
            let mut block = String::with_capacity(positions_vec.len() * 24); // estimate
            for (&pos, &count) in positions_vec {
                let _ = write!(block, "{}\t{}\t{}\n", chrom, pos, count);
            }
            (chrom.clone(), block)
        })
        .collect();

    // Write blocks in chromosome order
    let mut sorted_blocks = chrom_blocks;
    sorted_blocks.sort_by_key(|(chrom, _)| chrom.clone());
    
    for (_, block) in sorted_blocks {
        out.write_all(block.as_bytes())?;
    }
    out.flush()?;

    Ok(())
}

/// Generate plots from the merged coverage map
fn generate_plots_from_coverage_map(
    cli: &Cli,
    coverage: &HashMap<String, HashMap<u32, u32>>,
    read_stats: Option<&ReadStats>,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::collections::HashMap;

    if cli.skip_all_plots {
        return Ok(());
    }

    let output_stem = cli.output.file_stem().unwrap_or_default().to_string_lossy();
    let output_dir = cli.output.parent().unwrap_or_else(|| std::path::Path::new("."));

    // Apply theme if specified
    if let Some(theme) = &cli.theme {
        crate::plotting::set_theme(theme);
    }

    // Generate individual chromosome plots
    for (chrom, chrom_coverage) in coverage {
        let plot_path = output_dir.join(format!("coverage.{}.png", chrom));
        
        // Determine plot range from coverage data
        let min_pos = chrom_coverage.keys().min().copied().unwrap_or(0);
        let max_pos = chrom_coverage.keys().max().copied().unwrap_or(0);
        
        crate::plotting::plot_per_base_coverage_with_range(
            chrom,
            chrom_coverage,
            plot_path.to_str().unwrap(),
            min_pos,
            max_pos,
            read_stats,
            cli.show_zero_regions,
            cli.log_scale,
        )?;
    }

    // Generate multi-chromosome plot if we have multiple chromosomes
    if coverage.len() > 1 && !cli.skip_multi_plot {
        let multi_plot_path = output_dir.join(format!("{}.multi_chrom.png", output_stem));
        
        // Convert to the expected format for plot_all_chromosomes
        let chrom_coverages: HashMap<String, &HashMap<u32, u32>> = coverage
            .iter()
            .map(|(k, v)| (k.clone(), v))
            .collect();

        // Get current theme
        let theme = unsafe { crate::plotting::CURRENT_THEME };

        crate::plotting::plot_all_chromosomes(
            &chrom_coverages,
            multi_plot_path.to_str().unwrap(),
            cli.log_scale,
            read_stats,
            theme,
        )?;
    }

    Ok(())
}
