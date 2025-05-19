mod cli;
mod io;
mod utils;
mod plotting;

use crate::cli::Cli;
use std::path::PathBuf;

use clap::Parser;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    // Prepare BAM and BAI file paths
    let bam_path = &cli.input;
    let bai_path = bam_path.with_extension("bam.bai");
    if !bai_path.exists() {
        eprintln!("BAM index not found at {:?}. Please run 'samtools index {:?}' to create it.", bai_path, bam_path);
        std::process::exit(1);
    }

    // Extract read stats
    let read_stats = crate::utils::extract_read_stats(bam_path)?;

    // Generate cramino output if requested
    if cli.cramino_output {
        let cramino_path = if let Some(p) = &cli.cramino_output_path {
            p.clone()
        } else {
            // Default: use input filename with .cramino extension
            let mut default_path = PathBuf::from(bam_path);
            default_path.set_extension("cramino");
            default_path
        };
        
        println!("Generating cramino-like output at {:?}", cramino_path);
        
        // We don't have coverage info yet, so use 0.0 for now
        // Use genome_size from CLI if provided
        io::cramino::generate_cramino_output(
            bam_path,
            &cramino_path,
            Some(&read_stats),
            0.0, // We don't have coverage info yet
            cli.genome_size.unwrap_or(0),
        )?;
    }

    // Move coverage calculation and output logic to io module, pass read_stats
    io::run_coverage(&cli, Some(read_stats))?;
    Ok(())
}
