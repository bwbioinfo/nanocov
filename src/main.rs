mod cli;
mod io;
mod utils;
mod plotting;

use crate::cli::Cli;

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

    // Move coverage calculation and output logic to io module
    io::run_coverage(&cli)?;
    Ok(())
}
