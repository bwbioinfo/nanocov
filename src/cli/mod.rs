// src/cli/mod.rs
// CLI argument parsing module for nanocov
// Contains the Cli struct and related logic

use std::path::PathBuf;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "bam-coverage")]
#[command(about = "Calculates per-base coverage from a BAM file", long_about = None)]
pub struct Cli {
    /// Input BAM file
    #[arg(short, long)]
    pub input: PathBuf,

    /// BED file with regions to include (chrom, start, end)
    #[arg(short = 'b', long = "bed")]
    pub bed: Option<PathBuf>,

    /// BED file with full chromosome ranges (chrom, start, end for each chromosome)
    #[arg(long = "chrom-bed")]
    pub chrom_bed: Option<PathBuf>,

    /// Number of threads to use (default: half of available cores)
    #[arg(short = 't', long = "threads")]
    pub threads: Option<usize>,

    /// Output file path (default: coverage.tsv)
    #[arg(short = 'o', long = "output", default_value = "coverage.tsv")]
    pub output: PathBuf,

    /// Chunk size for parallel processing (default: 10,000)
    #[arg(short = 'c', long = "chunk-size", default_value_t = 10_000)]
    pub chunk_size: usize,
    
    /// Use SVG output format for plots instead of PNG
    #[arg(long = "svg")]
    pub svg_output: bool,
    
    /// Color theme for plots [latte, frappe, nord, gruvbox]
    #[arg(long = "theme")]
    pub theme: Option<String>,
    
    /// Show regions with zero coverage in plots
    #[arg(long = "show-zeros")]
    pub show_zero_regions: bool,
    
    /// Skip plotting (generate only TSV output)
    #[arg(long = "no-plot")]
    pub skip_plotting: bool,
    
    /// Skip multi-chromosome summary plot
    #[arg(long = "no-multi-plot")]
    pub skip_multi_plot: bool,
    
    /// Use logarithmic scale for the multi-chromosome plot
    #[arg(long = "log-scale")]
    pub log_scale: bool,
}
