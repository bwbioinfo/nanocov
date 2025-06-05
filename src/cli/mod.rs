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
    
    /// Generate cramino-like output (in addition to regular output)
    #[arg(long = "cramino")]
    pub cramino_output: bool,
    
    /// Path for cramino-like output file (default: uses input filename with .cramino extension)
    #[arg(long = "cramino-output")]
    pub cramino_output_path: Option<PathBuf>,
    
    /// Genome size in base pairs (used for coverage calculation in cramino output)
    #[arg(long = "genome-size")]
    pub genome_size: Option<u64>,
    
    /// Force streaming mode for memory efficiency (processes chromosomes one at a time)
    #[arg(long = "streaming")]
    pub force_streaming: bool,
    
    /// Memory limit in MB before switching to streaming mode (default: 500MB file size threshold)
    #[arg(long = "memory-limit")]
    pub memory_limit_mb: Option<u64>,
    
    /// Skip generating plots to save memory (output TSV only)
    #[arg(long = "no-plots")]
    pub skip_all_plots: bool,
}
