// src/io/cramino.rs
// Module for generating cramino-like output from BAM files

use crate::utils::ReadStats;
use noodles_bam as bam;
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::SystemTime;
use std::fmt::Write as FmtWrite;

/// Format for cramino-like output
pub struct CraminoOutput {
    pub file_name: String,
    pub num_alignments: u64,
    pub percent_from_total: f64,
    pub num_reads: u64,
    pub yield_gb: f64,
    pub mean_coverage: f64,
    pub yield_gb_greater_than_25kb: f64,
    pub n50: u32,
    pub n75: u32,
    pub median_length: f64,
    pub mean_length: f64,
    pub path: PathBuf,
    pub creation_time: String,
}

impl CraminoOutput {
    /// Create a new CraminoOutput struct with all zeros/empty values
    pub fn new_empty(path: PathBuf) -> Self {
        let file_name = path
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_default();
            
        // Format creation time as dd/mm/yyyy hh:mm:ss
        let creation_time = match path.metadata() {
            Ok(metadata) => {
                if let Ok(modified_time) = metadata.modified() {
                    if let Ok(duration) = modified_time.duration_since(SystemTime::UNIX_EPOCH) {
                        let secs = duration.as_secs();
                        let naive_time = chrono::DateTime::<chrono::Utc>::from_timestamp(secs as i64, 0)
                            .unwrap_or_else(|| chrono::DateTime::<chrono::Utc>::from_timestamp(0, 0).unwrap())
                            .naive_utc();
                        naive_time.format("%d/%m/%Y %H:%M:%S").to_string()
                    } else {
                        String::from("01/01/2000 00:00:00") // Fallback if time calculation fails
                    }
                } else {
                    String::from("01/01/2000 00:00:00") // Fallback if modified time can't be read
                }
            }
            Err(_) => String::from("01/01/2000 00:00:00"), // Fallback if metadata can't be read
        };

        Self {
            file_name,
            num_alignments: 0,
            percent_from_total: 0.0,
            num_reads: 0,
            yield_gb: 0.0,
            mean_coverage: 0.0,
            yield_gb_greater_than_25kb: 0.0,
            n50: 0,
            n75: 0,
            median_length: 0.0,
            mean_length: 0.0,
            path,
            creation_time,
        }
    }

    /// Create a new CraminoOutput struct from read statistics and BAM path
    pub fn from_read_stats(path: PathBuf, read_stats: &ReadStats, total_coverage: f64, genome_size: u64) -> Self {
        let mut result = Self::new_empty(path.clone());
        
        // Calculate N75
        // (We assume the ReadStats already has calculated N50 and sorted lengths)
        let n75 = if let Some(lengths) = &read_stats.lengths {
            if lengths.is_empty() {
                0
            } else {
                let total: u64 = lengths.iter().map(|&l| l as u64).sum();
                let mut acc = 0;
                let mut n75 = 0;
                for &l in lengths {
                    acc += l as u64;
                    if acc >= (total * 3) / 4 { // 75% of total length
                        n75 = l;
                        break;
                    }
                }
                n75
            }
        } else {
            0
        };
        
        // Calculate yield_gb and yield_gb_greater_than_25kb
        let total_bases = read_stats.num_bases;
        let yield_gb = total_bases as f64 / 1_000_000_000.0;
        
        let total_bases_gt_25kb = if let Some(lengths) = &read_stats.lengths {
            lengths.iter().filter(|&&l| l > 25_000).map(|&l| l as u64).sum::<u64>()
        } else {
            0
        };
        let yield_gb_gt_25kb = total_bases_gt_25kb as f64 / 1_000_000_000.0;
        
        // Calculate mean coverage
        let mean_coverage = if genome_size > 0 {
            total_bases as f64 / genome_size as f64
        } else {
            total_coverage // Use the provided total coverage if available
        };
        
        result.num_alignments = read_stats.num_reads;
        result.percent_from_total = 100.0; // Assume all alignments are counted
        result.num_reads = read_stats.num_reads;
        result.yield_gb = yield_gb;
        result.mean_coverage = mean_coverage;
        result.yield_gb_greater_than_25kb = yield_gb_gt_25kb;
        result.n50 = read_stats.n50;
        result.n75 = n75;
        result.median_length = read_stats.median_len;
        result.mean_length = read_stats.mean_len;
        
        result
    }
    
    /// Format as a string in cramino-like format
    pub fn format(&self) -> String {
        let mut output = String::new();
        
        writeln!(&mut output, "File name\t{}", self.file_name).unwrap();
        writeln!(&mut output, "Number of alignments\t{}", self.num_alignments).unwrap();
        writeln!(&mut output, "% from total alignments\t{:.2}", self.percent_from_total).unwrap();
        writeln!(&mut output, "Number of reads\t{}", self.num_reads).unwrap();
        writeln!(&mut output, "Yield [Gb]\t{:.2}", self.yield_gb).unwrap();
        writeln!(&mut output, "Mean coverage\t{:.2}", self.mean_coverage).unwrap();
        writeln!(&mut output, "Yield [Gb] (>25kb)\t{:.2}", self.yield_gb_greater_than_25kb).unwrap();
        writeln!(&mut output, "N50\t{}", self.n50).unwrap();
        writeln!(&mut output, "N75\t{}", self.n75).unwrap();
        writeln!(&mut output, "Median length\t{:.2}", self.median_length).unwrap();
        writeln!(&mut output, "Mean length\t{:.2}", self.mean_length).unwrap();
        writeln!(&mut output, "").unwrap(); // Empty line
        writeln!(&mut output, "Path\t{}", self.path.display()).unwrap();
        writeln!(&mut output, "Creation time\t{}", self.creation_time).unwrap();
        
        output
    }
    
    /// Write to a file
    pub fn write_to_file(&self, output_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(output_path)?;
        file.write_all(self.format().as_bytes())?;
        Ok(())
    }
}

/// Enhanced version of ReadStats to include more information needed for cramino output
pub struct EnhancedReadStats {
    pub n50: u32,
    pub n75: u32, 
    pub mean_len: f64,
    pub median_len: f64,
    pub _mean_qual: f64,
    pub _median_qual: f64,
    pub num_reads: u64,
    pub num_bases: u64,
    pub lengths: Option<Vec<u32>>, // Store lengths for additional calculations
}

impl From<&ReadStats> for EnhancedReadStats {
    fn from(stats: &ReadStats) -> Self {
        // Convert basic ReadStats to EnhancedReadStats
        // Note: Some fields like num_reads, num_bases, and lengths 
        // need to be populated separately
        Self {
            n50: stats.n50,
            n75: 0, // Will be calculated later
            mean_len: stats.mean_len,
            median_len: stats.median_len,
            _mean_qual: stats.mean_qual,
            _median_qual: stats.median_qual,
            num_reads: 0, // Will be calculated later
            num_bases: 0, // Will be calculated later
            lengths: None, // Will be populated later
        }
    }
}

/// Extract enhanced read statistics from a BAM file
pub fn extract_enhanced_read_stats(bam_path: &Path) -> Result<EnhancedReadStats, Box<dyn std::error::Error>> {
    let mut reader = bam::io::Reader::new(BufReader::new(File::open(bam_path)?));
    let _header = reader.read_header()?;
    let mut lengths = Vec::new();
    let mut quals = Vec::new();
    let mut num_reads: u64 = 0;
    let mut num_bases: u64 = 0;

    for result in reader.records() {
        let record = result?;
        let len = record.sequence().len() as u32;
        lengths.push(len);
        num_reads += 1;
        num_bases += len as u64;

        // Collect mean quality per read (if available)
        let qual = record.quality_scores();
        let qual_slice = qual.as_ref();
        if !qual_slice.is_empty() {
            let q: f64 = qual_slice.iter().map(|&q| q as u32).sum::<u32>() as f64 / qual_slice.len() as f64;
            quals.push(q);
        }
    }

    // N50 calculation
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let total: u64 = lengths.iter().map(|&l| l as u64).sum();
    let mut acc = 0;
    let mut n50 = 0;
    for &l in &lengths {
        acc += l as u64;
        if acc >= total / 2 {
            n50 = l;
            break;
        }
    }

    // N75 calculation
    let mut acc = 0;
    let mut n75 = 0;
    for &l in &lengths {
        acc += l as u64;
        if acc >= (total * 3) / 4 {
            n75 = l;
            break;
        }
    }

    // Mean/median length
    let mean_len = if lengths.is_empty() { 0.0 } else { lengths.iter().sum::<u32>() as f64 / lengths.len() as f64 };
    let median_len = if lengths.is_empty() {
        0.0
    } else {
        let mid = lengths.len() / 2;
        if lengths.len() % 2 == 0 {
            (lengths[mid - 1] + lengths[mid]) as f64 / 2.0
        } else {
            lengths[mid] as f64
        }
    };

    // Mean/median quality
    let mean_qual = if quals.is_empty() { 0.0 } else { quals.iter().sum::<f64>() / quals.len() as f64 };
    let median_qual = if quals.is_empty() {
        0.0
    } else {
        let mut sorted = quals.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mid = sorted.len() / 2;
        if sorted.len() % 2 == 0 {
            (sorted[mid - 1] + sorted[mid]) / 2.0
        } else {
            sorted[mid]
        }
    };

    Ok(EnhancedReadStats {
        n50,
        n75,
        mean_len,
        median_len,
        _mean_qual: mean_qual,
        _median_qual: median_qual,
        num_reads,
        num_bases,
        lengths: Some(lengths), // Store lengths for additional calculations
    })
}

/// Generate cramino-like output for a BAM file
pub fn generate_cramino_output(
    bam_path: &Path, 
    output_path: &Path, 
    read_stats: Option<&ReadStats>, 
    total_coverage: f64,
    genome_size: u64,
) -> Result<(), Box<dyn std::error::Error>> {
    let cramino_output = if let Some(stats) = read_stats {
        // Use provided stats if available
        let _enhanced_stats = EnhancedReadStats::from(stats);
        CraminoOutput::from_read_stats(bam_path.to_path_buf(), stats, total_coverage, genome_size)
    } else {
        // Otherwise extract stats from BAM
        let enhanced_stats = extract_enhanced_read_stats(bam_path)?;
        let path_buf = bam_path.to_path_buf();
        if enhanced_stats.num_reads == 0 {
            // Empty BAM file
            CraminoOutput::new_empty(path_buf)
        } else {
            CraminoOutput {
                file_name: path_buf.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default(),
                num_alignments: enhanced_stats.num_reads,
                percent_from_total: 100.0,
                num_reads: enhanced_stats.num_reads,
                yield_gb: enhanced_stats.num_bases as f64 / 1_000_000_000.0,
                mean_coverage: if genome_size > 0 { enhanced_stats.num_bases as f64 / genome_size as f64 } else { total_coverage },
                yield_gb_greater_than_25kb: if let Some(ref lengths) = enhanced_stats.lengths {
                    let total_gt_25kb: u64 = lengths.iter().filter(|&&l| l > 25_000).map(|&l| l as u64).sum();
                    total_gt_25kb as f64 / 1_000_000_000.0
                } else {
                    0.0
                },
                n50: enhanced_stats.n50,
                n75: enhanced_stats.n75,
                median_length: enhanced_stats.median_len,
                mean_length: enhanced_stats.mean_len,
                path: path_buf,
                creation_time: match bam_path.metadata() {
                    Ok(metadata) => {
                        if let Ok(modified_time) = metadata.modified() {
                            if let Ok(duration) = modified_time.duration_since(SystemTime::UNIX_EPOCH) {
                                let secs = duration.as_secs();
                                let naive_time = chrono::DateTime::<chrono::Utc>::from_timestamp(secs as i64, 0)
                                    .unwrap_or_else(|| chrono::DateTime::<chrono::Utc>::from_timestamp(0, 0).unwrap())
                                    .naive_utc();
                                naive_time.format("%d/%m/%Y %H:%M:%S").to_string()
                            } else {
                                String::from("01/01/2000 00:00:00")
                            }
                        } else {
                            String::from("01/01/2000 00:00:00")
                        }
                    }
                    Err(_) => String::from("01/01/2000 00:00:00"),
                },
            }
        }
    };

    cramino_output.write_to_file(output_path)?;
    Ok(())
}