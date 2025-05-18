// src/plotting/multi_chrom.rs
//
// Module for creating multi-chromosome overview plots

use crate::utils::ReadStats;
use crate::plotting::themes::ColorTheme;
use plotters::prelude::*;
use std::collections::HashMap;

/// Plot all chromosomes on a single chart for comparison
///
/// Creates a bar chart showing the mean coverage for each chromosome,
/// optionally using logarithmic scale for the y-axis.
///
/// # Arguments
/// * `chrom_coverages` - Map of chromosome names to coverage data
/// * `output_path` - Path to save the output plot
/// * `use_log_scale` - Whether to use log scale for y-axis
/// * `read_stats` - Optional read statistics for display
/// * `theme` - Color theme to use for plotting
pub fn plot_all_chromosomes(
    chrom_coverages: &HashMap<String, &HashMap<u32, u32>>,
    output_path: &str, 
    use_log_scale: bool,
    read_stats: Option<&ReadStats>,
    theme: &ColorTheme,
) -> Result<(), Box<dyn std::error::Error>> {
    // Filter to canonical chromosomes
    let canonical_chroms = [
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
        "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT", "M",
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
        "chr20", "chr21", "chr22", "chrX", "chrY", "chrM", "chrMT",
    ];
    
    // Calculate mean coverage for each chromosome
    let mut chrom_means = Vec::new();
    
    for (chrom, coverage) in chrom_coverages {
        // Skip non-canonical chromosomes or empty data
        let chrom_base = chrom.trim_start_matches("chr");
        if !canonical_chroms.contains(&chrom.as_str()) && !canonical_chroms.contains(&chrom_base) {
            continue;
        }
        
        if coverage.is_empty() {
            continue;
        }
        
        // Calculate mean coverage
        let total: u64 = coverage.values().map(|&v| v as u64).sum();
        let count = coverage.len();
        let mean = total as f64 / count as f64;
        
        // Format chromosome name consistently for display
        let display_name = if chrom.starts_with("chr") {
            chrom.trim_start_matches("chr").to_string()
        } else {
            chrom.to_string()
        };
        
        chrom_means.push((display_name, mean));
    }
    
    if chrom_means.is_empty() {
        eprintln!("No coverage data found for canonical chromosomes");
        return Ok(());
    }
    
    // Sort chromosomes naturally (1,2,3...10,11... X,Y,MT)
    chrom_means.sort_by(|(a, _), (b, _)| {
        // Check if both are numeric chromosomes
        if let (Ok(a_num), Ok(b_num)) = (a.parse::<u32>(), b.parse::<u32>()) {
            return a_num.cmp(&b_num);
        }
        
        // Special ordering for non-numeric chromosomes
        match (a.as_str(), b.as_str()) {
            // Numbers come before letters
            (a_str, b_str) if a_str.chars().next().unwrap().is_numeric() && 
                             !b_str.chars().next().unwrap().is_numeric() => {
                return std::cmp::Ordering::Less;
            }
            (a_str, b_str) if !a_str.chars().next().unwrap().is_numeric() && 
                             b_str.chars().next().unwrap().is_numeric() => {
                return std::cmp::Ordering::Greater;
            }
            
            // X comes before Y
            ("X", "Y") => return std::cmp::Ordering::Less,
            ("Y", "X") => return std::cmp::Ordering::Greater,
            
            // MT/M comes last
            ("MT", _) | ("M", _) => return std::cmp::Ordering::Greater,
            (_, "MT") | (_, "M") => return std::cmp::Ordering::Less,
            
            // Natural string comparison for anything else
            _ => return a.cmp(b),
        }
    });
    
    // Calculate global mean and find max/min for the plot
    let global_mean = chrom_means.iter().map(|(_, cov)| *cov).sum::<f64>() / chrom_means.len() as f64;
    let max_coverage = chrom_means.iter().map(|(_, cov)| *cov).fold(0.0, f64::max);
    
    // Setup drawing area
    let root = BitMapBackend::new(output_path, (1200, 800)).into_drawing_area();
    root.fill(&theme.base)?;
    
    // Split for title and stats
    let (title_area, chart_area) = root.split_vertically(80);
    
    // Draw title
    title_area.draw_text(
        &format!("Chromosome Coverage Overview ({})", 
                 if use_log_scale { "Log Scale" } else { "Linear Scale" }),
        &("sans-serif", 30).into_font().color(&theme.text),
        (600, 40),
    )?;
    
    // Get chromosome names for x-axis labels
    let x_labels: Vec<String> = chrom_means.iter()
        .map(|(name, _)| name.clone())
        .collect();
    
    // Setup the chart - different setup for log vs linear scale
    if use_log_scale {
        // Find minimum non-zero value for log scale
        let min_coverage = chrom_means.iter()
            .map(|(_, cov)| *cov)
            .filter(|&v| v > 0.0)
            .fold(max_coverage, f64::min)
            .max(0.1); // Minimum of 0.1
        
        let mut chart = ChartBuilder::on(&chart_area)
            .margin(10)
            .set_all_label_area_size(40)
            .build_cartesian_2d(
                0..chrom_means.len(),
                (min_coverage..max_coverage * 1.1).log_scale(),
            )?;
            
        chart.configure_mesh()
            .disable_x_mesh()
            .y_desc("Coverage (log scale)")
            .x_desc("Chromosome")
            .x_labels(chrom_means.len())
            .x_label_formatter(&|idx| {
                if *idx < x_labels.len() {
                    x_labels[*idx].clone()
                } else {
                    String::new()
                }
            })
            .draw()?;
            
        // Draw bars
        for (i, (_, coverage)) in chrom_means.iter().enumerate() {
            let color = get_color_for_coverage(*coverage, max_coverage, theme);
            chart.draw_series(std::iter::once(
                Rectangle::new([(i, min_coverage), (i + 1, *coverage)], color.filled())
            ))?;
        }
        
        // Draw global mean line if it's in range
        if global_mean >= min_coverage {
            chart.draw_series(LineSeries::new(
                vec![(0, global_mean), (chrom_means.len() - 1, global_mean)],
                theme.accent.stroke_width(2),
            ))?;
            
            // Add annotation
            chart.draw_series(std::iter::once(
                Text::new(
                    format!("Global Mean: {:.2}", global_mean),
                    (chrom_means.len() / 2, global_mean * 1.2),
                    ("sans-serif", 18).into_font().color(&theme.accent),
                )
            ))?;
        }
    } else {
        // Linear scale chart
        let mut chart = ChartBuilder::on(&chart_area)
            .margin(10)
            .set_all_label_area_size(40)
            .build_cartesian_2d(
                0..chrom_means.len(),
                0.0..max_coverage * 1.1,
            )?;
            
        chart.configure_mesh()
            .disable_x_mesh()
            .y_desc("Mean Coverage")
            .x_desc("Chromosome")
            .x_labels(chrom_means.len())
            .x_label_formatter(&|idx| {
                if *idx < x_labels.len() {
                    x_labels[*idx].clone()
                } else {
                    String::new()
                }
            })
            .draw()?;
            
        // Draw bars
        for (i, (_, coverage)) in chrom_means.iter().enumerate() {
            let color = get_color_for_coverage(*coverage, max_coverage, theme);
            chart.draw_series(std::iter::once(
                Rectangle::new([(i, 0.0), (i + 1, *coverage)], color.filled())
            ))?;
        }
        
        // Draw global mean line
        chart.draw_series(LineSeries::new(
            vec![(0, global_mean), (chrom_means.len() - 1, global_mean)],
            theme.accent.stroke_width(2),
        ))?;
        
        // Add annotation
        chart.draw_series(std::iter::once(
            Text::new(
                format!("Global Mean: {:.2}", global_mean),
                (chrom_means.len() / 2, global_mean * 1.2),
                ("sans-serif", 18).into_font().color(&theme.accent),
            )
        ))?;
    }
    
    // Add stats at the bottom
    let (max_chrom, max_cov) = chrom_means.iter()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(name, cov)| (name.as_str(), *cov))
        .unwrap_or(("", 0.0));
        
    let (min_chrom, min_cov) = chrom_means.iter()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(name, cov)| (name.as_str(), *cov))
        .unwrap_or(("", 0.0));
    
    let stats_text = format!(
        "Chromosomes: {}   |   Global Mean: {:.2}   |   Max: {:.2} ({})   |   Min: {:.2} ({})",
        chrom_means.len(), global_mean, max_cov, max_chrom, min_cov, min_chrom
    );
    
    chart_area.draw_text(
        &stats_text,
        &("sans-serif", 18).into_font().color(&theme.text),
        (600, 760),
    )?;
    
    // Add read stats if available
    if let Some(stats) = read_stats {
        let read_stats_text = format!(
            "N50: {}   |   Mean Read Length: {:.0}   |   Mean Quality: {:.2}",
            stats.n50, stats.mean_len, stats.mean_qual
        );
        
        chart_area.draw_text(
            &read_stats_text,
            &("sans-serif", 18).into_font().color(&theme.text),
            (600, 730),
        )?;
    }
    
    Ok(())
}

// Helper function to determine color based on relative coverage
fn get_color_for_coverage(coverage: f64, max_coverage: f64, theme: &ColorTheme) -> RGBColor {
    let rel_coverage = coverage / max_coverage;
    
    if rel_coverage < 0.3 {
        // Low coverage - blend from low to primary
        let blend = rel_coverage / 0.3;
        RGBColor(
            ((theme.low.0 as f64) * (1.0 - blend) + (theme.primary.0 as f64) * blend) as u8,
            ((theme.low.1 as f64) * (1.0 - blend) + (theme.primary.1 as f64) * blend) as u8,
            ((theme.low.2 as f64) * (1.0 - blend) + (theme.primary.2 as f64) * blend) as u8,
        )
    } else if rel_coverage > 0.7 {
        // High coverage - blend from primary to high
        let blend = (rel_coverage - 0.7) / 0.3;
        RGBColor(
            ((theme.primary.0 as f64) * (1.0 - blend) + (theme.high.0 as f64) * blend) as u8,
            ((theme.primary.1 as f64) * (1.0 - blend) + (theme.high.1 as f64) * blend) as u8,
            ((theme.primary.2 as f64) * (1.0 - blend) + (theme.high.2 as f64) * blend) as u8,
        )
    } else {
        // Medium coverage - use primary color
        theme.primary
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plotting::themes::CATPPUCCIN_LATTE;
    
    #[test]
    fn test_multi_chromosome_plot() {
        // Create test coverage data
        let mut chr1_coverage = HashMap::new();
        let mut chr2_coverage = HashMap::new();
        let mut chr_x_coverage = HashMap::new();
        
        for i in 0..100u32 {
            chr1_coverage.insert(i, 10 + (i % 5));  // Higher coverage
            chr2_coverage.insert(i, 5 + (i % 3));   // Medium coverage
            chr_x_coverage.insert(i, 2 + (i % 2));   // Lower coverage
        }
        
        let mut chrom_coverages = HashMap::new();
        chrom_coverages.insert("chr1".to_string(), &chr1_coverage);
        chrom_coverages.insert("chr2".to_string(), &chr2_coverage);
        chrom_coverages.insert("chrX".to_string(), &chr_x_coverage);
        
        // Test linear scale
        let out_path = "test-out/multi_chrom_test.png";
        let _ = std::fs::remove_file(out_path);
        plot_all_chromosomes(&chrom_coverages, out_path, false, None, &CATPPUCCIN_LATTE)
            .expect("Multi-chromosome plotting should succeed");
        
        assert!(std::fs::metadata(out_path).is_ok(), "Output file should exist");
    }
}