// src/plotting/multi_chrom.rs
//
// Module for creating multi-chromosome overview plots

use crate::utils::ReadStats;
use crate::plotting::themes::ColorTheme;
use plotters::prelude::*;
use std::collections::HashMap;

/// Plot all chromosomes on a single chart for comparison
///
/// Creates a continuous bar chart showing 100kb-binned coverage across all chromosomes,
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
    
    // Bin size for multi-chromosome plot
    const BIN_SIZE: u32 = 100_000; // 100kb bins
    
    // Collect and bin coverage data for each chromosome
    let mut chrom_data = Vec::new();
    
    for (chrom, coverage) in chrom_coverages {
        // Skip non-canonical chromosomes or empty data
        let chrom_base = chrom.trim_start_matches("chr");
        if !canonical_chroms.contains(&chrom.as_str()) && !canonical_chroms.contains(&chrom_base) {
            continue;
        }
        
        if coverage.is_empty() {
            continue;
        }
        
        // Format chromosome name consistently for display
        let display_name = if chrom.starts_with("chr") {
            chrom.trim_start_matches("chr").to_string()
        } else {
            chrom.to_string()
        };
        
        // Bin the coverage data
        let mut binned: std::collections::BTreeMap<u32, (u64, u32)> = std::collections::BTreeMap::new();
        for (&pos, &count) in coverage.iter() {
            let bin = (pos / BIN_SIZE) * BIN_SIZE;
            let entry = binned.entry(bin).or_insert((0, 0));
            entry.0 += count as u64;
            entry.1 += 1;
        }
        
        // Convert to average coverage per bin
        let binned_coverage: Vec<(u32, f64)> = binned
            .iter()
            .map(|(&bin, &(sum, n))| (bin, if n > 0 { sum as f64 / n as f64 } else { 0.0 }))
            .collect();
        
        if !binned_coverage.is_empty() {
            // Calculate chromosome length (approximate)
            let chrom_length = binned_coverage.iter().map(|(pos, _)| *pos).max().unwrap_or(0) + BIN_SIZE;
            chrom_data.push((display_name, binned_coverage, chrom_length));
        }
    }
    
    if chrom_data.is_empty() {
        eprintln!("No coverage data found for canonical chromosomes");
        return Ok(());
    }
    
    // Sort chromosomes naturally (1,2,3...10,11... X,Y,MT)
    chrom_data.sort_by(|(a, _, _), (b, _, _)| {
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
    
    // Create a continuous plot data with chromosome boundaries
    let mut plot_data: Vec<(i64, f64)> = Vec::new();
    let mut chrom_boundaries: Vec<(String, i64, i64)> = Vec::new(); // (name, start_x, end_x)
    let mut current_x = 0i64;
    
    for (chrom_name, binned_coverage, chrom_length) in &chrom_data {
        let start_x = current_x;
        
        // Add binned coverage data
        for (bin_pos, coverage) in binned_coverage {
            let x_pos = current_x + (*bin_pos as i64 / BIN_SIZE as i64);
            plot_data.push((x_pos, *coverage));
        }
        
        // Calculate chromosome span in terms of bins
        let chrom_span = (*chrom_length as i64 / BIN_SIZE as i64) + 1;
        let end_x = current_x + chrom_span;
        
        chrom_boundaries.push((chrom_name.clone(), start_x, end_x));
        current_x = end_x + 1; // Add gap between chromosomes
    }
    
    if plot_data.is_empty() {
        eprintln!("No binned coverage data available for plotting");
        return Ok(());
    }
    
    // Calculate max coverage for scaling
    let max_coverage = plot_data.iter().map(|(_, cov)| *cov).fold(0.0, f64::max);
    let global_mean = plot_data.iter().map(|(_, cov)| *cov).sum::<f64>() / plot_data.len() as f64;
    
    // Setup drawing area
    let root = BitMapBackend::new(output_path, (1600, 800)).into_drawing_area();
    root.fill(&theme.base)?;
    
    // Split for title
    let (title_area, mut chart_area) = root.split_vertically(80);
    
    // Draw title
    title_area.draw_text(
        &format!("Multi-Chromosome Coverage Overview - 100kb bins ({})", 
                if use_log_scale { "Log Scale" } else { "Linear Scale" }),
        &("sans-serif", 30).into_font().color(&theme.text),
        (800, 40),
    )?;
    
    let total_x_span = current_x - 1;
    
    // Create chart with appropriate scale and draw chromosome labels/boundaries first
    let chart_plotting_area = chart_area.margin(10, 10, 80, 60);
    
    // Draw chromosome labels and boundaries first (before creating the chart)
    for (chrom_name, start_x, end_x) in &chrom_boundaries {
        let center_x = (start_x + end_x) / 2;
        let x_pixel = 10.0 + (center_x as f64 / total_x_span as f64) * (1600.0 - 80.0 - 20.0); // Approximate chart width
        
        // Draw chromosome label
        chart_area.draw_text(
            chrom_name,
            &("sans-serif", 12).into_font().color(&theme.text),
            (x_pixel as i32, 720),
        )?;
        
        // Draw vertical boundary lines (except for the first chromosome)
        if *start_x > 0 {
            let boundary_x_pixel = 10.0 + (*start_x as f64 / total_x_span as f64) * (1600.0 - 80.0 - 20.0);
            
            chart_area.draw(&PathElement::new(
                vec![
                    (boundary_x_pixel as i32, 80),
                    (boundary_x_pixel as i32, 640)
                ],
                RGBColor(100, 100, 100).stroke_width(1),
            ))?;
        }
    }
    
    if use_log_scale {
        // Find minimum non-zero value for log scale
        let min_coverage = plot_data.iter()
            .map(|(_, cov)| *cov)
            .filter(|&v| v > 0.0)
            .fold(max_coverage, f64::min)
            .max(0.1); // Minimum of 0.1
        
        let mut chart = ChartBuilder::on(&chart_area)
            .margin(10)
            .set_label_area_size(LabelAreaPosition::Left, 60)
            .set_label_area_size(LabelAreaPosition::Bottom, 80)
            .build_cartesian_2d(
                0..total_x_span,
                (min_coverage..max_coverage * 1.1).log_scale(),
            )?;
            
        chart.configure_mesh()
            .disable_x_mesh()
            .y_desc("Coverage (log scale)")
            .x_desc("Chromosome")
            .x_labels(0)
            .draw()?;
            
        // Draw bars for each bin
        for window in plot_data.windows(2) {
            let (x1, y1) = window[0];
            let (x2, _) = window[1];
            
            if y1 > 0.0 {
                let color = get_color_for_coverage(y1, max_coverage, theme);
                let bar_width = (x2 - x1).max(1);
                chart.draw_series(std::iter::once(
                    Rectangle::new([(x1, min_coverage), (x1 + bar_width, y1)], color.filled())
                ))?;
            }
        }
        
        // Draw the last bar
        if let Some(&(x, y)) = plot_data.last() {
            if y > 0.0 {
                let color = get_color_for_coverage(y, max_coverage, theme);
                chart.draw_series(std::iter::once(
                    Rectangle::new([(x, min_coverage), (x + 1, y)], color.filled())
                ))?;
            }
        }
        
        // Draw global mean line if it's in range
        if global_mean >= min_coverage {
            chart.draw_series(LineSeries::new(
                vec![(0, global_mean), (total_x_span, global_mean)],
                theme.accent.stroke_width(2),
            ))?;
        }
    } else {
        // Linear scale chart
        let mut chart = ChartBuilder::on(&chart_area)
            .margin(10)
            .set_label_area_size(LabelAreaPosition::Left, 60)
            .set_label_area_size(LabelAreaPosition::Bottom, 80)
            .build_cartesian_2d(
                0..total_x_span,
                0.0..max_coverage * 1.1,
            )?;
            
        chart.configure_mesh()
            .disable_x_mesh()
            .y_desc("Coverage")
            .x_desc("Chromosome")
            .x_labels(0)
            .draw()?;
            
        // Draw bars for each bin
        for window in plot_data.windows(2) {
            let (x1, y1) = window[0];
            let (x2, _) = window[1];
            
            let color = get_color_for_coverage(y1, max_coverage, theme);
            let bar_width = (x2 - x1).max(1);
            chart.draw_series(std::iter::once(
                Rectangle::new([(x1, 0.0), (x1 + bar_width, y1)], color.filled())
            ))?;
        }
        
        // Draw the last bar
        if let Some(&(x, y)) = plot_data.last() {
            let color = get_color_for_coverage(y, max_coverage, theme);
            chart.draw_series(std::iter::once(
                Rectangle::new([(x, 0.0), (x + 1, y)], color.filled())
            ))?;
        }
        
        // Draw global mean line
        chart.draw_series(LineSeries::new(
            vec![(0, global_mean), (total_x_span, global_mean)],
            theme.accent.stroke_width(2),
        ))?;
    }
    
    // Add stats at the bottom
    let (max_chrom, max_cov) = chrom_data.iter()
        .flat_map(|(name, bins, _)| bins.iter().map(|(_, cov)| (name.as_str(), *cov)))
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or(("", 0.0));
        
    let (min_chrom, min_cov) = chrom_data.iter()
        .flat_map(|(name, bins, _)| bins.iter().map(|(_, cov)| (name.as_str(), *cov)))
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or(("", 0.0));
    
    let stats_text = format!(
        "Chromosomes: {}   |   Global Mean: {:.2}   |   Max: {:.2} ({})   |   Min: {:.2} ({})",
        chrom_data.len(), global_mean, max_cov, max_chrom, min_cov, min_chrom
    );
    
    chart_area.draw_text(
        &stats_text,
        &("sans-serif", 16).into_font().color(&theme.text),
        (800, 780),
    )?;
    
    // Add read stats if available
    if let Some(stats) = read_stats {
        let read_stats_text = format!(
            "N50: {}   |   Mean Read Length: {:.0}   |   Mean Quality: {:.2}",
            stats.n50, stats.mean_len, stats.mean_qual
        );
        
        chart_area.draw_text(
            &read_stats_text,
            &("sans-serif", 16).into_font().color(&theme.text),
            (800, 760),
        )?;
    }
    
    // Present the chart to finalize it (must be at the very end)
    root.present()?;

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
        
        for i in 0..100000u32 {
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