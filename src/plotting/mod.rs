// src/plotting/mod.rs
// Plotting module for nanocov - handles coverage visualization and statistics display

mod multi_chrom;
mod stats;
mod themes;
mod utils;

use stats::{calculate_coverage_stats, calculate_per_base_stats};
use themes::{CATPPUCCIN_FRAPPE, CATPPUCCIN_LATTE, ColorTheme, GRUVBOX_LIGHT, NORD};
use utils::format_number;

// Re-export multi-chromosome plotting functionality
pub use multi_chrom::plot_all_chromosomes;

// Default theme (can be overridden via CLI)
pub static mut CURRENT_THEME: &'static ColorTheme = &CATPPUCCIN_LATTE;

// Shorthand for accessing current theme colors
#[inline]
fn theme() -> &'static ColorTheme {
    unsafe { CURRENT_THEME }
}

use plotters::prelude::*;
use std::collections::HashMap;

/// Plots per-base coverage as a bar graph for a single chromosome.
///
/// This function takes a chromosome name, a map of base positions to coverage counts,
/// and an output path for the PNG file. It sorts the positions, determines the y-axis
/// scaling based on the data, and draws a bar for each covered base. The plot is saved
/// as a PNG file using the plotters crate.
///
/// # Arguments
/// * `chrom` - Chromosome name (for labeling the plot)
/// * `coverage` - Map from base position (u32, starting at 0) to coverage count (u32)
/// * `output_path` - Path to save the output PNG file
/// * `read_stats` - Optional read statistics to display in the sidebar
/// * `show_zero_regions` - Whether to show regions with zero coverage
///
/// # Details
/// - The y-axis is automatically scaled based on the data range
/// - Each bar is drawn from the current position to the next, with color gradient based on coverage
/// - The plot includes detailed statistics panels and legends
/// - The function is robust to empty or sparse data
use crate::utils::ReadStats;

/// Set the global color theme for all plots
///
/// # Arguments
/// * `theme_name` - One of: "latte", "frappe", "nord", "gruvbox"
#[inline]
pub fn set_theme(theme_name: &str) {
    unsafe {
        CURRENT_THEME = match theme_name.to_lowercase().as_str() {
            "frappe" => &CATPPUCCIN_FRAPPE,
            "nord" => &NORD,
            "gruvbox" => &GRUVBOX_LIGHT,
            _ => &CATPPUCCIN_LATTE, // Default to latte
        };
    }
}

/// Simple entry point for plotting coverage with automatic range determination
/// This is externally used by the io module
#[inline]
#[allow(dead_code)] // Used by io module
pub fn plot_per_base_coverage(
    chrom: &str,
    coverage: &HashMap<u32, u32>,
    output_path: &str,
    read_stats: Option<&ReadStats>,
    show_zero_regions: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    // Use the min/max of coverage as the range
    let mut positions: Vec<u32> = coverage.keys().copied().collect();
    positions.sort_unstable();

    if positions.is_empty() {
        eprintln!("Warning: No coverage data for chromosome {}", chrom);
        return Ok(());
    }

    let min_x = positions.iter().min().copied().unwrap_or(0);
    let max_x = positions.iter().max().copied().unwrap_or(0);
    plot_per_base_coverage_with_range(
        chrom,
        coverage,
        output_path,
        min_x,
        max_x,
        read_stats,
        show_zero_regions,
    )
}

pub fn plot_per_base_coverage_with_range(
    chrom: &str,
    coverage: &HashMap<u32, u32>,
    output_path: &str,
    plot_start: u32,
    plot_end: u32,
    read_stats: Option<&ReadStats>,
    show_zero_regions: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    // Adaptive binning logic for large ranges
    let range = plot_end.saturating_sub(plot_start);
    let bin_size = if range > 100_000_000 {
        1_000_000 // 1Mb bins for very large regions (>100Mb)
    } else if range > 10_000_000 {
        100_000 // 100kb bins for large regions (>10Mb)
    } else if range > 1_000_000 {
        10_000 // 10kb bins for medium regions (>1Mb)
    } else if range > 100_000 {
        1_000 // 1kb bins for smaller regions (>100kb)
    } else if range > 10_000 {
        100 // 100bp bins for tiny regions (>10kb)
    } else {
        1 // No binning for very small regions
    };

    let mut binned: std::collections::BTreeMap<u32, (u64, u32)> = std::collections::BTreeMap::new();
    for (&pos, &count) in coverage.iter() {
        if pos < plot_start || pos > plot_end {
            continue;
        }
        let bin = ((pos - plot_start) / bin_size) * bin_size + plot_start;
        let entry = binned.entry(bin).or_insert((0, 0));
        entry.0 += count as u64;
        entry.1 += 1;
    }
    let chrom_points: Vec<(i64, f64)> = binned
        .iter()
        .map(|(&bin, &(sum, n))| (bin as i64, if n > 0 { sum as f64 / n as f64 } else { 0.0 }))
        .collect();

    // Calculate and add zero coverage regions if requested
    let chrom_points = if show_zero_regions && !chrom_points.is_empty() {
        // Create a more complete dataset with zero regions
        let mut full_points = Vec::new();
        let mut last_pos = plot_start;

        for &(pos, coverage) in &chrom_points {
            // Add zeros for any gaps
            if pos as u32 > last_pos + bin_size {
                // Add zero points for missing regions
                let mut gap_pos = last_pos + bin_size;
                while gap_pos < pos as u32 {
                    full_points.push((gap_pos as i64, 0.0));
                    gap_pos += bin_size;
                }
            }
            full_points.push((pos, coverage));
            last_pos = pos as u32;
        }
        full_points
    } else {
        chrom_points
    };

    // Determine the maximum y (coverage) value with smarter scaling
    let y_max = chrom_points.iter().map(|&(_, y)| y).fold(0.0, f64::max);
    let y_max = if y_max < 3.0 {
        3.0 // Minimum scale for very low coverage
    } else if y_max < 10.0 {
        (y_max * 1.2).ceil() // Slightly larger for clearer view of small values
    } else if y_max < 100.0 {
        (y_max * 1.1).ceil() // 10% headroom for medium values
    } else {
        (y_max * 1.05).ceil() // 5% headroom for large values
    };

    // Debug output
    eprintln!("[DEBUG] Max coverage for {}: {}", chrom, y_max);

    // Set up the drawing area with BitMapBackend (PNG output)
    let root = BitMapBackend::new(output_path, (2200, 1000)).into_drawing_area();

    root.fill(&theme().base)?;

    // Split horizontally: left panel (stats), middle panel (main plot), right panel (scale bars)
    let (left_panel, remaining) = root.split_horizontally(400);
    let (plot_area, right_panel) = remaining.split_horizontally(1400);

    // Format bin size for display (kb or Mb for large bins)
    let bin_size_label = if bin_size >= 1_000_000 {
        format!("{} Mb", bin_size / 1_000_000)
    } else if bin_size >= 1_000 {
        format!("{} kb", bin_size / 1_000)
    } else {
        format!("{} bp", bin_size)
    };

    // Build the chart with large margins and axis labels for clarity
    let mut chart = ChartBuilder::on(&plot_area)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .set_label_area_size(LabelAreaPosition::Left, 75)
        .set_label_area_size(LabelAreaPosition::Bottom, 50)
        .margin(50)
        .caption(
            format!("Chromosome {} Coverage (bin: {})", chrom, bin_size_label),
            ("sans-serif", 40).into_font().color(&theme().text),
        )
        .build_cartesian_2d(plot_start as i64..plot_end as i64, 0f64..y_max)?;

    // Configure the mesh (axes, grid, and labels)
    chart
        .configure_mesh()
        .x_desc("Chromosome Position (Mb)")
        .y_desc("Coverage")
        .axis_desc_style(("sans-serif", 25).into_font().color(&theme().text))
        .x_label_formatter(&|x| format!("{:.2}", (*x as f64) / 1_000_000.0)) // Format x axis in Mb
        .x_labels(20) // Increase number of x-axis ticks
        .x_label_style(("sans-serif", 18).into_font().color(&theme().text))
        .y_label_style(("sans-serif", 18).into_font().color(&theme().text))
        .y_label_formatter(&|y| format!("{:.1}", y)) // Clean y-axis formatting
        .light_line_style(RGBAColor(100, 100, 100, 0.3)) // Subtle grid lines
        .bold_line_style(RGBAColor(100, 100, 100, 0.5)) // Bolder major grid lines
        .axis_style(&theme().text)
        .draw()?;

    // Draw each bar as a filled rectangle with gradient color based on coverage
    chart.draw_series(chrom_points.windows(2).map(|w| {
        let (x, y) = w[0];
        let (next_x, _) = w[1];
        let bar_end = if next_x > x {
            next_x
        } else {
            x + bin_size as i64
        };

        // Color gradient based on coverage level
        let fill_color = if y <= y_max * 0.3 {
            // Low coverage: blend from low coverage color to main color
            let blend_factor = (y / (y_max * 0.3)) as f64;
            let blend_factor = blend_factor.max(0.0).min(1.0);
            RGBColor(
                ((theme().low.0 as f64) * (1.0 - blend_factor)
                    + (theme().primary.0 as f64) * blend_factor) as u8,
                ((theme().low.1 as f64) * (1.0 - blend_factor)
                    + (theme().primary.1 as f64) * blend_factor) as u8,
                ((theme().low.2 as f64) * (1.0 - blend_factor)
                    + (theme().primary.2 as f64) * blend_factor) as u8,
            )
        } else if y >= y_max * 0.7 {
            // High coverage: blend from main color to high coverage color
            let blend_factor = ((y - y_max * 0.7) / (y_max * 0.3)) as f64;
            let blend_factor = blend_factor.max(0.0).min(1.0);
            RGBColor(
                ((theme().primary.0 as f64) * (1.0 - blend_factor)
                    + (theme().high.0 as f64) * blend_factor) as u8,
                ((theme().primary.1 as f64) * (1.0 - blend_factor)
                    + (theme().high.1 as f64) * blend_factor) as u8,
                ((theme().primary.2 as f64) * (1.0 - blend_factor)
                    + (theme().high.2 as f64) * blend_factor) as u8,
            )
        } else {
            // Medium coverage: use main color
            theme().primary
        };

        Rectangle::new([(x, 0.0), (bar_end, y)], fill_color.filled())
    }))?;

    // Draw the last bar if only one point or for the last position
    if let Some(&(x, y)) = chrom_points.last() {
        // Apply same color logic for last bar
        let fill_color = if y <= y_max * 0.3 {
            let blend_factor = (y / (y_max * 0.3)) as f64;
            let blend_factor = blend_factor.max(0.0).min(1.0);
            RGBColor(
                ((theme().low.0 as f64) * (1.0 - blend_factor)
                    + (theme().primary.0 as f64) * blend_factor) as u8,
                ((theme().low.1 as f64) * (1.0 - blend_factor)
                    + (theme().primary.1 as f64) * blend_factor) as u8,
                ((theme().low.2 as f64) * (1.0 - blend_factor)
                    + (theme().primary.2 as f64) * blend_factor) as u8,
            )
        } else if y >= y_max * 0.7 {
            let blend_factor = ((y - y_max * 0.7) / (y_max * 0.3)) as f64;
            let blend_factor = blend_factor.max(0.0).min(1.0);
            RGBColor(
                ((theme().primary.0 as f64) * (1.0 - blend_factor)
                    + (theme().high.0 as f64) * blend_factor) as u8,
                ((theme().primary.1 as f64) * (1.0 - blend_factor)
                    + (theme().high.1 as f64) * blend_factor) as u8,
                ((theme().primary.2 as f64) * (1.0 - blend_factor)
                    + (theme().high.2 as f64) * blend_factor) as u8,
            )
        } else {
            theme().primary
        };

        chart.draw_series(std::iter::once(Rectangle::new(
            [(x, 0.0), (x + bin_size as i64, y)],
            fill_color.filled(),
        )))?;
    }

    // Add threshold line if average coverage is available
    let y_vals: Vec<f64> = chrom_points.iter().map(|&(_, y)| y).collect();
    if !y_vals.is_empty() {
        let mean = y_vals.iter().sum::<f64>() / y_vals.len() as f64;
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(plot_start as i64, mean), (plot_end as i64, mean)],
            theme().accent.stroke_width(2),
        )))?;
    }

    // --- Apply annotations and legends in the right panel ---

    // Fill the right panel with the base color
    right_panel.fill(&theme().base)?;

    // Define scale bar dimensions
    let scale_box_width = 80;
    let scale_box_height = 30;
    let scale_spacing = 15;
    let padding = 30;
    let text_offset = 100;

    // Add color gradient box
    let gradient_height = 200;
    let gradient_start_y = padding + 80 + 3 * scale_box_height + 3 * scale_spacing + 20;

    right_panel.draw_text(
        "Coverage Scale",
        &("sans-serif", 24)
            .into_font()
            .style(FontStyle::Bold)
            .color(&theme().text),
        (padding, gradient_start_y),
    )?;

    // Draw color gradient
    for i in 0..gradient_height {
        let factor = (gradient_height - i) as f64 / gradient_height as f64;
        let color = if factor < 0.33 {
            // Low to medium transition
            let blend = factor / 0.33;
            RGBColor(
                ((theme().low.0 as f64) * (1.0 - blend) + (theme().primary.0 as f64) * blend) as u8,
                ((theme().low.1 as f64) * (1.0 - blend) + (theme().primary.1 as f64) * blend) as u8,
                ((theme().low.2 as f64) * (1.0 - blend) + (theme().primary.2 as f64) * blend) as u8,
            )
        } else if factor > 0.66 {
            // Medium to high transition
            let blend = (factor - 0.66) / 0.34;
            RGBColor(
                ((theme().primary.0 as f64) * (1.0 - blend) + (theme().high.0 as f64) * blend)
                    as u8,
                ((theme().primary.1 as f64) * (1.0 - blend) + (theme().high.1 as f64) * blend)
                    as u8,
                ((theme().primary.2 as f64) * (1.0 - blend) + (theme().high.2 as f64) * blend)
                    as u8,
            )
        } else {
            // Medium range
            theme().primary
        };

        right_panel.draw(&PathElement::new(
            vec![
                (padding, gradient_start_y + 40 + i),
                (padding + scale_box_width, gradient_start_y + 40 + i),
            ],
            color.stroke_width(1),
        ))?;
    }

    // Add labels for gradient
    right_panel.draw_text(
        "High",
        &("sans-serif", 16).into_font().color(&theme().text),
        (padding + text_offset, gradient_start_y + 40),
    )?;

    right_panel.draw_text(
        "Medium",
        &("sans-serif", 16).into_font().color(&theme().text),
        (
            padding + text_offset,
            gradient_start_y + 40 + gradient_height / 2,
        ),
    )?;

    right_panel.draw_text(
        "Low",
        &("sans-serif", 16).into_font().color(&theme().text),
        (
            padding + text_offset,
            gradient_start_y + 40 + gradient_height - 20,
        ),
    )?;

    // --- Draw read stats and coverage stats panels on the left ---
    let font = ("sans-serif", 24)
        .into_font()
        .style(FontStyle::Bold)
        .color(&theme().text);
    let padding = 30;
    let line_height = 32;
    let box_width = (400.0 * 1.1) as i32;
    let box_height = (6 * line_height + padding * 2) as i32;
    let num_boxes = 3;
    let spacing = 40;
    let total_boxes_height = num_boxes * box_height + (num_boxes - 1) * spacing;
    let available_height = 1000;
    let box_x = 30;
    let box_y_top = (available_height - total_boxes_height) / 2;
    let box_y_bottom = box_y_top + box_height + spacing;
    let per_base_box_y = box_y_bottom + box_height + spacing;
    // --- Top box: Read stats ---
    if let Some(stats) = read_stats {
        let stats_labels = [
            "N50:",
            "Mean Qual:",
            "Median Qual:",
            "Mean Length:",
            "Median Len:",
        ];
        let stats_values = [
            format_number(stats.n50 as u64),
            format!("{:.2}", stats.mean_qual),
            format!("{:.2}", stats.median_qual),
            format_number(stats.mean_len as u64),
            format_number(stats.median_len as u64),
        ];

        // Draw box background and border
        left_panel.draw(&Rectangle::new(
            [
                (box_x, box_y_top),
                (box_x + box_width, box_y_top + box_height),
            ],
            ShapeStyle {
                color: theme().overlay.to_rgba(),
                filled: true,
                stroke_width: 0,
            },
        ))?;
        left_panel.draw(&Rectangle::new(
            [
                (box_x, box_y_top),
                (box_x + box_width, box_y_top + box_height),
            ],
            ShapeStyle {
                color: theme().accent.to_rgba(),
                filled: false,
                stroke_width: 3,
            },
        ))?;

        // Draw title and stats
        left_panel.draw_text("Read Stats", &font, (box_x + padding, box_y_top + padding))?;
        for (i, (label, value)) in stats_labels.iter().zip(stats_values.iter()).enumerate() {
            left_panel.draw_text(
                label,
                &font,
                (
                    box_x + padding,
                    box_y_top + padding + ((i as i32 + 1) * line_height),
                ),
            )?;
            left_panel.draw_text(
                value,
                &font,
                (
                    box_x + box_width - padding - 160,
                    box_y_top + padding + ((i as i32 + 1) * line_height),
                ),
            )?;
        }
    }

    // Calculate statistics for display in the stats boxes
    let coverage_stats = calculate_coverage_stats(&chrom_points);
    let mean = coverage_stats.mean;
    let median = coverage_stats.median;
    let min = coverage_stats.min;
    let max = coverage_stats.max;

    // --- Bottom box: Coverage stats ---
    let stats_labels = ["Mean:", "Median:", "Min:", "Max:", "Bin:"];
    let stats_values = [
        format!("{:.2}", mean),
        format!("{:.2}", median),
        format!("{:.2}", min),
        format!("{:.2}", max),
        bin_size_label.clone(),
    ];
    left_panel.draw(&Rectangle::new(
        [
            (box_x, box_y_bottom),
            (box_x + box_width, box_y_bottom + box_height),
        ],
        ShapeStyle {
            color: theme().overlay.to_rgba(),
            filled: true,
            stroke_width: 0,
        },
    ))?;
    left_panel.draw(&Rectangle::new(
        [
            (box_x, box_y_bottom),
            (box_x + box_width, box_y_bottom + box_height),
        ],
        ShapeStyle {
            color: theme().accent.to_rgba(),
            filled: false,
            stroke_width: 3,
        },
    ))?;
    left_panel.draw_text(
        "Coverage Stats",
        &font,
        (box_x + padding, box_y_bottom + padding),
    )?;
    for (i, (label, value)) in stats_labels.iter().zip(stats_values.iter()).enumerate() {
        left_panel.draw_text(
            label,
            &font,
            (
                box_x + padding,
                box_y_bottom + padding + ((i as i32 + 1) * line_height),
            ),
        )?;
        left_panel.draw_text(
            value,
            &font,
            (
                box_x + box_width - padding - 160,
                box_y_bottom + padding + ((i as i32 + 1) * line_height),
            ),
        )?;
    }

    // --- Per-base coverage stats box ---
    let per_base_labels = [
        "Per-base Mean:",
        "Per-base Median:",
        "Per-base Min:",
        "Per-base Max:",
        "Per-base Stddev:",
    ];
    // Calculate per-base statistics using the helper function
    let per_base_stats = calculate_per_base_stats(coverage);
    let per_base_values = [
        format!("{:.2}", per_base_stats.mean),
        format!("{:.2}", per_base_stats.median),
        format!("{:.2}", per_base_stats.min),
        format!("{:.2}", per_base_stats.max),
        format!("{:.2}", per_base_stats.stddev),
    ];

    left_panel.draw(&Rectangle::new(
        [
            (box_x, per_base_box_y),
            (box_x + box_width, per_base_box_y + box_height),
        ],
        ShapeStyle {
            color: theme().overlay.to_rgba(),
            filled: true,
            stroke_width: 0,
        },
    ))?;
    left_panel.draw(&Rectangle::new(
        [
            (box_x, per_base_box_y),
            (box_x + box_width, per_base_box_y + box_height),
        ],
        ShapeStyle {
            color: theme().accent.to_rgba(),
            filled: false,
            stroke_width: 3,
        },
    ))?;
    left_panel.draw_text(
        "Per-base Coverage",
        &font,
        (box_x + padding, per_base_box_y + padding),
    )?;
    for (i, (label, value)) in per_base_labels
        .iter()
        .zip(per_base_values.iter())
        .enumerate()
    {
        left_panel.draw_text(
            label,
            &font,
            (
                box_x + padding,
                per_base_box_y + padding + ((i as i32 + 1) * line_height),
            ),
        )?;
        left_panel.draw_text(
            value,
            &font,
            (
                box_x + box_width - padding - 160,
                per_base_box_y + padding + ((i as i32 + 1) * line_height),
            ),
        )?;
    }

    // This repeated computation was unnecessary and can be removed

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_plot_output() {
        let mut coverage = HashMap::new();
        for i in 0..50u32 {
            coverage.insert(i, (i % 5) + 1);
        }

        let out_path = "test-out/coverage.test.png";
        let _ = fs::remove_file(out_path);
        plot_per_base_coverage("chrTest", &coverage, out_path, None, true)
            .expect("PNG plotting should succeed");

        assert!(fs::metadata(out_path).is_ok(), "Output PNG should exist");
        let meta = fs::metadata(out_path).unwrap();
        assert!(meta.len() > 0, "Output PNG should not be empty");
    }

    #[test]
    fn test_plot_per_base_coverage_with_different_themes() {
        let mut coverage = HashMap::new();
        // Simulate a small region with variable coverage
        for i in 0..100u32 {
            coverage.insert(i, (i % 10) + 1);
        }

        // Test with Nord theme
        set_theme("nord");
        let out_path = "test-out/coverage.test.nord.png";
        let _ = fs::remove_file(out_path);
        plot_per_base_coverage("chrTest", &coverage, out_path, None, false)
            .expect("plotting with Nord theme should succeed");

        // Test with Frappe theme
        set_theme("frappe");
        let out_path = "test-out/coverage.test.frappe.png";
        let _ = fs::remove_file(out_path);
        plot_per_base_coverage("chrTest", &coverage, out_path, None, false)
            .expect("plotting with Frappe theme should succeed");

        // Test with Gruvbox theme
        set_theme("gruvbox");
        let out_path = "test-out/coverage.test.gruvbox.png";
        let _ = fs::remove_file(out_path);
        plot_per_base_coverage("chrTest", &coverage, out_path, None, false)
            .expect("plotting with Frappe theme should succeed");

        // Test with Latte theme
        set_theme("latte");
        let out_path = "test-out/coverage.test.latte.png";
        let _ = fs::remove_file(out_path);
        plot_per_base_coverage("chrTest", &coverage, out_path, None, false)
            .expect("plotting with Frappe theme should succeed");

        // Reset to default theme
        set_theme("latte");
    }
}
