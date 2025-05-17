// src/plotting/mod.rs
// Plotting module for nanocov (currently empty, placeholder for future plotting logic)

// Catppuccin-inspired colors (Latte palette, light theme)
const CAT_BLUE: RGBColor = RGBColor(30, 102, 245); // Blue
const CAT_MAUVE: RGBColor = RGBColor(136, 57, 239); // Mauve
const CAT_GREEN: RGBColor = RGBColor(64, 160, 43); // Green
const CAT_BASE: RGBColor = RGBColor(239, 241, 245); // Base (background)
const CAT_OVERLAY: RGBColor = RGBColor(220, 224, 232); // Overlay (box)
const CAT_TEXT: RGBColor = RGBColor(76, 79, 105); // Text

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
///
/// # Details
/// - The y-axis is automatically scaled: for high coverage (>1000), it uses the next power of 10;
///   for low coverage (<5), it sets a minimum of 5; otherwise, it uses the next integer above the max.
/// - Each bar is drawn from the current position to the next, or width 1 for the last point.
/// - The plot is 3000x500 pixels, with large labels and margins for readability.
/// - The function is robust to empty or sparse data.
use crate::utils::ReadStats;

pub fn plot_per_base_coverage(
    chrom: &str,
    coverage: &HashMap<u32, u32>,
    output_path: &str,
    read_stats: Option<&ReadStats>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Use the min/max of coverage as the range
    let mut positions: Vec<u32> = coverage.keys().copied().collect();
    positions.sort_unstable();
    let min_x = positions.iter().min().copied().unwrap_or(0);
    let max_x = positions.iter().max().copied().unwrap_or(0);
    plot_per_base_coverage_with_range(chrom, coverage, output_path, min_x, max_x, read_stats)
}

pub fn plot_per_base_coverage_with_range(
    chrom: &str,
    coverage: &HashMap<u32, u32>,
    output_path: &str,
    plot_start: u32,
    plot_end: u32,
    read_stats: Option<&ReadStats>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Binning logic for large ranges
    let range = plot_end.saturating_sub(plot_start);
    let bin_size = if range > 10_000_000 {
        100_000
    } else if range > 1_000_000 {
        10_000
    } else if range > 100_000 {
        1_000
    } else if range > 10_000 {
        100
    } else {
        1
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

    // Determine the maximum y (coverage) value
    let y_max = chrom_points.iter().map(|&(_, y)| y).fold(0.0, f64::max);
    eprintln!("[DEBUG] Y-axis max for {}: {}", chrom, y_max);
    let y_max = if y_max < 3.0 { 3.0 } else { y_max };
    eprintln!("[DEBUG] Max coverage for {}: {}", chrom, y_max);

    // Set up the drawing area: 2000x1000 pixels for higher DPI, Catppuccin background
    let root = BitMapBackend::new(output_path, (2000, 1000)).into_drawing_area();
    root.fill(&CAT_BASE)?;

    // Split horizontally: left panel (stats), right panel (main plot)
    let (left_panel, plot_area) = root.split_horizontally(400);

    // Build the chart with large margins and axis labels for clarity
    let mut chart = ChartBuilder::on(&plot_area)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .set_label_area_size(LabelAreaPosition::Left, 75)
        .set_label_area_size(LabelAreaPosition::Bottom, 50)
        .margin(50)
        .caption(
            format!("Chromosome {} Coverage (bin: {} bp)", chrom, bin_size),
            ("sans-serif", 40).into_font().color(&CAT_TEXT),
        )
        .build_cartesian_2d(plot_start as i64..plot_end as i64, 0f64..y_max)?;

    // Configure the mesh (axes, grid, and labels)
    chart
        .configure_mesh()
        .x_desc("Chromosome Position (Mb)")
        .y_desc("Coverage")
        .axis_desc_style(("sans-serif", 25).into_font().color(&CAT_TEXT))
        .x_label_formatter(&|x| format!("{:.2}", (*x as f64) / 1_000_000.0)) // Format x axis in Mb
        .x_labels(20) // Increase number of x-axis ticks
        .x_label_style(("sans-serif", 18).into_font().color(&CAT_TEXT))
        .y_label_style(("sans-serif", 18).into_font().color(&CAT_TEXT))
        .axis_style(&CAT_TEXT)
        .draw()?;

    // Draw each bar as a filled rectangle from the current position to the next (or width 1 for the last)
    chart.draw_series(chrom_points.windows(2).map(|w| {
        let (x, y) = w[0];
        let (next_x, _) = w[1];
        let bar_end = if next_x > x {
            next_x
        } else {
            x + bin_size as i64
        };
        Rectangle::new([(x, 0.0), (bar_end, y)], CAT_BLUE.filled())
    }))?;
    // Draw the last bar if only one point or for the last position
    if let Some(&(x, y)) = chrom_points.last() {
        chart.draw_series(std::iter::once(Rectangle::new(
            [(x, 0.0), (x + bin_size as i64, y)],
            CAT_BLUE.filled(),
        )))?;
    }

    // --- Stats annotation in main plot area ---

    // --- Draw read stats and coverage stats panels on the left ---
    let font = ("sans-serif", 24)
        .into_font()
        .style(FontStyle::Bold)
        .color(&CAT_TEXT);
    let padding = 30;
    let line_height = 32;
    let box_width = 400;
    let box_height = ((6 * line_height + padding * 2) as f32 * 1.1) as i32; // 1 title + 5 stats, 10% more space
    let box_x = 30;
    let box_y_top = 40;
    let box_y_bottom = box_y_top + box_height + 40;

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
            format!("{}", stats.n50),
            format!("{:.2}", stats.mean_qual),
            format!("{:.2}", stats.median_qual),
            format!("{:.0}", stats.mean_len),
            format!("{:.0}", stats.median_len),
        ];
        left_panel.draw(&Rectangle::new(
            [
                (box_x, box_y_top),
                (box_x + box_width, box_y_top + box_height),
            ],
            ShapeStyle {
                color: CAT_OVERLAY.to_rgba(),
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
                color: CAT_MAUVE.to_rgba(),
                filled: false,
                stroke_width: 3,
            },
        ))?;
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
                    box_x + box_width - padding - 60,
                    box_y_top + padding + ((i as i32 + 1) * line_height),
                ),
            )?;
        }
    }

    // Compute stats for coverage box
    let mut y_vals: Vec<f64> = chrom_points.iter().map(|&(_, y)| y).collect();
    y_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mean = if !y_vals.is_empty() {
        y_vals.iter().sum::<f64>() / y_vals.len() as f64
    } else {
        0.0
    };
    let median = if !y_vals.is_empty() {
        let mid = y_vals.len() / 2;
        if y_vals.len() % 2 == 0 {
            (y_vals[mid - 1] + y_vals[mid]) / 2.0
        } else {
            y_vals[mid]
        }
    } else {
        0.0
    };
    let min = y_vals.first().copied().unwrap_or(0.0);
    let max = y_vals.last().copied().unwrap_or(0.0);

    // --- Bottom box: Coverage stats ---
    let stats_labels = ["Mean:", "Median:", "Min:", "Max:", "Bin:"];
    let stats_values = [
        format!("{:.2}", mean),
        format!("{:.2}", median),
        format!("{:.2}", min),
        format!("{:.2}", max),
        format!("{} bp", bin_size),
    ];
    left_panel.draw(&Rectangle::new(
        [
            (box_x, box_y_bottom),
            (box_x + box_width, box_y_bottom + box_height),
        ],
        ShapeStyle {
            color: CAT_OVERLAY.to_rgba(),
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
            color: CAT_MAUVE.to_rgba(),
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
                box_x + box_width - padding - 60,
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
    let mut per_base_vals: Vec<u32> = coverage.values().copied().collect();
    per_base_vals.sort_unstable();
    let per_base_mean = if per_base_vals.is_empty() {
        0.0
    } else {
        per_base_vals.iter().sum::<u32>() as f64 / per_base_vals.len() as f64
    };
    let per_base_median = if per_base_vals.is_empty() {
        0.0
    } else {
        let mid = per_base_vals.len() / 2;
        if per_base_vals.len() % 2 == 0 {
            (per_base_vals[mid - 1] + per_base_vals[mid]) as f64 / 2.0
        } else {
            per_base_vals[mid] as f64
        }
    };
    let per_base_min = per_base_vals.first().copied().unwrap_or(0) as f64;
    let per_base_max = per_base_vals.last().copied().unwrap_or(0) as f64;
    let per_base_stddev = if per_base_vals.len() > 1 {
        let mean = per_base_mean;
        let var = per_base_vals
            .iter()
            .map(|&v| (v as f64 - mean).powi(2))
            .sum::<f64>()
            / (per_base_vals.len() as f64 - 1.0);
        var.sqrt()
    } else {
        0.0
    };
    let per_base_values = [
        format!("{:.2}", per_base_mean),
        format!("{:.2}", per_base_median),
        format!("{:.2}", per_base_min),
        format!("{:.2}", per_base_max),
        format!("{:.2}", per_base_stddev),
    ];

    let per_base_box_y = box_y_bottom + box_height + 40;
    left_panel.draw(&Rectangle::new(
        [
            (box_x, per_base_box_y),
            (box_x + box_width, per_base_box_y + box_height),
        ],
        ShapeStyle {
            color: CAT_OVERLAY.to_rgba(),
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
            color: CAT_MAUVE.to_rgba(),
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
                box_x + box_width - padding - 60,
                per_base_box_y + padding + ((i as i32 + 1) * line_height),
            ),
        )?;
    }

    // Compute stats
    let mut y_vals: Vec<f64> = chrom_points.iter().map(|&(_, y)| y).collect();
    y_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mean = if !y_vals.is_empty() {
        y_vals.iter().sum::<f64>() / y_vals.len() as f64
    } else {
        0.0
    };
    let median = if !y_vals.is_empty() {
        let mid = y_vals.len() / 2;
        if y_vals.len() % 2 == 0 {
            (y_vals[mid - 1] + y_vals[mid]) / 2.0
        } else {
            y_vals[mid]
        }
    } else {
        0.0
    };
    let min = y_vals.first().copied().unwrap_or(0.0);
    let max = y_vals.last().copied().unwrap_or(0.0);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_plot_per_base_coverage_creates_png() {
        let mut coverage = HashMap::new();
        // Simulate a small region with variable coverage
        for i in 0..100u32 {
            coverage.insert(i, (i % 10) + 1);
        }
        let out_path = "test-out/coverage.test.png";
        let _ = fs::remove_file(out_path); // Clean up before test
        // Use None for read_stats in test
        plot_per_base_coverage("chrTest", &coverage, out_path, None)
            .expect("plotting should succeed");
        assert!(fs::metadata(out_path).is_ok(), "Output PNG should exist");
        // Optionally, check file size is nonzero
        let meta = fs::metadata(out_path).unwrap();
        assert!(meta.len() > 0, "Output PNG should not be empty");
    }
}
