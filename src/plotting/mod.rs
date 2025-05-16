// src/plotting/mod.rs
// Plotting module for nanocov (currently empty, placeholder for future plotting logic)

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
pub fn plot_per_base_coverage(
    chrom: &str,
    coverage: &HashMap<u32, u32>,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // Collect and sort all covered positions for plotting
    let mut positions: Vec<u32> = coverage.keys().copied().collect();
    positions.sort_unstable();
    // Determine the maximum x (position) and y (coverage) values
    let max_x = positions.iter().max().copied().unwrap_or(0) as i64;
    let y_max = coverage.values().copied().max().unwrap_or(1) as f64;
    // Adjust y_max for better visualization:
    //   - For very high coverage, use a log10 step above the max
    //   - For low coverage, set a minimum of 5
    //   - Otherwise, use the max value itself (tight axis)
    let y_max = if y_max > 1000.0 {
        10f64.powf((y_max.log10()).ceil())
    } else if y_max < 5.0 {
        5.0
    } else {
        y_max
    };
    // Print the max coverage value for debug
    eprintln!("[DEBUG] Max coverage for {}: {}", chrom, y_max);
    // Convert positions and coverage to a vector of (x, y) points for plotting
    let chrom_points: Vec<(i64, f64)> = positions
        .iter()
        .map(|&x| (x as i64, *coverage.get(&x).unwrap_or(&0) as f64))
        .collect();

    // Set up the drawing area: 3000x500 pixels, white background
    let root = BitMapBackend::new(output_path, (3000, 500)).into_drawing_area();
    root.fill(&WHITE)?;

    // Build the chart with large margins and axis labels for clarity
    let mut chart = ChartBuilder::on(&root)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .set_label_area_size(LabelAreaPosition::Left, 75)
        .set_label_area_size(LabelAreaPosition::Bottom, 50)
        .margin(50)
        .caption(
            format!("Chromosome {} Per-base Coverage", chrom),
            ("sans-serif", 40),
        )
        .build_cartesian_2d(0i64..max_x, 0f64..y_max)?;

    // Configure the mesh (axes, grid, and labels)
    chart
        .configure_mesh()
        .x_desc("Chromosome Position")
        .y_desc("Coverage")
        .axis_desc_style(("sans-serif", 25))
        .draw()?;

    // Draw each bar as a filled rectangle from the current position to the next (or width 1 for the last)
    // This ensures that even sparse data is visible as bars
    chart.draw_series(
        chrom_points.windows(2).map(|w| {
            let (x, y) = w[0];
            let (next_x, _) = w[1];
            // Draw a bar from x to next_x (or width 1 for last point)
            let bar_end = if next_x > x { next_x } else { x + 1 };
            Rectangle::new(
                [(x, 0.0), (bar_end, y)],
                BLUE.filled(),
            )
        })
    )?;
    // Handle the last point: always draw a bar of width 1 for the last covered base
    if let Some(&(x, y)) = chrom_points.last() {
        chart.draw_series(std::iter::once(
            Rectangle::new(
                [(x, 0.0), (x + 1, y)],
                BLUE.filled(),
            )
        ))?;
    }

    // Save the plot to disk
    root.present()?;
    println!("Plot saved to {}", output_path);
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
        plot_per_base_coverage("chrTest", &coverage, out_path).expect("plotting should succeed");
        assert!(fs::metadata(out_path).is_ok(), "Output PNG should exist");
        // Optionally, check file size is nonzero
        let meta = fs::metadata(out_path).unwrap();
        assert!(meta.len() > 0, "Output PNG should not be empty");
    }
}
