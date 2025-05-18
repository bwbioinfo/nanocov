// src/plotting/utils.rs
//
// Utility functions for plotting module

use plotters::style::RGBColor;

/// Format a number with thousands separators
///
/// # Arguments
/// * `num` - Number to format
///
/// # Returns
/// * String with formatted number (e.g., "1,234,567")
pub fn format_number(num: u64) -> String {
    let num_str = num.to_string();
    let mut result = String::new();
    let len = num_str.len();

    for (i, c) in num_str.chars().enumerate() {
        if i > 0 && (len - i) % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }

    result
}

/// Format a file size with appropriate units (B, KB, MB, GB)
///
/// # Arguments
/// * `size_bytes` - Size in bytes
///
/// # Returns
/// * String with formatted size (e.g., "1.23 MB")
#[allow(dead_code)] // Reserved for future use
pub fn format_file_size(size_bytes: u64) -> String {
    if size_bytes < 1024 {
        format!("{} B", size_bytes)
    } else if size_bytes < 1024 * 1024 {
        format!("{:.2} KB", size_bytes as f64 / 1024.0)
    } else if size_bytes < 1024 * 1024 * 1024 {
        format!("{:.2} MB", size_bytes as f64 / (1024.0 * 1024.0))
    } else {
        format!("{:.2} GB", size_bytes as f64 / (1024.0 * 1024.0 * 1024.0))
    }
}

/// Blend two colors based on a factor (0.0 to 1.0)
///
/// # Arguments
/// * `color1` - First color
/// * `color2` - Second color
/// * `factor` - Blend factor (0.0 = all color1, 1.0 = all color2)
///
/// # Returns
/// * Blended RGBColor
#[allow(dead_code)] // Used in debug builds
pub fn blend_colors(color1: &RGBColor, color2: &RGBColor, factor: f64) -> RGBColor {
    let factor = factor.max(0.0).min(1.0);
    
    RGBColor(
        ((color1.0 as f64) * (1.0 - factor) + (color2.0 as f64) * factor) as u8,
        ((color1.1 as f64) * (1.0 - factor) + (color2.1 as f64) * factor) as u8,
        ((color1.2 as f64) * (1.0 - factor) + (color2.2 as f64) * factor) as u8,
    )
}

/// Calculate appropriate bin size for a genomic range
///
/// # Arguments
/// * `range_size` - Size of the range in base pairs
///
/// # Returns
/// * Appropriate bin size in base pairs
#[allow(dead_code)] // Used in debug builds
pub fn calculate_bin_size(range_size: u32) -> u32 {
    if range_size > 100_000_000 {
        1_000_000     // 1Mb bins for very large regions (>100Mb)
    } else if range_size > 10_000_000 {
        100_000       // 100kb bins for large regions (>10Mb)
    } else if range_size > 1_000_000 {
        10_000        // 10kb bins for medium regions (>1Mb)
    } else if range_size > 100_000 {
        1_000         // 1kb bins for smaller regions (>100kb)
    } else if range_size > 10_000 {
        100           // 100bp bins for tiny regions (>10kb)
    } else {
        1             // No binning for very small regions
    }
}

/// Format bin size for display (kb or Mb for large bins)
///
/// # Arguments
/// * `bin_size` - Bin size in base pairs
///
/// # Returns
/// * Formatted string (e.g., "10 kb", "1 Mb")
#[allow(dead_code)] // Used in debug builds
pub fn format_bin_size(bin_size: u32) -> String {
    if bin_size >= 1_000_000 {
        format!("{} Mb", bin_size / 1_000_000)
    } else if bin_size >= 1_000 {
        format!("{} kb", bin_size / 1_000)
    } else {
        format!("{} bp", bin_size)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_format_number() {
        assert_eq!(format_number(1234), "1,234");
        assert_eq!(format_number(1000000), "1,000,000");
        assert_eq!(format_number(42), "42");
    }
    
    #[test]
    fn test_format_file_size() {
        assert_eq!(format_file_size(500), "500 B");
        assert_eq!(format_file_size(1500), "1.46 KB");
        assert_eq!(format_file_size(1500000), "1.43 MB");
    }
    
    #[test]
    fn test_blend_colors() {
        let color1 = RGBColor(0, 0, 0);
        let color2 = RGBColor(255, 255, 255);
        
        let blend_half = blend_colors(&color1, &color2, 0.5);
        assert_eq!(blend_half.0, 127);
        assert_eq!(blend_half.1, 127);
        assert_eq!(blend_half.2, 127);
    }
    
    #[test]
    fn test_calculate_bin_size() {
        assert_eq!(calculate_bin_size(5000), 1);
        assert_eq!(calculate_bin_size(50000), 100);
        assert_eq!(calculate_bin_size(500000), 1000);
        assert_eq!(calculate_bin_size(5000000), 10000);
        assert_eq!(calculate_bin_size(50000000), 100000);
        assert_eq!(calculate_bin_size(500000000), 1000000);
    }
}