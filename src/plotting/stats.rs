// src/plotting/stats.rs
//
// Coverage statistics calculation for plotting

use std::collections::HashMap;

/// Structure for storing coverage statistics
pub struct CoverageStats {
    pub mean: f64,
    pub median: f64,
    pub min: f64,
    pub max: f64,
    pub stddev: f64,
}

impl Default for CoverageStats {
    fn default() -> Self {
        Self {
            mean: 0.0,
            median: 0.0,
            min: 0.0,
            max: 0.0,
            stddev: 0.0,
        }
    }
}

/// Calculate statistics from binned coverage points
///
/// # Arguments
/// * `points` - Vector of (position, coverage) points
///
/// # Returns
/// * `CoverageStats` - Computed statistics
pub fn calculate_coverage_stats(points: &[(i64, f64)]) -> CoverageStats {
    if points.is_empty() {
        return CoverageStats::default();
    }

    // Extract coverage values
    let mut y_vals: Vec<f64> = points.iter().map(|&(_, y)| y).collect();
    y_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Calculate mean
    let mean = y_vals.iter().sum::<f64>() / y_vals.len() as f64;

    // Calculate median
    let median = if y_vals.len() % 2 == 0 {
        let mid = y_vals.len() / 2;
        (y_vals[mid - 1] + y_vals[mid]) / 2.0
    } else {
        y_vals[y_vals.len() / 2]
    };

    // Get min and max
    let min = *y_vals.first().unwrap_or(&0.0);
    let max = *y_vals.last().unwrap_or(&0.0);

    // Calculate standard deviation
    let variance = if y_vals.len() > 1 {
        y_vals.iter()
            .map(|&v| (v - mean).powi(2))
            .sum::<f64>() / (y_vals.len() - 1) as f64
    } else {
        0.0
    };
    let stddev = variance.sqrt();

    CoverageStats {
        mean,
        median,
        min,
        max,
        stddev,
    }
}

/// Calculate statistics from per-base coverage data
///
/// # Arguments
/// * `coverage` - HashMap mapping position to coverage count
///
/// # Returns
/// * `CoverageStats` - Computed statistics
pub fn calculate_per_base_stats(coverage: &HashMap<u32, u32>) -> CoverageStats {
    if coverage.is_empty() {
        return CoverageStats::default();
    }

    // Extract coverage values
    let mut values: Vec<u32> = coverage.values().copied().collect();
    values.sort_unstable();

    // Calculate mean
    let mean = values.iter().map(|&v| v as f64).sum::<f64>() / values.len() as f64;

    // Calculate median
    let median = if values.len() % 2 == 0 {
        let mid = values.len() / 2;
        (values[mid - 1] + values[mid]) as f64 / 2.0
    } else {
        values[values.len() / 2] as f64
    };

    // Get min and max
    let min = *values.first().unwrap_or(&0) as f64;
    let max = *values.last().unwrap_or(&0) as f64;

    // Calculate standard deviation
    let variance = if values.len() > 1 {
        values.iter()
            .map(|&v| (v as f64 - mean).powi(2))
            .sum::<f64>() / (values.len() - 1) as f64
    } else {
        0.0
    };
    let stddev = variance.sqrt();

    CoverageStats {
        mean,
        median,
        min,
        max,
        stddev,
    }
}

/// Calculate chromosome-wide coverage statistics
/// 
/// # Arguments
/// * `chrom_coverage` - Map of position to coverage for a chromosome
/// * `region_start` - Start position to consider (inclusive)
/// * `region_end` - End position to consider (inclusive)
/// 
/// # Returns
/// * Tuple containing (mean, median, coverage fraction)
#[allow(dead_code)] // Used in future extensions
pub fn calculate_chrom_stats(
    chrom_coverage: &HashMap<u32, u32>,
    region_start: u32,
    region_end: u32,
) -> (f64, f64, f64) {
    // Count of bases with coverage
    let covered_bases = chrom_coverage.keys()
        .filter(|&&pos| pos >= region_start && pos <= region_end)
        .count();
    
    // Total region size
    let region_size = (region_end - region_start + 1) as usize;
    
    // Coverage fraction (0.0 to 1.0)
    let coverage_fraction = if region_size > 0 {
        covered_bases as f64 / region_size as f64
    } else {
        0.0
    };
    
    // Extract statistics for covered bases
    let stats = calculate_per_base_stats(chrom_coverage);
    
    (stats.mean, stats.median, coverage_fraction)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_coverage_stats_calculation() {
        let points = vec![(1, 5.0), (2, 10.0), (3, 15.0)];
        let stats = calculate_coverage_stats(&points);
        
        assert_eq!(stats.mean, 10.0);
        assert_eq!(stats.median, 10.0);
        assert_eq!(stats.min, 5.0);
        assert_eq!(stats.max, 15.0);
        assert!((stats.stddev - 5.0).abs() < 0.001);
    }
    
    #[test]
    fn test_per_base_stats_calculation() {
        let mut coverage = HashMap::new();
        coverage.insert(1, 5);
        coverage.insert(2, 10);
        coverage.insert(3, 15);
        
        let stats = calculate_per_base_stats(&coverage);
        
        assert_eq!(stats.mean, 10.0);
        assert_eq!(stats.median, 10.0);
        assert_eq!(stats.min, 5.0);
        assert_eq!(stats.max, 15.0);
        assert!((stats.stddev - 5.0).abs() < 0.001);
    }
    
    #[test]
    fn test_empty_input() {
        let empty_points: Vec<(i64, f64)> = vec![];
        let empty_coverage = HashMap::new();
        
        let point_stats = calculate_coverage_stats(&empty_points);
        let coverage_stats = calculate_per_base_stats(&empty_coverage);
        
        assert_eq!(point_stats.mean, 0.0);
        assert_eq!(coverage_stats.mean, 0.0);
    }
}