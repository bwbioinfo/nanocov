// Test script to verify log scale plotting functionality
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create test coverage data with a wide range of values (good for log scale)
    let mut coverage: HashMap<u32, u32> = HashMap::new();
    
    // Add some very low coverage areas
    for i in 1000..1100 {
        coverage.insert(i, 1);
    }
    
    // Add some medium coverage areas
    for i in 2000..2100 {
        coverage.insert(i, 50);
    }
    
    // Add some high coverage areas
    for i in 3000..3100 {
        coverage.insert(i, 1000);
    }
    
    // Add some very high coverage areas  
    for i in 4000..4100 {
        coverage.insert(i, 10000);
    }
    
    println!("Testing log scale plotting with coverage range: 1 to 10,000");
    
    // Test linear scale
    nanocov::plotting::plot_per_base_coverage(
        "test_chr",
        &coverage,
        "test_linear_scale.png",
        None,
        false,
        false, // linear scale
    )?;
    
    println!("Generated linear scale plot: test_linear_scale.png");
    
    // Test log scale  
    nanocov::plotting::plot_per_base_coverage(
        "test_chr", 
        &coverage,
        "test_log_scale.png",
        None,
        false,
        true, // log scale
    )?;
    
    println!("Generated log scale plot: test_log_scale.png");
    println!("Compare the two plots to verify log scale is working correctly");
    
    Ok(())
}
