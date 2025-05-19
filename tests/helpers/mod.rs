// tests/helpers/mod.rs
// Helper functions for testing nanocov

use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::io::{BufRead, BufReader};

/// Create a directory if it doesn't exist
pub fn ensure_dir_exists(path: &str) -> std::io::Result<()> {
    if !Path::new(path).exists() {
        fs::create_dir_all(path)
    } else {
        Ok(())
    }
}

/// Clean up output files before testing
pub fn clean_test_output(prefix: &str) {
    // Find and remove files with the given prefix
    if let Ok(entries) = fs::read_dir("test-out") {
        for entry in entries.flatten() {
            let path = entry.path();
            if let Some(filename) = path.file_name() {
                if let Some(filename_str) = filename.to_str() {
                    if filename_str.starts_with(prefix) {
                        let _ = fs::remove_file(path);
                    }
                }
            }
        }
    }
}

/// Create test coverage data for multiple chromosomes
pub fn create_test_coverage_data() -> HashMap<String, HashMap<u32, u32>> {
    let mut data = HashMap::new();
    
    // Create chromosome 1 with higher coverage
    let mut chr1 = HashMap::new();
    for i in 0..1000u32 {
        chr1.insert(i, 20 + (i % 10));
    }
    data.insert("chr1".to_string(), chr1);
    
    // Create chromosome 2 with medium coverage
    let mut chr2 = HashMap::new();
    for i in 0..800u32 {
        chr2.insert(i, 10 + (i % 5));
    }
    data.insert("chr2".to_string(), chr2);
    
    // Create X chromosome with lower coverage
    let mut chrx = HashMap::new();
    for i in 0..600u32 {
        chrx.insert(i, 5 + (i % 3));
    }
    data.insert("chrX".to_string(), chrx);
    
    // Add a small chromosome with variable coverage
    let mut chr22 = HashMap::new();
    for i in 0..400u32 {
        let pos = i * 10; // Sparse coverage
        chr22.insert(pos, (i % 30) + 1); // More variable coverage
    }
    data.insert("chr22".to_string(), chr22);
    
    data
}

/// Get a list of all available themes
pub fn get_themes() -> Vec<&'static str> {
    vec!["latte", "frappe", "nord", "gruvbox"]
}

/// Verify that PNG files were created and have non-zero size
pub fn verify_png_files(prefix: &str) -> bool {
    let mut found = false;
    
    if let Ok(entries) = fs::read_dir("test-out") {
        for entry in entries.flatten() {
            let path = entry.path();
            if let Some(ext) = path.extension() {
                if ext == "png" {
                    if let Some(filename) = path.file_name() {
                        if let Some(filename_str) = filename.to_str() {
                            if filename_str.starts_with(prefix) {
                                if let Ok(metadata) = fs::metadata(&path) {
                                    if metadata.len() > 0 {
                                        found = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    found
}

/// Parse and validate a cramino output file, returning HashMap of fields to values
pub fn parse_cramino_file(filepath: &str) -> Result<HashMap<String, String>, String> {
    if !Path::new(filepath).exists() {
        return Err(format!("Cramino file does not exist: {}", filepath));
    }
    
    let file = match fs::File::open(filepath) {
        Ok(f) => f,
        Err(e) => return Err(format!("Could not open cramino file: {}", e)),
    };
    
    let reader = BufReader::new(file);
    let mut field_map = HashMap::new();
    
    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => return Err(format!("Error reading line: {}", e)),
        };
        
        // Skip empty lines
        if line.trim().is_empty() {
            continue;
        }
        
        // Parse field and value
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let field = parts[0].trim().to_string();
            let value = parts[1].trim().to_string();
            field_map.insert(field, value);
        }
    }
    
    // Validate required fields
    let required_fields = vec![
        "File name", 
        "Number of alignments", 
        "Number of reads", 
        "Yield [Gb]", 
        "Mean coverage", 
        "N50", 
        "N75", 
        "Mean length", 
        "Median length",
        "Path",
        "Creation time"
    ];
    
    for field in required_fields {
        if !field_map.contains_key(field) {
            return Err(format!("Required field missing in cramino output: {}", field));
        }
    }
    
    Ok(field_map)
}

/// Verify cramino output contains expected numeric values
pub fn verify_cramino_numeric_fields(field_map: &HashMap<String, String>) -> Result<(), String> {
    let numeric_fields = vec![
        "Number of alignments",
        "Number of reads",
        "Yield [Gb]",
        "Mean coverage", 
        "N50", 
        "N75", 
        "Mean length", 
        "Median length"
    ];
    
    for field in numeric_fields {
        if let Some(value) = field_map.get(field) {
            // Try to parse as f64 (handles both integer and decimal cases)
            if value.parse::<f64>().is_err() {
                return Err(format!("Field '{}' has non-numeric value: {}", field, value));
            }
        }
    }
    
    Ok(())
}

/// Create an empty cramino output structure for testing
pub fn create_empty_cramino_output() -> String {
    let timestamp = chrono::Local::now().format("%d/%m/%Y %H:%M:%S").to_string();
    
    format!(
        "File name\tempty.bam\n\
        Number of alignments\t0\n\
        % from total alignments\t0.00\n\
        Number of reads\t0\n\
        Yield [Gb]\t0.00\n\
        Mean coverage\t0.00\n\
        Yield [Gb] (>25kb)\t0.00\n\
        N50\t0\n\
        N75\t0\n\
        Median length\t0.00\n\
        Mean length\t0.00\n\
        \n\
        Path\ttest-data/empty.bam\n\
        Creation time\t{}", 
        timestamp
    )
}