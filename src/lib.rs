use std::collections::HashMap;
use std::path::PathBuf;
use std::io::BufRead;

pub fn parse_bed(path: &PathBuf) -> Result<HashMap<String, Vec<(u32, u32)>>, Box<dyn std::error::Error>> {
    let mut regions: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() < 3 {
            continue;
        }
        let chrom = fields[0].to_string();
        let start: u32 = fields[1].parse()?;
        let end: u32 = fields[2].parse()?;
        regions.entry(chrom).or_default().push((start, end));
    }
    Ok(regions)
}
