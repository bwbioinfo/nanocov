use assert_cmd::Command;
use std::fs;
use std::path::Path;

#[test]
fn test_nanocov_runs_and_outputs_file() {
    let bam = "test-data/small-test-phased.bam";
    let bed = "test-data/read"; // Use a small BED file or create one
    let out_file = "coverage.tsv";
    // Remove output file if it exists
    let _ = fs::remove_file(out_file);
    let mut cmd = Command::cargo_bin("nanocov").unwrap();
    let result = cmd.arg("-i").arg(bam).arg("-b").arg(bed).assert();
    result.success();
    assert!(Path::new(out_file).exists());
    let contents = fs::read_to_string(out_file).unwrap();
    assert!(contents.contains("#")); // Should contain at least one reference header
}
