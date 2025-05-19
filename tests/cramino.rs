use assert_cmd::Command;
use std::fs;
use std::path::Path;

#[test]
fn test_cramino_output_basic() {
    // Test basic cramino output generation
    let bam = "test-data/small-test-phased.bam";
    let cramino_output = "test-out/test_cramino_basic.txt";

    // Remove output file if it exists
    let _ = fs::remove_file(cramino_output);

    // Run command with cramino flag and explicit output path
    let mut cmd = Command::cargo_bin("nanocov").unwrap();
    let result = cmd
        .arg("-i")
        .arg(bam)
        .arg("--cramino")
        .arg("--cramino-output")
        .arg(cramino_output)
        .assert();

    // Check command succeeded
    result.success();

    // Verify output file exists
    assert!(
        Path::new(cramino_output).exists(),
        "Cramino output file was not created"
    );

    // Read file contents
    let contents = fs::read_to_string(cramino_output).unwrap();

    // Check for expected fields in the output
    assert!(contents.contains("File name"), "Missing 'File name' field");
    assert!(
        contents.contains("Number of alignments"),
        "Missing 'Number of alignments' field"
    );
    assert!(
        contents.contains("Number of reads"),
        "Missing 'Number of reads' field"
    );
    assert!(
        contents.contains("Yield [Gb]"),
        "Missing 'Yield [Gb]' field"
    );
    assert!(contents.contains("N50"), "Missing 'N50' field");
    assert!(contents.contains("N75"), "Missing 'N75' field");
    assert!(
        contents.contains("Mean length"),
        "Missing 'Mean length' field"
    );
    assert!(
        contents.contains("Median length"),
        "Missing 'Median length' field"
    );
    assert!(contents.contains("Path"), "Missing 'Path' field");
    assert!(
        contents.contains("Creation time"),
        "Missing 'Creation time' field"
    );
}

#[test]
fn test_cramino_output_with_genome_size() {
    // Test cramino output with genome size specified
    let bam = "test-data/small-test-phased.bam";
    let cramino_output = "test-out/test_cramino_genome_size.txt";

    // Remove output file if it exists
    let _ = fs::remove_file(cramino_output);

    // Run command with cramino flag, genome size, and explicit output path
    let mut cmd = Command::cargo_bin("nanocov").unwrap();
    let result = cmd
        .arg("-i")
        .arg(bam)
        .arg("--cramino")
        .arg("--genome-size")
        .arg("3000000") // 3Mb genome size
        .arg("--cramino-output")
        .arg(cramino_output)
        .assert();

    // Check command succeeded
    result.success();

    // Verify output file exists
    assert!(
        Path::new(cramino_output).exists(),
        "Cramino output file was not created"
    );

    // Read file contents
    let contents = fs::read_to_string(cramino_output).unwrap();

    // Check for expected fields and non-zero values
    assert!(
        contents.contains("Mean coverage"),
        "Missing 'Mean coverage' field"
    );

    // Extract the mean coverage value
    let lines: Vec<&str> = contents.lines().collect();
    for line in lines {
        if line.starts_with("Mean coverage") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                let coverage_str = parts[parts.len() - 1];
                let coverage: f64 = coverage_str.parse().unwrap_or(0.0);
                // Coverage should be non-zero when genome size is provided
                assert!(
                    coverage > 0.0,
                    "Mean coverage should be > 0 with genome size specified"
                );
            }
        }
    }
}

#[test]
fn test_cramino_format_consistency() {
    // Test that the cramino output format is consistent with the expected format
    let bam = "test-data/small-test-phased.bam";
    let cramino_output = "test-out/test_cramino_format.txt";

    // Remove output file if it exists
    let _ = fs::remove_file(cramino_output);

    // Run command
    let mut cmd = Command::cargo_bin("nanocov").unwrap();
    let result = cmd
        .arg("-i")
        .arg(bam)
        .arg("--cramino")
        .arg("--cramino-output")
        .arg(cramino_output)
        .assert();

    result.success();

    // Verify output file exists
    assert!(
        Path::new(cramino_output).exists(),
        "Cramino output file was not created"
    );

    // Read file contents
    let contents = fs::read_to_string(cramino_output).unwrap();

    // Check the exact format of the output
    let lines: Vec<&str> = contents.lines().collect();

    // These fields should appear in this specific order
    let expected_fields = [
        "File name",
        "Number of alignments",
        "% from total alignments",
        "Number of reads",
        "Yield [Gb]",
        "Mean coverage",
        "Yield [Gb] (>25kb)",
        "N50",
        "N75",
        "Median length",
        "Mean length",
    ];

    // Check field order
    for (i, &field) in expected_fields.iter().enumerate() {
        if i < lines.len() {
            assert!(
                lines[i].starts_with(field),
                "Line {} should start with '{}', but was '{}'",
                i + 1,
                field,
                lines[i]
            );
        } else {
            panic!("Not enough lines in output, missing field: {}", field);
        }
    }

    // Verify there's an empty line after the main fields
    let empty_line_index = expected_fields.len();
    if lines.len() > empty_line_index {
        assert_eq!(
            lines[empty_line_index], "",
            "Expected empty line after main fields"
        );
    }

    // And the next line should start with "Path"
    if lines.len() > empty_line_index + 1 {
        assert!(
            lines[empty_line_index + 1].starts_with("Path"),
            "Expected 'Path' after empty line, got '{}'",
            lines[empty_line_index + 1]
        );
    }

    // And the next line should start with "Creation time"
    if lines.len() > empty_line_index + 2 {
        assert!(
            lines[empty_line_index + 2].starts_with("Creation time"),
            "Expected 'Creation time' after Path, got '{}'",
            lines[empty_line_index + 2]
        );
    }
}
