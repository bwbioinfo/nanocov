use std::collections::HashMap;
use std::path::PathBuf;

use nanocov::{parse_bed};

#[test]
fn test_parse_bed_basic() {
    let bed_content = "chr1\t10\t20\nchr1\t30\t40\nchr2\t5\t15\n";
    let tmpfile = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(tmpfile.path(), bed_content).unwrap();
    let regions = parse_bed(&PathBuf::from(tmpfile.path())).unwrap();
    assert_eq!(regions["chr1"], vec![(10, 20), (30, 40)]);
    assert_eq!(regions["chr2"], vec![(5, 15)]);
}

#[test]
fn test_parse_bed_ignores_comments_and_blank() {
    let bed_content = "# comment\nchr1\t1\t2\n\nchr2\t3\t4\n";
    let tmpfile = tempfile::NamedTempFile::new().unwrap();
    std::fs::write(tmpfile.path(), bed_content).unwrap();
    let regions = parse_bed(&PathBuf::from(tmpfile.path())).unwrap();
    assert_eq!(regions["chr1"], vec![(1, 2)]);
    assert_eq!(regions["chr2"], vec![(3, 4)]);
}
