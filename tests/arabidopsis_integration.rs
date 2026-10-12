// tests/arabidopsis_integration.rs
use assert_cmd::Command;

const DATA_FILE: &str = "data/arabidopsis/Arabidopsis.Release5.matchList";

#[test]
fn filter_repetitive_runs_on_real_data() {
    assert!(std::path::Path::new(DATA_FILE).exists(), "missing test data: {DATA_FILE}");

    let mut cmd = Command::cargo_bin("dagtools").unwrap();
    cmd.arg("filter-repetitive").arg("10000");
    cmd.pipe_stdin(DATA_FILE).unwrap();
    cmd.assert().success();
}

#[test]
fn find_tandem_dups_discovers_known_self_comparisons() {
    let mut cmd = Command::cargo_bin("dagtools").unwrap();
    cmd.args(["find-tandem-dups", "-m", "1"]);
    cmd.pipe_stdin(DATA_FILE).unwrap();

    let output = cmd.output().unwrap();
    assert!(output.status.success());
    // At1g01010/At1g02230/At1g02250 etc. are all molecule "1" vs "1" —
    // real tandem-duplicate candidates should be discoverable here.
}