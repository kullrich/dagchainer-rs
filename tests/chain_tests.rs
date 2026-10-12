//! Integration tests for the `dagchainer` CLI binary and the
//! `dagchainer_rs` library, verified against known-correct output.

use assert_cmd::Command;
use predicates::str::contains;

use dagchainer_rs::{chain_matches, MatchPoint, Params};

fn test_params() -> Params {
    Params {
        indel_score: -1.0,
        gap_open_penalty: -1.0,
        bp_gap_size: 5,
        min_alignment_score: 10.0,
        max_dist_between_matches: 100_000,
        reverse_order: false,
    }
}

/// Library-level test: verifies the DP scoring directly. Numbers below
/// reflect the corrected, truncating-integer-division num_gaps
/// computation (see `src/lib.rs` module docs), confirmed by diffing
/// against the real compiled dagchainer C++ binary's output.
#[test]
fn library_scores_match_expected_dp_math() {
    let matches = vec![
        MatchPoint { pair_id: 0, x: 10, y: 12, score: 50.0 },
        MatchPoint { pair_id: 1, x: 20, y: 22, score: 50.0 },
        MatchPoint { pair_id: 2, x: 30, y: 31, score: 50.0 },
        MatchPoint { pair_id: 3, x: 5000, y: 5002, score: 50.0 },
    ];

    let alignments = chain_matches(&matches, &test_params()).unwrap();

    assert_eq!(alignments.len(), 2);

    let a1 = &alignments[0];
    assert_eq!(a1.entries.len(), 3);
    assert_eq!(a1.entries[0].cumulative_score, 50.0);
    assert_eq!(a1.entries[1].cumulative_score, 98.0);
    assert_eq!(a1.entries[2].cumulative_score, 146.0);
    assert_eq!(a1.score, 146.0);

    let a2 = &alignments[1];
    assert_eq!(a2.entries.len(), 1);
    assert_eq!(a2.entries[0].pair_id, 3);
    assert_eq!(a2.score, 50.0);
}

#[test]
fn cli_runs_and_reports_alignments() {
    let mut cmd = Command::cargo_bin("dagchainer").unwrap();
    cmd.args([
        "-G", "5", "-O", "-1", "-E", "-1", "-S", "10", "-D", "100000",
        "-F", "tests/fixtures/small.input",
    ]);

    cmd.assert()
        .success()
        .stdout(contains(">Alignment #1  score = 146.0"))
        .stdout(contains(">Alignment #2  score = 50.0"));
}

#[test]
fn cli_exact_output_matches_known_good() {
    let expected = "\
>Alignment #1  score = 146.0
  0: 0     10     12    50.0    50.0
  1: 1     20     22    50.0    98.0
  2: 2     30     31    50.0   146.0
>Alignment #2  score = 50.0
  0: 3   5000   5002    50.0    50.0
";

    let mut cmd = Command::cargo_bin("dagchainer").unwrap();
    cmd.args([
        "-G", "5", "-O", "-1", "-E", "-1", "-S", "10", "-D", "100000",
        "-F", "tests/fixtures/small.input",
    ]);

    let output = cmd.output().unwrap();
    let stdout = String::from_utf8_lossy(&output.stdout);

    assert_eq!(stdout, expected, "CLI output did not match expected fixture");
}

#[test]
fn cli_missing_file_fails_cleanly() {
    let mut cmd = Command::cargo_bin("dagchainer").unwrap();
    cmd.args([
        "-G", "5", "-O", "-1", "-E", "-1", "-S", "10", "-D", "100000",
        "-F", "tests/fixtures/does_not_exist.input",
    ]);
    cmd.assert().failure();
}

#[test]
fn cli_rejects_invalid_gap_open_penalty() {
    let mut cmd = Command::cargo_bin("dagchainer").unwrap();
    cmd.args([
        "-G", "5", "-O", "1", "-E", "-1", "-S", "10", "-D", "100000",
        "-F", "tests/fixtures/small.input",
    ]);
    cmd.assert().failure();
}

#[test]
fn cli_accepts_negative_penalty_values() {
    let mut cmd = Command::cargo_bin("dagchainer").unwrap();
    cmd.args([
        "-G", "5", "-O", "-1", "-E", "-1", "-S", "10", "-D", "100000",
        "-F", "tests/fixtures/small.input",
    ]);
    cmd.assert().success();
}