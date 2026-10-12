//! Port of `sum_genome_chain_nr_content.pl` — sum the non-redundant
//! (overlap-collapsed) base-pair content covered by DAGchainer chain
//! spans, separately for each side (genome A / genome B).
//!
//! ## Input format
//!
//! One line per chain, tab-separated, matching
//! [`crate::tools::to_chain_spans`]'s output:
//! `molA  lendA  rendA  molB  lendB  rendB  orient  num_genes`.
//!
//! Interval collapsing is delegated to [`crate::overlap_piler::collapse_coordsets`],
//! a faithful port of the original's `Overlap_piler.pm` dependency.
//!
//! ## Preserved behavior: malformed lines are fatal
//!
//! Any line starting with `#`, or containing no word character
//! (letter/digit/underscore), causes the original Perl to `die`
//! immediately. This is preserved: [`sum_nr_content`] returns an error
//! in that case rather than skipping the line, since this script's
//! input format (chain-span output) is expected to be strictly
//! well-formed with no headers.
//!
//! ## Deviation from the original Perl (deliberate improvements)
//!
//! - Data lines with the wrong field count, or non-numeric coordinate
//!   fields, are treated as hard errors (rather than Perl's silent
//!   `undef`-in-numeric-context corruption under `use warnings`),
//!   consistent with this script's own "fail loud on malformed input"
//!   philosophy.
//! - Molecule iteration order is sorted ascending (via `BTreeMap`),
//!   rather than Perl's randomized hash key order (Perl 5.18+). This
//!   only affects the order of `--verbose` pile lines; computed sums are
//!   unaffected by iteration order.
//! - The CLI's optional positional `[verbose]` argument became a proper
//!   `-v`/`--verbose` flag, for idiomatic `clap` ergonomics — a CLI
//!   interface change, not a logic change.

use std::collections::BTreeMap;
use std::io::{BufRead, Write};
use std::path::PathBuf;

use clap::Args as ClapArgs;
use thiserror::Error;

use crate::overlap_piler::collapse_coordsets;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Chain spans file (output of `to-chain-spans`).
    #[arg(value_name = "CHAIN_SPANS_FILE")]
    pub chain_spans_file: PathBuf,

    /// Print each non-redundant coordinate "pile" as it's computed.
    #[arg(short = 'v', long = "verbose")]
    pub verbose: bool,
}

#[derive(Debug, Error)]
pub enum ChainSpansError {
    #[error("cannot parse chain span line: {line:?}")]
    InvalidLine { line: String },
    #[error("malformed chain span line (expected 8 tab-separated fields): {line:?}")]
    MalformedDataLine { line: String },
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

/// Returns true if `line` contains at least one ASCII letter, digit, or
/// underscore — equivalent to Perl's `/\w/` test.
fn has_word_char(line: &str) -> bool {
    line.chars().any(|c| c.is_ascii_alphanumeric() || c == '_')
}

/// One parsed chain-span data line's fields actually used by this
/// script (orientation and gene count are parsed-but-unused upstream,
/// and are not needed here at all).
struct Row {
    mol_a: String,
    lend_a: i64,
    rend_a: i64,
    mol_b: String,
    lend_b: i64,
    rend_b: i64,
}

fn parse_line(line: &str) -> Result<Row, ChainSpansError> {
    if line.starts_with('#') || !has_word_char(line) {
        return Err(ChainSpansError::InvalidLine { line: line.to_string() });
    }

    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return Err(ChainSpansError::MalformedDataLine { line: line.to_string() });
    }

    let parse_coord = |s: &str| s.parse::<i64>().ok();
    let (Some(lend_a), Some(rend_a), Some(lend_b), Some(rend_b)) = (
        parse_coord(fields[1]),
        parse_coord(fields[2]),
        parse_coord(fields[4]),
        parse_coord(fields[5]),
    ) else {
        return Err(ChainSpansError::MalformedDataLine { line: line.to_string() });
    };

    Ok(Row {
        mol_a: fields[0].to_string(),
        lend_a,
        rend_a,
        mol_b: fields[3].to_string(),
        lend_b,
        rend_b,
    })
}

/// Compute the non-redundant summed length across all molecules'
/// coordinate piles for one side (genome A or B), optionally writing
/// each pile's details if `verbose` is set. Mirrors
/// `sum_coord_piles($org_token, \%data)` in the original Perl.
pub fn sum_coord_piles<W: Write>(
    org_token: &str,
    data: &BTreeMap<String, Vec<(i64, i64)>>,
    verbose: bool,
    writer: &mut W,
) -> std::io::Result<i64> {
    let mut sum_bp: i64 = 0;

    for (mol, coordpairs) in data {
        let piles = collapse_coordsets(coordpairs);
        for (lend, rend) in piles {
            let pile_len = rend - lend + 1;
            sum_bp += pile_len;
            if verbose {
                writeln!(writer, "{org_token}\t{mol}\t{lend}\t{rend}")?;
            }
        }
    }

    Ok(sum_bp)
}

/// Core logic: parse all chain-span lines, then compute and return
/// `(sum_bp_genome_a, sum_bp_genome_b)`, writing any verbose pile output
/// to `writer` as it's computed (A's piles fully before B's, matching
/// the original's call order).
pub fn sum_nr_content<R: BufRead, W: Write>(
    reader: R,
    verbose: bool,
    writer: &mut W,
) -> Result<(i64, i64), ChainSpansError> {
    let mut genome_a_data: BTreeMap<String, Vec<(i64, i64)>> = BTreeMap::new();
    let mut genome_b_data: BTreeMap<String, Vec<(i64, i64)>> = BTreeMap::new();

    for line in reader.lines() {
        let line = line?;
        let row = parse_line(&line)?;
        genome_a_data.entry(row.mol_a).or_default().push((row.lend_a, row.rend_a));
        genome_b_data.entry(row.mol_b).or_default().push((row.lend_b, row.rend_b));
    }

    let sum_bp_a = sum_coord_piles("A", &genome_a_data, verbose, writer)?;
    let sum_bp_b = sum_coord_piles("B", &genome_b_data, verbose, writer)?;

    Ok((sum_bp_a, sum_bp_b))
}

/// CLI entry point — a drop-in replacement for
/// `sum_genome_chain_nr_content.pl chain_spans_file [verbose]`
/// (with `[verbose]` now expressed as a `-v`/`--verbose` flag).
pub fn run(args: Args) -> anyhow::Result<()> {
    let file = std::fs::File::open(&args.chain_spans_file).map_err(|e| {
        anyhow::anyhow!("Error, cannot open file {}: {e}", args.chain_spans_file.display())
    })?;

    let stdout = std::io::stdout();
    let mut out = stdout.lock();

    let (sum_bp_a, sum_bp_b) = sum_nr_content(std::io::BufReader::new(file), args.verbose, &mut out)?;

    writeln!(out, "GenomeA: {sum_bp_a}")?;
    writeln!(out, "GenomeB: {sum_bp_b}")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_sum_without_overlap() {
        let input = "chr1\t100\t199\tchr2\t1\t50\t+\t3\nchr1\t300\t399\tchr2\t200\t249\t+\t2\n";
        let mut out = Vec::new();
        let (sum_a, sum_b) = sum_nr_content(input.as_bytes(), false, &mut out).unwrap();
        assert_eq!(sum_a, 200);
        assert_eq!(sum_b, 100);
        assert!(out.is_empty());
    }

    #[test]
    fn overlapping_chains_are_not_double_counted() {
        let input = "chr1\t100\t300\tchr2\t1\t100\t+\t5\nchr1\t200\t400\tchr2\t50\t150\t+\t4\n";
        let mut out = Vec::new();
        let (sum_a, sum_b) = sum_nr_content(input.as_bytes(), false, &mut out).unwrap();
        assert_eq!(sum_a, 301);
        assert_eq!(sum_b, 150);
    }

    #[test]
    fn adjacent_non_overlapping_chains_sum_correctly_despite_separate_piles() {
        // 199 and 200 share no coordinate under the real algorithm, so
        // these form two separate piles, but the sum is identical to
        // treating them as one contiguous span (100-299), since there's
        // no gap and no overlap to double-count.
        let input = "chr1\t100\t199\tchr2\t1\t1\t+\t1\nchr1\t200\t299\tchr2\t2\t2\t+\t1\n";
        let mut out = Vec::new();
        let (sum_a, _) = sum_nr_content(input.as_bytes(), false, &mut out).unwrap();
        assert_eq!(sum_a, 200); // 100bp + 100bp
    }

    #[test]
    fn verbose_output_prints_genome_a_piles_before_genome_b() {
        let input = "chr1\t100\t199\tchr2\t1\t50\t+\t3\n";
        let mut out = Vec::new();
        sum_nr_content(input.as_bytes(), true, &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();

        let a_idx = out.find("A\tchr1").unwrap();
        let b_idx = out.find("B\tchr2").unwrap();
        assert!(a_idx < b_idx);
    }

    #[test]
    fn verbose_pile_lines_have_expected_format() {
        let input = "chr1\t100\t199\tchr2\t1\t50\t+\t3\n";
        let mut out = Vec::new();
        sum_nr_content(input.as_bytes(), true, &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.contains("A\tchr1\t100\t199\n"));
        assert!(out.contains("B\tchr2\t1\t50\n"));
    }

    #[test]
    fn verbose_output_merges_overlapping_piles_into_one_line() {
        let input = "chr1\t100\t300\tchr2\t1\t10\t+\t1\nchr1\t200\t400\tchr2\t20\t30\t+\t1\n";
        let mut out = Vec::new();
        sum_nr_content(input.as_bytes(), true, &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.contains("A\tchr1\t100\t400\n"));
        assert_eq!(out.matches("A\tchr1").count(), 1); // merged into a single pile line
    }

    #[test]
    fn multiple_molecules_sorted_deterministically_in_verbose_output() {
        let input = "chrZ\t1\t10\tchrY\t1\t10\t+\t1\nchrA\t1\t10\tchrB\t1\t10\t+\t1\n";
        let mut out = Vec::new();
        sum_nr_content(input.as_bytes(), true, &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();

        let idx_a = out.find("A\tchrA").unwrap();
        let idx_z = out.find("A\tchrZ").unwrap();
        assert!(idx_a < idx_z);
    }

    #[test]
    fn header_line_is_a_hard_error() {
        let input = "# this is a header\nchr1\t1\t10\tchr2\t1\t10\t+\t1\n";
        let mut out = Vec::new();
        let result = sum_nr_content(input.as_bytes(), false, &mut out);
        assert!(matches!(result, Err(ChainSpansError::InvalidLine { .. })));
    }

    #[test]
    fn blank_line_is_a_hard_error() {
        let input = "chr1\t1\t10\tchr2\t1\t10\t+\t1\n\n";
        let mut out = Vec::new();
        let result = sum_nr_content(input.as_bytes(), false, &mut out);
        assert!(matches!(result, Err(ChainSpansError::InvalidLine { .. })));
    }

    #[test]
    fn malformed_field_count_is_a_hard_error() {
        let input = "chr1\t1\t10\tchr2\n";
        let mut out = Vec::new();
        let result = sum_nr_content(input.as_bytes(), false, &mut out);
        assert!(matches!(result, Err(ChainSpansError::MalformedDataLine { .. })));
    }

    #[test]
    fn non_numeric_coordinate_is_a_hard_error() {
        let input = "chr1\tnot_a_number\t10\tchr2\t1\t10\t+\t1\n";
        let mut out = Vec::new();
        let result = sum_nr_content(input.as_bytes(), false, &mut out);
        assert!(matches!(result, Err(ChainSpansError::MalformedDataLine { .. })));
    }

    #[test]
    fn single_interval_pile_length_is_inclusive() {
        let input = "chr1\t100\t100\tchr2\t1\t1\t+\t1\n";
        let mut out = Vec::new();
        let (sum_a, sum_b) = sum_nr_content(input.as_bytes(), false, &mut out).unwrap();
        assert_eq!(sum_a, 1);
        assert_eq!(sum_b, 1);
    }
}