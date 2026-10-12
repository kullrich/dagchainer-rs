//! Port of `min_chain_len_filter.pl` — filter DAGchainer `.aligncoords`
//! output, keeping only alignment ("chain") blocks whose header reports
//! at least `min_chain_len` aligned pairs.
//!
//! ## Input format
//!
//! Expects DAGchainer alignment output: a sequence of blocks, each
//! starting with a header line beginning with `#` that contains the
//! substring `num aligned pairs: <digits>`, e.g.:
//!
//! ```text
//! ## alignment chr1 vs. chr2 (num aligned pairs: 15):
//! chr1  geneA  100  200  chr2  geneB  150  250  1.2e-30  50
//! chr1  geneC  300  400  chr2  geneD  350  450  3.4e-20  48
//! ```
//!
//! ## Deviation from the original Perl (deliberate fix)
//!
//! `min_chain_len = 0` is accepted literally. The original Perl's
//! `my $min_chain_len = $ARGV[0] or die $usage;` treated an explicit `0`
//! as falsy and aborted with the usage message. Here, `0` simply means
//! "keep every block" (every chain has `>= 0` aligned pairs).
//!
//! ## Preserved behavior: blank/junk lines are stripped
//!
//! Any input line containing no "word" character (letter, digit, or
//! underscore) is dropped before block-grouping even begins — matching
//! the original Perl's `unless (/\w/) { next; }` guard.
//!
//! ## Preserved behavior: malformed headers abort processing
//!
//! If a block's header line does not contain a parseable
//! `num aligned pairs: <digits>` substring, the original Perl's
//! `confess` aborts the whole script with an error. This port returns a
//! [`MinChainLenError::MissingAlignedPairsCount`] in that case, which the
//! CLI entry point surfaces as a process failure, preserving the
//! "fail loud on malformed input" behavior rather than silently skipping.

use std::io::{BufRead, Write};

use clap::Args as ClapArgs;
use thiserror::Error;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Minimum number of aligned pairs a chain must have to be kept.
    /// `0` is valid (keeps every chain).
    #[arg(value_name = "MIN_CHAIN_LENGTH")]
    pub min_chain_len: usize,
}

#[derive(Debug, Error)]
pub enum MinChainLenError {
    #[error(
        "cannot determine num aligned pairs from header: {header:?}\n\nfull block:\n{block}"
    )]
    MissingAlignedPairsCount { header: String, block: String },
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

/// Returns true if `line` contains at least one ASCII letter, digit, or
/// underscore — equivalent to Perl's `/\w/` test used to strip
/// blank/junk lines.
fn has_word_char(line: &str) -> bool {
    line.chars().any(|c| c.is_ascii_alphanumeric() || c == '_')
}

/// Extract the aligned-pair count from a header line containing the
/// substring `num aligned pairs: <digits>`.
fn parse_num_pairs(header: &str) -> Option<usize> {
    const MARKER: &str = "num aligned pairs: ";
    let idx = header.find(MARKER)?;
    let after = &header[idx + MARKER.len()..];
    let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
    if digits.is_empty() {
        None
    } else {
        digits.parse().ok()
    }
}

/// Render a chain block (header + data lines) back into text, with the
/// same trailing blank-line separator the original Perl produces via
/// `print $text . "\n"` (where `$text` already ends in a newline from
/// the raw, un-chomped input lines).
fn render_block(lines: &[String]) -> String {
    let mut out = String::new();
    for line in lines {
        out.push_str(line);
        out.push('\n');
    }
    out.push('\n'); // the extra separator blank line.
    out
}

/// Decide whether to keep a chain block, returning its rendered text if
/// its header's aligned-pair count is `>= min_chain_len`.
fn process_chain_block(
    lines: &[String],
    min_chain_len: usize,
) -> Result<Option<String>, MinChainLenError> {
    let header = lines
        .first()
        .expect("chain block must have at least a header line");

    let num_pairs = parse_num_pairs(header).ok_or_else(|| {
        MinChainLenError::MissingAlignedPairsCount {
            header: header.clone(),
            block: lines.join("\n"),
        }
    })?;

    if num_pairs < min_chain_len {
        return Ok(None);
    }

    Ok(Some(render_block(lines)))
}

/// Core logic, generic over any `BufRead`/`Write`, mirroring the Perl
/// script's stdin -> stdout filter behaviour.
pub fn filter_min_chain_len<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
    min_chain_len: usize,
) -> Result<(), MinChainLenError> {
    let mut current_block: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if !has_word_char(&line) {
            continue;
        }

        if line.starts_with('#') {
            if !current_block.is_empty() {
                if let Some(text) = process_chain_block(&current_block, min_chain_len)? {
                    write!(writer, "{text}")?;
                }
            }
            current_block = vec![line];
        } else {
            current_block.push(line);
        }
    }

    if !current_block.is_empty() {
        if let Some(text) = process_chain_block(&current_block, min_chain_len)? {
            write!(writer, "{text}")?;
        }
    }

    Ok(())
}

/// CLI entry point, reading from stdin and writing to stdout — a drop-in
/// replacement for `min_chain_len_filter.pl min_chain_length < file.aligncoords`.
pub fn run(args: Args) -> anyhow::Result<()> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    filter_min_chain_len(stdin.lock(), stdout.lock(), args.min_chain_len).map_err(Into::into)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_blocks_meeting_threshold_drops_others() {
        let input = "\
## alignment chr1 vs. chr2 (num aligned pairs: 2):
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
chr1\tgeneC\t300\t400\tchr2\tgeneD\t350\t450\t1e-20\t48
## alignment chr1 vs. chr3 (num aligned pairs: 1):
chr1\tgeneE\t500\t600\tchr3\tgeneF\t550\t650\t1e-10\t40
";
        let mut out = Vec::new();
        filter_min_chain_len(input.as_bytes(), &mut out, 2).unwrap();
        let out = String::from_utf8(out).unwrap();

        let expected = "\
## alignment chr1 vs. chr2 (num aligned pairs: 2):
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
chr1\tgeneC\t300\t400\tchr2\tgeneD\t350\t450\t1e-20\t48

";
        assert_eq!(out, expected);
    }

    #[test]
    fn min_chain_len_zero_keeps_every_block() {
        let input = "\
## alignment a vs. b (num aligned pairs: 0):
## alignment c vs. d (num aligned pairs: 1):
chr1\tgeneA\t1\t2\tchr2\tgeneB\t1\t2\t1e-5\t10
";
        let mut out = Vec::new();
        filter_min_chain_len(input.as_bytes(), &mut out, 0).unwrap();
        let out = String::from_utf8(out).unwrap();

        // Both blocks survive, including the header-only one with 0 pairs.
        assert!(out.contains("num aligned pairs: 0"));
        assert!(out.contains("num aligned pairs: 1"));
    }

    #[test]
    fn blank_and_junk_lines_are_stripped() {
        let input = "\
## alignment a vs. b (num aligned pairs: 1):

chr1\tgeneA\t1\t2\tchr2\tgeneB\t1\t2\t1e-5\t10
----
";
        let mut out = Vec::new();
        filter_min_chain_len(input.as_bytes(), &mut out, 1).unwrap();
        let out = String::from_utf8(out).unwrap();

        // The blank line and the "----" line (no word chars) must not
        // appear in the output, and must not break block grouping.
        let expected = "\
## alignment a vs. b (num aligned pairs: 1):
chr1\tgeneA\t1\t2\tchr2\tgeneB\t1\t2\t1e-5\t10

";
        assert_eq!(out, expected);
    }

    #[test]
    fn header_only_block_with_no_data_lines_is_handled() {
        let input = "\
## alignment a vs. b (num aligned pairs: 0):
## alignment c vs. d (num aligned pairs: 5):
chr1\tgeneA\t1\t2\tchr2\tgeneB\t1\t2\t1e-5\t10
";
        let mut out = Vec::new();
        filter_min_chain_len(input.as_bytes(), &mut out, 1).unwrap();
        let out = String::from_utf8(out).unwrap();

        // First block (0 pairs) dropped; second (5 pairs) kept.
        assert!(!out.contains("num aligned pairs: 0"));
        assert!(out.contains("num aligned pairs: 5"));
    }

    #[test]
    fn malformed_header_returns_error() {
        let input = "## this header has no pair count\nchr1\tgeneA\t1\t2\tchr2\tgeneB\t1\t2\t1e-5\t10\n";
        let mut out = Vec::new();
        let result = filter_min_chain_len(input.as_bytes(), &mut out, 1);
        assert!(matches!(
            result,
            Err(MinChainLenError::MissingAlignedPairsCount { .. })
        ));
    }

    #[test]
    fn empty_input_produces_empty_output() {
        let mut out = Vec::new();
        filter_min_chain_len("".as_bytes(), &mut out, 1).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn parse_num_pairs_extracts_digits_correctly() {
        assert_eq!(
            parse_num_pairs("## alignment x vs. y (num aligned pairs: 123):"),
            Some(123)
        );
        assert_eq!(parse_num_pairs("## no count here"), None);
    }
}