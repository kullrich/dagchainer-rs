//! Port of `dagchainer_to_chain_spans.pl` — collapse each DAGchainer
//! alignment chain down to its overall coordinate span on each side.
//!
//! ## Input format
//!
//! A sequence of blocks, each starting with a `#`-prefixed header line
//! containing the two molecule/chromosome names (whitespace-delimited
//! tokens 3 and 5, i.e. index 2 and 4), optionally the literal substring
//! `(reverse)`, and `num aligned pairs: <digits>`, e.g.:
//!
//! ```text
//! ## alignment chr1 vs. chr2 (reverse) >Alignment #1  score = 100.0 (num aligned pairs: 3):
//! chr1  geneA  100  200  chr2  geneB  150  250  1e-30  50
//! chr1  geneC  300  400  chr2  geneD  350  450  1e-20  98
//! ```
//!
//! followed by tab-separated data lines:
//! `contigA geneA lendA rendA contigB geneB lendB rendB evalue dagchainscore`
//! — only `lendA, rendA, lendB, rendB` are actually used.
//!
//! ## Output
//!
//! One line per chain: `accA  lendA  rendA  accB  lendB  rendB  orient  num_gene_pairs`,
//! where `lendX`/`rendX` are the min/max of all collected coordinates on
//! that side across every data line in the chain.
//!
//! ## Deviations from the original Perl (deliberate fixes)
//!
//! 1. **Blank lines are skipped**, not treated as data. The original
//!    Perl has no blank-line guard at all; a blank line falls into the
//!    "data line" branch, where `split(/\t/, "")` yields no fields, so
//!    all of `lendA/rendA/lendB/rendB` become Perl `undef`, which is
//!    then pushed into the coordinate list and silently treated as `0`
//!    in the later numeric sort — corrupting that chain's computed
//!    span. This matters in practice because
//!    [`crate::tools::min_chain_len_filter`] (a natural upstream
//!    pipeline partner, given the naming convention) inserts exactly
//!    such a blank line after every surviving chain block. This port
//!    skips blank lines entirely, making the two tools safely
//!    composable.
//! 2. **Malformed data lines** (fewer than 8 tab-separated fields, or
//!    non-numeric coordinate fields) are skipped with a warning, rather
//!    than silently injecting `0` coordinates via Perl's
//!    `undef`-in-numeric-context coercion.
//! 3. **A data line before any header** is a hard
//!    [`ToChainSpansError::DataBeforeHeader`] error, rather than relying
//!    on Perl's autovivification of `$curr_chain` from `undef` into a
//!    malformed in-progress hash.
//!
//! ## Preserved behavior: empty-coordinate spans print as blank fields
//!
//! A header with zero following data lines produces empty `lendA/rendA/
//! lendB/rendB` output fields (mirroring Perl's `shift`/`pop` on an
//! empty array both yielding `undef`, which stringifies to `""`).

use std::io::{BufRead, Read, Write};
use std::path::PathBuf;

use clap::Args as ClapArgs;
use thiserror::Error;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Input file(s). If omitted, reads from stdin. Multiple files are
    /// treated as one continuous logical stream (mirroring Perl's `<>`
    /// diamond operator), which matters here since chain state persists
    /// across the whole input.
    #[arg(value_name = "FILE")]
    pub files: Vec<PathBuf>,
}

#[derive(Debug, Error)]
pub enum ToChainSpansError {
    #[error("malformed header line (fewer than 5 whitespace-separated tokens): {header:?}")]
    MalformedHeader { header: String },
    #[error("cannot determine num aligned pairs from header: {header:?}")]
    MissingAlignedPairsCount { header: String },
    #[error("data line encountered before any header line: {line:?}")]
    DataBeforeHeader { line: String },
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

/// One in-progress or completed alignment chain.
#[derive(Debug, Clone)]
struct Chain {
    acc_a: String,
    acc_b: String,
    orient: char,
    num_gene_pairs: u64,
    coords_a: Vec<i64>,
    coords_b: Vec<i64>,
}

/// Extract the aligned-pair count from a header line containing the
/// substring `num aligned pairs: <digits>`.
fn parse_num_pairs(header: &str) -> Option<u64> {
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

/// Parse a `#`-prefixed header line into a new, empty [`Chain`].
fn parse_header(line: &str) -> Result<Chain, ToChainSpansError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 5 {
        return Err(ToChainSpansError::MalformedHeader {
            header: line.to_string(),
        });
    }

    let acc_a = tokens[2].to_string();
    let acc_b = tokens[4].to_string();
    let orient = if line.contains("(reverse)") { '-' } else { '+' };
    let num_gene_pairs = parse_num_pairs(line).ok_or_else(|| {
        ToChainSpansError::MissingAlignedPairsCount {
            header: line.to_string(),
        }
    })?;

    Ok(Chain {
        acc_a,
        acc_b,
        orient,
        num_gene_pairs,
        coords_a: Vec::new(),
        coords_b: Vec::new(),
    })
}

/// Parse one tab-separated data line, appending its coordinates to
/// `chain` if well-formed. Returns `false` (and logs a warning) for
/// malformed lines, which are then skipped by the caller.
fn parse_data_line(line: &str, chain: &mut Chain) -> bool {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        log::warn!("Skipping data line with fewer than 8 tab-separated fields: {line:?}");
        return false;
    }

    let parse_coord = |s: &str| -> Option<i64> { s.parse().ok() };

    let (Some(lend_a), Some(rend_a), Some(lend_b), Some(rend_b)) = (
        parse_coord(fields[2]),
        parse_coord(fields[3]),
        parse_coord(fields[6]),
        parse_coord(fields[7]),
    ) else {
        log::warn!("Skipping data line with non-numeric coordinate field: {line:?}");
        return false;
    };

    chain.coords_a.push(lend_a);
    chain.coords_a.push(rend_a);
    chain.coords_b.push(lend_b);
    chain.coords_b.push(rend_b);
    true
}

/// Compute the (min, max) span of a coordinate list, or `(None, None)`
/// if empty — matching Perl's `shift`/`pop`-on-empty-array yielding
/// `undef` for both ends.
fn span(coords: &[i64]) -> (Option<i64>, Option<i64>) {
    if coords.is_empty() {
        return (None, None);
    }
    let min = *coords.iter().min().unwrap();
    let max = *coords.iter().max().unwrap();
    (Some(min), Some(max))
}

/// Format an optional coordinate for output: `""` for `None`, matching
/// Perl's stringification of `undef`.
fn fmt_opt(v: Option<i64>) -> String {
    v.map(|x| x.to_string()).unwrap_or_default()
}

/// Core logic, generic over any `BufRead`/`Write`, mirroring the Perl
/// script's stdin -> stdout filter behaviour (with the fixes above
/// applied).
pub fn to_chain_spans<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
) -> Result<(), ToChainSpansError> {
    let mut chains: Vec<Chain> = Vec::new();

    for line in reader.lines() {
        let line = line?;

        if line.trim().is_empty() {
            continue;
        }

        if line.starts_with('#') {
            chains.push(parse_header(&line)?);
        } else {
            let chain = chains
                .last_mut()
                .ok_or_else(|| ToChainSpansError::DataBeforeHeader { line: line.clone() })?;
            parse_data_line(&line, chain);
        }
    }

    for chain in &chains {
        let (lend_a, rend_a) = span(&chain.coords_a);
        let (lend_b, rend_b) = span(&chain.coords_b);

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chain.acc_a,
            fmt_opt(lend_a),
            fmt_opt(rend_a),
            chain.acc_b,
            fmt_opt(lend_b),
            fmt_opt(rend_b),
            chain.orient,
            chain.num_gene_pairs,
        )?;
    }

    Ok(())
}

/// CLI entry point — a drop-in replacement for
/// `dagchainer_to_chain_spans.pl < file.aligncoords` or
/// `dagchainer_to_chain_spans.pl file1 file2 ...`.
pub fn run(args: Args) -> anyhow::Result<()> {
    let stdout = std::io::stdout();
    let out = stdout.lock();

    if args.files.is_empty() {
        let stdin = std::io::stdin();
        to_chain_spans(stdin.lock(), out)?;
    } else {
        // Concatenate all files into one logical stream, matching
        // Perl's `<>` diamond-operator semantics — important since
        // chain-accumulation state persists across the whole input.
        let mut combined: Box<dyn Read> = Box::new(std::io::empty());
        for path in &args.files {
            let file = std::fs::File::open(path)
                .map_err(|e| anyhow::anyhow!("Could not open file {}: {e}", path.display()))?;
            combined = Box::new(combined.chain(file));
        }
        to_chain_spans(std::io::BufReader::new(combined), out)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_single_chain_forward_orientation() {
        let input = "\
## alignment chr1 vs. chr2 (num aligned pairs: 2):
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
chr1\tgeneC\t300\t400\tchr2\tgeneD\t350\t450\t1e-20\t98
";
        let mut out = Vec::new();
        to_chain_spans(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();

        assert_eq!(out, "chr1\t100\t400\tchr2\t150\t450\t+\t2\n");
    }

    #[test]
    fn reverse_orientation_is_detected() {
        let input = "\
## alignment chr1 vs. chr2 (reverse) >Alignment #1  score = 100.0 (num aligned pairs: 1):
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
";
        let mut out = Vec::new();
        to_chain_spans(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();

        assert_eq!(out, "chr1\t100\t200\tchr2\t150\t250\t-\t1\n");
    }

    #[test]
    fn multiple_chains_each_get_their_own_span() {
        let input = "\
## alignment chr1 vs. chr2 (num aligned pairs: 1):
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
## alignment chr1 vs. chr3 (num aligned pairs: 1):
chr1\tgeneC\t500\t600\tchr3\tgeneD\t700\t800\t1e-20\t40
";
        let mut out = Vec::new();
        to_chain_spans(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();

        let expected = "\
chr1\t100\t200\tchr2\t150\t250\t+\t1
chr1\t500\t600\tchr3\t700\t800\t+\t1
";
        assert_eq!(out, expected);
    }

    #[test]
    fn header_with_no_data_lines_produces_empty_span_fields() {
        let input = "## alignment chr1 vs. chr2 (num aligned pairs: 0):\n";
        let mut out = Vec::new();
        to_chain_spans(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();

        assert_eq!(out, "chr1\t\t\tchr2\t\t\t+\t0\n");
    }

    #[test]
    fn blank_separator_line_is_skipped_and_does_not_corrupt_span() {
        // Simulates the output of `min_chain_len_filter.pl`, which
        // inserts a trailing blank line after each surviving chain
        // block. Without the fix, this would inject spurious 0
        // coordinates into the preceding chain's span.
        let input = "\
## alignment chr1 vs. chr2 (num aligned pairs: 1):
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50

## alignment chr1 vs. chr3 (num aligned pairs: 1):
chr1\tgeneC\t500\t600\tchr3\tgeneD\t700\t800\t1e-20\t40
";
        let mut out = Vec::new();
        to_chain_spans(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();

        let expected = "\
chr1\t100\t200\tchr2\t150\t250\t+\t1
chr1\t500\t600\tchr3\t700\t800\t+\t1
";
        assert_eq!(out, expected);
    }

    #[test]
    fn malformed_data_line_is_skipped() {
        let input = "\
## alignment chr1 vs. chr2 (num aligned pairs: 1):
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
not_enough_fields
";
        let mut out = Vec::new();
        to_chain_spans(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();

        assert_eq!(out, "chr1\t100\t200\tchr2\t150\t250\t+\t1\n");
    }

    #[test]
    fn missing_num_aligned_pairs_is_an_error() {
        let input = "## alignment chr1 vs. chr2 (no count here):\n";
        let mut out = Vec::new();
        let result = to_chain_spans(input.as_bytes(), &mut out);
        assert!(matches!(
            result,
            Err(ToChainSpansError::MissingAlignedPairsCount { .. })
        ));
    }

    #[test]
    fn data_line_before_any_header_is_an_error() {
        let input = "chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50\n";
        let mut out = Vec::new();
        let result = to_chain_spans(input.as_bytes(), &mut out);
        assert!(matches!(result, Err(ToChainSpansError::DataBeforeHeader { .. })));
    }

    #[test]
    fn empty_input_produces_empty_output() {
        let mut out = Vec::new();
        to_chain_spans("".as_bytes(), &mut out).unwrap();
        assert!(out.is_empty());
    }
}