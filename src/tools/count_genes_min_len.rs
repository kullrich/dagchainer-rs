//! Port of `CountGenesMinLen.pl` — collect the unique set of gene models
//! (molecule + accession) participating in DAGchainer `.aligncoords`
//! chains that meet a minimum aligned-pair-count threshold.
//!
//! ## Input format
//!
//! Same as [`crate::tools::min_chain_len_filter`]: a sequence of blocks,
//! each starting with a `#`-prefixed header containing
//! `num aligned pairs: <digits>`, followed by tab-separated data lines
//! whose first, second, fifth, and sixth columns are
//! `chr_A, gene_A, _, _, chr_B, gene_B, ...`.
//!
//! ## Output
//!
//! One `chr; gene` line per unique gene model found in any surviving
//! chain (both sides of every pair, across all kept chains).
//!
//! ## Output order differs from the original Perl (deliberate, see below)
//!
//! The original Perl collects gene models into a hash and iterates
//! `keys %gene_models` directly. Since Perl 5.18, hash key iteration
//! order is randomized per-process specifically to prevent code from
//! depending on it — so the original script's output order is **not
//! even deterministic across repeated runs on identical input**. This
//! port collects into a `HashSet` (equivalent dedup semantics) but
//! sorts the final output lexicographically before printing, making
//! results reproducible. The *set* of lines produced is equivalent to
//! the original; only the order is (deliberately, and necessarily)
//! different.
//!
//! ## Deviation from the original Perl (deliberate fixes)
//!
//! - `min_chain_len = 0` is accepted literally (same fix as in
//!   [`crate::tools::min_chain_len_filter`]).
//! - Data lines with fewer than 6 tab-separated fields are skipped with
//!   a warning, rather than silently building a malformed gene-model key
//!   from Perl's `undef`-stringifies-to-empty-string behavior.

use std::collections::BTreeSet;
use std::io::{BufRead, Write};

use clap::Args as ClapArgs;
use thiserror::Error;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Minimum number of aligned pairs a chain must have for its gene
    /// models to be counted. `0` is valid (counts every chain).
    #[arg(value_name = "MIN_CHAIN_LENGTH")]
    pub min_chain_len: usize,
}

#[derive(Debug, Error)]
pub enum CountGenesMinLenError {
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

/// Parse one data line's gene-model keys (`chr; gene` for each side),
/// matching Perl's `($x[0], $x[1], $x[4], $x[5])` extraction.
fn parse_gene_models(line: &str) -> Option<(String, String)> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 6 {
        log::warn!("Skipping data line with fewer than 6 tab-separated fields: {line:?}");
        return None;
    }
    let key_a = format!("{}; {}", fields[0], fields[1]);
    let key_b = format!("{}; {}", fields[4], fields[5]);
    Some((key_a, key_b))
}

/// Process one chain block: if its header's aligned-pair count meets
/// `min_chain_len`, insert both gene-model keys from every data line
/// into `gene_models`.
fn process_chain_block(
    lines: &[String],
    min_chain_len: usize,
    gene_models: &mut BTreeSet<String>,
) -> Result<(), CountGenesMinLenError> {
    let header = lines
        .first()
        .expect("chain block must have at least a header line");

    let num_pairs = parse_num_pairs(header).ok_or_else(|| {
        CountGenesMinLenError::MissingAlignedPairsCount {
            header: header.clone(),
            block: lines.join("\n"),
        }
    })?;

    if num_pairs < min_chain_len {
        return Ok(());
    }

    for line in &lines[1..] {
        if let Some((key_a, key_b)) = parse_gene_models(line) {
            gene_models.insert(key_a);
            gene_models.insert(key_b);
        }
    }

    Ok(())
}

/// Core logic, generic over any `BufRead`/`Write`, mirroring the Perl
/// script's stdin -> stdout filter behaviour. Output lines are sorted
/// for determinism (see module docs).
pub fn count_genes_min_len<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
    min_chain_len: usize,
) -> Result<(), CountGenesMinLenError> {
    let mut current_block: Vec<String> = Vec::new();
    let mut gene_models: BTreeSet<String> = BTreeSet::new();

    for line in reader.lines() {
        let line = line?;
        if !has_word_char(&line) {
            continue;
        }

        if line.starts_with('#') {
            if !current_block.is_empty() {
                process_chain_block(&current_block, min_chain_len, &mut gene_models)?;
            }
            current_block = vec![line];
        } else {
            current_block.push(line);
        }
    }

    if !current_block.is_empty() {
        process_chain_block(&current_block, min_chain_len, &mut gene_models)?;
    }

    for gene_model in &gene_models {
        writeln!(writer, "{gene_model}")?;
    }

    Ok(())
}

/// CLI entry point, reading from stdin and writing to stdout — a drop-in
/// replacement for `CountGenesMinLen.pl min_chain_length < file.aligncoords`
/// (modulo output ordering — see module docs).
pub fn run(args: Args) -> anyhow::Result<()> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    count_genes_min_len(stdin.lock(), stdout.lock(), args.min_chain_len).map_err(Into::into)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collects_unique_gene_models_from_qualifying_chains() {
        let input = "\
## alignment chr1 vs. chr2 (num aligned pairs: 2):
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
chr1\tgeneC\t300\t400\tchr2\tgeneD\t350\t450\t1e-20\t48
## alignment chr1 vs. chr3 (num aligned pairs: 1):
chr1\tgeneE\t500\t600\tchr3\tgeneF\t550\t650\t1e-10\t40
";
        let mut out = Vec::new();
        count_genes_min_len(input.as_bytes(), &mut out, 2).unwrap();
        let out = String::from_utf8(out).unwrap();

        // Only the first block (2 pairs) qualifies; its 4 gene models
        // (sorted) should appear, none from the second (1-pair) block.
        let expected = "\
chr1; geneA
chr1; geneC
chr2; geneB
chr2; geneD
";
        assert_eq!(out, expected);
    }

    #[test]
    fn duplicate_gene_models_across_chains_are_deduplicated() {
        let input = "\
## alignment chr1 vs. chr2 (num aligned pairs: 1):
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
## alignment chr1 vs. chr3 (num aligned pairs: 1):
chr1\tgeneA\t100\t200\tchr3\tgeneC\t150\t250\t1e-20\t48
";
        let mut out = Vec::new();
        count_genes_min_len(input.as_bytes(), &mut out, 1).unwrap();
        let out = String::from_utf8(out).unwrap();

        // "chr1; geneA" appears in both chains but must only be listed once.
        let expected = "\
chr1; geneA
chr2; geneB
chr3; geneC
";
        assert_eq!(out, expected);
    }

    #[test]
    fn min_chain_len_zero_includes_every_chain() {
        let input = "\
## alignment a vs. b (num aligned pairs: 0):
## alignment c vs. d (num aligned pairs: 1):
chr1\tgeneA\t1\t2\tchr2\tgeneB\t1\t2\t1e-5\t10
";
        let mut out = Vec::new();
        count_genes_min_len(input.as_bytes(), &mut out, 0).unwrap();
        let out = String::from_utf8(out).unwrap();

        let expected = "\
chr1; geneA
chr2; geneB
";
        assert_eq!(out, expected);
    }

    #[test]
    fn blank_and_junk_lines_are_stripped() {
        let input = "\
## alignment a vs. b (num aligned pairs: 1):

chr1\tgeneA\t1\t2\tchr2\tgeneB\t1\t2\t1e-5\t10
----
";
        let mut out = Vec::new();
        count_genes_min_len(input.as_bytes(), &mut out, 1).unwrap();
        let out = String::from_utf8(out).unwrap();

        let expected = "\
chr1; geneA
chr2; geneB
";
        assert_eq!(out, expected);
    }

    #[test]
    fn malformed_header_returns_error() {
        let input = "## this header has no pair count\nchr1\tgeneA\t1\t2\tchr2\tgeneB\t1\t2\t1e-5\t10\n";
        let mut out = Vec::new();
        let result = count_genes_min_len(input.as_bytes(), &mut out, 1);
        assert!(matches!(
            result,
            Err(CountGenesMinLenError::MissingAlignedPairsCount { .. })
        ));
    }

    #[test]
    fn data_line_with_too_few_fields_is_skipped() {
        let input = "\
## alignment a vs. b (num aligned pairs: 1):
chr1\tgeneA\tonly_three_fields
";
        let mut out = Vec::new();
        count_genes_min_len(input.as_bytes(), &mut out, 1).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn empty_input_produces_empty_output() {
        let mut out = Vec::new();
        count_genes_min_len("".as_bytes(), &mut out, 1).unwrap();
        assert!(out.is_empty());
    }
}