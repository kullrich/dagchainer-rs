//! Port of `report_dist_between_chain_segments.pl` — for each DAGchainer
//! chain block, report the absolute midpoint distance between
//! consecutive gene-pair entries, on both sides (A and B).
//!
//! ## Input format
//!
//! A sequence of blocks, each starting with a `#`-prefixed header line
//! (content ignored beyond triggering a state reset), followed by
//! tab-separated data lines:
//! `contigA accA lendA rendA contigB accB lendB rendB [evalue] [chain_score]`.
//!
//! ## Output
//!
//! For every data line **after the first** within a chain block:
//! `contigA  delta_A  contigB  delta_B`, where `delta_X` is the absolute
//! difference between the current and previous line's midpoint on side
//! X. The first data line following any `#` header (or the start of
//! input) produces no output, since there is nothing yet to compare
//! against.
//!
//! ## Preserved behavior: numeric formatting matches Perl's default
//!
//! Midpoints and deltas are plain floating-point values (not rounded to
//! integers), and are formatted using Perl's default `%.15g`-equivalent
//! stringification via [`crate::tools::chain_length_ratio::perl_format_number`]
//! (see that module's docs for details/rationale).
//!
//! ## Deviation from the original Perl (deliberate improvement)
//!
//! Malformed data lines (fewer than 8 tab-separated fields, or
//! non-numeric coordinate fields) are skipped with a `log::warn!`,
//! rather than silently computing a corrupted midpoint via Perl's
//! `undef`-in-numeric-context coercion. Skipping a line also means it
//! does **not** update the "previous" comparison point — equivalent to
//! treating the malformed line as if it were never present, rather than
//! injecting a bogus `(0, 0)` midpoint into the comparison chain.

use std::io::{BufRead, Read, Write};
use std::path::PathBuf;

use clap::Args as ClapArgs;

use crate::tools::chain_length_ratio::perl_format_number;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Input file(s). If omitted, reads from stdin. Multiple files are
    /// treated as one continuous logical stream (mirroring Perl's `<>`
    /// diamond operator), since the "previous entry" comparison state
    /// persists across the whole input.
    #[arg(value_name = "FILE")]
    pub files: Vec<PathBuf>,
}

/// A single data line's midpoints on each side.
#[derive(Debug, Clone, Copy)]
struct Midpoints {
    mid_a: f64,
    mid_b: f64,
}

/// One parsed data line's fields actually used by this script.
struct Row<'a> {
    contig_a: &'a str,
    contig_b: &'a str,
    mid: Midpoints,
}

fn parse_data_line(line: &str) -> Option<Row<'_>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        log::warn!("Skipping line with fewer than 8 tab-separated fields: {line:?}");
        return None;
    }

    let parse = |s: &str| s.parse::<f64>().ok();
    let (Some(lend_a), Some(rend_a), Some(lend_b), Some(rend_b)) = (
        parse(fields[2]),
        parse(fields[3]),
        parse(fields[6]),
        parse(fields[7]),
    ) else {
        log::warn!("Skipping line with non-numeric coordinate field: {line:?}");
        return None;
    };

    Some(Row {
        contig_a: fields[0],
        contig_b: fields[4],
        mid: Midpoints {
            mid_a: (lend_a + rend_a) / 2.0,
            mid_b: (lend_b + rend_b) / 2.0,
        },
    })
}

/// Core logic, generic over any `BufRead`/`Write`, mirroring the Perl
/// script's stdin -> stdout filter behaviour.
pub fn report_segment_dist<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
) -> std::io::Result<()> {
    let mut prev: Option<Midpoints> = None;

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('#') {
            prev = None;
            continue;
        }

        let Some(row) = parse_data_line(&line) else {
            continue;
        };

        if let Some(p) = prev {
            let delta_a = (row.mid.mid_a - p.mid_a).abs();
            let delta_b = (row.mid.mid_b - p.mid_b).abs();
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                row.contig_a,
                perl_format_number(delta_a),
                row.contig_b,
                perl_format_number(delta_b),
            )?;
        }

        prev = Some(row.mid);
    }

    Ok(())
}

/// CLI entry point — a drop-in replacement for
/// `report_dist_between_chain_segments.pl < file.aligncoords` or
/// `report_dist_between_chain_segments.pl file1 file2 ...`.
pub fn run(args: Args) -> anyhow::Result<()> {
    let stdout = std::io::stdout();
    let out = stdout.lock();

    if args.files.is_empty() {
        let stdin = std::io::stdin();
        report_segment_dist(stdin.lock(), out)?;
    } else {
        // Concatenate all files into one logical stream, matching
        // Perl's `<>` diamond-operator semantics — important since
        // "previous entry" state persists across the whole input.
        let mut combined: Box<dyn Read> = Box::new(std::io::empty());
        for path in &args.files {
            let file = std::fs::File::open(path)
                .map_err(|e| anyhow::anyhow!("Could not open file {}: {e}", path.display()))?;
            combined = Box::new(combined.chain(file));
        }
        report_segment_dist(std::io::BufReader::new(combined), out)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_line_in_block_produces_no_output() {
        let input = "chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50\n";
        let mut out = Vec::new();
        report_segment_dist(input.as_bytes(), &mut out).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn computes_delta_between_consecutive_lines() {
        let input = "\
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
chr1\tgeneC\t300\t400\tchr2\tgeneD\t350\t450\t1e-20\t98
";
        // midA_1 = 150, midA_2 = 350 -> delta_A = 200
        // midB_1 = 200, midB_2 = 400 -> delta_B = 200
        let mut out = Vec::new();
        report_segment_dist(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert_eq!(out, "chr1\t200\tchr2\t200\n");
    }

    #[test]
    fn header_line_resets_comparison_state() {
        let input = "\
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
chr1\tgeneC\t300\t400\tchr2\tgeneD\t350\t450\t1e-20\t98
## new chain
chr1\tgeneE\t500\t600\tchr2\tgeneF\t550\t650\t1e-10\t40
chr1\tgeneG\t700\t800\tchr2\tgeneH\t750\t850\t1e-5\t30
";
        let mut out = Vec::new();
        report_segment_dist(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();

        // Two separate deltas, one per block; no cross-block delta
        // computed despite the 5th line following immediately after the
        // header.
        let expected = "\
chr1\t200\tchr2\t200
chr1\t200\tchr2\t200
";
        assert_eq!(out, expected);
    }

    #[test]
    fn uses_current_lines_contig_names_not_previous() {
        let input = "\
chrX\tgeneA\t100\t200\tchrY\tgeneB\t150\t250\t1e-30\t50
chrZ\tgeneC\t300\t400\tchrW\tgeneD\t350\t450\t1e-20\t98
";
        let mut out = Vec::new();
        report_segment_dist(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.starts_with("chrZ\t")); // current line's contig_A
        assert!(out.contains("chrW\t"));     // current line's contig_B
    }

    #[test]
    fn fractional_midpoints_format_like_perl() {
        // lend=100, rend=201 -> midpoint = 150.5
        // next: lend=100, rend=200 -> midpoint = 150.0, delta = 0.5
        let input = "\
chr1\tgeneA\t100\t201\tchr2\tgeneB\t100\t201\t1e-30\t50
chr1\tgeneC\t100\t200\tchr2\tgeneD\t100\t200\t1e-20\t98
";
        let mut out = Vec::new();
        report_segment_dist(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert_eq!(out, "chr1\t0.5\tchr2\t0.5\n");
    }

    #[test]
    fn malformed_line_is_skipped_and_does_not_update_state() {
        let input = "\
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
too\tfew\tfields
chr1\tgeneC\t300\t400\tchr2\tgeneD\t350\t450\t1e-20\t98
";
        let mut out = Vec::new();
        report_segment_dist(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        // The malformed middle line is skipped entirely; the delta is
        // still computed correctly between lines 1 and 3.
        assert_eq!(out, "chr1\t200\tchr2\t200\n");
    }

    #[test]
    fn empty_input_produces_empty_output() {
        let mut out = Vec::new();
        report_segment_dist("".as_bytes(), &mut out).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn multiple_headers_in_a_row_keep_state_reset() {
        let input = "\
## header 1
## header 2
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30\t50
";
        let mut out = Vec::new();
        report_segment_dist(input.as_bytes(), &mut out).unwrap();
        assert!(out.is_empty()); // only one data line total, no pair to compare
    }
}