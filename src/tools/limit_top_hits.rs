//! Port of `limit_to_top_N_hits.pl` — keep only the top N entries per
//! accession (2nd tab-separated column) from a BLAST-like tabular stream,
//! "top" meaning lowest (most significant) e-value.
//!
//! ## Deviations from the original Perl (deliberate fixes)
//!
//! 1. **Sort order fixed to ascending e-value.** The original Perl sorted
//!    secondarily by *descending* e-value (`$b->[1] <=> $a->[1]`), which
//!    despite the script's name actually kept the *worst* hits per
//!    accession and discarded the best ones. This port sorts ascending,
//!    so the N entries with the lowest (best) e-value per accession are
//!    kept.
//! 2. **`num = 0` is accepted**, not rejected. The original Perl's
//!    `my $num = $ARGV[0] or die $usage;` treated `0` as falsy and
//!    aborted. Here, `num` is a plain `usize` with no special-casing.
//!
//! ## Behavior note: the first entry per accession is always kept
//!
//! Due to the underlying control flow (inherited faithfully from the
//! original Perl), the *first* entry encountered for any accession
//! (after sorting) is always kept unconditionally — the `count > num`
//! check only applies to the *second and later* entries for that
//! accession. Consequently, `num = 0` and `num = 1` behave identically
//! (both keep exactly 1 entry per accession, the best one). `num = 2`
//! keeps up to 2 entries, `num = 3` keeps up to 3, and so on.

use std::cmp::Ordering;
use std::io::{BufRead, Write};

use clap::Args as ClapArgs;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Maximum number of entries to keep per accession (2nd tab column),
    /// "top" meaning lowest (most significant) e-value. Note: the first
    /// entry per accession is always kept regardless of this value (see
    /// module docs); `0` and `1` are therefore equivalent. `0` is valid.
    #[arg(value_name = "NUM")]
    pub num: usize,
}

/// One parsed input row: accession, e-value, and the original line text
/// (kept verbatim for output, exactly like the Perl `$line`).
#[derive(Debug, Clone)]
struct Entry {
    acc: String,
    evalue: f64,
    line: String,
}

/// Parse a tab-separated line, extracting column 2 (`acc`, 0-indexed 1)
/// and column 9 (`evalue`, 0-indexed 8). Lines with fewer than 9 columns
/// or a non-numeric e-value are skipped with a warning.
fn parse_line(line: &str) -> Option<Entry> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        log::warn!("Skipping line with fewer than 9 tab-separated fields: {line:?}");
        return None;
    }
    let acc = fields[1].to_string();
    let evalue: f64 = match fields[8].parse() {
        Ok(v) => v,
        Err(_) => {
            log::warn!("Skipping line with non-numeric e-value field: {line:?}");
            return None;
        }
    };
    Some(Entry { acc, evalue, line: line.to_string() })
}

/// Sort entries ascending by accession, then ascending by e-value (best
/// hit first) — the corrected ordering vs. the original Perl.
fn sort_entries(entries: &mut [Entry]) {
    entries.sort_by(|a, b| {
        a.acc
            .cmp(&b.acc)
            .then_with(|| a.evalue.partial_cmp(&b.evalue).unwrap_or(Ordering::Equal))
    });
}

/// Core logic, generic over any `BufRead`/`Write`, mirroring the Perl
/// script's stdin -> stdout filter behaviour (with the fixes above
/// applied).
pub fn limit_top_n<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
    num: usize,
) -> std::io::Result<()> {
    let mut entries: Vec<Entry> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(entry) = parse_line(&line) {
            entries.push(entry);
        }
    }

    sort_entries(&mut entries);

    let mut curr_acc: Option<String> = None;
    let mut count = 0usize;

    for entry in entries {
        if curr_acc.as_deref() == Some(entry.acc.as_str()) {
            count += 1;
            if count > num {
                continue;
            }
        } else {
            curr_acc = Some(entry.acc.clone());
            count = 1;
        }
        writeln!(writer, "{}", entry.line)?;
    }

    Ok(())
}

/// CLI entry point, reading from stdin and writing to stdout — a drop-in
/// replacement for `limit_to_top_N_hits.pl num < dataFile`.
pub fn run(args: Args) -> anyhow::Result<()> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    limit_top_n(stdin.lock(), stdout.lock(), args.num).map_err(Into::into)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_lowest_evalues_up_to_num_per_accession() {
        let input = "\
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.001
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.05
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.2
x\tgeneB\tx\tx\tx\tx\tx\tx\t1e-10
";

        // num=1: only the single best (lowest-evalue) entry per
        // accession survives, since the first entry is kept
        // unconditionally and the second one already exceeds the
        // count=1 budget.
        let mut out = Vec::new();
        limit_top_n(input.as_bytes(), &mut out, 1).unwrap();
        let out = String::from_utf8(out).unwrap();
        let expected = "\
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.001
x\tgeneB\tx\tx\tx\tx\tx\tx\t1e-10
";
        assert_eq!(out, expected);

        // num=2: the two best entries per accession survive.
        let mut out = Vec::new();
        limit_top_n(input.as_bytes(), &mut out, 2).unwrap();
        let out = String::from_utf8(out).unwrap();
        let expected = "\
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.001
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.05
x\tgeneB\tx\tx\tx\tx\tx\tx\t1e-10
";
        assert_eq!(out, expected);
    }

    #[test]
    fn num_zero_keeps_exactly_one_best_hit_per_accession() {
        let input = "\
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.001
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.05
x\tgeneB\tx\tx\tx\tx\tx\tx\t1e-10
";
        let mut out = Vec::new();
        limit_top_n(input.as_bytes(), &mut out, 0).unwrap();
        let out = String::from_utf8(out).unwrap();

        let expected = "\
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.001
x\tgeneB\tx\tx\tx\tx\tx\tx\t1e-10
";
        assert_eq!(out, expected);
    }

    #[test]
    fn num_zero_and_num_one_are_equivalent() {
        let input = "\
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.001
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.05
x\tgeneB\tx\tx\tx\tx\tx\tx\t1e-10
";
        let mut out0 = Vec::new();
        limit_top_n(input.as_bytes(), &mut out0, 0).unwrap();

        let mut out1 = Vec::new();
        limit_top_n(input.as_bytes(), &mut out1, 1).unwrap();

        assert_eq!(out0, out1);
    }

    #[test]
    fn malformed_lines_are_skipped() {
        let input = "too\tfew\tfields\n";
        let mut out = Vec::new();
        limit_top_n(input.as_bytes(), &mut out, 5).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn non_numeric_evalue_is_skipped() {
        let input = "x\tgeneA\tx\tx\tx\tx\tx\tx\tnot_a_number\n";
        let mut out = Vec::new();
        limit_top_n(input.as_bytes(), &mut out, 5).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn different_accessions_each_get_their_own_budget() {
        let input = "\
x\tgeneA\tx\tx\tx\tx\tx\tx\t0.01
x\tgeneB\tx\tx\tx\tx\tx\tx\t0.02
x\tgeneC\tx\tx\tx\tx\tx\tx\t0.03
";
        let mut out = Vec::new();
        limit_top_n(input.as_bytes(), &mut out, 0).unwrap();
        let out = String::from_utf8(out).unwrap();

        // Each accession appears only once in input, so all survive
        // regardless of num (each gets its unconditional first-entry).
        assert_eq!(out.lines().count(), 3);
    }
}