//! Port of `filter_repetitive_matches.pl` — collapse clusters of closely
//! spaced matches sharing an accession into a single representative
//! (lowest e-value) match, applied independently in both pairing
//! directions (A-grouped-by-B, then B-grouped-by-A).
//!
//! ## Deviation from the original Perl (deliberate fix)
//!
//! `window_link_size = 0` is accepted, not rejected. The original Perl's
//! `my $windowLinkSize = $ARGV[0] or die $usage;` treated `0` as falsy
//! and aborted. Here it's a plain `i64`; `0` simply means "only link
//! matches whose coordinate midpoints are exactly identical," which is a
//! well-defined (if narrow) window.

use std::collections::BTreeMap;
use std::io::{BufRead, Write};

use clap::Args as ClapArgs;
use thiserror::Error;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Window size (bp). Matches sharing an accession, whose paired
    /// coordinate midpoints fall within this distance of the previous
    /// one (sorted by midpoint), are linked into one cluster; only the
    /// lowest-e-value match per cluster is kept. `0` is valid.
    #[arg(value_name = "WINDOW_LINK_SIZE")]
    pub window_link_size: i64,
}

#[derive(Debug, Error)]
pub enum FilterRepetitiveError {
    #[error("WINDOW_LINK_SIZE must be >= 0, got {0}")]
    NegativeNotAllowed(i64),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

#[derive(Debug, Clone)]
struct Feature {
    acc: String,
    mid: i64,
}

#[derive(Debug, Clone)]
struct MatchRecord {
    chr_a: String,
    chr_b: String,
    feat_a: Feature,
    feat_b: Feature,
    evalue: f64,
    input: String,
}

/// Rounds like Perl's `int(x + 0.5)` for non-negative x (genomic
/// coordinates are always >= 0 in practice).
fn round_half_up(x: f64) -> i64 {
    (x + 0.5).floor() as i64
}

/// Parse one input line into a `MatchRecord`, applying the same
/// chromosome-order normalization (and corresponding feature swap) as
/// the original Perl.
fn parse_line(line: &str) -> Option<MatchRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        log::warn!("Skipping malformed line (fewer than 9 fields): {line:?}");
        return None;
    }

    let chr1 = fields[0];
    let acc1 = fields[1].to_string();
    let pos1_end5: f64 = fields[2].parse().ok()?;
    let pos1_end3: f64 = fields[3].parse().ok()?;

    let chr2 = fields[4];
    let acc2 = fields[5].to_string();
    let pos2_end5: f64 = fields[6].parse().ok()?;
    let pos2_end3: f64 = fields[7].parse().ok()?;

    let evalue: f64 = fields[8].parse().ok()?;

    let mut feat_a = Feature { acc: acc1, mid: round_half_up((pos1_end5 + pos1_end3) / 2.0) };
    let mut feat_b = Feature { acc: acc2, mid: round_half_up((pos2_end5 + pos2_end3) / 2.0) };
    let (mut chr_a, mut chr_b) = (chr1.to_string(), chr2.to_string());

    if chr1 > chr2 {
        std::mem::swap(&mut feat_a, &mut feat_b);
        std::mem::swap(&mut chr_a, &mut chr_b);
    }

    Some(MatchRecord {
        chr_a,
        chr_b,
        feat_a,
        feat_b,
        evalue,
        input: line.to_string(),
    })
}

/// Which side of the match pair to treat as "accession key" vs
/// "coordinate key" during one filtering pass. Mirrors the Perl
/// function's `($accKey, $coordKey)` string dispatch, type-safe.
#[derive(Clone, Copy)]
enum Side {
    A,
    B,
}

fn acc_of(m: &MatchRecord, side: Side) -> &str {
    match side {
        Side::A => &m.feat_a.acc,
        Side::B => &m.feat_b.acc,
    }
}

fn mid_of(m: &MatchRecord, side: Side) -> i64 {
    match side {
        Side::A => m.feat_a.mid,
        Side::B => m.feat_b.mid,
    }
}

/// Equivalent to Perl's `filter_matches($accKey, $coordKey, $list)`:
/// sort by (acc_side.acc, coord_side.mid), bin consecutive entries that
/// share acc_side.acc and are within `window_link_size` of the previous
/// entry's coord_side.mid, then keep only the lowest-e-value match per
/// bin.
fn filter_matches(
    list: Vec<MatchRecord>,
    acc_side: Side,
    coord_side: Side,
    window_link_size: i64,
) -> Vec<MatchRecord> {
    if list.is_empty() {
        return list;
    }

    let mut list = list;
    list.sort_by(|a, b| {
        acc_of(a, acc_side)
            .cmp(acc_of(b, acc_side))
            .then_with(|| mid_of(a, coord_side).cmp(&mid_of(b, coord_side)))
    });

    let mut bins: Vec<Vec<MatchRecord>> = Vec::new();
    for m in list {
        let needs_new_bin = match bins.last() {
            None => true,
            Some(prev_bin) => {
                let prev = prev_bin.last().unwrap();
                acc_of(prev, acc_side) != acc_of(&m, acc_side)
                    || (mid_of(&m, coord_side) - mid_of(prev, coord_side) > window_link_size)
            }
        };

        if needs_new_bin {
            bins.push(vec![m]);
        } else {
            bins.last_mut().unwrap().push(m);
        }
    }

    bins.into_iter()
        .map(|mut bin| {
            if bin.len() == 1 {
                bin.pop().unwrap()
            } else {
                bin.sort_by(|a, b| {
                    a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal)
                });
                bin.into_iter().next().unwrap()
            }
        })
        .collect()
}

/// Core logic: read matches, group by chromosome pair, filter twice
/// (once keyed on feature A / linked via feature B's midpoint, once
/// keyed on feature B / linked via feature A's midpoint), write
/// survivors — output order matches the original Perl (chromosome pairs
/// visited in sorted key order).
pub fn filter_repetitive<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
    window_link_size: i64,
) -> Result<(), FilterRepetitiveError> {
    if window_link_size < 0 {
        return Err(FilterRepetitiveError::NegativeNotAllowed(window_link_size));
    }

    let mut grouped: BTreeMap<String, Vec<MatchRecord>> = BTreeMap::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(m) = parse_line(&line) {
            let key = format!("{},{}", m.chr_a, m.chr_b);
            grouped.entry(key).or_default().push(m);
        }
    }

    for (_chromo_pair, match_list) in grouped {
        let filtered = filter_matches(match_list, Side::A, Side::B, window_link_size);
        let filtered = filter_matches(filtered, Side::B, Side::A, window_link_size);

        for m in filtered {
            write!(writer, "{}", m.input)?;
            if !m.input.ends_with('\n') {
                writeln!(writer)?;
            }
        }
    }

    Ok(())
}

/// CLI entry point, reading from stdin and writing to stdout — a drop-in
/// replacement for `filter_repetitive_matches.pl windowLinkSize < inputFile`.
pub fn run(args: Args) -> anyhow::Result<()> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    filter_repetitive(stdin.lock(), stdout.lock(), args.window_link_size).map_err(Into::into)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(chr_a: &str, acc_a: &str, mid_a: i64, chr_b: &str, acc_b: &str, mid_b: i64, evalue: f64) -> String {
        // end5/end3 both set to `mid` for simplicity -> midpoint == mid.
        format!(
            "{chr_a}\t{acc_a}\t{mid_a}\t{mid_a}\t{chr_b}\t{acc_b}\t{mid_b}\t{mid_b}\t{evalue}\n"
        )
    }

    #[test]
    fn collapses_cluster_to_lowest_evalue() {
        // Three matches on the same chromosome pair, all involving acc_a,
        // with featB midpoints close together (within window=50) -> should
        // collapse to the single lowest-evalue entry.
        let input = format!(
            "{}{}{}",
            line("chr1", "accA", 100, "chr2", "accB1", 1000, 0.01),
            line("chr1", "accA", 100, "chr2", "accB2", 1020, 0.5),
            line("chr1", "accA", 100, "chr2", "accB3", 1040, 0.2),
        );

        let mut out = Vec::new();
        filter_repetitive(input.as_bytes(), &mut out, 50).unwrap();
        let out = String::from_utf8(out).unwrap();

        // Only the lowest-evalue (0.01) entry should survive.
        assert_eq!(out.lines().count(), 1);
        assert!(out.contains("0.01"));
    }

    #[test]
    fn distinct_matches_outside_window_are_kept_separately() {
        let input = format!(
            "{}{}",
            line("chr1", "accA", 100, "chr2", "accB1", 1000, 0.01),
            line("chr1", "accA", 100, "chr2", "accB2", 5000, 0.02), // far outside window
        );

        let mut out = Vec::new();
        filter_repetitive(input.as_bytes(), &mut out, 50).unwrap();
        let out = String::from_utf8(out).unwrap();

        assert_eq!(out.lines().count(), 2);
    }

    #[test]
    fn window_zero_only_links_identical_midpoints() {
        let input = format!(
            "{}{}",
            line("chr1", "accA", 100, "chr2", "accB1", 1000, 0.01),
            line("chr1", "accA", 100, "chr2", "accB2", 1001, 0.02), // 1bp apart
        );

        let mut out = Vec::new();
        filter_repetitive(input.as_bytes(), &mut out, 0).unwrap();
        let out = String::from_utf8(out).unwrap();

        // With window=0, these are NOT linked (1bp apart > 0), so both survive.
        assert_eq!(out.lines().count(), 2);
    }

    #[test]
    fn negative_window_is_rejected() {
        let input = line("chr1", "accA", 100, "chr2", "accB", 1000, 0.01);
        let mut out = Vec::new();
        let result = filter_repetitive(input.as_bytes(), &mut out, -1);
        assert!(matches!(result, Err(FilterRepetitiveError::NegativeNotAllowed(-1))));
    }

    #[test]
    fn malformed_lines_are_skipped() {
        let input = "too\tfew\tfields\n";
        let mut out = Vec::new();
        filter_repetitive(input.as_bytes(), &mut out, 10).unwrap();
        assert!(out.is_empty());
    }
}