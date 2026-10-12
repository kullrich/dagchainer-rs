//! Port of `extract_only_chromo_pairs_with_dups.pl` — filter a match file
//! down to only those lines whose chromosome/molecule pair also appears
//! somewhere in a separate "diagonal" file (typically DAGchainer chain
//! output), i.e. only pairs already known to have at least one syntenic
//! chain / duplication block.
//!
//! ## Input formats
//!
//! Both files are tab-separated, with the chromosome/molecule name in
//! column 1 and column 5 (0-indexed 0 and 4) — matching the standard
//! DAGchainer match-file layout:
//! `mol_A  acc_A  end5_A  end3_A  mol_B  acc_B  end5_B  end3_B  [...]`.
//! Only columns 1 and 5 are actually inspected; all other columns
//! (including the exact file "type" — raw matches vs. chain output) are
//! irrelevant to this script's logic.
//!
//! ## Output
//!
//! Every line from `match_file` whose (sorted) chromosome pair appears
//! anywhere in `diag_file`, printed verbatim, in original order,
//! without deduplication.
//!
//! ## Deviation from the original Perl (deliberate improvement)
//!
//! Lines with fewer than 5 tab-separated fields (in either input file)
//! are skipped with a `log::warn!`, rather than silently building a
//! malformed pair key from Perl's `undef`-stringifies-to-empty-string
//! behavior (e.g. a line with only 2 fields would key as `",<field0>"`
//! or similar nonsense in the original).
//!
//! Both file arguments are plain required positional CLI arguments, so
//! there is no Perl `$ARGV[0] or die` truthiness trap to replicate here
//! (unlike scripts with a numeric threshold argument, where an explicit
//! `0` could be misinterpreted as "missing").

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

use clap::Args as ClapArgs;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Match file to filter.
    #[arg(value_name = "MATCH_FILE")]
    pub match_file: PathBuf,

    /// Diagonal/chain file whose chromosome pairs define the keep-set.
    #[arg(value_name = "DIAG_FILE")]
    pub diag_file: PathBuf,
}

/// Build the sorted, comma-joined chromosome-pair key from columns 1 and
/// 5 (0-indexed 0 and 4) of a tab-separated line. Returns `None` (with a
/// warning logged) if the line has fewer than 5 fields.
fn chromo_pair_key(line: &str) -> Option<String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 5 {
        log::warn!("Skipping line with fewer than 5 tab-separated fields: {line:?}");
        return None;
    }
    let (chr_a, chr_b) = if fields[0] <= fields[4] {
        (fields[0], fields[4])
    } else {
        (fields[4], fields[0])
    };
    Some(format!("{chr_a},{chr_b}"))
}

/// Parse the diagonal/chain file, building the set of chromosome pairs
/// it mentions.
pub fn build_chromo_pair_set<R: BufRead>(reader: R) -> std::io::Result<BTreeSet<String>> {
    let mut pairs = BTreeSet::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(key) = chromo_pair_key(&line) {
            pairs.insert(key);
        }
    }
    Ok(pairs)
}

/// Core logic: filter `match_reader`'s lines, writing (verbatim, in
/// order) every line whose chromosome pair is present in `chromo_pairs`.
pub fn filter_matches_by_chromo_pairs<R: BufRead, W: Write>(
    match_reader: R,
    chromo_pairs: &BTreeSet<String>,
    mut writer: W,
) -> std::io::Result<()> {
    for line in match_reader.lines() {
        let line = line?;
        if let Some(key) = chromo_pair_key(&line) {
            if chromo_pairs.contains(&key) {
                writeln!(writer, "{line}")?;
            }
        }
    }
    Ok(())
}

/// CLI entry point — a drop-in replacement for
/// `extract_only_chromo_pairs_with_dups.pl matchFile diagFile`.
pub fn run(args: Args) -> anyhow::Result<()> {
    let diag_file = File::open(&args.diag_file)
        .map_err(|e| anyhow::anyhow!("Cannot open {}: {e}", args.diag_file.display()))?;
    let chromo_pairs = build_chromo_pair_set(BufReader::new(diag_file))?;

    let match_file = File::open(&args.match_file)
        .map_err(|e| anyhow::anyhow!("Cannot open {}: {e}", args.match_file.display()))?;

    let stdout = std::io::stdout();
    filter_matches_by_chromo_pairs(BufReader::new(match_file), &chromo_pairs, stdout.lock())?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_only_lines_with_chromo_pair_present_in_diag_file() {
        let diag = "chr1\tx\tx\tx\tchr2\tx\tx\tx\n"; // registers pair "chr1,chr2"
        let matches = "\
chr1\tgeneA\t1\t2\tchr2\tgeneB\t3\t4\n
chr1\tgeneC\t1\t2\tchr3\tgeneD\t3\t4\n
chr2\tgeneE\t1\t2\tchr1\tgeneF\t3\t4\n
";
        let pairs = build_chromo_pair_set(diag.as_bytes()).unwrap();
        let mut out = Vec::new();
        filter_matches_by_chromo_pairs(matches.as_bytes(), &pairs, &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();

        // Line 1 (chr1,chr2) and line 3 (chr2,chr1 -> same sorted pair)
        // should survive; line 2 (chr1,chr3) should not.
        assert!(out.contains("geneA"));
        assert!(out.contains("geneE"));
        assert!(!out.contains("geneC"));
    }

    #[test]
    fn pair_order_is_irrelevant_due_to_sorting() {
        let diag = "chrB\tx\tx\tx\tchrA\tx\tx\tx\n"; // sorted key: "chrA,chrB"
        let matches = "chrA\tg1\t1\t2\tchrB\tg2\t3\t4\n"; // same sorted key
        let pairs = build_chromo_pair_set(diag.as_bytes()).unwrap();
        let mut out = Vec::new();
        filter_matches_by_chromo_pairs(matches.as_bytes(), &pairs, &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.contains("g1"));
    }

    #[test]
    fn self_chromosome_pairs_work() {
        let diag = "chr1\tx\tx\tx\tchr1\tx\tx\tx\n"; // key: "chr1,chr1"
        let matches = "chr1\tgA\t1\t2\tchr1\tgB\t3\t4\n";
        let pairs = build_chromo_pair_set(diag.as_bytes()).unwrap();
        let mut out = Vec::new();
        filter_matches_by_chromo_pairs(matches.as_bytes(), &pairs, &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.contains("gA"));
    }

    #[test]
    fn lines_with_too_few_fields_are_skipped_in_both_files() {
        let diag = "only\tthree\tfields\n";
        let matches = "chr1\tg1\t1\t2\tchr2\tg2\t3\t4\n";
        let pairs = build_chromo_pair_set(diag.as_bytes()).unwrap();
        assert!(pairs.is_empty());

        let mut out = Vec::new();
        filter_matches_by_chromo_pairs(matches.as_bytes(), &pairs, &mut out).unwrap();
        assert!(out.is_empty()); // pairs set is empty, so nothing matches
    }

    #[test]
    fn duplicate_matching_lines_are_all_kept_unfiltered() {
        let diag = "chr1\tx\tx\tx\tchr2\tx\tx\tx\n";
        let matches = "\
chr1\tgA\t1\t2\tchr2\tgB\t3\t4\n
chr1\tgA\t1\t2\tchr2\tgB\t3\t4\n
";
        let pairs = build_chromo_pair_set(diag.as_bytes()).unwrap();
        let mut out = Vec::new();
        filter_matches_by_chromo_pairs(matches.as_bytes(), &pairs, &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert_eq!(out.lines().count(), 2); // no dedup, both kept
    }

    #[test]
    fn empty_diag_file_yields_no_matches() {
        let pairs = build_chromo_pair_set("".as_bytes()).unwrap();
        let matches = "chr1\tg1\t1\t2\tchr2\tg2\t3\t4\n";
        let mut out = Vec::new();
        filter_matches_by_chromo_pairs(matches.as_bytes(), &pairs, &mut out).unwrap();
        assert!(out.is_empty());
    }
}