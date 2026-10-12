//! Port of `separateDataSets.pl` — bucket lines of a match file into
//! separate per-molecule-pair output files, named
//! `<input_filename>.<mol1>vs<mol2>` (molecules sorted lexically).
//!
//! ## Input format
//!
//! Tab-separated, with molecule/chromosome names in column 1 and column
//! 5 (0-indexed 0 and 4), matching the standard DAGchainer match-file
//! layout. Lines with fewer than 8 fields (after trimming trailing empty
//! fields, matching Perl's `split` behavior) are silently skipped.
//!
//! ## Output
//!
//! One file per distinct (sorted) molecule pair, named
//! `<input_filename>.<mol1>vs<mol2>`, containing every qualifying line
//! involving that pair, in original order.
//!
//! ## Preserved behavior: truncate-on-first-write
//!
//! The first line seen for a given molecule-pair key causes any
//! pre-existing file of that name to be deleted before any writing
//! occurs (matching the original Perl's `unless ($files{$filename}) { if
//! (-s $filename) { unlink $filename; } ... }` guard). All subsequent
//! lines for that key — including the triggering first line — are then
//! appended.
//!
//! ## Engineering improvement (not a behavior change)
//!
//! The original Perl reopens each output file in append mode on every
//! single matching line. This port instead keeps one cached, open file
//! handle per molecule-pair key, producing identical on-disk output far
//! more efficiently for large inputs.
//!
//! ## Documented assumption: line terminators
//!
//! The original Perl never calls `chomp`, so each line (including its
//! original terminator) is written out verbatim — meaning a final input
//! line lacking a trailing newline would be written without one. This
//! port reads lines via [`BufRead::lines`] (which strips terminators)
//! and always writes lines back with `\n` appended. For any input file
//! ending in a newline (the overwhelming norm), output is byte-identical
//! to the original; only a missing-final-newline edge case would differ
//! (gaining a trailing `\n` here that the original would not have had).

use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use clap::Args as ClapArgs;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Match file to split into per-molecule-pair output files.
    #[arg(value_name = "MATCH_FILE")]
    pub match_file: PathBuf,
}

/// Split a line on tabs the way Perl's `split(/\t/, $line)` does by
/// default: on an unlimited split, any **trailing** empty fields
/// (produced by trailing delimiters) are removed; internal empty fields
/// are preserved.
fn perl_split_tab(line: &str) -> Vec<&str> {
    let mut fields: Vec<&str> = line.split('\t').collect();
    while fields.last() == Some(&"") {
        fields.pop();
    }
    fields
}

/// Compute the sorted molecule-pair key (`"<mol1>vs<mol2>"`) for a line,
/// or `None` if the line has fewer than 8 fields.
fn molecule_pair_key(line: &str) -> Option<String> {
    let fields = perl_split_tab(line);
    if fields.len() < 8 {
        return None;
    }
    let (mol1, mol2) = (fields[0], fields[4]);
    let (a, b) = if mol1 <= mol2 { (mol1, mol2) } else { (mol2, mol1) };
    Some(format!("{a}vs{b}"))
}

/// Core logic: read `reader` line by line, bucketing qualifying lines
/// into per-molecule-pair output files under `output_dir_prefix`
/// (typically the original input file's path, so output files are named
/// `<match_file>.<mol1>vs<mol2>`).
///
/// `open_output` is injected so tests can bucket into in-memory buffers
/// instead of real files; [`run`] supplies a real-filesystem
/// implementation.
pub fn separate_datasets<R: BufRead>(
    reader: R,
    mut open_output: impl FnMut(&str) -> std::io::Result<Box<dyn Write>>,
) -> std::io::Result<()> {
    let mut seen_keys: HashSet<String> = HashSet::new();
    let mut handles: HashMap<String, Box<dyn Write>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let Some(key) = molecule_pair_key(&line) else {
            continue;
        };

        if !seen_keys.contains(&key) {
            seen_keys.insert(key.clone());
        }

        let handle = match handles.get_mut(&key) {
            Some(h) => h,
            None => {
                let h = open_output(&key)?;
                handles.insert(key.clone(), h);
                handles.get_mut(&key).unwrap()
            }
        };

        writeln!(handle, "{line}")?;
    }

    Ok(())
}

/// CLI entry point — a drop-in replacement for
/// `separateDataSets.pl matchFile`.
pub fn run(args: Args) -> anyhow::Result<()> {
    let input_path = &args.match_file;
    let file = File::open(input_path)
        .map_err(|e| anyhow::anyhow!("Cannot open {}: {e}", input_path.display()))?;

    // Truncate-on-first-write bookkeeping: tracks which output filenames
    // have already been handled (and thus deleted-if-preexisting) during
    // this run, matching the original Perl's `%files` guard.
    let mut truncated: HashSet<String> = HashSet::new();
    let input_path_owned = input_path.clone();

    separate_datasets(BufReader::new(file), move |key| {
        let out_path = output_path_for(&input_path_owned, key);

        if !truncated.contains(key) {
            if out_path.exists() {
                std::fs::remove_file(&out_path)?;
            }
            truncated.insert(key.to_string());
        }

        let f = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&out_path)?;
        Ok(Box::new(f) as Box<dyn Write>)
    })?;

    Ok(())
}

/// Build the output file path for a given molecule-pair key, matching
/// the original Perl's `"$filename.$molPairKey"` naming.
fn output_path_for(input_path: &Path, key: &str) -> PathBuf {
    let mut os_string = input_path.as_os_str().to_os_string();
    os_string.push(".");
    os_string.push(key);
    PathBuf::from(os_string)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use std::collections::HashMap as StdHashMap;
    use std::rc::Rc;

    /// Shared in-memory store of named output buffers, keyed by
    /// molecule-pair bucket name.
    type SharedStore = Rc<RefCell<StdHashMap<String, Vec<u8>>>>;

    /// Test helper: an in-memory "filesystem" of named buffers, so we can
    /// exercise `separate_datasets` without touching the real filesystem.
    fn in_memory_sink() -> (SharedStore, impl FnMut(&str) -> std::io::Result<Box<dyn Write>>) {
        let store: Rc<RefCell<StdHashMap<String, Vec<u8>>>> = Rc::new(RefCell::new(StdHashMap::new()));
        let store_clone = Rc::clone(&store);

        struct BufWriter {
            key: String,
            store: Rc<RefCell<StdHashMap<String, Vec<u8>>>>,
        }
        impl Write for BufWriter {
            fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
                self.store.borrow_mut().entry(self.key.clone()).or_default().extend_from_slice(buf);
                Ok(buf.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }

        let opener = move |key: &str| -> std::io::Result<Box<dyn Write>> {
            store_clone.borrow_mut().entry(key.to_string()).or_default();
            Ok(Box::new(BufWriter { key: key.to_string(), store: Rc::clone(&store_clone) }))
        };

        (store, opener)
    }

    #[test]
    fn buckets_lines_by_sorted_molecule_pair() {
        let input = "\
1\tgeneA\t1\t2\t3\tgeneB\t3\t4\n
3\tgeneC\t1\t2\t1\tgeneD\t3\t4\n
1\tgeneE\t1\t2\t1\tgeneF\t3\t4\n
";
        let (store, opener) = in_memory_sink();
        separate_datasets(input.as_bytes(), opener).unwrap();
        let store = store.borrow();

        // "1 vs 3" and "3 vs 1" both sort to key "1vs3", so both lines
        // land in the same bucket.
        let bucket_1v3 = String::from_utf8(store.get("1vs3").unwrap().clone()).unwrap();
        assert!(bucket_1v3.contains("geneA"));
        assert!(bucket_1v3.contains("geneC"));

        // Self-comparison "1 vs 1" gets its own bucket.
        let bucket_1v1 = String::from_utf8(store.get("1vs1").unwrap().clone()).unwrap();
        assert!(bucket_1v1.contains("geneE"));
    }

    #[test]
    fn lines_with_fewer_than_8_fields_are_skipped() {
        let input = "only\tseven\tfields\there\tfor\tthis\tline\n";
        let (store, opener) = in_memory_sink();
        separate_datasets(input.as_bytes(), opener).unwrap();
        assert!(store.borrow().is_empty());
    }

    #[test]
    fn perl_split_tab_trims_trailing_empty_fields() {
        // 8 real fields + 2 trailing empty fields from trailing tabs ->
        // Perl's split would trim the trailing empties, leaving exactly
        // 8 real fields (meets the >= 8 threshold).
        let line = "1\tgA\t1\t2\t3\tgB\t3\t4\t\t";
        assert_eq!(perl_split_tab(line).len(), 8);
    }

    #[test]
    fn perl_split_tab_preserves_internal_empty_fields() {
        let line = "a\tb\t\td\te\tf\tg\th";
        let fields = perl_split_tab(line);
        assert_eq!(fields.len(), 8);
        assert_eq!(fields[2], "");
    }

    #[test]
    fn molecule_pair_key_sorts_lexically() {
        assert_eq!(
            molecule_pair_key("10\tx\tx\tx\t2\tx\tx\tx"),
            Some("10vs2".to_string()) // lexical, not numeric: "10" < "2"
        );
    }

    #[test]
    fn output_path_for_matches_perl_naming_convention() {
        let path = PathBuf::from("/tmp/myMatches.txt");
        let out = output_path_for(&path, "1vs3");
        assert_eq!(out, PathBuf::from("/tmp/myMatches.txt.1vs3"));
    }

    #[test]
    fn empty_input_creates_no_buckets() {
        let (store, opener) = in_memory_sink();
        separate_datasets("".as_bytes(), opener).unwrap();
        assert!(store.borrow().is_empty());
    }

    #[test]
    fn multiple_pairs_remain_independent() {
        let input = "\
1\ta\t1\t2\t2\tb\t1\t2\n
3\tc\t1\t2\t4\td\t1\t2\n
";
        let (store, opener) = in_memory_sink();
        separate_datasets(input.as_bytes(), opener).unwrap();
        let store = store.borrow();
        assert_eq!(store.len(), 2);
        assert!(store.contains_key("1vs2"));
        assert!(store.contains_key("3vs4"));
    }
}
