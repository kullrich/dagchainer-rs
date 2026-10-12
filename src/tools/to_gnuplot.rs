//! Port of `dagchainer_to_gnuplot.pl` — render DAGchainer `.aligncoords`
//! chains as a gnuplot dot-plot, with contigs/scaffolds laid out along
//! each axis in descending-length order.
//!
//! See module-level deviations/quirks documented in the project's PR
//! description / commit message; summarized briefly on each function
//! below where most relevant.
//! ## Deviations from the original Perl (deliberate fixes)
//!
//! - `--accession` filtering now matches a chain if **either** `acc_a`
//!   or `acc_b` equals the requested value, rather than only checking
//!   `acc_a` as the original Perl did (despite its own help text
//!   implying both sides should be considered).
//! - Scaffolding tie-breaks for equal-length contigs are deterministic
//!   (descending length, then ascending name) rather than relying on
//!   Perl's non-reproducible hash key iteration order.
//! - A chain with zero data lines, or a data line before any header, is
//!   now a hard error rather than silently producing corrupted/malformed
//!   state.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use clap::Args as ClapArgs;
use regex::Regex;
use thiserror::Error;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// FASTA file for genome A.
    #[arg(long = "genomeA", value_name = "FILE")]
    pub genome_a: PathBuf,

    /// FASTA file for genome B.
    #[arg(long = "genomeB", value_name = "FILE")]
    pub genome_b: PathBuf,

    /// DAGchainer `.aligncoords` output file.
    #[arg(long = "aligncoords", value_name = "FILE")]
    pub aligncoords: PathBuf,

    /// Prefix for generated plotting files.
    #[arg(long = "prefix", default_value = "gp")]
    pub prefix: String,

    /// Minimum number of genes a chain ("diagonal") must have to be
    /// plotted.
    #[arg(long = "min_diag_genes", default_value_t = 2)]
    pub min_diag_genes: usize,

    /// Plot only chains whose `accA` equals this value (see module docs
    /// re: this only checking `accA`, matching the original Perl).
    #[arg(long = "accession")]
    pub accession: Option<String>,

    /// New, optional addition (not present in the original Perl): skip
    /// invoking the external `gnuplot` binary after writing the script
    /// and data files. Useful for testing/CI/headless environments.
    #[arg(long = "no-exec")]
    pub no_exec: bool,
}

#[derive(Debug, Error)]
pub enum ToGnuplotError {
    #[error("malformed header line (no 'X vs. Y' pattern found): {header:?}")]
    MalformedHeader { header: String },
    #[error("data line encountered before any header line: {line:?}")]
    DataBeforeHeader { line: String },
    #[error("chain {acc_a} vs {acc_b} has a header but zero data lines")]
    EmptyChain { acc_a: String, acc_b: String },
    #[error("no length found for contig {contig:?} (missing or zero in genome sequence)")]
    MissingContigLength { contig: String },
    #[error("no scaffold position for contig {contig:?} (side {side})")]
    MissingScaffoldPosition { side: char, contig: String },
    #[error("no data points to plot (all four output streams are empty)")]
    NoDataToPlot,
    #[error("malformed seq_lengths cache line: {line:?}")]
    MalformedCacheLine { line: String },
    #[error("FASTA file has sequence data before any '>' header: {line:?}")]
    FastaDataBeforeHeader { line: String },
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

// ---------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------

/// Raw (unsorted) gene boundary coordinates for one match, on each side.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MatchPair {
    /// (end5, end3) on molecule A — not necessarily `end5 <= end3`.
    pub a: (i64, i64),
    /// (end5, end3) on molecule B.
    pub b: (i64, i64),
}

/// One parsed DAGchainer alignment chain.
#[derive(Debug, Clone, PartialEq)]
pub struct DagChain {
    pub acc_a: String,
    pub acc_b: String,
    /// `'+'` or `'-'`; `'-'` indicates molecule B is on the reverse strand.
    pub orient: char,
    pub match_pairs: Vec<MatchPair>,
    /// (min, max) span across all match coordinates on side A.
    pub range_a: (i64, i64),
    /// (min, max) span across all match coordinates on side B.
    pub range_b: (i64, i64),
}

/// The four gnuplot data streams produced by [`render_plot_data`].
#[derive(Debug, Default, Clone, PartialEq)]
pub struct PlotData {
    pub chains_forward: String,
    pub chains_reverse: String,
    pub dots_forward: String,
    pub dots_reverse: String,
}

// ---------------------------------------------------------------------
// Sequence length parsing (FASTA + cache)
// ---------------------------------------------------------------------

/// Parse a FASTA file's sequence lengths, keyed by the first
/// whitespace-delimited token after `>` on each header line.
///
/// Matches the original Perl's own explicit error: sequence data
/// appearing before any header is a hard error, not silently ignored.
pub fn parse_seqlengths_from_fasta_text<R: BufRead>(
    reader: R,
) -> Result<BTreeMap<String, i64>, ToGnuplotError> {
    let mut lengths: BTreeMap<String, i64> = BTreeMap::new();
    let mut acc: Option<String> = None;

    for line in reader.lines() {
        let line = line?;
        if let Some(stripped) = line.strip_prefix('>') {
            let acc_name = stripped.split_whitespace().next().unwrap_or("").to_string();
            acc = Some(acc_name);
        } else {
            let seq_len = line.chars().filter(|c| !c.is_whitespace()).count() as i64;
            if seq_len == 0 {
                continue;
            }
            match &acc {
                Some(a) => {
                    *lengths.entry(a.clone()).or_insert(0) += seq_len;
                }
                None => {
                    return Err(ToGnuplotError::FastaDataBeforeHeader { line });
                }
            }
        }
    }

    Ok(lengths)
}

/// Parse a previously written `<genome>.seq_lengths` cache file:
/// tab-separated `accession<TAB>length` per line.
pub fn parse_seqlengths_cache<R: BufRead>(
    reader: R,
) -> Result<BTreeMap<String, i64>, ToGnuplotError> {
    let mut lengths = BTreeMap::new();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let mut parts = line.splitn(2, '\t');
        let acc = parts.next().unwrap_or("").to_string();
        let len: i64 = parts
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| ToGnuplotError::MalformedCacheLine { line: line.clone() })?;
        lengths.insert(acc, len);
    }
    Ok(lengths)
}

/// Load sequence lengths for `genome_path`, using (and creating, if
/// absent) a `<genome_path>.seq_lengths` cache file alongside it. This
/// is the only I/O-performing (non-pure) function in sequence-length
/// handling; the parsing itself is delegated to the pure functions
/// above.
pub fn load_seq_lengths(genome_path: &Path) -> Result<BTreeMap<String, i64>, ToGnuplotError> {
    let mut cache_path = genome_path.as_os_str().to_os_string();
    cache_path.push(".seq_lengths");
    let cache_path = PathBuf::from(cache_path);

    if cache_path.exists() {
        let f = File::open(&cache_path)?;
        parse_seqlengths_cache(BufReader::new(f))
    } else {
        let f = File::open(genome_path)?;
        let lengths = parse_seqlengths_from_fasta_text(BufReader::new(f))?;

        let mut out = File::create(&cache_path)?;
        for (acc, len) in &lengths {
            writeln!(out, "{acc}\t{len}")?;
        }
        Ok(lengths)
    }
}

// ---------------------------------------------------------------------
// DAGchainer output parsing
// ---------------------------------------------------------------------

/// Parse DAGchainer `.aligncoords`-style output into a list of
/// [`DagChain`]s, including computed per-side coordinate spans.
///
/// ## Deviations from the original Perl (documented, deliberate)
///
/// - A data line before any header is a hard
///   [`ToGnuplotError::DataBeforeHeader`], rather than relying on
///   Perl's autovivification of an undef scalar into a malformed chain.
/// - A chain whose header is followed by zero data lines is a hard
///   [`ToGnuplotError::EmptyChain`], rather than silently producing a
///   corrupted `(undef, undef)` range (which Perl would later treat as
///   `(0, 0)` in numeric context).
/// - Malformed data lines (fewer than 8 tab fields, or non-numeric
///   coordinates) are skipped with a `log::warn!`, rather than silently
///   parsing as `0`.
pub fn parse_dagchainer_output<R: BufRead>(reader: R) -> Result<Vec<DagChain>, ToGnuplotError> {
    let header_re = Regex::new(r"(\S+) vs\. (\S+)").unwrap();
    let mut chains: Vec<DagChain> = Vec::new();

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('#') {
            let caps = header_re
                .captures(&line)
                .ok_or_else(|| ToGnuplotError::MalformedHeader { header: line.clone() })?;
            let acc_a = caps[1].to_string();
            let acc_b = caps[2].to_string();
            let orient = if line.contains("(reverse)") { '-' } else { '+' };

            chains.push(DagChain {
                acc_a,
                acc_b,
                orient,
                match_pairs: Vec::new(),
                range_a: (0, 0),
                range_b: (0, 0),
            });
        } else if line.chars().any(|c| c.is_ascii_alphanumeric() || c == '_') {
            let chain = chains
                .last_mut()
                .ok_or_else(|| ToGnuplotError::DataBeforeHeader { line: line.clone() })?;

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 8 {
                log::warn!("Skipping malformed data line (fewer than 8 fields): {line:?}");
                continue;
            }
            let parse = |s: &str| s.parse::<i64>().ok();
            match (parse(fields[2]), parse(fields[3]), parse(fields[6]), parse(fields[7])) {
                (Some(e5a), Some(e3a), Some(e5b), Some(e3b)) => {
                    chain.match_pairs.push(MatchPair { a: (e5a, e3a), b: (e5b, e3b) });
                }
                _ => {
                    log::warn!("Skipping data line with non-numeric coordinates: {line:?}");
                }
            }
        }
        // else: blank or punctuation-only line; silently skipped,
        // matching the original Perl's `elsif (/\w/)` guard.
    }

    for chain in &mut chains {
        if chain.match_pairs.is_empty() {
            return Err(ToGnuplotError::EmptyChain {
                acc_a: chain.acc_a.clone(),
                acc_b: chain.acc_b.clone(),
            });
        }

        let mut coords_a = Vec::with_capacity(chain.match_pairs.len() * 2);
        let mut coords_b = Vec::with_capacity(chain.match_pairs.len() * 2);
        for mp in &chain.match_pairs {
            coords_a.push(mp.a.0);
            coords_a.push(mp.a.1);
            coords_b.push(mp.b.0);
            coords_b.push(mp.b.1);
        }

        chain.range_a = (
            *coords_a.iter().min().unwrap(),
            *coords_a.iter().max().unwrap(),
        );
        chain.range_b = (
            *coords_b.iter().min().unwrap(),
            *coords_b.iter().max().unwrap(),
        );
    }

    Ok(chains)
}

// ---------------------------------------------------------------------
// Filtering, scaffolding, reciprocation
// ---------------------------------------------------------------------

/// Filter chains by optional accession and minimum diagonal length.
///
/// ## Deviation from the original Perl (deliberate fix)
///
/// The original Perl's `if ($accession && $dagchain->{accA} ne
/// $accession) { next; }` only ever checked `accA`, despite the
/// script's help text describing `--accession` as filtering "chains
/// involving molecule accession" — which implies either side should
/// match. This port checks **both** `acc_a` and `acc_b`, matching a
/// chain if either side equals the requested accession.
///
/// ## Preserved quirk: contig sets are populated before the length filter
///
/// `contigs_a`/`contigs_b` are populated for every chain that passes the
/// accession filter, **regardless** of whether it also passes
/// `min_diag_genes` — so a too-short chain can still cause its molecule
/// to appear in the scaffolding/axis labels.
pub fn filter_chains(
    chains: Vec<DagChain>,
    accession_filter: Option<&str>,
    min_diag_genes: usize,
) -> (Vec<DagChain>, BTreeSet<String>, BTreeSet<String>) {
    let mut contigs_a = BTreeSet::new();
    let mut contigs_b = BTreeSet::new();
    let mut kept = Vec::new();

    for chain in chains {
        if let Some(acc) = accession_filter {
            if chain.acc_a != acc && chain.acc_b != acc {
                continue;
            }
        }

        contigs_a.insert(chain.acc_a.clone());
        contigs_b.insert(chain.acc_b.clone());

        if chain.match_pairs.len() >= min_diag_genes {
            kept.push(chain);
        }
    }

    (kept, contigs_a, contigs_b)
}

/// Assign each contig a scaffold (axis) offset, laid out in descending
/// length order, longest first, so long contigs occupy the
/// lower-position end of the axis.
///
/// ## Deviation from the original Perl (deliberate fix)
///
/// Ties (contigs of equal length) are broken by ascending contig name,
/// for deterministic, reproducible output. The original Perl's
/// `reverse sort {...}  keys %$contigs_href` relies on Perl's hash key
/// iteration order for tie-breaking (via `sort`'s documented
/// stability), which has been randomized per-process since Perl 5.18 —
/// so the original's tie order was already non-reproducible across
/// separate runs on identical input.
pub fn create_scaffolding(
    seq_lengths: &BTreeMap<String, i64>,
    contigs: &BTreeSet<String>,
) -> Result<BTreeMap<String, i64>, ToGnuplotError> {
    let mut entries: Vec<(&String, i64)> = Vec::with_capacity(contigs.len());
    for contig in contigs {
        let len = seq_lengths.get(contig).copied().unwrap_or(0);
        if len == 0 {
            return Err(ToGnuplotError::MissingContigLength { contig: contig.clone() });
        }
        entries.push((contig, len));
    }

    entries.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));

    let mut scaffolding = BTreeMap::new();
    let mut pos: i64 = 1;
    for (contig, len) in entries {
        scaffolding.insert(contig.clone(), pos);
        pos += len;
    }
    Ok(scaffolding)
}

/// For self-genome comparisons, duplicate each chain into both
/// orientations: the original (A vs B) and its reciprocal (B vs A, with
/// every match pair's sides swapped), interleaved in input order.
pub fn reciprocate_chains(chains: &[DagChain]) -> Vec<DagChain> {
    let mut result = Vec::with_capacity(chains.len() * 2);
    for chain in chains {
        let reciprocal = DagChain {
            acc_a: chain.acc_b.clone(),
            acc_b: chain.acc_a.clone(),
            orient: chain.orient,
            range_a: chain.range_b,
            range_b: chain.range_a,
            match_pairs: chain
                .match_pairs
                .iter()
                .map(|mp| MatchPair { a: mp.b, b: mp.a })
                .collect(),
        };
        result.push(chain.clone());
        result.push(reciprocal);
    }
    result
}

// ---------------------------------------------------------------------
// Rendering: plot data streams + gnuplot script
// ---------------------------------------------------------------------

fn two_point_block(p1: (i64, i64), p2: (i64, i64)) -> String {
    format!("{}\t{}\n{}\t{}\n\n", p1.0, p1.1, p2.0, p2.1)
}

/// Render the four gnuplot data streams (chain-bounds and per-gene dots,
/// split by forward/reverse orientation) for a set of chains, given
/// each side's scaffold offsets.
pub fn render_plot_data(
    chains: &[DagChain],
    scaffolds_a: &BTreeMap<String, i64>,
    scaffolds_b: &BTreeMap<String, i64>,
) -> Result<PlotData, ToGnuplotError> {
    let mut data = PlotData::default();

    for chain in chains {
        let scaff_pos_a = *scaffolds_a.get(&chain.acc_a).ok_or_else(|| {
            ToGnuplotError::MissingScaffoldPosition { side: 'A', contig: chain.acc_a.clone() }
        })?;
        let scaff_pos_b = *scaffolds_b.get(&chain.acc_b).ok_or_else(|| {
            ToGnuplotError::MissingScaffoldPosition { side: 'B', contig: chain.acc_b.clone() }
        })?;

        let range_a = (chain.range_a.0 + scaff_pos_a, chain.range_a.1 + scaff_pos_a);
        let mut range_b = (chain.range_b.0 + scaff_pos_b, chain.range_b.1 + scaff_pos_b);

        let reverse = chain.orient == '-';
        if reverse {
            range_b = (range_b.1, range_b.0);
        }

        let (chain_buf, dots_buf) = if reverse {
            (&mut data.chains_reverse, &mut data.dots_reverse)
        } else {
            (&mut data.chains_forward, &mut data.dots_forward)
        };

        chain_buf.push_str(&two_point_block(
            (range_a.0, range_b.0),
            (range_a.1, range_b.1),
        ));

        for mp in &chain.match_pairs {
            let gene_a = (mp.a.0 + scaff_pos_a, mp.a.1 + scaff_pos_a);
            let gene_b = (mp.b.0 + scaff_pos_b, mp.b.1 + scaff_pos_b);
            dots_buf.push_str(&two_point_block((gene_a.0, gene_b.0), (gene_a.1, gene_b.1)));
        }
    }

    Ok(data)
}

/// Build the list of gnuplot `xtics`/`ytics` entries from a scaffold
/// offset map, sorted ascending by offset (equivalently: longest-contig
/// first, since offsets are assigned in that order and are strictly
/// increasing).
pub fn build_tics(scaffolds: &BTreeMap<String, i64>) -> Vec<String> {
    let mut entries: Vec<(&String, &i64)> = scaffolds.iter().collect();
    entries.sort_by_key(|(_, pos)| **pos);
    entries.into_iter().map(|(acc, pos)| format!("\"{acc}\" {pos}")).collect()
}

/// Build the full gnuplot script text.
pub fn build_gnuplot_script(
    prefix: &str,
    xtics: &[String],
    ytics: &[String],
    plot_data: &PlotData,
) -> Result<String, ToGnuplotError> {
    let mut script = String::new();
    script.push_str(
        "set size 1,1\n\
         set grid\n\
         unset key\n\
         set border 10\n\
         set ticscale 0 0\n\
         set xrange [0:]\n\
         set yrange [0:]\n\
         set xtics rotate\n",
    );
    script.push_str(&format!("set xtics (\\\n{})\n", xtics.join(", \\\n")));
    script.push_str(&format!("set ytics (\\\n{})\n", ytics.join(", \\\n")));
    script.push_str(
        "set style line 1  lt 1 lw 2 pt 6 ps 1\n\
         set style line 2  lt 2 lw 2 pt 6 ps 1\n\
         set style line 3  lt 3 lw 2 pt 6 ps 1\n\
         set style line 4  lt 4 lw 2 pt 6 ps 1\n",
    );

    let mut files_to_plot = Vec::new();
    if !plot_data.chains_forward.is_empty() {
        files_to_plot.push(format!(" \"{prefix}.chains_F\" w lp ls 1 "));
    }
    if !plot_data.chains_reverse.is_empty() {
        files_to_plot.push(format!(" \"{prefix}.chains_R\" w lp ls 2 "));
    }
    if !plot_data.dots_forward.is_empty() {
        files_to_plot.push(format!(" \"{prefix}.dots_F\" w lp ls 3 "));
    }
    if !plot_data.dots_reverse.is_empty() {
        files_to_plot.push(format!(" \"{prefix}.dots_R\" w lp ls 4 "));
    }

    if files_to_plot.is_empty() {
        return Err(ToGnuplotError::NoDataToPlot);
    }

    script.push_str(&format!("plot {}\n\n", files_to_plot.join(", \\\n")));
    Ok(script)
}

// ---------------------------------------------------------------------
// CLI entry point (I/O orchestration)
// ---------------------------------------------------------------------

/// CLI entry point — a drop-in replacement for
/// `dagchainer_to_gnuplot.pl --genomeA ... --genomeB ... --aligncoords ...`.
pub fn run(args: Args) -> anyhow::Result<()> {
    let seq_lengths_a = load_seq_lengths(&args.genome_a)?;
    let seq_lengths_b = load_seq_lengths(&args.genome_b)?;

    let aligncoords_file = File::open(&args.aligncoords)?;
    let mut chains = parse_dagchainer_output(BufReader::new(aligncoords_file))?;

    // Preserved quirk: self-comparison is plain string/path equality of
    // the two CLI arguments, not file-content comparison.
    if args.genome_a == args.genome_b {
        chains = reciprocate_chains(&chains);
    }

    let (chains, contigs_a, contigs_b) =
        filter_chains(chains, args.accession.as_deref(), args.min_diag_genes);

    let scaffolds_a = create_scaffolding(&seq_lengths_a, &contigs_a)?;
    let scaffolds_b = create_scaffolding(&seq_lengths_b, &contigs_b)?;

    let plot_data = render_plot_data(&chains, &scaffolds_a, &scaffolds_b)?;

    std::fs::write(format!("{}.chains_F", args.prefix), &plot_data.chains_forward)?;
    std::fs::write(format!("{}.chains_R", args.prefix), &plot_data.chains_reverse)?;
    std::fs::write(format!("{}.dots_F", args.prefix), &plot_data.dots_forward)?;
    std::fs::write(format!("{}.dots_R", args.prefix), &plot_data.dots_reverse)?;

    let xtics = build_tics(&scaffolds_a);
    let ytics = build_tics(&scaffolds_b);
    let script = build_gnuplot_script(&args.prefix, &xtics, &ytics, &plot_data)?;

    let script_path = format!("{}.gnuplot_script", args.prefix);
    std::fs::write(&script_path, &script)?;

    if !args.no_exec {
        let status = std::process::Command::new("gnuplot")
            .arg(&script_path)
            .arg("-")
            .status()
            .map_err(|e| anyhow::anyhow!("Failed to invoke gnuplot: {e}"))?;
        if !status.success() {
            anyhow::bail!("gnuplot exited with status: {status}");
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_input() -> &'static str {
        "\
## alignment chr1 vs. chr2 (num aligned pairs: 2):
chr1\tgeneA\t100\t200\tchr2\tgeneB\t150\t250\t1e-30
chr1\tgeneC\t300\t400\tchr2\tgeneD\t350\t450\t1e-20
## alignment chr1 vs. chr3 (reverse) (num aligned pairs: 1):
chr1\tgeneE\t500\t600\tchr3\tgeneF\t550\t650\t1e-10
"
    }

    #[test]
    fn parses_chains_with_correct_ranges_and_orientation() {
        let chains = parse_dagchainer_output(sample_input().as_bytes()).unwrap();
        assert_eq!(chains.len(), 2);

        assert_eq!(chains[0].acc_a, "chr1");
        assert_eq!(chains[0].acc_b, "chr2");
        assert_eq!(chains[0].orient, '+');
        assert_eq!(chains[0].range_a, (100, 400));
        assert_eq!(chains[0].range_b, (150, 450));

        assert_eq!(chains[1].acc_a, "chr1");
        assert_eq!(chains[1].acc_b, "chr3");
        assert_eq!(chains[1].orient, '-');
        assert_eq!(chains[1].range_a, (500, 600));
        assert_eq!(chains[1].range_b, (550, 650));
    }

    #[test]
    fn data_before_header_is_an_error() {
        let input = "chr1\tgeneA\t1\t2\tchr2\tgeneB\t1\t2\t1e-5\n";
        let result = parse_dagchainer_output(input.as_bytes());
        assert!(matches!(result, Err(ToGnuplotError::DataBeforeHeader { .. })));
    }

    #[test]
    fn empty_chain_is_an_error() {
        let input = "## alignment chr1 vs. chr2 (num aligned pairs: 0):\n";
        let result = parse_dagchainer_output(input.as_bytes());
        assert!(matches!(result, Err(ToGnuplotError::EmptyChain { .. })));
    }

    #[test]
    fn malformed_header_is_an_error() {
        let input = "## no vs pattern here\nchr1\tg1\t1\t2\tchr2\tg2\t1\t2\t1e-5\n";
        let result = parse_dagchainer_output(input.as_bytes());
        assert!(matches!(result, Err(ToGnuplotError::MalformedHeader { .. })));
    }

    #[test]
    fn filter_chains_respects_accession_and_min_diag_genes() {
        let chains = parse_dagchainer_output(sample_input().as_bytes()).unwrap();
        let (kept, contigs_a, contigs_b) = filter_chains(chains, None, 2);

        assert_eq!(kept.len(), 1); // only the chr1-vs-chr2 chain has 2 pairs
        assert_eq!(kept[0].acc_b, "chr2");

        // Both chains passed the (absent) accession filter, so both
        // contribute to the contig sets, even though one was filtered
        // out by min_diag_genes.
        assert!(contigs_a.contains("chr1"));
        assert!(contigs_b.contains("chr2"));
        assert!(contigs_b.contains("chr3"));
    }

    #[test]
    fn filter_chains_accession_filter_matches_either_side() {
        let chains = parse_dagchainer_output(sample_input().as_bytes()).unwrap();

        // "chr2" only ever appears as acc_b in this input. With the fix,
        // filtering by accession="chr2" should now match the chr1-vs-chr2
        // chain (previously, with the acc_a-only quirk, this would have
        // returned nothing).
        let (kept, _, _) = filter_chains(chains.clone(), Some("chr2"), 1);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].acc_b, "chr2");

        // "chr1" appears as acc_a on both chains, so filtering by "chr1"
        // should still match both.
        let (kept, _, _) = filter_chains(chains, Some("chr1"), 1);
        assert_eq!(kept.len(), 2);
    }

    #[test]
    fn create_scaffolding_orders_longest_first_with_deterministic_ties() {
        let mut lengths = BTreeMap::new();
        lengths.insert("short".to_string(), 100);
        lengths.insert("long".to_string(), 1000);
        lengths.insert("tie_b".to_string(), 500);
        lengths.insert("tie_a".to_string(), 500);

        let contigs: BTreeSet<String> =
            ["short", "long", "tie_b", "tie_a"].iter().map(|s| s.to_string()).collect();

        let scaffolding = create_scaffolding(&lengths, &contigs).unwrap();

        // "long" (1000) gets the first slot.
        assert_eq!(scaffolding["long"], 1);
        // Then the two 500-length contigs, tie broken alphabetically:
        // "tie_a" before "tie_b".
        assert_eq!(scaffolding["tie_a"], 1001);
        assert_eq!(scaffolding["tie_b"], 1501);
        // "short" comes last.
        assert_eq!(scaffolding["short"], 2001);
    }

    #[test]
    fn create_scaffolding_errors_on_missing_length() {
        let lengths = BTreeMap::new();
        let contigs: BTreeSet<String> = ["chrX".to_string()].into_iter().collect();
        let result = create_scaffolding(&lengths, &contigs);
        assert!(matches!(result, Err(ToGnuplotError::MissingContigLength { .. })));
    }

    #[test]
    fn reciprocate_chains_interleaves_original_and_swapped() {
        let chains = parse_dagchainer_output(
            "## alignment chr1 vs. chr2 (num aligned pairs: 1):\nchr1\tg1\t10\t20\tchr2\tg2\t30\t40\t1e-5\n"
                .as_bytes(),
        )
            .unwrap();

        let recip = reciprocate_chains(&chains);
        assert_eq!(recip.len(), 2);
        assert_eq!(recip[0].acc_a, "chr1");
        assert_eq!(recip[0].acc_b, "chr2");
        assert_eq!(recip[1].acc_a, "chr2");
        assert_eq!(recip[1].acc_b, "chr1");
        assert_eq!(recip[1].range_a, recip[0].range_b);
        assert_eq!(recip[1].match_pairs[0].a, recip[0].match_pairs[0].b);
    }

    #[test]
    fn render_plot_data_places_forward_and_reverse_chains_correctly() {
        let chains = parse_dagchainer_output(sample_input().as_bytes()).unwrap();
        let (chains, contigs_a, contigs_b) = filter_chains(chains, None, 1);

        let mut lengths_a = BTreeMap::new();
        lengths_a.insert("chr1".to_string(), 10_000);
        let mut lengths_b = BTreeMap::new();
        lengths_b.insert("chr2".to_string(), 10_000);
        lengths_b.insert("chr3".to_string(), 10_000);

        let scaffolds_a = create_scaffolding(&lengths_a, &contigs_a).unwrap();
        let scaffolds_b = create_scaffolding(&lengths_b, &contigs_b).unwrap();

        let plot_data = render_plot_data(&chains, &scaffolds_a, &scaffolds_b).unwrap();

        assert!(!plot_data.chains_forward.is_empty());
        assert!(!plot_data.chains_reverse.is_empty());
        assert!(!plot_data.dots_forward.is_empty());
        assert!(!plot_data.dots_reverse.is_empty());
    }

    #[test]
    fn build_gnuplot_script_errors_when_nothing_to_plot() {
        let empty = PlotData::default();
        let result = build_gnuplot_script("gp", &[], &[], &empty);
        assert!(matches!(result, Err(ToGnuplotError::NoDataToPlot)));
    }

    #[test]
    fn build_gnuplot_script_includes_only_nonempty_streams() {
        let data = PlotData { chains_forward: "1\t2\n".to_string(), ..Default::default() };

        let script = build_gnuplot_script("myprefix", &["\"chr1\" 1".into()], &["\"chr2\" 1".into()], &data).unwrap();
        assert!(script.contains("myprefix.chains_F"));
        assert!(!script.contains("myprefix.chains_R"));
        assert!(!script.contains("myprefix.dots_F"));
    }

    #[test]
    fn build_tics_sorts_by_position() {
        let mut scaffolds = BTreeMap::new();
        scaffolds.insert("b".to_string(), 500);
        scaffolds.insert("a".to_string(), 1);

        let tics = build_tics(&scaffolds);
        assert_eq!(tics, vec!["\"a\" 1".to_string(), "\"b\" 500".to_string()]);
    }

    #[test]
    fn parse_fasta_lengths_sums_multiline_sequences() {
        let fasta = ">chr1 some description\nACGT\nACGT\n>chr2\nAC\n";
        let lengths = parse_seqlengths_from_fasta_text(fasta.as_bytes()).unwrap();
        assert_eq!(lengths["chr1"], 8);
        assert_eq!(lengths["chr2"], 2);
    }

    #[test]
    fn parse_fasta_lengths_errors_on_data_before_header() {
        let fasta = "ACGT\n>chr1\nACGT\n";
        let result = parse_seqlengths_from_fasta_text(fasta.as_bytes());
        assert!(matches!(result, Err(ToGnuplotError::FastaDataBeforeHeader { .. })));
    }

    #[test]
    fn parse_seqlengths_cache_roundtrip() {
        let cache = "chr1\t1000\nchr2\t2000\n";
        let lengths = parse_seqlengths_cache(cache.as_bytes()).unwrap();
        assert_eq!(lengths["chr1"], 1000);
        assert_eq!(lengths["chr2"], 2000);
    }
}
