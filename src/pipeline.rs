//! Port of `run_DAG_chainer.pl` — the orchestration layer that turns a
//! raw tabular match file into DAGchainer `.aligncoords` chain output,
//! by bucketing matches per molecule pair and chaining each bucket
//! (forward and reverse-complement) via the core [`crate::chain_matches`]
//! engine.
//!
//! ## Architectural simplification vs. the original Perl
//!
//! The original Perl writes a temporary "Delcher format" file per
//! molecule pair, shells out to the compiled `dagchainer` C++ binary
//! **twice** (forward, then `-r` reverse), captures its stdout text, and
//! re-parses that text back into Perl structs. Since this Rust port has
//! [`crate::chain_matches`] available as a native library function, none
//! of that is necessary: this module calls it directly, in-process, with
//! zero file I/O, zero subprocess spawning, and zero text round-tripping
//! per molecule pair. Consequently, the original's `-K`
//! (`KEEP_DELCHER_FILES`) troubleshooting flag has no equivalent here —
//! there is no transient file to keep.
//!
//! ## Input format
//!
//! Tab-separated, 9 columns:
//! `mol_1  acc_1  end5_1  end3_1  mol_2  acc_2  end5_2  end3_2  e_value`.
//!
//! ## Output format
//!
//! DAGchainer `.aligncoords` text: for each surviving chain, a header
//! line `## alignment {mol_A} vs. {mol_B} {match_header} (num aligned
//! pairs: {N}):`, followed by one tab-separated data line per aligned
//! pair: `mol_A acc_A end5_A end3_A mol_B acc_B end5_B end3_B e_value
//! dag_position_score`.
//!
//! ## Example
//! ```
//! use dagchainer_rs::pipeline::{parse_input_lines, run_pipeline, PipelineParams};
//!
//! let input = "\
//! chrA\tg1\t100\t200\tchrB\th1\t150\t250\t1e-50
//! chrA\tg2\t300\t400\tchrB\th2\t350\t450\t1e-40
//! chrA\tg3\t500\t600\tchrB\th3\t550\t650\t1e-30
//! ";
//!
//! let raw_matches = parse_input_lines(input.as_bytes()).unwrap();
//! let params = PipelineParams { min_num_aligned_pairs: 2, min_alignment_score: Some(10.0), ..Default::default() };
//!
//! let mut output = Vec::new();
//! run_pipeline(&raw_matches, &params, &mut output).unwrap();
//! assert!(String::from_utf8(output).unwrap().contains("num aligned pairs: 3"));
//! ```

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, Write};
use std::path::PathBuf;

use clap::Args as ClapArgs;
use thiserror::Error;

use crate::tools::find_tandem_dups::format_c_style_e;
use crate::{chain_matches, Alignment, DagChainerError, MatchPoint, Params as ChainParams};

#[derive(Debug, Error)]
pub enum PipelineError {
    #[error("cannot set both --ignore-tandem and --tandem-only (they are mutually exclusive)")]
    ConflictingTandemOptions,
    #[error(transparent)]
    Chain(#[from] DagChainerError),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

/// One parsed input line's raw fields, prior to any filtering.
#[derive(Debug, Clone)]
pub struct RawMatch {
    pub mol_1: String,
    pub acc_1: String,
    pub end5_1: i64,
    pub end3_1: i64,
    pub mol_2: String,
    pub acc_2: String,
    pub end5_2: i64,
    pub end3_2: i64,
    pub e_value: f64,
}

/// Parameters controlling the chaining pipeline, equivalent to
/// `run_DAG_chainer.pl`'s CLI options.
#[derive(Debug, Clone, PartialEq)]
pub struct PipelineParams {
    /// Length (bp) of a single gap (`-g`/`-G`).
    pub gap_length: i64,
    /// Gap open penalty, must be `<= 0` (`-o`).
    pub gap_open_penalty: f64,
    /// Gap extension penalty, must be `<= 0` (`-e`).
    pub gap_extension_penalty: f64,
    /// Maximum match score derived from e-value (`-M`).
    pub max_match_score: f64,
    /// If set, use this fixed score for every match instead of deriving
    /// one from e-value (`-Z`).
    pub constant_match_score: Option<f64>,
    /// Maximum distance (bp) allowed between two chained matches (`-D`).
    pub max_dist_between_matches: i64,
    /// Include self-molecule (tandem) comparisons (`-s`).
    pub include_self: bool,
    /// Ignore tandem/overlapping self-chain alignments (`-I`).
    pub ignore_tandem: bool,
    /// Report only tandem alignments; implies `include_self` (`-T`).
    pub tandem_only: bool,
    /// Minimum number of aligned pairs required to report a chain (`-A`).
    pub min_num_aligned_pairs: usize,
    /// Maximum e-value for an input match to be considered (`-E`).
    pub max_e_value: f64,
    /// Minimum alignment (chain) score to report (`-x`). If `None`, computed
    /// as `trunc(min_num_aligned_pairs * 0.5 * max_match_score)`, matching
    /// the original Perl's default-computation formula — but, unlike the
    /// original, an explicit `Some(0.0)` here is honored literally rather
    /// than silently overridden (see module docs on Perl truthiness traps).
    pub min_alignment_score: Option<f64>,
}

impl Default for PipelineParams {
    fn default() -> Self {
        PipelineParams {
            gap_length: 10_000,
            gap_open_penalty: 0.0,
            gap_extension_penalty: -3.0,
            max_match_score: 50.0,
            constant_match_score: None,
            max_dist_between_matches: 200_000,
            include_self: false,
            ignore_tandem: false,
            tandem_only: false,
            min_num_aligned_pairs: 6,
            max_e_value: 1e-5,
            min_alignment_score: None,
        }
    }
}

impl PipelineParams {
    pub fn validate(&self) -> Result<(), PipelineError> {
        if self.ignore_tandem && self.tandem_only {
            return Err(PipelineError::ConflictingTandemOptions);
        }
        Ok(())
    }

    /// Resolve the effective minimum alignment score: the explicit value
    /// if set, else computed from `min_num_aligned_pairs` and
    /// `max_match_score`, truncated toward zero (matching the original
    /// Perl's `int(...)`).
    pub fn effective_min_alignment_score(&self) -> f64 {
        match self.min_alignment_score {
            Some(v) => v,
            None => (self.min_num_aligned_pairs as f64 * 0.5 * self.max_match_score).trunc(),
        }
    }

    /// Whether self-molecule comparisons should be retained during
    /// parsing — `include_self`, or forced on by `tandem_only` (matching
    /// the original Perl's `if ($TANDEM_ONLY) { $INCLUDE_SELF = 1; }`).
    fn effective_include_self(&self) -> bool {
        self.include_self || self.tandem_only
    }
}

/// Rounds like Perl's `int(x + 0.5)` for non-negative x.
fn round_half_up(x: f64) -> i64 {
    (x + 0.5).floor() as i64
}

/// Returns true if `line` contains at least one ASCII letter, digit, or
/// underscore — equivalent to Perl's `/\w/` test, used to skip blank
/// lines.
fn has_word_char(line: &str) -> bool {
    line.chars().any(|c| c.is_ascii_alphanumeric() || c == '_')
}

/// Parse raw input lines into [`RawMatch`] records. Blank lines are
/// skipped; lines with fewer than 9 tab-separated fields, or non-numeric
/// coordinate/e-value fields, are skipped with a warning (a deliberate
/// improvement over the original Perl, which would silently process such
/// lines using `undef`-as-zero coercion under `use warnings`).
pub fn parse_input_lines<R: BufRead>(reader: R) -> Result<Vec<RawMatch>, PipelineError> {
    let mut out = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if !has_word_char(&line) {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            log::warn!("Skipping malformed line (fewer than 9 fields): {line:?}");
            continue;
        }

        let parse_i64 = |s: &str| s.parse::<i64>().ok();
        let parse_f64 = |s: &str| s.parse::<f64>().ok();

        let (Some(end5_1), Some(end3_1), Some(end5_2), Some(end3_2), Some(e_value)) = (
            parse_i64(fields[2]),
            parse_i64(fields[3]),
            parse_i64(fields[6]),
            parse_i64(fields[7]),
            parse_f64(fields[8]),
        ) else {
            log::warn!("Skipping line with non-numeric field: {line:?}");
            continue;
        };

        out.push(RawMatch {
            mol_1: fields[0].to_string(),
            acc_1: fields[1].to_string(),
            end5_1,
            end3_1,
            mol_2: fields[4].to_string(),
            acc_2: fields[5].to_string(),
            end5_2,
            end3_2,
            e_value,
        });
    }

    Ok(out)
}

/// One accession's genomic feature: molecule, coordinates, and midpoint.
#[derive(Debug, Clone)]
struct Feature {
    mol: String,
    acc: String,
    end5: i64,
    end3: i64,
    mid: i64,
}

/// A deduplicated accession-pair match: the two features involved, and
/// the (possibly updated-to-minimum) e-value between them.
#[derive(Debug, Clone)]
struct MatchRecord {
    feature_a: Feature,
    feature_b: Feature,
    e_value: f64,
}

/// All matches, plus the sorted-molecule-pair -> match-index grouping
/// used to drive per-molecule-pair chaining.
struct ParsedData {
    matches: Vec<MatchRecord>,
    mol_pair_indices: BTreeMap<(String, String), Vec<usize>>,
}

/// Build deduplicated per-accession features and per-molecule-pair match
/// groupings from raw input lines, applying all of
/// `run_DAG_chainer.pl`'s `parse_input_file` filtering and bookkeeping
/// rules.
///
/// ## Preserved behavior: shared mutable match state
///
/// The original Perl stores matches as hash references, so updating an
/// accession pair's e-value (when a lower-e-value duplicate is seen)
/// transparently updates the *same* object already referenced from the
/// per-molecule-pair list. This port replicates that by storing matches
/// in one `Vec`, with both `acc_pair_index` and `mol_pair_indices`
/// holding indices into it — updating `matches[idx].e_value` is visible
/// through both lookup paths, exactly mirroring Perl's shared-reference
/// semantics.
fn build_match_data(raw_matches: &[RawMatch], params: &PipelineParams) -> ParsedData {
    let mut acc_info: HashMap<String, Feature> = HashMap::new();
    let mut acc_pair_index: HashMap<(String, String), usize> = HashMap::new();
    let mut matches: Vec<MatchRecord> = Vec::new();
    let mut mol_pair_indices: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();

    let include_self = params.effective_include_self();

    for raw in raw_matches {
        let mut e_value = raw.e_value;
        if e_value < 1e-250 {
            e_value = 1e-250; // can't take log of 0
        }

        if raw.acc_1 == raw.acc_2 {
            continue;
        }
        if e_value > params.max_e_value {
            continue;
        }
        if raw.mol_1 == raw.mol_2 && !include_self {
            continue;
        }
        if raw.mol_1 != raw.mol_2 && params.tandem_only {
            continue;
        }

        let feature_1 = acc_info
            .entry(raw.acc_1.clone())
            .or_insert_with(|| Feature {
                mol: raw.mol_1.clone(),
                acc: raw.acc_1.clone(),
                end5: raw.end5_1,
                end3: raw.end3_1,
                mid: round_half_up((raw.end5_1 + raw.end3_1) as f64 / 2.0),
            })
            .clone();

        let feature_2 = acc_info
            .entry(raw.acc_2.clone())
            .or_insert_with(|| Feature {
                mol: raw.mol_2.clone(),
                acc: raw.acc_2.clone(),
                end5: raw.end5_2,
                end3: raw.end3_2,
                mid: round_half_up((raw.end5_2 + raw.end3_2) as f64 / 2.0),
            })
            .clone();

        let (mut feature_a, mut feature_b, mol_a, mol_b) = if raw.mol_1 < raw.mol_2 {
            (feature_1, feature_2, raw.mol_1.clone(), raw.mol_2.clone())
        } else {
            (feature_2, feature_1, raw.mol_2.clone(), raw.mol_1.clone())
        };

        if mol_a == mol_b && feature_a.mid > feature_b.mid {
            std::mem::swap(&mut feature_a, &mut feature_b);
        }

        let acc_pair_key = if raw.acc_1 <= raw.acc_2 {
            (raw.acc_1.clone(), raw.acc_2.clone())
        } else {
            (raw.acc_2.clone(), raw.acc_1.clone())
        };

        match acc_pair_index.get(&acc_pair_key) {
            Some(&idx) => {
                if matches[idx].e_value > e_value {
                    matches[idx].e_value = e_value;
                }
            }
            None => {
                let idx = matches.len();
                matches.push(MatchRecord { feature_a, feature_b, e_value });
                acc_pair_index.insert(acc_pair_key, idx);
                mol_pair_indices.entry((mol_a, mol_b)).or_default().push(idx);
            }
        }
    }

    ParsedData { matches, mol_pair_indices }
}

/// Compute a match's chaining score from its e-value, matching
/// `run_DAG_chainer.pl`'s `scoringF`: `-log10(e) * 10`, rounded to one
/// decimal place, capped at `max_match_score` — or the fixed
/// `constant_match_score`, if set.
fn scoring_f(e_value: f64, constant_match_score: Option<f64>, max_match_score: f64) -> f64 {
    if let Some(c) = constant_match_score {
        return c;
    }
    let mut score = -e_value.log10();
    score *= 10.0;
    score = (score + 0.5).floor();
    score /= 10.0;
    if score > max_match_score {
        max_match_score
    } else {
        score
    }
}

/// Run both forward and reverse-orientation chaining for one molecule
/// pair's matches, returning `(is_reverse, Alignment)` pairs in
/// forward-then-reverse order (matching the original's call sequence:
/// forward `run_DAG_chainer` invocation, then the `-r` reverse one).
fn chain_one_molecule_pair(
    indices: &[usize],
    matches: &[MatchRecord],
    params: &PipelineParams,
) -> Result<Vec<(bool, Alignment)>, PipelineError> {
    let points: Vec<MatchPoint> = indices
        .iter()
        .enumerate()
        .map(|(pair_id, &idx)| {
            let m = &matches[idx];
            MatchPoint {
                pair_id: pair_id as i64,
                x: m.feature_a.mid,
                y: m.feature_b.mid,
                score: scoring_f(m.e_value, params.constant_match_score, params.max_match_score),
            }
        })
        .collect();

    let chain_params_fwd = ChainParams {
        indel_score: params.gap_extension_penalty,
        gap_open_penalty: params.gap_open_penalty,
        bp_gap_size: params.gap_length,
        min_alignment_score: params.effective_min_alignment_score(),
        max_dist_between_matches: params.max_dist_between_matches,
        reverse_order: false,
    };
    let mut chain_params_rev = chain_params_fwd;
    chain_params_rev.reverse_order = true;

    let fwd = chain_matches(&points, &chain_params_fwd)?;
    let rev = chain_matches(&points, &chain_params_rev)?;

    let mut result = Vec::with_capacity(fwd.len() + rev.len());
    result.extend(fwd.into_iter().map(|a| (false, a)));
    result.extend(rev.into_iter().map(|a| (true, a)));
    Ok(result)
}

/// Rendered output for one alignment: an always-present informational
/// message (mirrors the original's `print "# $mol1 vs. $mol2 ..."`
/// console line), and the `.aligncoords`-format content to write, or
/// `None` if the alignment should be suppressed (too few pairs, or
/// filtered as tandem/overlapping/non-tandem per `-I`/`-T`).
struct AlignmentOutput {
    info_message: String,
    content: Option<String>,
}

/// Mirrors the original pipeline's double-rounding of the cumulative
/// chain score before it reaches the final `.aligncoords` file: the
/// compiled `dagchainer` C++ binary prints `path_score` via
/// `printf("%7.1f", ...)`, which rounds to 1 decimal place;
/// `run_DAG_chainer.pl` then re-parses that already-rounded text and
/// truncates it to an integer via `sprintf("%d", ...)`. This two-step
/// "round-to-1-decimal, then truncate" is **not** equivalent to directly
/// truncating the raw score — any true value whose second decimal digit
/// is `>= 5` gets rounded up at the printf step before truncation ever
/// sees it, which plain `.trunc()` on the raw value misses entirely.
fn dag_position_score_as_printed(raw_score: f64) -> i64 {
    let rounded_to_1_decimal = (raw_score * 10.0).round() / 10.0;
    rounded_to_1_decimal.trunc() as i64
}

/// Build the output for one alignment, replicating
/// `run_DAG_chainer.pl`'s `print_alignment` (minus its direct
/// `ALIGNCOORDS` file-handle writes, which the caller handles instead).
///
/// ## Note on tandem/overlap control-flow fidelity
///
/// The original Perl computes `IS_TANDEM`/`IS_OVERLAPPING` inside an
/// `if ((!$ignore_alignment) && ($mol1 eq $mol2) && ($IGNORE_TANDEM ||
/// $TANDEM_ONLY))` block, then applies a **separate, unconditionally
/// executed** `if`/`elsif` afterward to decide final suppression. This
/// port nests that final decision inside the same gated block instead.
/// This is provably equivalent in every reachable case: whenever the
/// outer gate is false, `IS_TANDEM`/`IS_OVERLAPPING` remain at their
/// initial `false` values in the original too (they are *only* ever set
/// inside that same gate), and plugging `false`/`false` into the
/// "separate" `if`/`elsif` either changes nothing (`ignore_tandem`
/// branch) or would only fire when `tandem_only` is set — but
/// `tandem_only` forces `mol1 == mol2` for every surviving match at
/// parse time (cross-molecule matches are discarded earlier whenever
/// `tandem_only` is set), so that branch can only ever actually fire
/// when the outer gate *was* true anyway.
fn build_alignment_output(
    mol1: &str,
    mol2: &str,
    is_reverse: bool,
    alignment: &Alignment,
    indices: &[usize],
    matches: &[MatchRecord],
    params: &PipelineParams,
) -> AlignmentOutput {
    let match_header = if is_reverse {
        format!("(reverse) Alignment #{}  score = {:.1}", alignment.id, alignment.score)
    } else {
        format!("Alignment #{}  score = {:.1}", alignment.id, alignment.score)
    };

    let num_aligned_pairs = alignment.entries.len();
    let mut ignore_alignment = num_aligned_pairs < params.min_num_aligned_pairs;

    let info_message = format!("# {mol1} vs. {mol2} {match_header} {num_aligned_pairs} aligned pairs.\n");

    if !ignore_alignment && mol1 == mol2 && (params.ignore_tandem || params.tandem_only) {
        let mut seen: HashSet<&str> = HashSet::new();
        let mut is_tandem = false;
        let mut acc1_coords: Vec<i64> = Vec::new();
        let mut acc2_coords: Vec<i64> = Vec::new();

        for entry in &alignment.entries {
            let idx = indices[entry.pair_id as usize];
            let m = &matches[idx];
            let (acc_a, acc_b) = (m.feature_a.acc.as_str(), m.feature_b.acc.as_str());

            if (acc_a != "-" && seen.contains(acc_a)) || (acc_b != "-" && seen.contains(acc_b)) {
                is_tandem = true;
                break;
            }
            seen.insert(acc_a);
            seen.insert(acc_b);

            acc1_coords.push(m.feature_a.end5);
            acc1_coords.push(m.feature_a.end3);
            acc2_coords.push(m.feature_b.end5);
            acc2_coords.push(m.feature_b.end3);
        }

        acc1_coords.sort_unstable();
        acc2_coords.sort_unstable();

        let is_overlapping = if acc1_coords.len() >= 2 && acc2_coords.len() >= 2 {
            let (lend_a, rend_a) = (acc1_coords[0], *acc1_coords.last().unwrap());
            let (lend_b, rend_b) = (acc2_coords[0], *acc2_coords.last().unwrap());
            lend_a <= rend_b && rend_a >= lend_b
        } else {
            false
        };

        if (params.ignore_tandem && (is_tandem || is_overlapping))
            || (params.tandem_only && !is_tandem)
        {
            ignore_alignment = true;
        }
    }

    if ignore_alignment {
        return AlignmentOutput { info_message, content: None };
    }

    let mut content = format!(
        "## alignment {mol1} vs. {mol2} {match_header} (num aligned pairs: {num_aligned_pairs}):\n"
    );

    for entry in &alignment.entries {
        let idx = indices[entry.pair_id as usize];
        let m = &matches[idx];
        content.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            m.feature_a.mol,
            m.feature_a.acc,
            m.feature_a.end5,
            m.feature_a.end3,
            m.feature_b.mol,
            m.feature_b.acc,
            m.feature_b.end5,
            m.feature_b.end3,
            format_c_style_e(m.e_value),
            dag_position_score_as_printed(entry.cumulative_score),
        ));
    }

    AlignmentOutput { info_message, content: Some(content) }
}

/// Run the full chaining pipeline over already-parsed matches, writing
/// `.aligncoords`-format output for every surviving chain, across every
/// molecule pair, in sorted-molecule-pair order (forward chains before
/// reverse chains within each pair) — matching
/// `run_DAG_chainer.pl`'s overall structure.
///
/// This function is pure with respect to its `writer`: it performs no
/// file I/O, logging, or console output of its own (unlike the original
/// Perl, which interleaves `print` progress messages with file writes).
/// Informational per-alignment messages are available via
/// [`run_pipeline_verbose`] for callers that want them.
pub fn run_pipeline<W: Write>(
    raw_matches: &[RawMatch],
    params: &PipelineParams,
    writer: &mut W,
) -> Result<(), PipelineError> {
    run_pipeline_verbose(raw_matches, params, writer, &mut |_msg| {})
}

/// As [`run_pipeline`], but additionally invokes `on_info` with each
/// alignment's informational message (mirroring the original Perl's
/// `print "# $mol1 vs. $mol2 ... aligned pairs.\n"` console line),
/// regardless of whether that alignment was ultimately suppressed.
pub fn run_pipeline_verbose<W: Write>(
    raw_matches: &[RawMatch],
    params: &PipelineParams,
    writer: &mut W,
    on_info: &mut dyn FnMut(&str),
) -> Result<(), PipelineError> {
    params.validate()?;

    let parsed = build_match_data(raw_matches, params);

    for ((mol_a, mol_b), indices) in &parsed.mol_pair_indices {
        let alignments = chain_one_molecule_pair(indices, &parsed.matches, params)?;
        for (is_reverse, alignment) in &alignments {
            let output = build_alignment_output(
                mol_a, mol_b, *is_reverse, alignment, indices, &parsed.matches, params,
            );
            on_info(&output.info_message);
            if let Some(content) = output.content {
                write!(writer, "{content}")?;
            }
        }
    }

    Ok(())
}

/// Convenience wrapper: parse input text, then run the full pipeline.
pub fn run_pipeline_from_reader<R: BufRead, W: Write>(
    reader: R,
    params: &PipelineParams,
    writer: &mut W,
) -> Result<(), PipelineError> {
    let raw_matches = parse_input_lines(reader)?;
    run_pipeline(&raw_matches, params, writer)
}

// ---------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Input match file (9-column tab format, see module docs).
    #[arg(short = 'i', long = "input", value_name = "FILE")]
    pub input: PathBuf,

    /// Length (bp) of a single gap.
    #[arg(short = 'g', long = "gap-length", default_value_t = 10_000)]
    pub gap_length: i64,

    /// Gap open penalty (must be <= 0).
    #[arg(short = 'o', long = "gap-open", default_value_t = 0.0, allow_hyphen_values = true)]
    pub gap_open_penalty: f64,

    /// Gap extension penalty (must be <= 0).
    #[arg(short = 'e', long = "gap-extension", default_value_t = -3.0, allow_hyphen_values = true)]
    pub gap_extension_penalty: f64,

    /// Maximum match score (used unless --constant-match-score is set).
    #[arg(short = 'M', long = "max-match-score", default_value_t = 50.0)]
    pub max_match_score: f64,

    /// Use a fixed score for every match instead of one derived from e-value.
    #[arg(short = 'Z', long = "constant-match-score")]
    pub constant_match_score: Option<f64>,

    /// Maximum distance (bp) allowed between two chained matches.
    #[arg(short = 'D', long = "max-dist", default_value_t = 200_000)]
    pub max_dist_between_matches: i64,

    /// Include self-molecule (tandem) comparisons.
    #[arg(short = 's', long = "include-self")]
    pub include_self: bool,

    /// Ignore tandem/overlapping self-chain alignments (requires -s).
    #[arg(short = 'I', long = "ignore-tandem")]
    pub ignore_tandem: bool,

    /// Report only tandem alignments (implies --include-self).
    #[arg(short = 'T', long = "tandem-only")]
    pub tandem_only: bool,

    /// Minimum number of aligned pairs required to report a chain.
    #[arg(short = 'A', long = "min-aligned-pairs", default_value_t = 6)]
    pub min_num_aligned_pairs: usize,

    /// Maximum e-value for an input match to be considered.
    #[arg(short = 'E', long = "max-evalue", default_value_t = 1e-5)]
    pub max_e_value: f64,

    /// Minimum alignment (chain) score to report (default: computed from
    /// --min-aligned-pairs and --max-match-score).
    #[arg(short = 'x', long = "min-alignment-score")]
    pub min_alignment_score: Option<f64>,
}

impl From<&Args> for PipelineParams {
    fn from(args: &Args) -> Self {
        PipelineParams {
            gap_length: args.gap_length,
            gap_open_penalty: args.gap_open_penalty,
            gap_extension_penalty: args.gap_extension_penalty,
            max_match_score: args.max_match_score,
            constant_match_score: args.constant_match_score,
            max_dist_between_matches: args.max_dist_between_matches,
            include_self: args.include_self,
            ignore_tandem: args.ignore_tandem,
            tandem_only: args.tandem_only,
            min_num_aligned_pairs: args.min_num_aligned_pairs,
            max_e_value: args.max_e_value,
            min_alignment_score: args.min_alignment_score,
        }
    }
}

/// CLI entry point — a drop-in (same-semantics) replacement for
/// `run_DAG_chainer.pl -i matchFile [options]`, writing output to
/// `<input>.aligncoords`, matching the original's naming convention so
/// downstream tools (`min-chain-len`, `to-chain-spans`, `to-gnuplot`,
/// etc.) work unchanged.
pub fn run(args: Args) -> anyhow::Result<()> {
    let params = PipelineParams::from(&args);
    params.validate()?;

    let input_file = std::fs::File::open(&args.input)
        .map_err(|e| anyhow::anyhow!("Cannot open {}: {e}", args.input.display()))?;
    let raw_matches = parse_input_lines(std::io::BufReader::new(input_file))?;
    println!("Parsed {} input match lines", raw_matches.len());

    let output_path = format!("{}.aligncoords", args.input.display());
    let mut output = std::io::BufWriter::new(std::fs::File::create(&output_path)?);

    let mut alignment_count = 0usize;
    run_pipeline_verbose(&raw_matches, &params, &mut output, &mut |msg| {
        print!("{msg}");
        alignment_count += 1;
    })?;
    output.flush()?;

    println!("Wrote {output_path} ({alignment_count} alignments considered)");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> PipelineParams {
        PipelineParams {
            min_num_aligned_pairs: 2,
            min_alignment_score: Some(10.0),
            max_e_value: 1.0,
            ..Default::default()
        }
    }

    #[test]
    fn scoring_f_computes_and_caps_correctly() {
        // -log10(1e-5) = 5.0 -> *10 = 50.0 -> round = 50.0 -> /10 = 5.0
        assert_eq!(scoring_f(1e-5, None, 50.0), 5.0);

        // -log10(1e-3) = 3.0 -> *10 = 30.0 -> round = 30.0 -> /10 = 3.0
        assert_eq!(scoring_f(1e-3, None, 50.0), 3.0);

        // -log10(1e-60) = 60.0 -> *10 = 600.0 -> round = 600.0 -> /10 = 60.0,
        // which exceeds max_match_score (50.0) and must be capped.
        assert_eq!(scoring_f(1e-60, None, 50.0), 50.0);

        // constant_match_score overrides everything, regardless of e-value.
        assert_eq!(scoring_f(0.5, Some(99.0), 50.0), 99.0);
    }

    #[test]
    fn build_match_data_dedups_accession_pairs_keeping_min_evalue() {
        let raw = vec![
            RawMatch {
                mol_1: "chrA".into(), acc_1: "g1".into(), end5_1: 100, end3_1: 200,
                mol_2: "chrB".into(), acc_2: "h1".into(), end5_2: 150, end3_2: 250,
                e_value: 1e-10,
            },
            RawMatch {
                mol_1: "chrA".into(), acc_1: "g1".into(), end5_1: 100, end3_1: 200,
                mol_2: "chrB".into(), acc_2: "h1".into(), end5_2: 150, end3_2: 250,
                e_value: 1e-20, // lower e-value for the same accession pair
            },
        ];
        let parsed = build_match_data(&raw, &default_params());
        assert_eq!(parsed.matches.len(), 1);
        assert_eq!(parsed.matches[0].e_value, 1e-20);
    }

    #[test]
    fn build_match_data_filters_self_accession() {
        let raw = vec![RawMatch {
            mol_1: "chrA".into(), acc_1: "g1".into(), end5_1: 100, end3_1: 200,
            mol_2: "chrA".into(), acc_2: "g1".into(), end5_2: 100, end3_2: 200,
            e_value: 1e-10,
        }];
        let parsed = build_match_data(&raw, &default_params());
        assert!(parsed.matches.is_empty());
    }

    #[test]
    fn build_match_data_filters_self_molecule_unless_include_self() {
        let raw = vec![RawMatch {
            mol_1: "chrA".into(), acc_1: "g1".into(), end5_1: 100, end3_1: 200,
            mol_2: "chrA".into(), acc_2: "g2".into(), end5_2: 300, end3_2: 400,
            e_value: 1e-10,
        }];
        let params_no_self = default_params();
        assert!(build_match_data(&raw, &params_no_self).matches.is_empty());

        let params_with_self = PipelineParams { include_self: true, ..default_params() };
        assert_eq!(build_match_data(&raw, &params_with_self).matches.len(), 1);
    }

    #[test]
    fn build_match_data_tandem_only_filters_cross_molecule_and_forces_include_self() {
        let raw = vec![
            RawMatch {
                mol_1: "chrA".into(), acc_1: "g1".into(), end5_1: 100, end3_1: 200,
                mol_2: "chrB".into(), acc_2: "h1".into(), end5_2: 150, end3_2: 250,
                e_value: 1e-10,
            },
            RawMatch {
                mol_1: "chrA".into(), acc_1: "g2".into(), end5_1: 300, end3_1: 400,
                mol_2: "chrA".into(), acc_2: "g3".into(), end5_2: 500, end3_2: 600,
                e_value: 1e-10,
            },
        ];
        let params = PipelineParams { tandem_only: true, ..default_params() };
        let parsed = build_match_data(&raw, &params);
        assert_eq!(parsed.matches.len(), 1); // only the self-molecule match survives
        assert_eq!(parsed.matches[0].feature_a.mol, "chrA");
        assert_eq!(parsed.matches[0].feature_b.mol, "chrA");
    }

    #[test]
    fn build_match_data_orders_features_by_lexical_molecule_then_midpoint() {
        // mol_1 ("chrB") > mol_2 ("chrA") lexically, so A/B should swap to
        // put "chrA" as feature_a.
        let raw = vec![RawMatch {
            mol_1: "chrB".into(), acc_1: "h1".into(), end5_1: 150, end3_1: 250,
            mol_2: "chrA".into(), acc_2: "g1".into(), end5_2: 100, end3_2: 200,
            e_value: 1e-10,
        }];
        let params = PipelineParams { include_self: false, ..default_params() };
        let parsed = build_match_data(&raw, &params);
        assert_eq!(parsed.matches[0].feature_a.mol, "chrA");
        assert_eq!(parsed.matches[0].feature_a.acc, "g1");
        assert_eq!(parsed.matches[0].feature_b.mol, "chrB");
        assert_eq!(parsed.matches[0].feature_b.acc, "h1");
    }

    #[test]
    fn build_match_data_same_molecule_orders_by_midpoint() {
        // Same molecule on both sides; feature with larger midpoint should
        // end up as feature_b regardless of input column order.
        let raw = vec![RawMatch {
            mol_1: "chrA".into(), acc_1: "g_high".into(), end5_1: 900, end3_1: 1000, // mid=950
            mol_2: "chrA".into(), acc_2: "g_low".into(), end5_2: 100, end3_2: 200,   // mid=150
            e_value: 1e-10,
        }];
        let params = PipelineParams { include_self: true, ..default_params() };
        let parsed = build_match_data(&raw, &params);
        assert_eq!(parsed.matches[0].feature_a.acc, "g_low");
        assert_eq!(parsed.matches[0].feature_b.acc, "g_high");
    }

    #[test]
    fn full_pipeline_chains_three_colinear_matches() {
        let input = "\
chrA\tg1\t100\t200\tchrB\th1\t150\t250\t1e-50
chrA\tg2\t300\t400\tchrB\th2\t350\t450\t1e-40
chrA\tg3\t500\t600\tchrB\th3\t550\t650\t1e-30
";
        let raw = parse_input_lines(input.as_bytes()).unwrap();
        let params = PipelineParams {
            min_num_aligned_pairs: 2,
            min_alignment_score: Some(10.0),
            max_e_value: 1.0,
            max_dist_between_matches: 100_000,
            gap_length: 100,
            gap_open_penalty: -1.0,
            gap_extension_penalty: -1.0,
            ..Default::default()
        };

        let mut output = Vec::new();
        run_pipeline(&raw, &params, &mut output).unwrap();
        let output = String::from_utf8(output).unwrap();

        assert!(output.contains("num aligned pairs: 3"));
        assert!(output.contains("g1"));
        assert!(output.contains("g2"));
        assert!(output.contains("g3"));
    }

    #[test]
    fn pipeline_produces_both_forward_and_reverse_chains() {
        let input = "\
chrA\tg1\t100\t200\tchrB\th1\t150\t250\t1e-50
chrA\tg2\t300\t400\tchrB\th2\t350\t450\t1e-40
chrA\tg3\t500\t600\tchrB\th3\t550\t650\t1e-30
";
        let raw = parse_input_lines(input.as_bytes()).unwrap();
        let params = PipelineParams {
            min_num_aligned_pairs: 1,
            min_alignment_score: Some(1.0),
            max_e_value: 1.0,
            max_dist_between_matches: 100_000,
            gap_length: 100,
            gap_open_penalty: -1.0,
            gap_extension_penalty: -1.0,
            ..Default::default()
        };

        let mut messages = Vec::new();
        let mut output = Vec::new();
        run_pipeline_verbose(&raw, &params, &mut output, &mut |msg| messages.push(msg.to_string())).unwrap();

        assert!(messages.iter().any(|m| !m.contains("(reverse)")), "expected at least one forward-chain message");
        assert!(messages.iter().any(|m| m.contains("(reverse)")), "expected at least one reverse-chain message");
    }

    #[test]
    fn min_e_value_filters_out_weak_matches() {
        let input = "chrA\tg1\t100\t200\tchrB\th1\t150\t250\t0.5\n"; // e-value above threshold
        let raw = parse_input_lines(input.as_bytes()).unwrap();
        let params = PipelineParams { max_e_value: 1e-5, ..default_params() };
        let parsed = build_match_data(&raw, &params);
        assert!(parsed.matches.is_empty());
    }

    #[test]
    fn conflicting_tandem_options_is_rejected() {
        let params = PipelineParams { ignore_tandem: true, tandem_only: true, ..Default::default() };
        assert!(matches!(params.validate(), Err(PipelineError::ConflictingTandemOptions)));
    }

    #[test]
    fn effective_min_alignment_score_computed_when_unset() {
        let params = PipelineParams {
            min_alignment_score: None,
            min_num_aligned_pairs: 6,
            max_match_score: 50.0,
            ..Default::default()
        };
        // trunc(6 * 0.5 * 50) = trunc(150.0) = 150.0
        assert_eq!(params.effective_min_alignment_score(), 150.0);
    }

    #[test]
    fn effective_min_alignment_score_explicit_zero_is_honored() {
        // This is the fixed Perl-truthiness-trap behavior: an explicit 0.0
        // must NOT be silently replaced by the computed default.
        let params = PipelineParams { min_alignment_score: Some(0.0), ..Default::default() };
        assert_eq!(params.effective_min_alignment_score(), 0.0);
    }

    #[test]
    fn parse_input_lines_skips_malformed_and_blank_lines() {
        let input = "\n   \ntoo\tfew\tfields\nchrA\tg1\t100\t200\tchrB\th1\t150\t250\t1e-10\n";
        let raw = parse_input_lines(input.as_bytes()).unwrap();
        assert_eq!(raw.len(), 1);
        assert_eq!(raw[0].acc_1, "g1");
    }

    #[test]
    fn empty_input_produces_empty_output() {
        let raw = parse_input_lines("".as_bytes()).unwrap();
        let mut output = Vec::new();
        run_pipeline(&raw, &default_params(), &mut output).unwrap();
        assert!(output.is_empty());
    }

    #[test]
    fn dag_position_score_matches_double_rounding_of_original_pipeline() {
        // True value rounds UP to the next integer at the 1-decimal step,
        // which plain truncation of the raw value would miss.
        assert_eq!(dag_position_score_as_printed(11631.96), 11632);
        // True value rounds DOWN (stays below the next integer) at the
        // 1-decimal step, matching plain truncation in this case.
        assert_eq!(dag_position_score_as_printed(11631.94), 11631);
        // Exactly on a whole number already.
        assert_eq!(dag_position_score_as_printed(100.0), 100);
    }
}
