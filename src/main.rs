//! CLI binary: drop-in replacement for the original `dagchainer` C++ tool.
//!
//! ```text
//! dagchainer -G <gap_bp> -O <gap_open> -E <gap_ext> -S <min_score>
//!            -D <max_dist> -F <file> [-r]
//! ```

use std::fs::File;
use std::io::{BufReader, Write};

use clap::Parser;
use dagchainer_rs::{chain_matches, read_matches, Params};

/// Rust port of DAGchainer - chain colinear pairwise matches into synteny blocks.
#[derive(Parser, Debug)]
#[command(
    name = "dagchainer",
    about = "Do DP on a DAG of matches to produce chains of colinear matches.",
    long_about = None
)]
struct Cli {
    /// Length (bp) of a single gap.
    #[arg(short = 'G', long = "gap-length", value_name = "BP")]
    bp_gap_size: i64,

    /// Gap open penalty (must be <= 0).
    #[arg(short = 'O', long = "gap-open", value_name = "PENALTY", allow_hyphen_values = true)]
    gap_open_penalty: f64,

    /// Gap extension penalty (must be <= 0).
    #[arg(short = 'E', long = "gap-extension", value_name = "PENALTY", allow_hyphen_values = true)]
    indel_score: f64,

    /// Minimum alignment (chain) score to report.
    #[arg(short = 'S', long = "min-score", value_name = "SCORE")]
    min_alignment_score: f64,

    /// Maximum distance (bp) allowed between two chained matches.
    #[arg(short = 'D', long = "max-dist", default_value_t = 100_000)]
    max_dist_between_matches: i64,

    /// Input match file: `pair_id  x  y  score` per line.
    #[arg(short = 'F', long = "file", value_name = "FILE")]
    filename: String,

    /// Reverse the second coordinate list (e.g. minus-strand synteny).
    #[arg(short = 'r', long = "reverse")]
    reverse_order: bool,
}

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let params = Params {
        indel_score: cli.indel_score,
        gap_open_penalty: cli.gap_open_penalty,
        bp_gap_size: cli.bp_gap_size,
        min_alignment_score: cli.min_alignment_score,
        max_dist_between_matches: cli.max_dist_between_matches,
        reverse_order: cli.reverse_order,
    };

    let file = File::open(&cli.filename)
        .map_err(|e| anyhow::anyhow!("Could not open file {}: {}", cli.filename, e))?;
    let matches = read_matches(BufReader::new(file))?;

    if matches.is_empty() {
        eprintln!("No scores read. Nothing to do");
        return Ok(());
    }

    let alignments = chain_matches(&matches, &params)?;

    let stdout = std::io::stdout();
    let mut out = stdout.lock();

    for alignment in &alignments {
        writeln!(
            out,
            ">Alignment #{}  score = {:.1}",
            alignment.id, alignment.score
        )?;
        for entry in &alignment.entries {
            writeln!(
                out,
                "{:3}: {} {:6} {:6} {:7.1} {:7.1}",
                entry.index,
                entry.pair_id,
                entry.x,
                entry.y,
                entry.match_score,
                entry.cumulative_score
            )?;
        }
    }

    Ok(())
}