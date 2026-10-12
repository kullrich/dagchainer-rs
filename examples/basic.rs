//! Minimal example showing how to use `dagchainer-rs` as a library,
//! without going through the CLI binary.
//!
//! Run with:
//! ```bash
//! cargo run --example basic
//! ```

use dagchainer_rs::{chain_matches, MatchPoint, Params};

fn main() -> anyhow::Result<()> {
    // Three colinear matches (increasing in both x and y) that should
    // chain together, plus one distant match that should form its own
    // separate alignment.
    let matches = vec![
        MatchPoint { pair_id: 0, x: 10, y: 12, score: 50.0 },
        MatchPoint { pair_id: 1, x: 20, y: 22, score: 50.0 },
        MatchPoint { pair_id: 2, x: 30, y: 31, score: 50.0 },
        MatchPoint { pair_id: 3, x: 5000, y: 5002, score: 50.0 },
    ];

    let params = Params {
        indel_score: -1.0,
        gap_open_penalty: -1.0,
        bp_gap_size: 5,
        min_alignment_score: 10.0,
        max_dist_between_matches: 100_000,
        reverse_order: false,
    };

    let alignments = chain_matches(&matches, &params)?;

    println!("Found {} alignment(s):\n", alignments.len());

    for alignment in &alignments {
        println!(">Alignment #{}  score = {:.1}", alignment.id, alignment.score);
        for entry in &alignment.entries {
            println!(
                "{:3}: {} {:6} {:6} {:7.1} {:7.1}",
                entry.index,
                entry.pair_id,
                entry.x,
                entry.y,
                entry.match_score,
                entry.cumulative_score
            );
        }
        println!();
    }

    Ok(())
}