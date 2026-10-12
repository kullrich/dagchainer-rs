//! # dagchainer-rs
//!
//! Rust port of the core DAG-chaining algorithm from DAGchainer
//! (A. L. Delcher, modified by B. Haas and K. K. Ullrich).
//!
//! Given a set of pairwise "match" points `(pair_id, x, y, score)`, this
//! crate finds chains of co-linear matches (increasing in both `x` and `y`)
//! that maximize cumulative score under an affine gap penalty model. This
//! is the core building block for synteny / collinear-orthology detection
//! used by CRBHits.
//!
//! ## Example
//! ```
//! use dagchainer_rs::{MatchPoint, Params, chain_matches};
//!
//! let matches = vec![
//!     MatchPoint { pair_id: 0, x: 10, y: 12, score: 50.0 },
//!     MatchPoint { pair_id: 1, x: 20, y: 22, score: 50.0 },
//!     MatchPoint { pair_id: 2, x: 30, y: 31, score: 50.0 },
//! ];
//!
//! let params = Params {
//!     indel_score: -1.0,
//!     gap_open_penalty: -1.0,
//!     bp_gap_size: 5,
//!     min_alignment_score: 10.0,
//!     max_dist_between_matches: 100_000,
//!     reverse_order: false,
//! };
//!
//! let alignments = chain_matches(&matches, &params).unwrap();
//! assert_eq!(alignments.len(), 1);
//! assert_eq!(alignments[0].entries.len(), 3);
//! ```

use std::io::BufRead;
use thiserror::Error;

/// A single pairwise match (hit) between two sequences/coordinate systems.
///
/// Mirrors the 4-column input format of the original `dagchainer` binary:
/// `pair_id  x  y  score`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MatchPoint {
    /// Identifier tying this point back to an original match/gene pair.
    pub pair_id: i64,
    /// Coordinate in the first sequence/molecule (e.g. gene midpoint).
    pub x: i64,
    /// Coordinate in the second sequence/molecule.
    pub y: i64,
    /// Match score (e.g. derived from -log10(e-value)).
    pub score: f64,
}

/// Chaining parameters, equivalent to the CLI flags of `dagchainer.cpp`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Params {
    /// Gap extension penalty (`-E`). Must be <= 0.
    pub indel_score: f64,
    /// Gap open penalty (`-O`). Must be <= 0.
    pub gap_open_penalty: f64,
    /// Length (bp) of a single gap (`-G`). Must be > 0.
    pub bp_gap_size: i64,
    /// Minimum alignment (chain) score to report (`-S`). Must be > 0.
    pub min_alignment_score: f64,
    /// Maximum allowed distance between two chained matches (`-D`).
    pub max_dist_between_matches: i64,
    /// Reverse the second coordinate (`-r`), e.g. for minus-strand synteny.
    pub reverse_order: bool,
}

/// One entry (matched pair) within a reported alignment/chain.
#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentEntry {
    /// Position within the chain (0-based, in chain order).
    pub index: usize,
    pub pair_id: i64,
    pub x: i64,
    /// `y`, corrected back to original orientation if `reverse_order` was set.
    pub y: i64,
    pub match_score: f64,
    /// Cumulative chain score up to and including this entry.
    pub cumulative_score: f64,
}

/// A single reported chain of co-linear matches.
#[derive(Debug, Clone, PartialEq)]
pub struct Alignment {
    /// 1-based alignment number, in the order they were discovered.
    pub id: usize,
    /// Total chain score.
    pub score: f64,
    pub entries: Vec<AlignmentEntry>,
}

#[derive(Debug, Error)]
pub enum DagChainerError {
    #[error("gap open penalty must be <= 0, got {0}")]
    InvalidGapOpenPenalty(f64),
    #[error("gap extension (indel) penalty must be <= 0, got {0}")]
    InvalidIndelScore(f64),
    #[error("bp_gap_size must be > 0, got {0}")]
    InvalidBpGapSize(i64),
    #[error("min_alignment_score must be > 0, got {0}")]
    InvalidMinAlignmentScore(f64),
    #[error("max_dist_between_matches must be >= 0, got {0}")]
    InvalidMaxDist(i64),
    #[error("failed to parse match line {line_no}: {line:?}")]
    ParseError { line_no: usize, line: String },
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

fn validate_params(p: &Params) -> Result<(), DagChainerError> {
    if p.gap_open_penalty > 0.0 {
        return Err(DagChainerError::InvalidGapOpenPenalty(p.gap_open_penalty));
    }
    if p.indel_score > 0.0 {
        return Err(DagChainerError::InvalidIndelScore(p.indel_score));
    }
    if p.bp_gap_size <= 0 {
        return Err(DagChainerError::InvalidBpGapSize(p.bp_gap_size));
    }
    if p.min_alignment_score <= 0.0 {
        return Err(DagChainerError::InvalidMinAlignmentScore(
            p.min_alignment_score,
        ));
    }
    if p.max_dist_between_matches < 0 {
        return Err(DagChainerError::InvalidMaxDist(p.max_dist_between_matches));
    }
    Ok(())
}

/// Parse match points from any `BufRead` source, one per line, whitespace
/// separated: `pair_id x y score`. Mirrors the original `fscanf("%d %d %d %f")`.
pub fn read_matches<R: BufRead>(reader: R) -> Result<Vec<MatchPoint>, DagChainerError> {
    let mut out = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let mut fields = trimmed.split_whitespace();
        let parse_err = || DagChainerError::ParseError {
            line_no: line_no + 1,
            line: line.clone(),
        };
        let pair_id: i64 = fields.next().ok_or_else(parse_err)?.parse().map_err(|_| parse_err())?;
        let x: i64 = fields.next().ok_or_else(parse_err)?.parse().map_err(|_| parse_err())?;
        let y: i64 = fields.next().ok_or_else(parse_err)?.parse().map_err(|_| parse_err())?;
        let score: f64 = fields.next().ok_or_else(parse_err)?.parse().map_err(|_| parse_err())?;
        out.push(MatchPoint { pair_id, x, y, score });
    }
    Ok(out)
}

/// Internal marker for the DP predecessor pointer, matching the original
/// `from` array semantics: `-1` = root (no predecessor), `-2` = consumed by
/// a previously emitted chain, `>=0` = index of predecessor.
#[derive(Debug, Clone, Copy, PartialEq)]
enum From {
    Root,
    Used,
    Index(usize),
}

/// Internal working representation (post reverse-transform, pre-sort).
#[derive(Debug, Clone, Copy)]
struct Node {
    pair_id: i64,
    x: i64,
    y: i64,
    /// Original (pre reverse-transform) y, used for reporting.
    display_y: i64,
    score: f64,
}

/// Run the DAG-chaining DP and return all reported alignments (chains),
/// in discovery order, matching the behaviour of `dagchainer.cpp`'s
/// `Print_Chains`.
///
/// This repeatedly:
/// 1. Runs an O(n^2) DP over the (sorted, deduplicated) match points to find
///    best cumulative chain scores under an affine gap penalty.
/// 2. Greedily emits chains in descending score order, as long as they don't
///    reuse an already-consumed match point.
/// 3. If a conflict is found (two high-scoring chains want to share a node),
///    stops early, removes all consumed points, and restarts the DP on the
///    remainder — exactly mirroring the original's `done` loop.
pub fn chain_matches(
    matches_in: &[MatchPoint],
    params: &Params,
) -> Result<Vec<Alignment>, DagChainerError> {
    validate_params(params)?;

    if matches_in.is_empty() {
        return Ok(Vec::new());
    }

    let max_y = matches_in.iter().map(|m| m.y).max().unwrap_or(0);

    let mut nodes: Vec<Node> = matches_in
        .iter()
        .map(|m| {
            let y = if params.reverse_order {
                max_y - m.y + 1
            } else {
                m.y
            };
            Node {
                pair_id: m.pair_id,
                x: m.x,
                y,
                display_y: m.y,
                score: m.score,
            }
        })
        .collect();

    nodes.sort_by(|a, b| a.x.cmp(&b.x).then(a.y.cmp(&b.y)));

    // Deduplicate: keep the first of any (x, y) duplicate, warn on the rest.
    let mut deduped: Vec<Node> = Vec::with_capacity(nodes.len());
    for node in nodes {
        match deduped.last() {
            Some(prev) if prev.x == node.x && prev.y == node.y => {
                log::warn!(
                    "Duplicate score entries: ({},{},{:.1}) and ({},{},{:.1}) - discarding latter",
                    prev.x, prev.y, prev.score, node.x, node.y, node.score
                );
            }
            _ => deduped.push(node),
        }
    }

    let mut alignments = Vec::new();
    let mut working = deduped;

    loop {
        let n = working.len();
        if n == 0 {
            break;
        }

        let mut path_score: Vec<f64> = working.iter().map(|m| m.score).collect();
        let mut from: Vec<From> = vec![From::Root; n];

        for j in 1..n {
            for i in (0..j).rev() {
                let del_x = working[j].x - working[i].x - 1;
                let del_y = working[j].y - working[i].y - 1;

                if del_x >= 0 && del_y >= 0 {
                    if del_x > params.max_dist_between_matches
                        && del_y > params.max_dist_between_matches
                    {
                        break;
                    }
                    if del_x > params.max_dist_between_matches
                        || del_y > params.max_dist_between_matches
                    {
                        continue;
                    }

                    // IMPORTANT: this replicates a quirk in the original dagchainer.cpp,
                    // not true round-half-up rounding. The C++ source computes:
                    //   (int)( ((del_x+del_y)+abs(del_x-del_y)) / (2*BP_GAP_SIZE) + 0.5 )
                    // Since del_x, del_y, and BP_GAP_SIZE are all `int` in C++, the
                    // division `/` is performed as INTEGER (truncating) division *before*
                    // the `+ 0.5` is ever applied in floating point. Because the division
                    // result is therefore already a whole number, adding 0.5 and
                    // truncating via `(int)` again is a no-op: the final value is simply
                    // `max(del_x, del_y) / BP_GAP_SIZE` using truncating integer division,
                    // with NO rounding to nearest ever actually occurring. This is almost
                    // certainly an unintentional bug in the original (the author likely
                    // intended float division with round-to-nearest), but it is preserved
                    // here deliberately for exact numerical compatibility with the
                    // original dagchainer binary and all existing literature/pipelines
                    // built on its output.
                    let max_del = del_x.max(del_y);
                    let num_gaps = max_del / params.bp_gap_size;

                    let mut x_score = path_score[i] + working[j].score;
                    if num_gaps > 0 {
                        x_score += params.gap_open_penalty + num_gaps as f64 * params.indel_score;
                    }

                    if x_score > path_score[j] {
                        path_score[j] = x_score;
                        from[j] = From::Index(i);
                    }
                }
            }
        }

        struct HighEntry {
            score: f64,
            rc: i64,
            sub: usize,
        }

        let mut high: Vec<HighEntry> = (0..n)
            .filter(|&i| path_score[i] >= params.min_alignment_score)
            .map(|i| HighEntry {
                score: path_score[i],
                rc: working[i].x + working[i].y,
                sub: i,
            })
            .collect();

        high.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(b.rc.cmp(&a.rc))
        });

        let mut done = true;

        for h in &high {
            if from[h.sub] == From::Used {
                continue;
            }

            let mut ans = Vec::new();
            let mut j = h.sub;
            loop {
                match from[j] {
                    From::Index(prev) => {
                        ans.push(j);
                        j = prev;
                    }
                    _ => {
                        ans.push(j);
                        break;
                    }
                }
            }

            if from[j] == From::Used {
                // Conflict: this chain's root was already consumed by a
                // previous chain in *this* pass. Stop, purge used nodes,
                // and restart the DP (matches original algorithm).
                done = false;
                break;
            }

            ans.reverse();
            let mut entries = Vec::with_capacity(ans.len());
            for (idx, &node_idx) in ans.iter().enumerate() {
                from[node_idx] = From::Used;
                let node = &working[node_idx];
                entries.push(AlignmentEntry {
                    index: idx,
                    pair_id: node.pair_id,
                    x: node.x,
                    y: node.display_y,
                    match_score: node.score,
                    cumulative_score: path_score[node_idx],
                });
            }

            alignments.push(Alignment {
                id: alignments.len() + 1,
                score: path_score[h.sub],
                entries,
            });
        }

        if done {
            break;
        }

        // Remove all consumed ("Used") nodes and redo the DP on the rest.
        let mut retained = Vec::with_capacity(n);
        for (i, node) in working.into_iter().enumerate() {
            if from[i] != From::Used {
                retained.push(node);
            }
        }
        working = retained;
    }

    Ok(alignments)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> Params {
        Params {
            indel_score: -1.0,
            gap_open_penalty: -1.0,
            bp_gap_size: 5,
            min_alignment_score: 10.0,
            max_dist_between_matches: 100_000,
            reverse_order: false,
        }
    }

    #[test]
    fn simple_colinear_chain() {
        let matches = vec![
            MatchPoint { pair_id: 0, x: 10, y: 12, score: 50.0 },
            MatchPoint { pair_id: 1, x: 20, y: 22, score: 50.0 },
            MatchPoint { pair_id: 2, x: 30, y: 31, score: 50.0 },
        ];
        let alignments = chain_matches(&matches, &params()).unwrap();
        assert_eq!(alignments.len(), 1);
        assert_eq!(alignments[0].entries.len(), 3);
        assert_eq!(alignments[0].entries[0].pair_id, 0);
        assert_eq!(alignments[0].entries[2].pair_id, 2);
    }

    #[test]
    fn duplicate_points_are_deduped() {
        let matches = vec![
            MatchPoint { pair_id: 0, x: 1, y: 1, score: 50.0 },
            MatchPoint { pair_id: 1, x: 1, y: 1, score: 99.0 },
        ];
        // Should not panic, and should produce a single chain of length 1
        // (the first-seen point after sorting is kept).
        let alignments = chain_matches(&matches, &params()).unwrap();
        assert_eq!(alignments.len(), 1);
        assert_eq!(alignments[0].entries.len(), 1);
    }

    #[test]
    fn invalid_params_rejected() {
        let mut p = params();
        p.gap_open_penalty = 1.0;
        assert!(chain_matches(&[], &p).is_err());
    }

    #[test]
    fn reverse_order_restores_original_y() {
        let matches = vec![
            MatchPoint { pair_id: 0, x: 1, y: 5, score: 50.0 },
            MatchPoint { pair_id: 1, x: 2, y: 3, score: 50.0 },
        ];
        let mut p = params();
        p.reverse_order = true;
        let alignments = chain_matches(&matches, &p).unwrap();
        // Regardless of internal transform, reported y must be original.
        let ys: Vec<i64> = alignments
            .iter()
            .flat_map(|a| a.entries.iter().map(|e| e.y))
            .collect();
        assert!(ys.contains(&5));
        assert!(ys.contains(&3));
    }

    #[test]
    fn num_gaps_uses_truncating_integer_division_not_rounding() {
        // Verified against real dagchainer C++ output: with bp_gap_size=10000,
        // a max_del of 6759 must yield num_gaps=0 (truncating division),
        // NOT num_gaps=1 (which true round-half-up would give, since
        // 6759/10000 + 0.5 = 1.1759 -> floor = 1). See module docs on this
        // deliberately-preserved C++ integer-division quirk.
        let matches = vec![
            MatchPoint { pair_id: 0, x: 0, y: 0, score: 50.0 },
            MatchPoint { pair_id: 1, x: 6760, y: 6760, score: 50.0 }, // del_x=del_y=6759
        ];
        let params = Params {
            indel_score: -3.0,
            gap_open_penalty: 0.0,
            bp_gap_size: 10_000,
            min_alignment_score: 1.0,
            max_dist_between_matches: 100_000,
            reverse_order: false,
        };
        let alignments = chain_matches(&matches, &params).unwrap();
        assert_eq!(alignments.len(), 1);
        // No gap penalty applied (num_gaps=0): cumulative score = 50+50 = 100,
        // NOT 97 (which would reflect one incorrectly-applied -3 penalty).
        assert_eq!(alignments[0].score, 100.0);
    }
}

pub mod tools;      // accessory utilities (ported Perl scripts)
pub mod overlap_piler;
pub mod pipeline;