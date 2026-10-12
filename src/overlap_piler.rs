//! Port of `Overlap_piler.pm` — collapse a set of (possibly unordered,
//! possibly overlapping) coordinate pairs into the minimal set of
//! non-overlapping "piles" (spans), each pile being the overall
//! `(min, max)` extent of every original coordinate pair merged into it.
//!
//! ## Algorithm
//!
//! 1. Normalize each pair so its first element is the smaller
//!    (`lend <= rend`), matching the original's own internal
//!    `sort {$a<=>$b} ($end5, $end3)`.
//! 2. Sort all pairs ascending by `lend`.
//! 3. Walk the sorted list, maintaining a running pile
//!    `(exp_left, exp_right)`. A pair merges into the current pile if
//!    `exp_right >= lend` — i.e. its start falls at or before the
//!    pile's current rightmost extent (including exactly at it: sharing
//!    a boundary coordinate counts as overlap). Otherwise, the current
//!    pile is closed and a new one starts at that pair.
//!
//! This is a direct simplification of the original Perl's merge test
//! `$exp_left <= $rend && $exp_right >= $lend`: once inputs are sorted
//! by `lend`, `exp_left` (the pile's minimum `lend`) is always exactly
//! the `lend` of the pair that started the pile — every later pair's
//! `lend` is `>= exp_left` by sort order, so `min(exp_left, lend)` never
//! changes it. That makes `exp_left <= lend <= rend` unconditionally
//! true for any pair under consideration, leaving `exp_right >= lend` as
//! the only condition that ever actually decides a merge.
//!
//! ## Important: merges only on genuine overlap or shared boundary
//!
//! `(100, 200)` and `(200, 300)` merge (position `200` is shared).
//! `(100, 200)` and `(201, 300)` do **not** merge — they share no
//! coordinate at all, and are returned as two separate piles. There is
//! no "touching"/adjacency-without-overlap merge rule.

/// Collapse a set of coordinate pairs into the minimal set of
/// non-overlapping spans, returned in ascending order by start
/// coordinate.
///
/// Each pair may be given in either order (`(end5, end3)` or
/// `(end3, end5)`); both are normalized internally. An empty input
/// returns an empty result.
pub fn collapse_coordsets(coordsets: &[(i64, i64)]) -> Vec<(i64, i64)> {
    if coordsets.is_empty() {
        return Vec::new();
    }

    let mut normalized: Vec<(i64, i64)> = coordsets
        .iter()
        .map(|&(a, b)| if a <= b { (a, b) } else { (b, a) })
        .collect();

    normalized.sort_unstable_by_key(|&(lend, _)| lend);

    let mut piles: Vec<(i64, i64)> = Vec::new();
    let mut exp_left = normalized[0].0;
    let mut exp_right = normalized[0].1;

    for &(lend, rend) in &normalized[1..] {
        if exp_right >= lend {
            if rend > exp_right {
                exp_right = rend;
            }
            // exp_left never needs updating post-sort (see module docs).
        } else {
            piles.push((exp_left, exp_right));
            exp_left = lend;
            exp_right = rend;
        }
    }
    piles.push((exp_left, exp_right));

    piles
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_interval_returns_itself() {
        assert_eq!(collapse_coordsets(&[(100, 200)]), vec![(100, 200)]);
    }

    #[test]
    fn overlapping_intervals_merge() {
        assert_eq!(collapse_coordsets(&[(100, 200), (150, 300)]), vec![(100, 300)]);
    }

    #[test]
    fn shared_boundary_coordinate_merges() {
        // 200 is shared by both intervals -> must merge.
        assert_eq!(collapse_coordsets(&[(100, 200), (200, 300)]), vec![(100, 300)]);
    }

    #[test]
    fn adjacent_but_non_overlapping_intervals_do_not_merge() {
        // 200 and 201 share no coordinate -> must remain separate piles,
        // matching the real Overlap_piler.pm exactly (no "+1" adjacency
        // merging, confirmed against the actual source).
        assert_eq!(
            collapse_coordsets(&[(100, 200), (201, 300)]),
            vec![(100, 200), (201, 300)]
        );
    }

    #[test]
    fn distant_intervals_remain_separate() {
        assert_eq!(
            collapse_coordsets(&[(100, 200), (500, 600)]),
            vec![(100, 200), (500, 600)]
        );
    }

    #[test]
    fn unsorted_input_is_handled() {
        assert_eq!(
            collapse_coordsets(&[(500, 600), (100, 200), (150, 550)]),
            vec![(100, 600)]
        );
    }

    #[test]
    fn reversed_pairs_are_normalized() {
        assert_eq!(collapse_coordsets(&[(200, 100)]), vec![(100, 200)]);
    }

    #[test]
    fn three_distinct_clusters() {
        let input = vec![(100, 200), (150, 250), (500, 600), (900, 1000), (950, 1100)];
        let result = collapse_coordsets(&input);
        assert_eq!(result, vec![(100, 250), (500, 600), (900, 1100)]);
    }

    #[test]
    fn empty_input_returns_empty() {
        assert_eq!(collapse_coordsets(&[]), Vec::<(i64, i64)>::new());
    }

    #[test]
    fn fully_contained_interval_is_absorbed() {
        assert_eq!(
            collapse_coordsets(&[(100, 200), (120, 180)]),
            vec![(100, 200)]
        );
    }

    #[test]
    fn mixed_merges_and_gaps() {
        let input = vec![(10, 20), (15, 25), (30, 40), (40, 50), (100, 110)];
        let result = collapse_coordsets(&input);
        // (10,20)+(15,25) overlap -> (10,25)
        // (30,40)+(40,50) share boundary 40 -> (30,50)
        // (100,110) stands alone
        assert_eq!(result, vec![(10, 25), (30, 50), (100, 110)]);
    }
}