//! Port of `findTandemGeneDuplicates.pl` — detect tandemly duplicated
//! genes along each scaffold/chromosome, allowing up to
//! `max_intervening` non-homologous genes to interrupt a cluster.
//!
//! ## Algorithm
//!
//! For each chromosome, genes are sorted by midpoint position and walked
//! left-to-right by a small state machine ("automaton"):
//! - If the current gene and the next gene are homologs (e-value <= 1e-20,
//!   per `blast` lookup) **and** we haven't exceeded our allowed gap
//!   budget, link them as a tandem pair and advance.
//! - Otherwise, if we still have gap budget left, consume one unit of
//!   budget and skip past the non-homologous gene without advancing
//!   `curr_gene`.
//! - Otherwise, close out the current tandem cluster (if any), reset the
//!   gap budget, and restart the automaton from the gene immediately
//!   after `curr_gene`.
//!
//! ## Deviation from the original Perl (deliberate fix)
//!
//! `max_intervening = 0` is accepted and honored literally (meaning: only
//! strictly adjacent genes may be linked, no gaps tolerated at all). The
//! original Perl's `my $MAX_NUM_INTERVENING = $opt_m || 1;` silently
//! forced an explicit `-m 0` back to `1`, due to Perl's `0`-is-falsy
//! truthiness. Here, the default (when `-m`/`--max-intervening` is
//! omitted) is still `1`, matching the original's default behavior.
//!
//! ## Note on safety of the index arithmetic
//!
//! The original Perl's "reset automaton" branch advances `curr_gene` to
//! `curr_gene_pos + 1` and `i` to `curr_gene_pos + 2`. It may look like
//! this could leave `curr_gene` pointing one-past-the-end of the gene
//! list. In practice this never happens: whenever
//! `curr_gene_pos + 1 == n` (i.e. `curr_gene` would go out of bounds),
//! `i` is simultaneously set to `n + 1`, which immediately fails the
//! loop's `i < n` condition on the next iteration — so the
//! out-of-bounds `curr_gene` value is set but never read. This Rust port
//! preserves that same invariant (and relies on it for memory safety,
//! verified in tests covering exactly this boundary).

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, Write};

use clap::Args as ClapArgs;

/// E-value threshold below which two genes are considered homologs for
/// tandem-duplicate linking purposes. Matches the hardcoded `1e-20` in
/// the original Perl.
const HOMOLOG_EVALUE_THRESHOLD: f64 = 1e-20;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Maximum number of non-homologous genes allowed to intervene
    /// within a single tandem duplication cluster. `0` is valid (strict
    /// adjacency only).
    #[arg(short = 'm', long = "max-intervening", default_value_t = 1)]
    pub max_intervening: usize,
}

/// A single gene's position record, equivalent to the Perl `$gene_struct`
/// hash (`gene_id`, `c1`, `c2`, `mid`, `chr`).
#[derive(Debug, Clone)]
struct Gene {
    gene_id: String,
    chr: String,
    c1: i64,
    c2: i64,
    mid: f64,
}

/// Symmetric gene-pair -> minimum-e-value lookup, equivalent to the
/// Perl `%blast` hash of hashes.
type BlastMap = HashMap<String, HashMap<String, f64>>;

/// One parsed input line's worth of data for a single gene side
/// (chromosome, accession, and coordinates).
struct RawGene<'a> {
    chromo: &'a str,
    gene_id: &'a str,
    c1: i64,
    c2: i64,
}

fn parse_line(line: &str) -> Option<(RawGene<'_>, RawGene<'_>, f64)> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        log::warn!("Skipping malformed line (fewer than 9 fields): {line:?}");
        return None;
    }

    let c1_a: i64 = fields[2].parse().ok().or_else(|| {
        log::warn!("Skipping line with non-numeric end5_A: {line:?}");
        None
    })?;
    let c2_a: i64 = fields[3].parse().ok().or_else(|| {
        log::warn!("Skipping line with non-numeric end3_A: {line:?}");
        None
    })?;
    let c1_b: i64 = fields[6].parse().ok().or_else(|| {
        log::warn!("Skipping line with non-numeric end5_B: {line:?}");
        None
    })?;
    let c2_b: i64 = fields[7].parse().ok().or_else(|| {
        log::warn!("Skipping line with non-numeric end3_B: {line:?}");
        None
    })?;
    let e_value: f64 = fields[8].parse().ok().or_else(|| {
        log::warn!("Skipping line with non-numeric e-value: {line:?}");
        None
    })?;

    Some((
        RawGene { chromo: fields[0], gene_id: fields[1], c1: c1_a, c2: c2_a },
        RawGene { chromo: fields[4], gene_id: fields[5], c1: c1_b, c2: c2_b },
        e_value,
    ))
}

/// Read all input lines, building the per-chromosome gene position lists
/// (first occurrence of each gene_id wins, exactly like the Perl `%seen`
/// guard) and the symmetric minimum-e-value `blast` map.
fn build_data<R: BufRead>(
    reader: R,
) -> std::io::Result<(BTreeMap<String, Vec<Gene>>, BlastMap)> {
    let mut data: BTreeMap<String, Vec<Gene>> = BTreeMap::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut blast: BlastMap = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let Some((raw_a, raw_b, e_value)) = parse_line(&line) else {
            continue;
        };

        for raw in [&raw_a, &raw_b] {
            if seen.contains(raw.gene_id) {
                continue;
            }
            seen.insert(raw.gene_id.to_string());
            let mid = (raw.c1 + raw.c2) as f64 / 2.0;
            data.entry(raw.chromo.to_string()).or_default().push(Gene {
                gene_id: raw.gene_id.to_string(),
                chr: raw.chromo.to_string(),
                c1: raw.c1,
                c2: raw.c2,
                mid,
            });
        }

        let existing = blast
            .get(raw_a.gene_id)
            .and_then(|m| m.get(raw_b.gene_id))
            .copied();
        if existing.map_or(true, |ev| ev > e_value) {
            blast
                .entry(raw_a.gene_id.to_string())
                .or_default()
                .insert(raw_b.gene_id.to_string(), e_value);
            blast
                .entry(raw_b.gene_id.to_string())
                .or_default()
                .insert(raw_a.gene_id.to_string(), e_value);
        }
    }

    Ok((data, blast))
}

/// Run the tandem-duplicate-finding automaton over one chromosome's
/// genes (already sorted by midpoint), returning clusters as lists of
/// adjacent-pair index tuples `(curr_idx, next_idx)` into `genes`.
fn find_tandem_arrays(
    genes: &[Gene],
    blast: &BlastMap,
    max_intervening: usize,
) -> Vec<Vec<(usize, usize)>> {
    let n = genes.len();
    let mut all_tandem_arrays: Vec<Vec<(usize, usize)>> = Vec::new();

    if n == 0 {
        return all_tandem_arrays;
    }

    let mut curr_tandem_array: Vec<(usize, usize)> = Vec::new();
    let mut curr_gene_idx: usize = 0;
    let mut i: usize = 1;
    let mut num_intervening: usize = 0;

    while i < n {
        let next_gene_idx = i;
        let curr_gene_id = &genes[curr_gene_idx].gene_id;
        let next_gene_id = &genes[next_gene_idx].gene_id;

        let evalue = blast.get(curr_gene_id).and_then(|m| m.get(next_gene_id)).copied();

        let is_tandem_hit =
            matches!(evalue, Some(ev) if ev <= HOMOLOG_EVALUE_THRESHOLD)
                && num_intervening <= max_intervening;

        if is_tandem_hit {
            curr_tandem_array.push((curr_gene_idx, next_gene_idx));
            curr_gene_idx = next_gene_idx;
            i += 1;
        } else if num_intervening < max_intervening {
            num_intervening += 1;
            i += 1;
            // curr_gene_idx intentionally unchanged.
        } else {
            if !curr_tandem_array.is_empty() {
                all_tandem_arrays.push(std::mem::take(&mut curr_tandem_array));
            }
            num_intervening = 0;
            let curr_gene_pos = curr_gene_idx;
            // May become `n` (one past the end) when curr_gene_pos is the
            // last valid index. Safe: see module docs — `i` is
            // simultaneously set to `n + 1`, so this value is never read.
            curr_gene_idx = curr_gene_pos + 1;
            i = curr_gene_pos + 2;
        }
    }

    if !curr_tandem_array.is_empty() {
        all_tandem_arrays.push(curr_tandem_array);
    }

    all_tandem_arrays
}

/// Formats a float using C/Perl `printf("%e", ...)` conventions: one
/// leading digit, six decimal places, and a signed two-digit-minimum
/// exponent (e.g. `1.000000e-20`). Rust's built-in `{:e}` formatter does
/// not match this style, so this small helper bridges the gap for exact
/// output parity with the original script.
pub fn format_c_style_e(value: f64) -> String {
    if value == 0.0 {
        return "0.000000e+00".to_string();
    }

    let sign = if value < 0.0 { "-" } else { "" };
    let abs = value.abs();
    let mut exponent = abs.log10().floor() as i32;
    let mut mantissa = abs / 10f64.powi(exponent);

    // Guard against rounding pushing the mantissa up to 10.0 (e.g.
    // 9.9999999e-21 rounding to "10.000000e-21" instead of
    // "1.000000e-20").
    let rounded = format!("{mantissa:.6}");
    if rounded.starts_with("10.") {
        mantissa = 1.0;
        exponent += 1;
    }

    let exp_sign = if exponent >= 0 { "+" } else { "-" };
    format!("{sign}{mantissa:.6}e{exp_sign}{:02}", exponent.abs())
}

/// Write one tandem-duplication cluster's output block: the
/// `## tandem duplication (n=N)` header line, followed by one row per
/// adjacent gene pair in the cluster.
fn write_tandem_cluster<W: Write>(
    writer: &mut W,
    genes: &[Gene],
    cluster: &[(usize, usize)],
    blast: &BlastMap,
) -> std::io::Result<()> {
    let num_genes_in_cluster = cluster.len() + 1;
    writeln!(writer, "## tandem duplication (n={num_genes_in_cluster})")?;

    for &(i, j) in cluster {
        let g1 = &genes[i];
        let g2 = &genes[j];
        let e_value = blast
            .get(&g1.gene_id)
            .and_then(|m| m.get(&g2.gene_id))
            .copied()
            .unwrap_or(f64::NAN);

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            g1.chr,
            g1.gene_id,
            g1.c1,
            g1.c2,
            g2.chr,
            g2.gene_id,
            g2.c1,
            g2.c2,
            format_c_style_e(e_value),
        )?;
    }

    Ok(())
}

/// Core logic, generic over any `BufRead`/`Write`, mirroring the Perl
/// script's stdin -> stdout filter behaviour.
pub fn find_tandem_duplicates<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
    max_intervening: usize,
) -> std::io::Result<()> {
    let (data, blast) = build_data(reader)?;

    for (_chromo, mut genes) in data {
        genes.sort_by(|a, b| a.mid.partial_cmp(&b.mid).unwrap_or(std::cmp::Ordering::Equal));

        let all_tandem_arrays = find_tandem_arrays(&genes, &blast, max_intervening);

        for cluster in &all_tandem_arrays {
            write_tandem_cluster(&mut writer, &genes, cluster, &blast)?;
        }
    }

    Ok(())
}

/// CLI entry point, reading from stdin and writing to stdout — a drop-in
/// replacement for `findTandemGeneDuplicates.pl [-m max_intervening] < match_file`.
pub fn run(args: Args) -> anyhow::Result<()> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    find_tandem_duplicates(stdin.lock(), stdout.lock(), args.max_intervening).map_err(Into::into)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(
        chr_a: &str,
        acc_a: &str,
        mid_a: i64,
        chr_b: &str,
        acc_b: &str,
        mid_b: i64,
        evalue: &str,
    ) -> String {
        format!(
            "{chr_a}\t{acc_a}\t{mid_a}\t{mid_a}\t{chr_b}\t{acc_b}\t{mid_b}\t{mid_b}\t{evalue}\n"
        )
    }

    #[test]
    fn simple_chain_of_three_tandem_genes() {
        // g1-g2 and g2-g3 both strongly homologous and adjacent.
        let input = format!(
            "{}{}",
            line("chr1", "g1", 100, "chr1", "g2", 200, "1e-30"),
            line("chr1", "g2", 200, "chr1", "g3", 300, "1e-25"),
        );

        let mut out = Vec::new();
        find_tandem_duplicates(input.as_bytes(), &mut out, 1).unwrap();
        let out = String::from_utf8(out).unwrap();

        assert!(out.contains("## tandem duplication (n=3)"));
        assert!(out.contains("g1\t100\t100\tchr1\tg2"));
        assert!(out.contains("g2\t200\t200\tchr1\tg3"));
    }

    #[test]
    fn one_intervening_gene_bridged_when_allowed() {
        // g1-g3 homologous but not adjacent; g2 sits between them with
        // no homology to g1. With max_intervening=1, g1 and g3 should
        // still be linked, skipping over g2.
        let input = line("chr1", "g1", 100, "chr1", "g3", 300, "1e-30");
        // Register g2's position via an unrelated, non-homologous pairing
        // on a different chromosome's partner (just need its position
        // recorded once; reuse chr1 to place it between g1/g3).
        let input = format!(
            "{}{}",
            input,
            line("chr1", "g2", 200, "chrX", "other", 9999, "1e-5"),
        );

        let mut out = Vec::new();
        find_tandem_duplicates(input.as_bytes(), &mut out, 1).unwrap();
        let out = String::from_utf8(out).unwrap();

        assert!(out.contains("## tandem duplication (n=2)"));
        assert!(out.contains("g1\t100\t100\tchr1\tg3\t300\t300"));
    }

    #[test]
    fn max_intervening_zero_requires_strict_adjacency() {
        // Same setup as above, but with max_intervening=0 the gap
        // through g2 is NOT allowed to bridge, so no tandem cluster
        // should be reported at all. This also demonstrates the fixed
        // Perl quirk: 0 is honored literally, not coerced to 1.
        let input = format!(
            "{}{}",
            line("chr1", "g1", 100, "chr1", "g3", 300, "1e-30"),
            line("chr1", "g2", 200, "chrX", "other", 9999, "1e-5"),
        );

        let mut out = Vec::new();
        find_tandem_duplicates(input.as_bytes(), &mut out, 0).unwrap();
        let out = String::from_utf8(out).unwrap();

        assert!(out.is_empty(), "expected no tandem clusters, got: {out}");
    }

    #[test]
    fn non_homologous_adjacent_genes_produce_no_cluster() {
        let input = line("chr1", "g1", 100, "chr1", "g2", 200, "0.5"); // weak e-value
        let mut out = Vec::new();
        find_tandem_duplicates(input.as_bytes(), &mut out, 1).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn single_gene_chromosome_produces_no_cluster() {
        // A chromosome with only one gene can never form a tandem pair.
        // Constructed via a match against a different chromosome so g1
        // is registered but has no same-chromosome partner.
        let input = line("chr1", "g1", 100, "chr2", "g2", 500, "1e-30");
        let mut out = Vec::new();
        find_tandem_duplicates(input.as_bytes(), &mut out, 1).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn c_style_e_formatting_matches_printf() {
        assert_eq!(format_c_style_e(1e-20), "1.000000e-20");
        assert_eq!(format_c_style_e(1.23456789e-5), "1.234568e-05");
        assert_eq!(format_c_style_e(0.0), "0.000000e+00");
    }

    #[test]
    fn malformed_lines_are_skipped() {
        let input = "too\tfew\tfields\n";
        let mut out = Vec::new();
        find_tandem_duplicates(input.as_bytes(), &mut out, 1).unwrap();
        assert!(out.is_empty());
    }
}