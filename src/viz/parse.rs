//! Parsers for DAGchainer-family file formats, replacing
//! `parseMatchFile` / `parseDiagonalFile` from `JsyntenyView.java`.

use super::types::MatchRecord;
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a raw match file: tab-separated
/// `mol1 acc1 end5_1 end3_1 mol2 acc2 end5_2 end3_2 [score]`.
///
/// Equivalent to `JsyntenyView.parseMatchFile`.
pub fn parse_match_file<R: BufRead>(reader: R) -> anyhow::Result<Vec<MatchRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            log::warn!("Skipping malformed line (fewer than 8 fields): {line}");
            continue;
        }

        let score = fields.get(8).and_then(|s| s.parse::<f64>().ok());

        records.push(MatchRecord {
            mol_a: fields[0].to_string(),
            acc_a: fields[1].to_string(),
            end5_a: fields[2].parse()?,
            end3_a: fields[3].parse()?,
            mol_b: fields[4].to_string(),
            acc_b: fields[5].to_string(),
            end5_b: fields[6].parse()?,
            end3_b: fields[7].parse()?,
            score,
            diagonal: None,
        });
    }

    Ok(records)
}

/// Parse DAGchainer's own alignment output (the format produced by
/// `run_DAG_chainer.pl` / the `.aligncoords` file, or raw `dagchainer`
/// stdout) and assign a `diagonal` (chain) ID to each match, by joining
/// back against an existing set of parsed records.
///
/// Equivalent to `JsyntenyView.parseDiagonalFile`, but instead of mutating
/// a shared hashtable in place, this returns a lookup you can apply to
/// your `MatchRecord`s via `acc_a`/`acc_b` pair keys.
pub fn parse_chain_assignments<R: BufRead>(
    reader: R,
) -> anyhow::Result<HashMap<(String, String), usize>> {
    let mut assignments = HashMap::new();
    let mut current_chain = 0usize;
    let mut in_header = false;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.starts_with('>') || trimmed.starts_with("##") {
            current_chain += 1;
            in_header = true;
            continue;
        }
        if trimmed.is_empty() {
            continue;
        }
        in_header = false;

        // Supports both raw dagchainer stdout ("N: pairID x y score cum")
        // and aligncoords format (mol acc end5 end3 mol acc end5 end3 ...).
        let fields: Vec<&str> = trimmed.split_whitespace().collect();

        // Heuristic: aligncoords lines have >= 9 tab/space-separated fields
        // and accession names in positions 1 and 5.
        if fields.len() >= 8 && !trimmed.contains(':') {
            let acc_a = fields[1].to_string();
            let acc_b = fields[5].to_string();
            let key = canonical_acc_key(&acc_a, &acc_b);
            assignments.insert(key, current_chain);
        }
    }

    let _ = in_header;
    Ok(assignments)
}

fn canonical_acc_key(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.to_string(), b.to_string())
    } else {
        (b.to_string(), a.to_string())
    }
}