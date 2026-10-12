//! Port of `chain_length_ratio.pl` — compute the ratio of chain lengths
//! (molecule A span vs. molecule B span) for each line of a DAGchainer
//! alignment/coordinate file.
//!
//! ## Input format
//!
//! Tab-separated: `contigA  lendA  rendA  contigB  lendB  rendB  [...]`
//! — any columns beyond the first six are read but ignored, matching the
//! original Perl's unused `@res`.
//!
//! ## Output format
//!
//! `contigA  lenA  contigB  lenB  ratio`, where `lenX = rendX - lendX + 1`
//! and `ratio = lenA / lenB`.
//!
//! ## Preserved behavior: division by zero is fatal
//!
//! If `lenB == 0`, the original Perl triggers a runtime
//! "Illegal division by zero" error, immediately aborting the entire
//! script. This port preserves that: [`chain_length_ratio`] returns an
//! error (and the CLI aborts) rather than skipping the offending line or
//! emitting `inf`/`NaN`.
//!
//! ## Preserved behavior: numeric formatting matches Perl's default
//!
//! Perl stringifies floating-point numbers using (on virtually all
//! modern builds) `%.15g`-equivalent formatting: 15 significant digits,
//! switching between fixed and scientific notation depending on
//! magnitude, with trailing zeros trimmed. Rust's default `f64`
//! formatting uses a different algorithm (shortest round-trippable
//! representation) and would **not** match Perl's output byte-for-byte
//! (e.g. Rust prints `2.3333333333333335`, Perl prints
//! `2.33333333333333`). [`perl_format_number`] replicates Perl's
//! behavior explicitly. This assumes a standard 64-bit `double` build of
//! Perl (`DBL_DIG == 15`), true for essentially all modern Perl
//! installations.
//!
//! ## Deviation from the original Perl (deliberate improvement)
//!
//! Lines with fewer than 6 tab-separated fields, or non-numeric
//! coordinate fields, are skipped with a warning rather than silently
//! computing with Perl's `undef`-coerced-to-`0` arithmetic (which, under
//! `use warnings`, would also print a runtime warning to stderr but
//! continue with incorrect `0` values).

use std::io::{BufRead, Write};
use std::path::PathBuf;

use clap::Args as ClapArgs;
use thiserror::Error;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Input file(s). If omitted, reads from stdin. Multiple files are
    /// read in sequence (mirroring Perl's `<>` diamond operator).
    #[arg(value_name = "FILE")]
    pub files: Vec<PathBuf>,
}

#[derive(Debug, Error)]
pub enum ChainLengthRatioError {
    #[error("Illegal division by zero: lenB == 0 for line {line:?}")]
    DivisionByZero { line: String },
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

/// One parsed input row's first six required fields; any further columns
/// are intentionally discarded, matching the original Perl's `@res`.
struct Row<'a> {
    contig_a: &'a str,
    lend_a: i64,
    rend_a: i64,
    contig_b: &'a str,
    lend_b: i64,
    rend_b: i64,
}

fn parse_line(line: &str) -> Option<Row<'_>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 6 {
        log::warn!("Skipping line with fewer than 6 tab-separated fields: {line:?}");
        return None;
    }

    let lend_a: i64 = fields[1].parse().ok().or_else(|| {
        log::warn!("Skipping line with non-numeric lendA: {line:?}");
        None
    })?;
    let rend_a: i64 = fields[2].parse().ok().or_else(|| {
        log::warn!("Skipping line with non-numeric rendA: {line:?}");
        None
    })?;
    let lend_b: i64 = fields[4].parse().ok().or_else(|| {
        log::warn!("Skipping line with non-numeric lendB: {line:?}");
        None
    })?;
    let rend_b: i64 = fields[5].parse().ok().or_else(|| {
        log::warn!("Skipping line with non-numeric rendB: {line:?}");
        None
    })?;

    Some(Row {
        contig_a: fields[0],
        lend_a,
        rend_a,
        contig_b: fields[3],
        lend_b,
        rend_b,
    })
}

/// Replicates Perl's default floating-point stringification: equivalent
/// to C's `printf("%.15g", value)`, with trailing zeros trimmed. See
/// module docs for rationale and assumptions.
pub fn perl_format_number(value: f64) -> String {
    const PRECISION: i32 = 15;

    if value == 0.0 {
        return "0".to_string();
    }
    if !value.is_finite() {
        // Not expected to be reached in normal operation (lenB == 0 is
        // caught earlier as a hard error), but handled defensively.
        return if value.is_nan() {
            "NaN".to_string()
        } else if value > 0.0 {
            "Inf".to_string()
        } else {
            "-Inf".to_string()
        };
    }

    let neg = value.is_sign_negative();
    let abs = value.abs();

    // Get a correctly-rounded scientific representation with
    // (PRECISION - 1) digits after the decimal point (PRECISION total
    // significant digits), using Rust's built-in (correctly-rounded)
    // exponential formatter as the source of truth for rounding.
    let sci = format!("{:.*e}", (PRECISION - 1) as usize, abs);
    let (mantissa_str, exp_str) = sci.split_once('e').expect("exponential format always has 'e'");
    let exponent: i32 = exp_str.parse().expect("exponent is always a valid integer");

    let digits: String = mantissa_str.chars().filter(|c| *c != '.').collect();

    let body = if exponent >= -4 && exponent < PRECISION {
        digits_to_fixed(&digits, exponent)
    } else {
        digits_to_sci(&digits, exponent)
    };

    if neg {
        format!("-{body}")
    } else {
        body
    }
}

/// Render a digit string (representing `0.<digits> * 10^(exponent+1)`,
/// i.e. `digits[0]` is the leading significant digit) as fixed-point
/// notation, trimming trailing zeros and a bare trailing decimal point.
fn digits_to_fixed(digits: &str, exponent: i32) -> String {
    let point_pos = exponent + 1; // number of digits before the decimal point
    let mut s = String::new();

    if point_pos <= 0 {
        s.push_str("0.");
        for _ in 0..(-point_pos) {
            s.push('0');
        }
        s.push_str(digits);
    } else if (point_pos as usize) >= digits.len() {
        s.push_str(digits);
        for _ in 0..(point_pos as usize - digits.len()) {
            s.push('0');
        }
    } else {
        let (int_part, frac_part) = digits.split_at(point_pos as usize);
        s.push_str(int_part);
        s.push('.');
        s.push_str(frac_part);
    }

    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    s
}

/// Render a digit string as scientific notation (`d.ddde±NN`), trimming
/// trailing zeros from the mantissa and using a minimum two-digit,
/// explicitly signed exponent — matching C/Perl `%g` conventions.
fn digits_to_sci(digits: &str, exponent: i32) -> String {
    let mut mantissa = String::new();
    mantissa.push_str(&digits[0..1]);
    let rest = digits[1..].trim_end_matches('0');
    if !rest.is_empty() {
        mantissa.push('.');
        mantissa.push_str(rest);
    }
    let sign = if exponent >= 0 { "+" } else { "-" };
    format!("{mantissa}e{sign}{:02}", exponent.abs())
}

/// Core per-line computation: returns the formatted output line, or an
/// error if `lenB == 0` (matching Perl's fatal division-by-zero).
fn compute_ratio_line(row: &Row<'_>, original_line: &str) -> Result<String, ChainLengthRatioError> {
    let len_a = row.rend_a - row.lend_a + 1;
    let len_b = row.rend_b - row.lend_b + 1;

    if len_b == 0 {
        return Err(ChainLengthRatioError::DivisionByZero {
            line: original_line.to_string(),
        });
    }

    let ratio = len_a as f64 / len_b as f64;

    Ok(format!(
        "{}\t{}\t{}\t{}\t{}",
        row.contig_a,
        len_a,
        row.contig_b,
        len_b,
        perl_format_number(ratio)
    ))
}

/// Core logic for a single input stream, generic over any
/// `BufRead`/`Write`.
pub fn chain_length_ratio<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
) -> Result<(), ChainLengthRatioError> {
    for line in reader.lines() {
        let line = line?;
        let Some(row) = parse_line(&line) else {
            continue;
        };
        let out_line = compute_ratio_line(&row, &line)?;
        writeln!(writer, "{out_line}")?;
    }
    Ok(())
}

/// CLI entry point — a drop-in replacement for
/// `chain_length_ratio.pl < file` or `chain_length_ratio.pl file1 file2 ...`.
pub fn run(args: Args) -> anyhow::Result<()> {
    let stdout = std::io::stdout();
    let mut out = stdout.lock();

    if args.files.is_empty() {
        let stdin = std::io::stdin();
        chain_length_ratio(stdin.lock(), &mut out)?;
    } else {
        for path in &args.files {
            let file = std::fs::File::open(path)
                .map_err(|e| anyhow::anyhow!("Could not open file {}: {e}", path.display()))?;
            chain_length_ratio(std::io::BufReader::new(file), &mut out)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_ratio_computation() {
        let input = "chrA\t100\t199\tchrB\t50\t99\n"; // lenA=100, lenB=50
        let mut out = Vec::new();
        chain_length_ratio(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert_eq!(out, "chrA\t100\tchrB\t50\t2\n");
    }

    #[test]
    fn repeating_decimal_matches_perl_formatting() {
        // lenA=7 (e.g. 1..7), lenB=3 (e.g. 1..3) -> ratio = 7/3
        let input = "chrA\t1\t7\tchrB\t1\t3\n";
        let mut out = Vec::new();
        chain_length_ratio(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert_eq!(out, "chrA\t7\tchrB\t3\t2.33333333333333\n");
    }

    #[test]
    fn whole_number_ratio_has_no_trailing_decimal() {
        let input = "chrA\t1\t4\tchrB\t1\t2\n"; // lenA=4, lenB=2, ratio=2
        let mut out = Vec::new();
        chain_length_ratio(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert_eq!(out, "chrA\t4\tchrB\t2\t2\n");
    }

    #[test]
    fn extra_trailing_columns_are_ignored() {
        let input = "chrA\t1\t4\tchrB\t1\t2\tscore=99\tjunk\n";
        let mut out = Vec::new();
        chain_length_ratio(input.as_bytes(), &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert_eq!(out, "chrA\t4\tchrB\t2\t2\n");
    }

    #[test]
    fn division_by_zero_is_fatal() {
        // lenB = rendB - lendB + 1 = 0 when rendB == lendB - 1.
        let input = "chrA\t1\t10\tchrB\t5\t4\n";
        let mut out = Vec::new();
        let result = chain_length_ratio(input.as_bytes(), &mut out);
        assert!(matches!(
            result,
            Err(ChainLengthRatioError::DivisionByZero { .. })
        ));
    }

    #[test]
    fn malformed_lines_are_skipped() {
        let input = "too\tfew\tfields\n";
        let mut out = Vec::new();
        chain_length_ratio(input.as_bytes(), &mut out).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn non_numeric_coordinates_are_skipped() {
        let input = "chrA\tnot_a_number\t199\tchrB\t50\t99\n";
        let mut out = Vec::new();
        chain_length_ratio(input.as_bytes(), &mut out).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn perl_format_number_matches_known_perl_outputs() {
        assert_eq!(perl_format_number(2.0), "2");
        assert_eq!(perl_format_number(7.0 / 3.0), "2.33333333333333");
        assert_eq!(perl_format_number(0.000001), "1e-06");
        assert_eq!(perl_format_number(1.5), "1.5");
        assert_eq!(perl_format_number(-2.5), "-2.5");
        assert_eq!(perl_format_number(0.0), "0");
        assert_eq!(perl_format_number(100.0), "100");
        assert_eq!(perl_format_number(0.1), "0.1");
    }

    #[test]
    fn perl_format_number_handles_very_small_and_large_magnitudes() {
        // exponent < -4 -> scientific notation
        assert_eq!(perl_format_number(0.00001234), "1.234e-05");
        // exponent >= 15 -> scientific notation
        assert_eq!(perl_format_number(1e20), "1e+20");
    }
}