//! Interactive HTML dot-plot generation — a modern, GUI-toolkit-free
//! replacement for `JsyntenyView.java`.
//!
//! Instead of opening Swing windows (one `XYplot` per chromosome pair),
//! this writes one self-contained, zoomable/pannable/hoverable HTML file
//! per molecule pair (or a single multi-panel HTML, your choice), viewable
//! in any browser — no installation, no GUI framework, cross-platform.

use super::types::MatchRecord;
use plotly::common::{Marker, Mode};
use plotly::layout::Axis;
use plotly::{Layout, Plot, Scatter};
use std::collections::HashMap;
use std::path::Path;

/// A small, readable color palette for distinguishing chains ("diagonals"),
/// replacing the Java original's `getRandomColor()` — deterministic and
/// colorblind-friendlier than pure random RGB.
const PALETTE: &[&str] = &[
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
];

fn color_for_diagonal(diag: Option<usize>) -> &'static str {
    match diag {
        Some(d) => PALETTE[d % PALETTE.len()],
        None => "#cccccc", // unassigned matches, shown in grey
    }
}

/// Group matches by molecule pair and write one interactive HTML dot plot
/// per pair into `out_dir`.
pub fn write_interactive_plots(
    records: &[MatchRecord],
    out_dir: &Path,
) -> anyhow::Result<Vec<std::path::PathBuf>> {
    std::fs::create_dir_all(out_dir)?;

    let mut grouped: HashMap<(String, String), Vec<&MatchRecord>> = HashMap::new();
    for r in records {
        grouped.entry(r.molecule_pair_key()).or_default().push(r);
    }

    let mut written = Vec::new();

    for ((mol_a, mol_b), matches) in grouped {
        let mut plot = Plot::new();

        // Group by diagonal so each chain gets its own legend entry/color,
        // mirroring the Java original's per-diagonal coloring.
        let mut by_diag: HashMap<Option<usize>, (Vec<f64>, Vec<f64>, Vec<String>)> =
            HashMap::new();

        for m in &matches {
            let entry = by_diag.entry(m.diagonal).or_insert_with(|| {
                (Vec::new(), Vec::new(), Vec::new())
            });
            entry.0.push(m.mid_a());
            entry.1.push(m.mid_b());
            entry.2.push(format!("{} vs {}", m.acc_a, m.acc_b)); // hover text
        }

        for (diag, (xs, ys, hover)) in by_diag {
            let name = match diag {
                Some(d) => format!("chain {d}"),
                None => "unassigned".to_string(),
            };
            let trace = Scatter::new(xs, ys)
                .mode(Mode::Markers)
                .name(&name)
                .text_array(hover)
                .marker(Marker::new().color(color_for_diagonal(diag)).size(6));
            plot.add_trace(trace);
        }

        let layout = Layout::new()
            .title(format!("{mol_a} vs {mol_b}").as_str())
            .x_axis(Axis::new().title(format!("chr {mol_b}").as_str()))
            .y_axis(Axis::new().title(format!("chr {mol_a}").as_str()));
        plot.set_layout(layout);

        let filename = out_dir.join(format!("{mol_a}_vs_{mol_b}.html"));
        plot.write_html(&filename);
        written.push(filename);
    }

    Ok(written)
}