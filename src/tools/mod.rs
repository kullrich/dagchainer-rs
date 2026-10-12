//! Rust ports of the DAGchainer accessory Perl scripts.
//!
//! Each submodule corresponds 1:1 to an original script from the
//! `dagchainer` distribution, exposed both as a library function
//! (for use from `CRBHits-rs` or other Rust code) and as a `dagtools`
//! CLI subcommand.

pub mod limit_top_hits;
pub mod filter_repetitive;
pub mod find_tandem_dups;
pub mod min_chain_len_filter;
pub mod to_gnuplot;
pub mod to_chain_spans;
pub mod chain_length_ratio;
pub mod count_genes_min_len;
pub mod extract_chromo_dups;
pub mod report_segment_dist;
pub mod separate_datasets;
pub mod sum_nr_content;