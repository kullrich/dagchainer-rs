//! `dagtools` — a multi-call CLI bundling DAGchainer's orchestration
//! pipeline and all accessory utilities (originally separate Perl
//! scripts) behind one binary, `samtools`-style.

use clap::{Parser, Subcommand};

/// DAGchainer orchestration + accessory toolkit.
#[derive(Parser, Debug)]
#[command(name = "dagtools", about = "DAGchainer pipeline and accessory tools")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Run the full chaining pipeline (run_DAG_chainer.pl port): match
    /// file -> .aligncoords chain output.
    Pipeline(dagchainer_rs::pipeline::Args),

    /// Keep only the top N hits per query (limit_to_top_N_hits.pl)
    LimitTopHits(dagchainer_rs::tools::limit_top_hits::Args),

    /// Filter overly repetitive matches (filter_repetitive_matches.pl)
    FilterRepetitive(dagchainer_rs::tools::filter_repetitive::Args),

    /// Identify tandem gene duplicates (findTandemGeneDuplicates.pl)
    FindTandemDups(dagchainer_rs::tools::find_tandem_dups::Args),

    /// Filter chains by minimum length (min_chain_len_filter.pl)
    MinChainLen(dagchainer_rs::tools::min_chain_len_filter::Args),

    /// Convert to gnuplot dotplot format (dagchainer_to_gnuplot.pl)
    ToGnuplot(dagchainer_rs::tools::to_gnuplot::Args),

    /// Convert to chain span coordinates (dagchainer_to_chain_spans.pl)
    ToChainSpans(dagchainer_rs::tools::to_chain_spans::Args),

    /// Compute chain length ratios (chain_to_length_ratio.pl)
    ChainLengthRatio(dagchainer_rs::tools::chain_length_ratio::Args),

    /// Count genes in chains above min length (count_genes_min_chain_len.pl)
    CountGenesMinLen(dagchainer_rs::tools::count_genes_min_len::Args),

    /// Extract chromosome pairs containing duplicates (extract_only_chromo_pairs_with_dups.pl)
    ExtractChromoDups(dagchainer_rs::tools::extract_chromo_dups::Args),

    /// Report distances between chain segments (report_dist_between_chain_segments.pl)
    ReportSegmentDist(dagchainer_rs::tools::report_segment_dist::Args),

    /// Separate combined datasets (separateDataSets.pl)
    SeparateDatasets(dagchainer_rs::tools::separate_datasets::Args),

    /// Sum non-redundant genome content covered by chains (sum_genome_chain_nr_content.pl)
    SumNrContent(dagchainer_rs::tools::sum_nr_content::Args),
}

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Pipeline(args) => dagchainer_rs::pipeline::run(args),
        Commands::LimitTopHits(args) => dagchainer_rs::tools::limit_top_hits::run(args),
        Commands::FilterRepetitive(args) => dagchainer_rs::tools::filter_repetitive::run(args),
        Commands::FindTandemDups(args) => dagchainer_rs::tools::find_tandem_dups::run(args),
        Commands::MinChainLen(args) => dagchainer_rs::tools::min_chain_len_filter::run(args),
        Commands::ToGnuplot(args) => dagchainer_rs::tools::to_gnuplot::run(args),
        Commands::ToChainSpans(args) => dagchainer_rs::tools::to_chain_spans::run(args),
        Commands::ChainLengthRatio(args) => dagchainer_rs::tools::chain_length_ratio::run(args),
        Commands::CountGenesMinLen(args) => dagchainer_rs::tools::count_genes_min_len::run(args),
        Commands::ExtractChromoDups(args) => dagchainer_rs::tools::extract_chromo_dups::run(args),
        Commands::ReportSegmentDist(args) => dagchainer_rs::tools::report_segment_dist::run(args),
        Commands::SeparateDatasets(args) => dagchainer_rs::tools::separate_datasets::run(args),
        Commands::SumNrContent(args) => dagchainer_rs::tools::sum_nr_content::run(args),
    }
}