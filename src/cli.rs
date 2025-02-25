// src/cli.rs

use clap::{Parser, ValueEnum};

/// Distribution options for sequence length.
#[derive(ValueEnum, Clone, Debug)]
pub enum SeqLenDistribution {
    /// Normal (Gaussian) distribution
    Norm,
    /// Uniform distribution
    Unif,
}

/// RNA Structure Triplet Dataset Generator Command-Line Arguments
#[derive(Parser, Debug)]
#[command(
    name = "rna_generator",
    version,
    about = "Generate RNA clusters (one anchor per cluster with multiple positives)."
)]
pub struct Cli {
    // ----------------------
    // Cluster Generation
    // ----------------------
    /// Number of clusters (each cluster has one anchor structure)
    #[arg(long)]
    pub n_clusters: usize,

    /// Number of positive structures per cluster
    #[arg(long)]
    pub structures_per_cluster: usize,

    // ----------------------
    // Sequence Generation
    // ----------------------
    /// Minimum sequence length
    #[arg(long, default_value_t = 50)]
    pub seq_min_len: usize,

    /// Maximum sequence length
    #[arg(long, default_value_t = 100)]
    pub seq_max_len: usize,

    /// Distribution of sequence lengths
    #[arg(long, value_enum, default_value_t = SeqLenDistribution::Unif)]
    pub seq_len_distribution: SeqLenDistribution,

    /// Mean for normal distribution of sequence length
    #[arg(long, default_value_t = 75.0)]
    pub seq_len_mean: f64,

    /// Standard deviation for normal distribution of sequence length
    #[arg(long, default_value_t = 10.0)]
    pub seq_len_sd: f64,

    // ----------------------
    // Stem Modifications
    // ----------------------
    /// Number of stem modification cycles
    #[arg(long, default_value_t = 1)]
    pub n_stem_indels: usize,

    /// Minimum stem size
    #[arg(long, default_value_t = 2)]
    pub stem_min_size: usize,

    /// Maximum stem size
    #[arg(long, default_value_t = 10)]
    pub stem_max_size: usize,

    /// Maximum modifications per stem
    #[arg(long, default_value_t = 1)]
    pub stem_max_n_modifications: usize,

    // ----------------------
    // Loop Modifications
    // ----------------------
    /// Number of hairpin loop modification cycles
    #[arg(long, default_value_t = 1)]
    pub n_hloop_indels: usize,

    /// Number of internal loop modification cycles
    #[arg(long, default_value_t = 1)]
    pub n_iloop_indels: usize,

    /// Number of bulge modification cycles
    #[arg(long, default_value_t = 1)]
    pub n_bulge_indels: usize,

    /// Number of multi loop modification cycles
    #[arg(long, default_value_t = 1)]
    pub n_mloop_indels: usize,

    // ----------------------
    // Loop Size Constraints
    // ----------------------
    /// Minimum hairpin loop size
    #[arg(long, default_value_t = 3)]
    pub hloop_min_size: usize,

    /// Maximum hairpin loop size
    #[arg(long, default_value_t = 10)]
    pub hloop_max_size: usize,

    /// Minimum internal loop size
    #[arg(long, default_value_t = 2)]
    pub iloop_min_size: usize,

    /// Maximum internal loop size
    #[arg(long, default_value_t = 10)]
    pub iloop_max_size: usize,

    /// Minimum bulge loop size
    #[arg(long, default_value_t = 1)]
    pub bulge_min_size: usize,

    /// Maximum bulge loop size
    #[arg(long, default_value_t = 1)]
    pub bulge_max_size: usize,

    /// Minimum multi loop size
    #[arg(long, default_value_t = 2)]
    pub mloop_min_size: usize,

    /// Maximum multi loop size
    #[arg(long, default_value_t = 15)]
    pub mloop_max_size: usize,

    // ----------------------
    // Loop Modification Limits
    // ----------------------
    /// Maximum modifications per hairpin loop
    #[arg(long, default_value_t = 1)]
    pub hloop_max_n_modifications: usize,

    /// Maximum modifications per internal loop
    #[arg(long, default_value_t = 1)]
    pub iloop_max_n_modifications: usize,

    /// Maximum modifications per bulge loop
    #[arg(long, default_value_t = 1)]
    pub bulge_max_n_modifications: usize,

    /// Maximum modifications per multi loop
    #[arg(long, default_value_t = 1)]
    pub mloop_max_n_modifications: usize,

    // ----------------------
    // Performance & Output
    // ----------------------
    /// Number of parallel workers
    #[arg(long, default_value_t = 4)]
    pub num_workers: usize,

    /// Directory to save output CSV, metadata, etc.
    #[arg(long, default_value_t = String::from("output"))]
    pub output_dir: String,

    /// Number of triplets per worker batch (now used for clusters)
    #[arg(long, default_value_t = 64)]
    pub batch_size: usize,

    // ----------------------
    // Visualization (optional - will ignore in final Rust)
    // ----------------------
    /// Generate structure plots (currently not implemented in Rust)
    #[arg(long, default_value_t = false)]
    pub plot: bool,

    /// Number of triplets to plot (no-op in Rust version)
    #[arg(long, default_value_t = 5)]
    pub num_plots: usize,

    // ----------------------
    // Misc / Debug
    // ----------------------
    /// Enable debug logging
    #[arg(long, default_value_t = false)]
    pub debug: bool,

    /// Enable detailed timing logs (optional use in Rust)
    #[arg(long, default_value_t = false)]
    pub timing_log: bool,
}

/// Helper to parse arguments from CLI.
pub fn parse_cli() -> Cli {
    Cli::parse()
}
