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
    about = "Generate RNA triplets (anchor/positive/negative) with modifications."
)]
pub struct Cli {
    // ----------------------
    // Sequence Generation
    // ----------------------
    /// Number of structures (anchors) to generate
    #[arg(long, default_value_t = 100)]
    pub num_structures: usize,

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

    /// Maximum length variation for negative structures
    #[arg(long, default_value_t = 0)]
    pub neg_len_variation: isize,

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
    // Appending Parameters
    // ----------------------
    /// Probability that an appending event will occur for a triplet
    #[arg(long, default_value_t = 0.3)]
    pub appending_event_probability: f64,

    /// Probability to append on both sides (otherwise left or right)
    #[arg(long, default_value_t = 0.33)]
    pub both_sides_appending_probability: f64,

    /// Minimum linker length (in bases) for appending event
    #[arg(long, default_value_t = 2)]
    pub linker_min: usize,

    /// Maximum linker length (in bases) for appending event
    #[arg(long, default_value_t = 8)]
    pub linker_max: usize,

    /// Factor to multiply the anchor length for appended RNA length sampling
    #[arg(long, default_value_t = 1.0)]
    pub appending_size_factor: f64,

    // ----------------------
    // Modification Normalization
    // ----------------------
    /// Enable normalization of modification counts based on anchor length
    #[arg(long, default_value_t = false)]
    pub mod_normalization: bool,

    /// Normalization length to scale modifications (default: 100)
    #[arg(long, default_value_t = 100)]
    pub normalization_len: usize,

    // ----------------------
    // Performance & Output
    // ----------------------
    /// Number of parallel workers
    #[arg(long, default_value_t = 4)]
    pub num_workers: usize,

    /// Directory to save output CSV, metadata, etc.
    #[arg(long, default_value_t = String::from("output"))]
    pub output_dir: String,

    /// Number of triplets to generate per worker batch
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
    // Dataset Splitting
    // ----------------------
    /// Enable splitting the dataset into train/val sets
    #[arg(long, default_value_t = false)]
    pub split: bool,

    /// Fraction of data for training
    #[arg(long)]
    pub train_fraction: Option<f64>,

    /// Fraction of data for validation
    #[arg(long)]
    pub val_fraction: Option<f64>,

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
