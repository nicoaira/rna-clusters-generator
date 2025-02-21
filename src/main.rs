// src/main.rs
//
// The main entry point of our RNA dataset generator.
// Usage example:
//    cargo run -- \
//      --num_structures 100 \
//      --seq_min_len 80 \
//      --seq_max_len 120 \
//      --seq_len_distribution norm \
//      --seq_len_mean 100 \
//      --seq_len_sd 30 \
//      --mod_normalization \
//      --normalization_len 20 \
//      --neg_len_variation 20 \
//      --n_stem_indels 2 \
//      --stem_max_size 14 \
//      --stem_min_size 3 \
//      --stem_max_n_modifications 1 \
//      --n_hloop_indels 2 \
//      --hloop_min_size 1 \
//      --hloop_max_size 10 \
//      --hloop_max_n_modifications 1 \
//      --n_iloop_indels 2 \
//      --iloop_min_size 2 \
//      --iloop_max_size 10 \
//      --iloop_max_n_modifications 1 \
//      --n_bulge_indels 2 \
//      --bulge_min_size 1 \
//      --bulge_max_size 8 \
//      --bulge_max_n_modifications 1 \
//      --n_mloop_indels 2 \
//      --mloop_min_size 2 \
//      --mloop_max_size 15 \
//      --mloop_max_n_modifications 1 \
//      --num_workers 1 \
//      --output_dir output_2102 \
//      --split \
//      --train_fraction 0.97 \
//      --val_fraction 0.03 \
//      --plot \
//      --num_plots 50 \
//      --batch_size 32 \
//      --debug
//
// This mirrors your original Python "generate_data.py" usage.

#![recursion_limit = "256"]

mod bulgegraph;        // Our custom bulge graph logic
mod cli;               // Command-line argument parser
mod data_generation;   // Data generation logic

use std::fs::{create_dir_all, File};
use std::path::Path;

use env_logger::Env;
use log::{debug, error, info};

use crate::cli::parse_cli;
use crate::data_generation::{
    build_metadata, generate_all_triplets, save_dataset_csv, split_dataset,
};

fn main() {
    // 1) Parse CLI
    let args = parse_cli();

    // 2) Initialize logging. If `--debug`, set RUST_LOG=debug, else default to error.
    let mut log_level = "error";
    if args.debug {
        log_level = "debug";
    }
    env_logger::Builder::from_env(Env::default().default_filter_or(log_level)).init();
    debug!("Debug mode enabled");

    // 3) Create output directory
    let out_dir = &args.output_dir;
    create_dir_all(out_dir).unwrap_or_else(|_| {
        panic!("Failed to create output directory: {}", out_dir);
    });
    info!("Created or found output directory: {}", out_dir);

    // 4) Build metadata
    let metadata = build_metadata(&args);
    // also save metadata.json
    let metadata_path = Path::new(out_dir).join("metadata.json");
    {
        let f = File::create(&metadata_path).unwrap_or_else(|_| {
            panic!("Failed to create metadata file: {:?}", metadata_path);
        });
        serde_json::to_writer_pretty(f, &metadata).expect("Failed to write metadata JSON");
        info!("Saved metadata to {:?}", metadata_path);
    }

    // 5) Generate triplets
    info!("Starting RNA triplet generation...");
    let mut triplets = generate_all_triplets(&args);
    info!(
        "Generated {} triplets. (Anchor/Positive/Negative)",
        triplets.len()
    );

    // 6) Save full dataset CSV
    let csv_path = Path::new(out_dir).join("rna_triplets.csv");
    {
        let path_str = csv_path.to_string_lossy().to_string();
        save_dataset_csv(&path_str, &triplets, &metadata)
            .unwrap_or_else(|err| error!("Error saving full CSV: {}", err));
        info!("Full dataset saved to {:?}", csv_path);
    }

    // 7) Optionally split dataset
    if args.split {
        let (t_frac, v_frac) = decide_split_fractions(args.train_fraction, args.val_fraction);
        info!("Splitting dataset: train_frac={}, val_frac={}", t_frac, v_frac);
        let (train_data, val_data) = split_dataset(&mut triplets, t_frac);

        let train_csv = Path::new(out_dir).join("rna_triplets_train.csv");
        {
            let path_str = train_csv.to_string_lossy().to_string();
            save_dataset_csv(&path_str, &train_data, &metadata)
                .unwrap_or_else(|err| error!("Error saving train CSV: {}", err));
            info!("Train split saved to {:?}", train_csv);
        }

        let val_csv = Path::new(out_dir).join("rna_triplets_val.csv");
        {
            let path_str = val_csv.to_string_lossy().to_string();
            save_dataset_csv(&path_str, &val_data, &metadata)
                .unwrap_or_else(|err| error!("Error saving val CSV: {}", err));
            info!("Val split saved to {:?}", val_csv);
        }
    }

    // 8) If `--plot`, we do not implement actual plotting in Rust, but we can mention it
    if args.plot {
        info!("Plotting was requested, but this Rust version does not implement it. Ignoring...");
    }

    info!("Data generation complete. Output in: {}", out_dir);
}

/// Decide the (train, val) fractions. If both are provided, ensure they sum to 1. Otherwise fill the missing.
fn decide_split_fractions(
    train_fraction: Option<f64>,
    val_fraction: Option<f64>,
) -> (f64, f64) {
    match (train_fraction, val_fraction) {
        (Some(tf), Some(vf)) => {
            if (tf + vf - 1.0).abs() > 1e-6 {
                error!("train_fraction + val_fraction != 1. Using defaults (0.8, 0.2).");
                (0.8, 0.2)
            } else {
                (tf, vf)
            }
        }
        (Some(tf), None) => {
            let vf = 1.0 - tf;
            (tf, vf)
        }
        (None, Some(vf)) => {
            let tf = 1.0 - vf;
            (tf, vf)
        }
        (None, None) => (0.8, 0.2),
    }
}
