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
    build_metadata, generate_all_clusters, save_dataset_csv,
};

fn main() {
    // 1) Parse CLI
    let args = parse_cli();

    // 2) Initialize logging. If `--debug`, set RUST_LOG=debug, else default to error.
    let log_level = if args.debug { "debug" } else { "error" };
    env_logger::Builder::from_env(Env::default().default_filter_or(log_level)).init();
    debug!("Debug mode enabled");

    // 3) Create output directory
    create_dir_all(&args.output_dir).unwrap_or_else(|_| {
        panic!("Failed to create output directory: {}", &args.output_dir);
    });
    info!("Output directory: {}", &args.output_dir);

    // 4) Build metadata (updated parameters will be captured)
    let metadata = build_metadata(&args);
    let metadata_path = Path::new(&args.output_dir).join("metadata.json");
    {
        let f = File::create(&metadata_path).unwrap_or_else(|_| {
            panic!("Failed to create metadata file: {:?}", metadata_path);
        });
        serde_json::to_writer_pretty(f, &metadata).expect("Failed to write metadata JSON");
        info!("Metadata saved to {:?}", metadata_path);
    }

    // 5) Generate clusters with a customized progress bar
    info!("Starting RNA cluster generation...");
    // Create progress bar based on total number of clusters
    let total_clusters = args.n_clusters;
    let pb = indicatif::ProgressBar::new(total_clusters as u64);
    pb.set_style(indicatif::ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} clusters generated")
        .progress_chars("##-"));
    
    let clusters = generate_all_clusters(&args); // This function updates pb internally
    pb.finish_with_message("Cluster generation complete");
    info!("Generated {} records across clusters", clusters.len());

    // 6) Save CSV dataset
    let csv_path = Path::new(&args.output_dir).join("rna_clusters.csv");
    {
        let path_str = csv_path.to_string_lossy().to_string();
        save_dataset_csv(&path_str, &clusters, &metadata)
            .unwrap_or_else(|err| error!("Error saving CSV: {}", err));
        info!("CSV saved to {:?}", csv_path);
    }

    // 7) (REMOVED) Splitting block removed since splitting is no longer supported.
    // if args.split {
    //     let (t_frac, v_frac) = { ... };
    //     info!("Splitting dataset: train_frac={}, val_frac={}", t_frac, v_frac);
    //     let mut cluster_data = clusters.clone();
    //     let (train_data, val_data) = split_dataset(&mut cluster_data, t_frac);
    //     // Save train and val CSV files ...
    // }

    // 8) Plotting block remains unchanged
    if args.plot {
        use std::process::Command;
        use rand::seq::SliceRandom;
        use std::fs;
        
        info!("Plotting requested, sampling {} records...", args.num_plots);
        let mut rng = rand::thread_rng();
        // Sample from clusters
        let sampled: Vec<_> = clusters
            .choose_multiple(&mut rng, args.num_plots)
            .cloned()
            .collect();

        let plot_dir = Path::new(&args.output_dir).join("plots");
        fs::create_dir_all(&plot_dir).expect("Failed to create plot directory");
        let temp_dir = Path::new(&args.output_dir).join("temp_plots");
        fs::create_dir_all(&temp_dir).expect("Failed to create temporary plot directory");
        let abs_temp_dir = temp_dir.canonicalize().unwrap().to_string_lossy().into_owned();

        let mut batch_script = String::new();
        for record in &sampled {
            let label = if record.anchor == 1 { "anchor" } else { "positive" };
            let block = format!(
r#"rnartist {{
  png {{
    path = "{}"
    width = 500.0
    height = 500.0
  }}
  ss {{
    bn {{
      value = "{}"
      seq = "{}"
      name = "cluster_{}_{}"
    }}
  }}
  theme {{
    details {{
      value = 5
    }}
    scheme {{
      value = "Pumpkin Vegas"
    }}
  }}
}}"#,
                abs_temp_dir,
                record.structure,
                record.seq,
                record.cluster,
                label
            );
            batch_script.push_str(&block);
            batch_script.push_str("\n");
        }

        let batch_script_path = temp_dir.join("batch_plot.kts");
        fs::write(&batch_script_path, batch_script).expect("Failed to write batch plot script");

        let status = Command::new("rnartistcore")
            .arg(batch_script_path.to_string_lossy().as_ref())
            .status()
            .expect("Failed to execute rnartistcore");
        if !status.success() {
            eprintln!("rnartistcore failed on batch plotting script");
        }
    }

    info!("Data generation complete. Output in: {}", &args.output_dir);
}
