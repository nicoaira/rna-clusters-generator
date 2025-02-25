// src/data_generation.rs
//
// This module replicates the core Python logic: random RNA, folding via RNAfold,
// identifying stems/loops, applying modifications, generating negative sample,
// returning the anchor/positive/negative triplets. It also includes routines
// for dataset splitting and CSV/JSON output.
//

use std::collections::HashMap;
use std::fs::File;
use std::io::{Write, BufRead, BufReader};
use std::process::{Command, Stdio};

use rand::distributions::{Distribution, Uniform};
use rand::{thread_rng, Rng};
use rand_distr::Normal;
use serde::Serialize;
use rand::seq::SliceRandom;


use crate::bulgegraph::BulgeGraph;
use crate::cli::{Cli, SeqLenDistribution};

/// New record type for cluster datasets (no negative sample)
#[derive(Debug, Serialize, Clone)]
pub struct RnaClusterRecord {
    pub cluster: usize,
    pub seq: String,
    pub structure: String,
    pub anchor: u8, // 1 for anchor, 0 for positive
}

/// Generate a random RNA sequence of the given length (uniform distribution of {A,C,G,U}).
pub fn generate_random_rna(length: usize) -> String {
    let mut rng = thread_rng();
    let alphabet = ['A', 'C', 'G', 'U'];
    let dist = Uniform::from(0..alphabet.len());
    (0..length)
        .map(|_| {
            let idx = dist.sample(&mut rng);
            alphabet[idx]
        })
        .collect()
}

/// Fold an RNA sequence using `RNAfold` in a child process, capturing the MFE structure.
/// We only return the dot-bracket; we do not parse free energy from the output.
#[allow(dead_code)]
pub fn fold_rna(seq: &str) -> Result<String, String> {
    // We'll run: echo "ACG..." | RNAfold --noPS
    // and parse the second line from stdout for the structure
    let mut child = Command::new("RNAfold")
        .arg("--noPS") // do not generate PostScript
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Could not spawn RNAfold: {}", e))?;

    {
        // Write the sequence into RNAfold's stdin
        let mut stdin = child.stdin.take().ok_or("Failed to open stdin")?;
        std::io::Write::write_all(&mut stdin, seq.as_bytes())
            .map_err(|e| format!("Error writing to RNAfold stdin: {}", e))?;
        // we can also append a newline if needed
        std::io::Write::write_all(&mut stdin, b"\n")
            .map_err(|e| format!("Error writing newline to RNAfold stdin: {}", e))?;
    }

    // Capture and read
    let output = child
        .wait_with_output()
        .map_err(|e| format!("Error waiting for RNAfold: {}", e))?;

    if !output.status.success() {
        return Err(format!("RNAfold exited with status {}", output.status));
    }

    // Example RNAfold output (2 lines):
    // ACGUACGU
    // ....()(). (-0.70)
    let stdout_str = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout_str.trim().lines().collect();
    if lines.len() < 2 {
        return Err(format!(
            "Unexpected RNAfold output (need at least 2 lines):\n{}",
            stdout_str
        ));
    }
    // the second line typically has something like "....()(). (-0.70)"
    // we'll parse everything up to the first space
    // e.g. "....()()." => structure
    let line2 = lines[1].trim();
    let space_idx = line2
        .find(' ')
        .unwrap_or_else(|| line2.len().saturating_sub(1));
    let structure = &line2[..space_idx];
    Ok(structure.to_string())
}

/// Helper: choose an integer length from either a uniform or normal distribution
pub fn choose_seq_length(
    distribution: &SeqLenDistribution,
    min_len: usize,
    max_len: usize,
    mean: f64,
    sd: f64,
) -> usize {
    match distribution {
        SeqLenDistribution::Unif => {
            let mut rng = thread_rng();
            rng.gen_range(min_len..=max_len)
        }
        SeqLenDistribution::Norm => {
            let normal_dist = Normal::new(mean, sd).unwrap();
            let mut rng = thread_rng();
            loop {
                // sample until we get within [min_len, max_len]
                let sample = normal_dist.sample(&mut rng).round();
                if sample >= min_len as f64 && sample <= max_len as f64 {
                    return sample as usize;
                }
            }
        }
    }
}

// ============ Negative sample via dinuc shuffle ============

/// Return a new sequence that preserves the original dinucleotide frequencies.
#[allow(dead_code)]
pub fn dinuc_shuffle(seq: &str) -> String {
    if seq.len() <= 1 {
        return seq.to_string();
    }
    // We'll build a graph of adjacency: for each base, list possible next bases
    let mut graph: HashMap<char, Vec<char>> = HashMap::new();
    let chars: Vec<char> = seq.chars().collect();
    for i in 0..(chars.len() - 1) {
        let src = chars[i];
        let dst = chars[i + 1];
        graph.entry(src).or_default().push(dst);
    }

    // Shuffle each adjacency list
    let mut rng = thread_rng();
    for val in graph.values_mut() {
        val.shuffle(&mut rng);
    }

    // Construct an Eulerian-like path from the first character
    let mut stack = vec![chars[0]];
    let mut trail = Vec::new();
    while let Some(&current) = stack.last() {
        if let Some(targets) = graph.get_mut(&current) {
            if !targets.is_empty() {
                let next_node = targets.pop().unwrap();
                stack.push(next_node);
            } else {
                trail.push(stack.pop().unwrap());
            }
        } else {
            // no adjacency
            trail.push(stack.pop().unwrap());
        }
    }
    trail.reverse();
    if trail.len() != seq.len() {
        // fallback: return original
        return seq.to_string();
    }
    trail.into_iter().collect()
}

/// Generate a negative sample by dinuc shuffling and optional length variation.
#[allow(dead_code)]
pub fn generate_negative_sample(seq: &str, allowed_variation: isize) -> (String, String) {
    let mut rng = thread_rng();
    let mut new_seq = dinuc_shuffle(seq);

    // optionally vary length
    if allowed_variation != 0 {
        let variation = rng.gen_range(-allowed_variation..=allowed_variation);
        if variation > 0 {
            // lengthen
            for _ in 0..variation {
                let pos = rng.gen_range(0..=new_seq.len());
                let base = *['A', 'C', 'G', 'U'].choose(&mut rng).unwrap();
                new_seq.insert(pos, base);
            }
        } else if variation < 0 {
            // shorten
            let absvar = variation.abs() as usize;
            if absvar < new_seq.len() {
                for _ in 0..absvar {
                    let pos = rng.gen_range(0..new_seq.len());
                    new_seq.remove(pos);
                }
            }
        }
    }

    // Return placeholder structure; will update it in batch.
    let placeholder = ".".repeat(new_seq.len());
    (new_seq, placeholder)
}

// ============ Node identification ============

/// Quick helper to classify a node by name:
///   - "s" => Stem
///   - "h" => Hairpin
///   - "i" => Internal
///   - "m" => Multi
///   - "b" => We usually interpret as a bulge if it has length=1, else internal
///   - "f0" => 5' unpaired
///   - "t0" => 3' unpaired
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeType {
    Stem,
    Hairpin,
    Internal,
    Multi,
    Bulge,
    FivePrime,
    ThreePrime,
    Unknown,
}

fn classify_node(node_name: &str, coords: &[usize]) -> NodeType {
    if node_name.starts_with('s') {
        NodeType::Stem
    } else if node_name.starts_with('h') {
        NodeType::Hairpin
    } else if node_name.starts_with('i') {
        NodeType::Internal
    } else if node_name.starts_with('m') {
        NodeType::Multi
    } else if node_name == "f0" {
        NodeType::FivePrime
    } else if node_name == "t0" {
        NodeType::ThreePrime
    } else if node_name.starts_with('b') {
        // Heuristic: if it has exactly 1 or 2 coords => might be bulge
        // (the Python code used some special logic, but let's keep it simple)
        if coords.len() == 1 || coords.len() == 2 {
            NodeType::Bulge
        } else {
            NodeType::Internal
        }
    } else {
        NodeType::Unknown
    }
}

// ============ Modification logic ============

/// For insertion or deletion
#[derive(Debug)]
enum ModAction {
    Insert,
    Delete,
}

/// Decide whether we do Insert or Delete based on current_size, min, max.
fn choose_action(current_size: usize, min_size: usize, max_size: usize) -> ModAction {
    let mut rng = thread_rng();
    if current_size <= min_size {
        ModAction::Insert
    } else if current_size >= max_size {
        ModAction::Delete
    } else {
        // random
        if rng.gen_bool(0.5) {
            ModAction::Insert
        } else {
            ModAction::Delete
        }
    }
}

/// Insert a base at a random position in the node (for loops).
/// We'll replace one '.' in the structure with '.' + new base. The structure is re-folded later anyway.
fn insert_base(seq: &mut String, structure: &mut String, insertion_idx: usize) {
    // pick a random base
    let base = *['A', 'C', 'G', 'U'].choose(&mut thread_rng()).unwrap();
    seq.insert(insertion_idx, base);
    structure.insert(insertion_idx, '.');
}

/// Delete a base at a given position in the node
fn delete_base(seq: &mut String, structure: &mut String, deletion_idx: usize) {
    if deletion_idx < seq.len() {
        seq.remove(deletion_idx);
        structure.remove(deletion_idx);
    }
}

// ============ Public facing "generate triplet" ============

/// Generate an anchor without folding now.
/// Later the batch function will update its structure.
pub fn generate_anchor(args: &Cli) -> (String, String) {
    let length = choose_seq_length(
        &args.seq_len_distribution,
        args.seq_min_len,
        args.seq_max_len,
        args.seq_len_mean,
        args.seq_len_sd,
    );
    let anchor_seq = generate_random_rna(length);
    // Fold the anchor to obtain an initial dot-bracket structure instead of a placeholder of dots.
    let anchor_structure = fold_rna(&anchor_seq).unwrap_or(".".repeat(length));
    (anchor_seq, anchor_structure)
}

/// Modify a stem node: choose a random eligible stem and perform an insertion or deletion.
fn modify_stem<F>(
    seq: String,
    structure: String,
    args: &Cli,
    mod_counts: &mut HashMap<String, usize>,
    rebuild_bg: F,
) -> (String, String, HashMap<String, usize>)
where
    F: Fn(&str, &mut String) -> BulgeGraph,
{
    let mut pos_seq = seq;
    let mut pos_struct = structure;

    // Rebuild the bulge graph to get node mapping.
    let bg = rebuild_bg(&pos_seq, &mut pos_struct);
    let node_map = bg.node_mapping();

    let eligible_nodes: Vec<String> = node_map
        .iter()
        .filter_map(|(node_name, coords)| {
            if let NodeType::Stem = classify_node(node_name, coords) {
                let count = mod_counts.get(node_name).copied().unwrap_or(0);
                if count < args.stem_max_n_modifications && coords.len() >= args.stem_min_size {
                    Some(node_name.clone())
                } else {
                    None
                }
            } else {
                None
            }
        })
        .collect();

    if eligible_nodes.is_empty() {
        return (pos_seq, pos_struct, mod_counts.clone());
    }

    // Pick a random eligible stem.
    let node_name = eligible_nodes[thread_rng().gen_range(0..eligible_nodes.len())].clone();
    let coords = node_map.get(&node_name).unwrap();
    let current_size = coords.len();
    let action = choose_action(current_size, args.stem_min_size, args.stem_max_size);
    match action {
        ModAction::Insert => {
            let mut sorted_coords = coords.clone();
            sorted_coords.sort();
            let n = sorted_coords.len();
            let left_half = &sorted_coords[..n/2];
            let right_half = &sorted_coords[n/2..];
            let left_last = left_half[left_half.len()-1].saturating_sub(1);
            let right_first = right_half[0].saturating_sub(1);
            let ins_left = left_last + 1; // <-- NEW: define insertion point for left side
            let complement_pairs = [('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G')];
            let (base_left, base_right) = complement_pairs[thread_rng().gen_range(0..complement_pairs.len())];
            if ins_left < right_first {
                pos_seq.insert(right_first, base_right);
                pos_struct.insert(right_first, ')');
                pos_seq.insert(ins_left, base_left);
                pos_struct.insert(ins_left, '(');
            } else {
                pos_seq.insert(ins_left, base_left);
                pos_struct.insert(ins_left, '(');
                pos_seq.insert(right_first, base_right);
                pos_struct.insert(right_first, ')');
            }
        }
        ModAction::Delete => {
            let mut sorted_coords = coords.clone();
            sorted_coords.sort();
            let n = sorted_coords.len();
            if n >= 2 {
                let left_idx = sorted_coords[n/2 - 1].saturating_sub(1);
                let right_idx = sorted_coords[n/2].saturating_sub(1);
                if right_idx > left_idx && right_idx < pos_seq.len() {
                    pos_seq.remove(right_idx);
                    pos_struct.remove(right_idx);
                    pos_seq.remove(left_idx);
                    pos_struct.remove(left_idx);
                }
            }
        }
    }
    *mod_counts.entry(node_name).or_default() += 1;
    (pos_seq, pos_struct, mod_counts.clone())
}

/// Modify a loop region (hairpin/internal/bulge/multi) by selecting an eligible node and performing an insertion or deletion.
fn modify_loop_region<F>(
    seq: String,
    structure: String,
    args: &Cli,
    loop_type: &str,
    mod_counts: &mut HashMap<String, usize>,
    rebuild_bg: F,
) -> (String, String, HashMap<String, usize>)
where
    F: Fn(&str, &mut String) -> BulgeGraph,
{
    let mut pos_seq = seq;
    let mut pos_struct = structure;
    let bg = rebuild_bg(&pos_seq, &mut pos_struct);
    let node_map = bg.node_mapping();

    // Determine constraints based on loop type.
    let (min_size, max_size, max_mods) = match loop_type {
        "hairpin" => (args.hloop_min_size, args.hloop_max_size, args.hloop_max_n_modifications),
        "internal" => (args.iloop_min_size, args.iloop_max_size, args.iloop_max_n_modifications),
        "bulge" => (args.bulge_min_size, args.bulge_max_size, args.bulge_max_n_modifications),
        "multi" => (args.mloop_min_size, args.mloop_max_size, args.mloop_max_n_modifications),
        _ => (1, 9999, 1),
    };

    let required_type = match loop_type {
        "hairpin" => NodeType::Hairpin,
        "internal" => NodeType::Internal,
        "bulge" => NodeType::Bulge,
        "multi" => NodeType::Multi,
        _ => NodeType::Unknown,
    };

    let eligible_nodes: Vec<String> = node_map
        .iter()
        .filter_map(|(node_name, coords)| {
            if classify_node(node_name, coords) == required_type {
                let cur_size = coords.len();
                let count = mod_counts.get(node_name).copied().unwrap_or(0);
                if count < max_mods && cur_size >= min_size && cur_size <= max_size {
                    Some(node_name.clone())
                } else {
                    None
                }
            } else {
                None
            }
        })
        .collect();

    if eligible_nodes.is_empty() {
        return (pos_seq, pos_struct, mod_counts.clone());
    }

    let node_name = eligible_nodes[thread_rng().gen_range(0..eligible_nodes.len())].clone();
    let coords = node_map.get(&node_name).unwrap();
    let current_size = coords.len();
    let action = choose_action(current_size, min_size, max_size);
    match action {
        ModAction::Insert => {
            let idx = coords[thread_rng().gen_range(0..coords.len())].saturating_sub(1);
            insert_base(&mut pos_seq, &mut pos_struct, idx);
        }
        ModAction::Delete => {
            if current_size > min_size {
                let idx = coords[thread_rng().gen_range(0..coords.len())].saturating_sub(1);
                if idx < pos_seq.len() {
                    delete_base(&mut pos_seq, &mut pos_struct, idx);
                }
            }
        }
    }
    *mod_counts.entry(node_name).or_default() += 1;
    (pos_seq, pos_struct, mod_counts.clone())
}

/// Generate a positive sample from an existing anchor by applying modifications.
pub fn generate_positive(
    anchor_seq: &str,
    anchor_structure: &str,
    args: &Cli,
) -> (String, String) {
    let mut pos_seq = anchor_seq.to_string();
    let mut pos_struct = anchor_structure.to_string();

    // Rename unused variable "length" to "_length"
    let _length = anchor_seq.len();
    // Remove unnecessary mut on local counters:
    let n_stem_indels = args.n_stem_indels;
    let n_hloop_indels = args.n_hloop_indels;
    let n_iloop_indels = args.n_iloop_indels;
    let n_bulge_indels = args.n_bulge_indels;
    let n_mloop_indels = args.n_mloop_indels;

    let mut mod_counts: HashMap<String, usize> = HashMap::new();
    // Define a rebuild closure that builds the bulge graph from the current structure.
    let rebuild_bg = |_: &str, st: &mut String| -> BulgeGraph {
        BulgeGraph::from_dot_bracket(st).unwrap_or_else(|_| BulgeGraph::new())
    };

    // Apply stem modifications.
    for _ in 0..n_stem_indels {
        let (new_seq, new_struct, new_mod_counts) =
            modify_stem(pos_seq, pos_struct, args, &mut mod_counts, &rebuild_bg);
        pos_seq = new_seq;
        pos_struct = new_struct;
        mod_counts = new_mod_counts;
    }
    // Apply hairpin loop modifications.
    for _ in 0..n_hloop_indels {
        let (new_seq, new_struct, new_mod_counts) =
            modify_loop_region(pos_seq, pos_struct, args, "hairpin", &mut mod_counts, &rebuild_bg);
        pos_seq = new_seq;
        pos_struct = new_struct;
        mod_counts = new_mod_counts;
    }
    // Apply internal loop modifications.
    for _ in 0..n_iloop_indels {
        let (new_seq, new_struct, new_mod_counts) =
            modify_loop_region(pos_seq, pos_struct, args, "internal", &mut mod_counts, &rebuild_bg);
        pos_seq = new_seq;
        pos_struct = new_struct;
        mod_counts = new_mod_counts;
    }
    // Apply bulge modifications.
    for _ in 0..n_bulge_indels {
        let (new_seq, new_struct, new_mod_counts) =
            modify_loop_region(pos_seq, pos_struct, args, "bulge", &mut mod_counts, &rebuild_bg);
        pos_seq = new_seq;
        pos_struct = new_struct;
        mod_counts = new_mod_counts;
    }
    // Apply multi loop modifications.
    for _ in 0..n_mloop_indels {
        let (new_seq, new_struct, new_mod_counts) =
            modify_loop_region(pos_seq, pos_struct, args, "multi", &mut mod_counts, &rebuild_bg);
        pos_seq = new_seq;
        pos_struct = new_struct;
        mod_counts = new_mod_counts;
    }

    (pos_seq, pos_struct)
}

/// Generate a cluster: one anchor and multiple positives generated from it.
pub fn generate_cluster(args: &Cli, cluster_id: usize) -> Vec<RnaClusterRecord> {
    let (anchor_seq, anchor_structure) = generate_anchor(args);
    let mut records = Vec::with_capacity(args.structures_per_cluster);
    // Add the anchor record (flag 1)
    records.push(RnaClusterRecord {
        cluster: cluster_id,
        seq: anchor_seq.clone(),
        structure: anchor_structure.clone(),
        anchor: 1,
    });
    // Generate (structures_per_cluster - 1) positive records (flag 0)
    for _ in 1..args.structures_per_cluster {
        let (pos_seq, pos_struct) = generate_positive(&anchor_seq, &anchor_structure, args);
        records.push(RnaClusterRecord {
            cluster: cluster_id,
            seq: pos_seq,
            structure: pos_struct,
            anchor: 0,
        });
    }
    records
}

/// Generate all clusters in batches.
pub fn generate_all_clusters(args: &Cli) -> Vec<RnaClusterRecord> {
    use rayon::prelude::*;
    let total_clusters = args.n_clusters;
    let batch_size = args.batch_size;
    let pb = indicatif::ProgressBar::new(total_clusters as u64);
    let chunk_count = (total_clusters + batch_size - 1) / batch_size;
    let chunks: Vec<usize> = (0..chunk_count).collect();

    let results: Vec<RnaClusterRecord> = rayon::ThreadPoolBuilder::new()
        .num_threads(args.num_workers)
        .build()
        .unwrap()
        .install(|| {
            chunks.par_iter().flat_map(|chunk_id| {
                let start = chunk_id * batch_size;
                let cur_batch = std::cmp::min(batch_size, total_clusters - start);
                let mut cluster_records = Vec::new();
                for i in 0..cur_batch {
                    let cluster_id = start + i + 1; // clusters numbered starting at 1.
                    let recs = generate_cluster(args, cluster_id);
                    cluster_records.extend(recs);
                }
                pb.inc(cur_batch as u64);
                cluster_records
            }).collect()
        });
    pb.finish();
    results
}


/// Change save_dataset_csv to be generic over any serializable type.
pub fn save_dataset_csv<T: Serialize>(
    file_path: &str,
    data: &[T],
    metadata: &serde_json::Value,
) -> Result<(), String> {
    let file = File::create(file_path)
        .map_err(|e| format!("Failed to create file {}: {}", file_path, e))?;
    let mut writer = std::io::BufWriter::new(file);

    let meta_str = serde_json::to_string(metadata)
        .map_err(|e| format!("Cannot serialize metadata: {}", e))?;
    writeln!(writer, "# Metadata: {}", meta_str)
        .map_err(|e| format!("Failed to write metadata comment: {}", e))?;

    let mut wtr = csv::Writer::from_writer(writer);
    for row in data {
        wtr.serialize(row).map_err(|e| e.to_string())?;
    }
    wtr.flush().map_err(|e| e.to_string())?;
    Ok(())
}

/// Update metadata to remove keys for removed options.
pub fn build_metadata(args: &Cli) -> serde_json::Value {
    let run_id = uuid::Uuid::new_v4();
    serde_json::json!({
        "run_id": run_id.to_string(),
        "timestamp": chrono::Local::now().to_rfc3339(),
        "parameters": {
            "n_clusters": args.n_clusters,
            "structures_per_cluster": args.structures_per_cluster,
            "seq_min_len": args.seq_min_len,
            "seq_max_len": args.seq_max_len,
            "seq_len_distribution": format!("{:?}", args.seq_len_distribution),
            "seq_len_mean": args.seq_len_mean,
            "seq_len_sd": args.seq_len_sd,
            "n_stem_indels": args.n_stem_indels,
            "stem_min_size": args.stem_min_size,
            "stem_max_size": args.stem_max_size,
            "stem_max_n_modifications": args.stem_max_n_modifications,
            "n_hloop_indels": args.n_hloop_indels,
            "hloop_min_size": args.hloop_min_size,
            "hloop_max_size": args.hloop_max_size,
            "hloop_max_n_modifications": args.hloop_max_n_modifications,
            "n_iloop_indels": args.n_iloop_indels,
            "iloop_min_size": args.iloop_min_size,
            "iloop_max_size": args.iloop_max_size,
            "iloop_max_n_modifications": args.iloop_max_n_modifications,
            "n_bulge_indels": args.n_bulge_indels,
            "bulge_min_size": args.bulge_min_size,
            "bulge_max_size": args.bulge_max_size,
            "bulge_max_n_modifications": args.bulge_max_n_modifications,
            "n_mloop_indels": args.n_mloop_indels,
            "mloop_min_size": args.mloop_min_size,
            "mloop_max_size": args.mloop_max_size,
            "mloop_max_n_modifications": args.mloop_max_n_modifications,
            "num_workers": args.num_workers,
            "batch_size": args.batch_size,
            "plot": args.plot,
            "num_plots": args.num_plots,
            "debug": args.debug,
            "timing_log": args.timing_log
        }
    })
}

// ----- New function: batched RNAfold -----
// Add an attribute to suppress the dead code warning.
#[allow(dead_code)]
pub fn fold_rna_batch(sequences: &[(String, String)]) -> Result<HashMap<String, String>, String> {
    // Build multi-FASTA input. Each entry: ">ID\nSEQUENCE"
    let input = sequences.iter()
        .map(|(id, seq)| format!(">{}\n{}", id, seq))
        .collect::<Vec<_>>()
        .join("\n");
    
    let mut child = std::process::Command::new("RNAfold")
        .arg("--noPS")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to spawn RNAfold: {}", e))?;
    
    {
        let mut stdin = child.stdin.take().ok_or("Failed to open RNAfold stdin")?;
        std::io::Write::write_all(&mut stdin, input.as_bytes())
            .map_err(|e| format!("Error writing to RNAfold stdin: {}", e))?;
    }
    let output = child.wait_with_output().map_err(|e| format!("RNAfold error: {}", e))?;
    if !output.status.success() {
        return Err(format!("RNAfold exited with status: {}", output.status));
    }
    let reader = BufReader::new(&output.stdout[..]);
    let mut results = HashMap::new();
    // let mut current_id = String::new();  // Removed unused assignment
    let _structure_line = String::new();
    
    // Parse FASTA output.
    // Assume the output format:
    // >ID
    // SEQUENCE (same as input)
    // dot-bracket and free energy info line; we take first token as structure.
    let mut lines = reader.lines();
    while let Some(line_res) = lines.next() {
        let line = line_res.map_err(|e| format!("Reading RNAfold output: {}", e))?;
        if line.starts_with('>') {
            let current_id = line[1..].trim().to_string();
            // skip the sequence line:
            lines.next();
            // next line is structure line.
            if let Some(sline) = lines.next() {
                let s = sline.map_err(|e| format!("Error reading structure: {}", e))?;
                // take everything up to first space.
                let structure = s.split_whitespace().next().unwrap_or("").to_string();
                results.insert(current_id.clone(), structure);
            }
        }
    }
    Ok(results)
}
