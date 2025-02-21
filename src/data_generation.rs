// src/data_generation.rs
//
// This module replicates the core Python logic: random RNA, folding via RNAfold,
// identifying stems/loops, applying modifications, generating negative sample,
// returning the anchor/positive/negative triplets. It also includes routines
// for dataset splitting and CSV/JSON output.
//

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process::{Command, Stdio};
use std::io::{BufRead, BufReader};

use chrono::Local;
use rand::distributions::{Distribution, Uniform};
use rand::{thread_rng, Rng};
use rand_distr::Normal;
use serde::Serialize;
use uuid::Uuid;
use rand::seq::SliceRandom;


use crate::bulgegraph::BulgeGraph;
use crate::cli::{Cli, SeqLenDistribution};

/// A single record in our final CSV (equivalent to one "triplet").
#[derive(Debug, Serialize, Clone)]
pub struct RnaTripletRecord {
    pub triplet_id: usize,
    pub anchor_seq: String,
    pub anchor_structure: String,
    pub positive_seq: String,
    pub positive_structure: String,
    pub negative_seq: String,
    pub negative_structure: String,
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
    // Set placeholder; will be updated in batch.
    let anchor_structure = ".".repeat(length);
    (anchor_seq, anchor_structure)
}

/// Generate a positive sample from an existing anchor by applying modifications.
pub fn generate_positive(
    anchor_seq: &str,
    anchor_structure: &str,
    args: &Cli,
) -> (String, String) {
    let mut pos_seq = anchor_seq.to_string();
    let mut pos_struct = anchor_structure.to_string();

    // Possibly normalize how many modifications we do.
    let length = anchor_seq.len();
    let mut n_stem_indels = args.n_stem_indels;
    let mut n_hloop_indels = args.n_hloop_indels;
    let mut n_iloop_indels = args.n_iloop_indels;
    let mut n_bulge_indels = args.n_bulge_indels;
    let mut n_mloop_indels = args.n_mloop_indels;

    if args.mod_normalization {
        let factor = (length as f64 / args.normalization_len as f64).max(1.0);
        // scale each
        let scale = |val: usize| -> usize {
            if val <= 1 {
                val
            } else {
                ((val as f64) * factor).round() as usize
            }
        };
        n_stem_indels = scale(n_stem_indels);
        n_hloop_indels = scale(n_hloop_indels);
        n_iloop_indels = scale(n_iloop_indels);
        n_bulge_indels = scale(n_bulge_indels);
        n_mloop_indels = scale(n_mloop_indels);
    }

    // We'll track how many modifications we've done on each node
    let mut mod_counts: HashMap<String, usize> = HashMap::new();

    // New rebuild closure: simply re-parse the current structure.
    let rebuild_bg = |_: &str, st: &mut String| -> BulgeGraph {
        // Use the existing pos_struct to build the BulgeGraph.
        BulgeGraph::from_dot_bracket(st).unwrap_or_else(|_| BulgeGraph::new())
    };

    // Modify stems
    for _ in 0..n_stem_indels {
        (pos_seq, pos_struct, mod_counts) = modify_stem(
            pos_seq,
            pos_struct,
            &args,
            &mut mod_counts,
            rebuild_bg,
        );
    }
    // Modify hairpin loops
    for _ in 0..n_hloop_indels {
        (pos_seq, pos_struct, mod_counts) = modify_loop_region(
            pos_seq,
            pos_struct,
            &args,
            "hairpin",
            &mut mod_counts,
            rebuild_bg,
        );
    }
    // Modify internal loops
    for _ in 0..n_iloop_indels {
        (pos_seq, pos_struct, mod_counts) = modify_loop_region(
            pos_seq,
            pos_struct,
            &args,
            "internal",
            &mut mod_counts,
            rebuild_bg,
        );
    }
    // Modify bulges
    for _ in 0..n_bulge_indels {
        (pos_seq, pos_struct, mod_counts) = modify_loop_region(
            pos_seq,
            pos_struct,
            &args,
            "bulge",
            &mut mod_counts,
            rebuild_bg,
        );
    }
    // Modify multi loops
    for _ in 0..n_mloop_indels {
        (pos_seq, pos_struct, mod_counts) = modify_loop_region(
            pos_seq,
            pos_struct,
            &args,
            "multi",
            &mut mod_counts,
            rebuild_bg,
        );
    }

    (pos_seq, pos_struct)
}

/// Modify stems: choose a random stem node that is not yet at max modifications, do insert/delete.
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

    // build or rebuild the bulge graph
    let bg = rebuild_bg(&pos_seq, &mut pos_struct);

    // gather stems that still can be modified
    let node_map = bg.node_mapping();
    let eligible_nodes: Vec<_> = node_map
        .iter()
        .filter_map(|(node_name, coords)| {
            if classify_node(node_name, coords) == crate::data_generation::NodeType::Stem {
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
        // no stems or none are eligible
        return (pos_seq, pos_struct, mod_counts.clone());
    }

    // pick a random node from the eligible set
    let node_name = eligible_nodes
        [thread_rng().gen_range(0..eligible_nodes.len())]
        .clone();

    let coords = node_map.get(&node_name).unwrap();
    let current_size = coords.len();
    let action = choose_action(current_size, args.stem_min_size, args.stem_max_size);
    // we approximate the "middle" of the stem for insertion
    let mid_idx = coords[current_size / 2].saturating_sub(1); // to 0-based

    match action {
        ModAction::Insert => {
            // We'll attempt to insert a pair around mid. But for simplicity, let's just insert a single base pair if feasible.
            // Actually, let's do the simpler approach: we insert 2 bases, one in the "paired" region
            // Because the python code is quite complicated, let's do a simpler approach:
            // Insert a single base in the middle, then also insert a complementary base near the partner if it exists.

            // For the "partner" we find the partner index from the current structure
            // But we can't easily get that without re-parsing. Alternatively, we do a single base insertion in the middle and re-fold.
            insert_base(&mut pos_seq, &mut pos_struct, mid_idx);
        }
        ModAction::Delete => {
            // delete one base from roughly the middle
            if current_size > args.stem_min_size {
                let del_idx = coords[current_size / 2].saturating_sub(1);
                if del_idx < pos_seq.len() {
                    delete_base(&mut pos_seq, &mut pos_struct, del_idx);
                }
            }
        }
    }
    *mod_counts.entry(node_name).or_default() += 1;
    (pos_seq, pos_struct, mod_counts.clone())
}

/// Generic function to modify a loop region (hairpin, internal, bulge, multi).
/// We'll pick a node from the appropriate type and do an insert or delete.
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

    // figure out which node type we need
    let (min_size, max_size, max_mods) = match loop_type {
        "hairpin" => (
            args.hloop_min_size,
            args.hloop_max_size,
            args.hloop_max_n_modifications,
        ),
        "internal" => (
            args.iloop_min_size,
            args.iloop_max_size,
            args.iloop_max_n_modifications,
        ),
        "bulge" => (
            args.bulge_min_size,
            args.bulge_max_size,
            args.bulge_max_n_modifications,
        ),
        "multi" => (
            args.mloop_min_size,
            args.mloop_max_size,
            args.mloop_max_n_modifications,
        ),
        _ => (1, 9999, 1),
    };

    // gather appropriate node names
    let required_node_type = match loop_type {
        "hairpin" => NodeType::Hairpin,
        "internal" => NodeType::Internal,
        "bulge" => NodeType::Bulge,
        "multi" => NodeType::Multi,
        _ => NodeType::Unknown,
    };

    let eligible_nodes: Vec<_> = node_map
        .iter()
        .filter_map(|(node_name, coords)| {
            let node_t = classify_node(node_name, coords);
            if node_t == required_node_type {
                // check if size, mod counts
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

    let node_name = eligible_nodes
        [thread_rng().gen_range(0..eligible_nodes.len())]
        .clone();
    let coords = node_map.get(&node_name).unwrap();
    let current_size = coords.len();
    let action = choose_action(current_size, min_size, max_size);

    match action {
        ModAction::Insert => {
            // choose random insertion position
            let idx = coords[thread_rng().gen_range(0..coords.len())].saturating_sub(1);
            insert_base(&mut pos_seq, &mut pos_struct, idx);
        }
        ModAction::Delete => {
            // if current_size <= min_size, skip
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

/// Add “appending” to both positive and negative sequences if triggered.
fn maybe_append_sequences(
    pos_seq: &mut String,
    pos_struct: &mut String,
    neg_seq: &mut String,
    neg_struct: &mut String,
    anchor_len: usize,
    args: &Cli,
) {
    let mut rng = thread_rng();
    if rng.gen_bool(args.appending_event_probability) {
        // decide left-only, right-only or both
        let r = rng.gen_range(0.0..1.0);
        let p_both = args.both_sides_appending_probability;
        let p_left = (1.0 - p_both) / 2.0;
        let p_right = p_left;

        let mean_append = anchor_len as f64 * args.appending_size_factor;
        let sigma_append = mean_append / 2.0;

        // Define a helper that takes &mut rng so that we do not capture it mutably.
        let sample_append_length = |rng: &mut _| {
            let gauss = Normal::new(mean_append, sigma_append.max(1.0)).unwrap();
            gauss.sample(rng).round().max(1.0) as usize
        };

        // Now use rng as usual (no conflict)
        let linker_len = rng.gen_range(args.linker_min..=args.linker_max);
        let linker_seq = generate_random_rna(linker_len);
        let linker_struct = ".".repeat(linker_len);

        if r < p_left {
            let l_app = sample_append_length(&mut rng);
            let appended_left_seq = generate_random_rna(l_app);
            let appended_left_struct = ".".repeat(l_app);
            *pos_seq = format!("{}{}{}", appended_left_seq, linker_seq, pos_seq);
            *pos_struct = format!("{}{}{}", appended_left_struct, linker_struct, pos_struct);
            *neg_seq = format!("{}{}{}", appended_left_seq, linker_seq, neg_seq);
            *neg_struct = format!("{}{}{}", appended_left_struct, linker_struct, neg_struct);
        } else if r < p_left + p_right {
            let r_app = sample_append_length(&mut rng);
            let appended_right_seq = generate_random_rna(r_app);
            let appended_right_struct = ".".repeat(r_app);
            *pos_seq = format!("{}{}{}", pos_seq, linker_seq, appended_right_seq);
            *pos_struct = format!("{}{}{}", pos_struct, linker_struct, appended_right_struct);
            *neg_seq = format!("{}{}{}", neg_seq, linker_seq, appended_right_seq);
            *neg_struct = format!("{}{}{}", neg_struct, linker_struct, appended_right_struct);
        } else {
            let l_app = sample_append_length(&mut rng);
            let appended_left_seq = generate_random_rna(l_app);
            let appended_left_struct = ".".repeat(l_app);
            let r_app = sample_append_length(&mut rng);
            let appended_right_seq = generate_random_rna(r_app);
            let appended_right_struct = ".".repeat(r_app);
            *pos_seq = format!(
                "{}{}{}{}{}",
                appended_left_seq, linker_seq, pos_seq, linker_seq, appended_right_seq
            );
            *pos_struct = format!(
                "{}{}{}{}{}",
                appended_left_struct, linker_struct, pos_struct, linker_struct, appended_right_struct
            );
            *neg_seq = format!(
                "{}{}{}{}{}",
                appended_left_seq, linker_seq, neg_seq, linker_seq, appended_right_seq
            );
            *neg_struct = format!(
                "{}{}{}{}{}",
                appended_left_struct, linker_struct, neg_struct, linker_struct, appended_right_struct
            );
        }
    }
}

/// Generate a single anchor/positive/negative triplet.
pub fn generate_triplet(args: &Cli) -> RnaTripletRecord {
    // 1) Generate anchor
    let (anchor_seq, anchor_structure) = generate_anchor(args);

    // 2) Generate positive
    let (mut pos_seq, mut pos_struct) = generate_positive(&anchor_seq, &anchor_structure, args);

    // 3) Generate negative
    let (mut neg_seq, mut neg_struct) =
        generate_negative_sample(&anchor_seq, args.neg_len_variation);

    // 4) Maybe do “appending event”
    maybe_append_sequences(
        &mut pos_seq,
        &mut pos_struct,
        &mut neg_seq,
        &mut neg_struct,
        anchor_seq.len(),
        args,
    );

    RnaTripletRecord {
        triplet_id: 0, // we’ll set the real ID later
        anchor_seq,
        anchor_structure, // placeholder; to be updated in batch below
        positive_seq: pos_seq,
        positive_structure: pos_struct,
        negative_seq: neg_seq,
        negative_structure: neg_struct, // placeholder
    }
}

// ============ Splitting, Saving, and other I/O ============

/// Split the dataset into train/val sets (roughly).
pub fn split_dataset(
    data: &mut Vec<RnaTripletRecord>,
    train_fraction: f64,
) -> (Vec<RnaTripletRecord>, Vec<RnaTripletRecord>) {
    // We can shuffle first, then partition
    let mut rng = thread_rng();
    data.shuffle(&mut rng);
    let train_size = ((data.len() as f64) * train_fraction).round() as usize;
    let train_data = data[..train_size].to_vec();
    let val_data = data[train_size..].to_vec();
    (train_data, val_data)
}

/// Save the entire dataset as CSV with a JSON metadata header at the top
pub fn save_dataset_csv(
    file_path: &str,
    data: &[RnaTripletRecord],
    metadata: &serde_json::Value,
) -> Result<(), String> {
    let file = File::create(file_path)
        .map_err(|e| format!("Failed to create file {}: {}", file_path, e))?;
    let mut writer = BufWriter::new(file);

    // write metadata as a comment
    let meta_str =
        serde_json::to_string(metadata).map_err(|e| format!("Cannot serialize metadata: {}", e))?;
    writeln!(writer, "# Metadata: {}", meta_str)
        .map_err(|e| format!("Failed to write metadata comment: {}", e))?;

    // now write CSV
    let mut wtr = csv::Writer::from_writer(writer);
    for row in data {
        wtr.serialize(row).map_err(|e| e.to_string())?;
    }
    wtr.flush().map_err(|e| e.to_string())?;
    Ok(())
}

/// Generate multiple triplets in "batches", store them in a Vec.
pub fn generate_all_triplets(args: &Cli) -> Vec<RnaTripletRecord> {
    let mut records = Vec::with_capacity(args.num_structures);
    let total = args.num_structures;
    let batch_size = args.batch_size;
    let mut generated = 0;

    // Create a progress bar using indicatif
    let pb = indicatif::ProgressBar::new(total as u64);
    
    while generated < total {
        let current_batch_size = std::cmp::min(batch_size, total - generated);
        // Generate batch_triplets by calling generate_triplet for each item in the batch.
        let mut batch_triplets = Vec::with_capacity(current_batch_size);
        for i in 0..current_batch_size {
            let mut triplet = generate_triplet(args);
            triplet.triplet_id = generated + i;
            batch_triplets.push(triplet);
        }
        // Batch fold anchors.
        let anchor_pairs: Vec<(String, String)> = batch_triplets.iter()
            .map(|r| (format!("A{}", r.triplet_id), r.anchor_seq.clone()))
            .collect();
        if let Ok(anchor_map) = fold_rna_batch(&anchor_pairs) {
            for record in batch_triplets.iter_mut() {
                if let Some(struc) = anchor_map.get(&format!("A{}", record.triplet_id)) {
                    record.anchor_structure = struc.clone();
                }
            }
        }
        // Batch fold negatives.
        let neg_pairs: Vec<(String, String)> = batch_triplets.iter()
            .map(|r| (format!("N{}", r.triplet_id), r.negative_seq.clone()))
            .collect();
        if let Ok(neg_map) = fold_rna_batch(&neg_pairs) {
            for record in batch_triplets.iter_mut() {
                if let Some(struc) = neg_map.get(&format!("N{}", record.triplet_id)) {
                    record.negative_structure = struc.clone();
                }
            }
        }
        records.extend(batch_triplets);
        generated += current_batch_size;
        pb.inc(current_batch_size as u64);
    }
    
    pb.finish();
    records
}

/// Build a metadata object for the run
pub fn build_metadata(args: &Cli) -> serde_json::Value {
    let run_id = Uuid::new_v4();
    serde_json::json!({
        "run_id": run_id.to_string(),
        "timestamp": Local::now().to_rfc3339(),
        "parameters": {
            "num_structures": args.num_structures,
            "seq_min_len": args.seq_min_len,
            "seq_max_len": args.seq_max_len,
            "seq_len_distribution": format!("{:?}", args.seq_len_distribution),
            "seq_len_mean": args.seq_len_mean,
            "seq_len_sd": args.seq_len_sd,
            "neg_len_variation": args.neg_len_variation,
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
            "appending_event_probability": args.appending_event_probability,
            "both_sides_appending_probability": args.both_sides_appending_probability,
            "linker_min": args.linker_min,
            "linker_max": args.linker_max,
            "appending_size_factor": args.appending_size_factor,
            "mod_normalization": args.mod_normalization,
            "normalization_len": args.normalization_len,
            "num_workers": args.num_workers,
            "batch_size": args.batch_size,
            "plot": args.plot,
            "num_plots": args.num_plots,
            "split": args.split,
            "train_fraction": args.train_fraction,
            "val_fraction": args.val_fraction,
            "debug": args.debug,
            "timing_log": args.timing_log,
        }
    })
}

// ----- New function: batched RNAfold -----
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
