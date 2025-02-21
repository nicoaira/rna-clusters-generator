// src/bulgegraph.rs
//
// This module provides functionality to parse a dot-bracket RNA secondary
// structure and build a "bulge graph" representation. It is a stand-in
// for the forgi.BulgeGraph methods in Python.
//

use std::collections::{HashMap, HashSet};

/// A representation of a BulgeGraph node mapping (like the Python "define <node> ...")
/// 
/// Keys are node labels (e.g., "s0", "i0", etc.) and values are the list of 1-based
/// indices belonging to that node.
pub type NodeMapping = HashMap<String, Vec<usize>>;

/// An internal struct to represent the bulge graph in memory.
#[derive(Debug)]
pub struct BulgeGraph {
    /// For each node name, the list of 1-based coordinates.
    defines: NodeMapping,

    /// Adjacency: for each node name, a set of connected node names.
    edges: HashMap<String, HashSet<String>>,

    /// "Weights" used in merges. Usually 1, but can sum up if merged.
    weights: HashMap<String, usize>,

    /// A counter to generate fresh node names if needed.
    name_counter: usize,
}

impl BulgeGraph {
    /// Create a new, empty `BulgeGraph`.
    pub fn new() -> Self {
        Self {
            defines: HashMap::new(),
            edges: HashMap::new(),
            weights: HashMap::new(),
            name_counter: 0,
        }
    }

    /// Main entry point: parse a dot-bracket structure into BulgeGraph.
    pub fn from_dot_bracket(db: &str) -> Result<Self, String> {
        // parse the bracket into (pos, partner)
        let pairs = parse_dot_bracket(db).map_err(|e| e.to_string())?;

        // convert to a sorted list of (pos, partner) for pairs
        let mut bg = BulgeGraph::new();
        bg.from_tuples(pairs);
        Ok(bg)
    }

    /// Access the final node mapping (“defines”) after building.
    pub fn node_mapping(&self) -> &NodeMapping {
        &self.defines
    }

    /// Internally build the graph from the (pos, partner) pairs.
    fn from_tuples(&mut self, mut tuples: Vec<(usize, usize)>) {
        // Sort by the position (i.e. the first element in each tuple)
        tuples.sort_by_key(|(a, _)| *a);

        let mut stems = Vec::new();
        let mut bulges = Vec::new();

        let mut iter = tuples.into_iter();
        let first = iter.next().expect("No tuples provided!");
        let (mut prev_from, mut prev_to) = first;

        let (mut start_from, mut start_to) = (prev_from, prev_to);
        let mut last_paired = prev_from;

        // Build up the list of stems and bulges
        for (from_bp, to_bp) in iter {
            if to_bp != 0 && prev_to != 0 && (to_bp as isize - prev_to as isize).abs() == 1 {
                let prev_diff = prev_to as isize - prev_from as isize;
                let curr_diff = to_bp as isize - from_bp as isize;
                let check1 = (prev_diff > 0 && curr_diff > 0) || (prev_diff < 0 && curr_diff < 0);
                let check2 = (to_bp as isize - prev_to as isize) == -(from_bp as isize - prev_from as isize);
                if check1 && check2 {
                    // extend the current stem segment
                    prev_from = from_bp;
                    prev_to = to_bp;
                    last_paired = from_bp;
                    continue;
                }
            }

            // if both unpaired => treat as continuous bulge (just continue, no new entry)
            if to_bp == 0 && prev_to == 0 {
                prev_from = from_bp;
                prev_to = to_bp;
                continue;
            } else {
                // if the last region was paired, finalize that stem
                if prev_to != 0 {
                    let new_stem = ((start_from - 1, start_to - 1), (prev_from - 1, prev_to - 1));
                    let sorted = sort_stem_tuple(new_stem);
                    if !stems.contains(&sorted) {
                        stems.push(sorted);
                    }
                    last_paired = from_bp;
                    start_from = from_bp;
                    start_to = to_bp;
                } else {
                    // finalize a bulge
                    let new_bulge = (last_paired - 1, prev_from - 1);
                    bulges.push(new_bulge);

                    start_from = from_bp;
                    start_to = to_bp;
                }
                prev_from = from_bp;
                prev_to = to_bp;
            }
        }

        // handle the last segment
        if prev_to != 0 {
            let new_stem = ((start_from - 1, start_to - 1), (prev_from - 1, prev_to - 1));
            let sorted = sort_stem_tuple(new_stem);
            if !stems.contains(&sorted) {
                stems.push(sorted);
            }
        } else {
            let new_bulge = (last_paired - 1, prev_from - 1);
            bulges.push(new_bulge);
        }

        // fill the BulgeGraph with stems and bulges
        self.from_stems_and_bulges(&stems, &bulges);
    }

    fn from_stems_and_bulges(
        &mut self,
        stems: &Vec<((usize, usize), (usize, usize))>,
        bulges: &Vec<(usize, usize)>,
    ) {
        // Insert stems as y# => [s1, s2, e1, e2]
        for (i, stem) in stems.iter().enumerate() {
            let (s1, e1) = stem.0;
            let (s2, e2) = stem.1;
            let node_id = format!("y{}", i);

            self.defines.insert(
                node_id.clone(),
                vec![
                    s1.min(s2) + 1,
                    s1.max(s2) + 1,
                    e1.min(e2) + 1,
                    e1.max(e2) + 1,
                ],
            );
            self.weights.insert(node_id.clone(), 1);
        }

        // Insert bulges as b# => 2 coords
        for (i, &(bstart, bend)) in bulges.iter().enumerate() {
            let node_id = format!("b{}", i);
            let sorted = if bstart <= bend {
                vec![bstart + 1, bend + 1]
            } else {
                vec![bend + 1, bstart + 1]
            };
            self.defines.insert(node_id.clone(), sorted);
            self.weights.insert(node_id.clone(), 1);
        }

        // create adjacency
        self.create_bulge_graph(stems, bulges);
        self.create_stem_graph(stems, bulges.len());
        self.collapse();
        self.sort_defines();
        self.relabel_nodes();
        self.remove_degenerate_nodes();
    }

    fn create_bulge_graph(
        &mut self,
        stems: &Vec<((usize, usize), (usize, usize))>,
        bulges: &Vec<(usize, usize)>,
    ) {
        for (i, stem) in stems.iter().enumerate() {
            for (j, &(bstart, bend)) in bulges.iter().enumerate() {
                if any_difference_of_one(*stem, (bstart, bend)) {
                    let stem_node = format!("y{}", i);
                    let bulge_node = format!("b{}", j);
                    self.add_edge(&stem_node, &bulge_node);
                }
            }
        }
    }

    fn create_stem_graph(
        &mut self,
        stems: &Vec<((usize, usize), (usize, usize))>,
        mut bulge_counter: usize,
    ) {
        let n = stems.len();

        // For each pair of stems i < j
        for i in 0..n {
            for j in (i + 1)..n {
                // check adjacency of each side
                for k1 in 0..2 {
                    for k2 in 0..2 {
                        for l1 in 0..2 {
                            for l2 in 0..2 {
                                // skip degenerate sides
                                let si0 = stems[i].0;
                                let si1 = stems[i].1;
                                let sj0 = stems[j].0;
                                let sj1 = stems[j].1;

                                let s_i_0_l1 = if l1 == 0 { si0.0 } else { si0.1 };
                                let s_i_1_l1 = if l1 == 0 { si1.0 } else { si1.1 };

                                if k1 == 1 && s_i_0_l1 == s_i_1_l1 {
                                    continue;
                                }

                                let s_j_0_l2 = if l2 == 0 { sj0.0 } else { sj0.1 };
                                let s_j_1_l2 = if l2 == 0 { sj1.0 } else { sj1.1 };

                                if k2 == 1 && s_j_0_l2 == s_j_1_l2 {
                                    continue;
                                }

                                // actual coords s1, s2
                                let s1 = if k1 == 0 {
                                    if l1 == 0 { si0.0 } else { si0.1 }
                                } else {
                                    if l1 == 0 { si1.0 } else { si1.1 }
                                };
                                let s2 = if k2 == 0 {
                                    if l2 == 0 { sj0.0 } else { sj0.1 }
                                } else {
                                    if l2 == 0 { sj1.0 } else { sj1.1 }
                                };

                                // if adjacent => new b# with empty coords
                                if (s1 as isize - s2 as isize).abs() == 1 {
                                    let bn = format!("b{}", bulge_counter);
                                    bulge_counter += 1;
                                    self.defines.insert(bn.clone(), vec![]);
                                    self.weights.insert(bn.clone(), 1);

                                    let y_i = format!("y{}", i);
                                    let y_j = format!("y{}", j);
                                    self.add_edge(&bn, &y_i);
                                    self.add_edge(&bn, &y_j);
                                }
                            }
                        }
                    }
                }
            }
        }

        // 0-nt hairpins => check each y#
        for i in 0..n {
            let y_name = format!("y{}", i);
            if let Some(d) = self.defines.get(&y_name) {
                if d.len() == 4 {
                    let (_s1, e1, s2, _e2) = (d[0], d[1], d[2], d[3]);
                    if (s2 as isize - e1 as isize).abs() == 1 {
                        let bn = format!("b{}", bulge_counter);
                        bulge_counter += 1;
                        self.defines.insert(bn.clone(), vec![]);
                        self.weights.insert(bn.clone(), 1);
                        self.add_edge(&bn, &y_name);
                    }
                }
            }
        }
    }

    fn collapse(&mut self) {
        let mut merged_once = true;
        while merged_once {
            merged_once = false;
            // gather non-stem nodes
            let all_nodes: Vec<String> = self
                .defines
                .keys()
                .filter(|k| !k.starts_with('y') && !k.starts_with('s'))
                .cloned()
                .collect();

            let len_nodes = all_nodes.len();
            'outer: for i in 0..len_nodes {
                for j in (i + 1)..len_nodes {
                    let b1 = &all_nodes[i];
                    let b2 = &all_nodes[j];

                    if !self.defines.contains_key(b1) || !self.defines.contains_key(b2) {
                        continue;
                    }
                    let e_b1 = match self.edges.get(b1) {
                        Some(s) => s,
                        None => continue,
                    };
                    let e_b2 = match self.edges.get(b2) {
                        Some(s) => s,
                        None => continue,
                    };
                    if e_b1 == e_b2 && e_b1.len() == 2 {
                        // merge
                        self.merge_vertices(&[b1.clone(), b2.clone()]);
                        merged_once = true;
                        break 'outer;
                    }
                }
            }
        }
    }

    fn sort_defines(&mut self) {
        let keys: Vec<String> = self.defines.keys().cloned().collect();
        for k in keys {
            if let Some(d) = self.defines.get_mut(&k) {
                if d.len() == 4 {
                    let (l0, l1, r0, r1) = (d[0], d[1], d[2], d[3]);
                    if l0 > r0 {
                        d[0] = r0;
                        d[1] = r1;
                        d[2] = l0;
                        d[3] = l1;
                    }
                }
            }
        }
    }

    fn relabel_nodes(&mut self) {
        let mut stems = vec![];
        let mut hairpins = vec![];
        let mut interior = vec![];
        let mut multi = vec![];
        let mut fivep = vec![];
        let mut threep = vec![];

        let seq_len = self.total_length_guess();

        // classify each node
        for node in self.defines.keys() {
            // treat 'y' or 's' => stems
            if node.starts_with('y') || node.starts_with('s') {
                stems.push(node.clone());
                continue;
            }
            let coords = self.defines[node].clone();
            let deg = self.edges.get(node).map_or(0, |s| s.len());
            let w = *self.weights.get(node).unwrap_or(&1);

            // Empty + 1 neighbor => hairpin
            if coords.is_empty() && deg == 1 {
                hairpins.push(node.clone());
                continue;
            }
            // Empty + 2 neighbors => multi
            if coords.is_empty() && deg == 2 {
                multi.push(node.clone());
                continue;
            }
            // 5' unpaired
            if deg == 1 && !coords.is_empty() && coords[0] == 1 {
                fivep.push(node.clone());
                continue;
            }
            // 3' unpaired
            if deg == 1 && !coords.is_empty() && coords.last() == Some(&seq_len) {
                threep.push(node.clone());
                continue;
            }
            // single-edge but not 5'/3' => hairpin
            if deg == 1 && !coords.is_empty() && coords[0] != 1 && coords.last() != Some(&seq_len) {
                hairpins.push(node.clone());
                continue;
            }
            // 2 edges, weight=1 => multi
            if deg == 2 && w == 1 && !coords.is_empty() {
                multi.push(node.clone());
                continue;
            }
            // else => interior
            interior.push(node.clone());
        }

        // Sort them for consistent labeling
        stems.sort_by_key(|n| self.defines[n][0]);
        hairpins.sort_by_key(|n| self.defines[n].first().cloned().unwrap_or(0));
        multi.sort_by_key(|n| self.define_a(n));
        interior.sort_by_key(|n| self.defines[n][0]);

        // rename single f0 or t0
        if fivep.len() == 1 {
            let old = &fivep[0];
            self.relabel_node(old, "f0");
        }
        if threep.len() == 1 {
            let old = &threep[0];
            self.relabel_node(old, "t0");
        }

        // stems => s#
        for (i, old) in stems.iter().enumerate() {
            let new_name = format!("s{}", i);
            self.relabel_node(old, &new_name);
        }
        // interior => i#
        for (i, old) in interior.iter().enumerate() {
            let new_name = format!("i{}", i);
            self.relabel_node(old, &new_name);
        }
        // multi => m#
        for (i, old) in multi.iter().enumerate() {
            let new_name = format!("m{}", i);
            self.relabel_node(old, &new_name);
        }
        // hairpin => h#
        for (i, old) in hairpins.iter().enumerate() {
            let new_name = format!("h{}", i);
            self.relabel_node(old, &new_name);
        }
    }

    fn remove_degenerate_nodes(&mut self) {
        let mut to_remove = vec![];
        for (k, coords) in &self.defines {
            if k.starts_with('h') && coords.is_empty() {
                to_remove.push(k.clone());
            }
        }
        for r in to_remove {
            self.remove_vertex(&r);
        }
    }

    fn total_length_guess(&self) -> usize {
        let mut max_val = 0;
        for coords in self.defines.values() {
            for &c in coords {
                if c > max_val {
                    max_val = c;
                }
            }
        }
        max_val
    }

    fn define_a(&self, node: &str) -> usize {
        self.defines[node].first().cloned().unwrap_or(0)
    }

    fn merge_vertices(&mut self, vs: &[String]) -> String {
        let new_vertex = self.get_vertex(None);
        self.weights.insert(new_vertex.clone(), 0);

        let mut all_connections = HashSet::new();
        let mut merged_coords = Vec::new();
        let mut total_w = 0;

        for v in vs {
            if !self.defines.contains_key(v) {
                continue;
            }
            if let Some(nb) = self.edges.get(v).cloned() {
                for n in nb {
                    all_connections.insert(n);
                }
            }
            let c = self.defines[v].clone();
            merged_coords.extend(c);
            total_w += self.weights[v];

            self.remove_vertex(v);
        }

        self.defines.insert(new_vertex.clone(), merged_coords);
        self.weights.insert(new_vertex.clone(), total_w);

        self.edges.insert(new_vertex.clone(), HashSet::new());
        for c in &all_connections {
            if c != &new_vertex {
                self.add_edge(&new_vertex, c);
            }
        }

        new_vertex
    }

    fn relabel_node(&mut self, old: &str, new: &str) {
        if self.defines.contains_key(old) {
            let coords = self.defines.remove(old).unwrap();
            self.defines.insert(new.to_string(), coords);
        }
        if let Some(adj) = self.edges.remove(old) {
            self.edges.insert(new.to_string(), adj);
        }
        if let Some(w) = self.weights.remove(old) {
            self.weights.insert(new.to_string(), w);
        }

        // fix adjacency references
        for (_k, neighbors) in self.edges.iter_mut() {
            if neighbors.remove(old) {
                neighbors.insert(new.to_string());
            }
        }
    }

    fn remove_vertex(&mut self, v: &str) {
        if let Some(nb) = self.edges.get(v).cloned() {
            for n in nb {
                if let Some(ns) = self.edges.get_mut(&n) {
                    ns.remove(v);
                }
            }
        }
        self.edges.remove(v);
        self.defines.remove(v);
        self.weights.remove(v);
    }

    fn add_edge(&mut self, a: &str, b: &str) {
        self.edges
            .entry(a.to_owned())
            .or_default()
            .insert(b.to_owned());
        self.edges
            .entry(b.to_owned())
            .or_default()
            .insert(a.to_owned());
    }

    fn get_vertex(&mut self, preferred: Option<String>) -> String {
        if let Some(x) = preferred {
            x
        } else {
            let name = format!("x{}", self.name_counter);
            self.name_counter += 1;
            name
        }
    }
}

/// Parse a dot-bracket into (position, partner). Partner=0 if unpaired.
/// Returns error string if unmatched parentheses are found.
fn parse_dot_bracket(db: &str) -> Result<Vec<(usize, usize)>, &'static str> {
    let mut stack = Vec::new();
    let mut pairs = vec![(0, 0); db.len()];

    for (i, ch) in db.chars().enumerate() {
        let pos = i + 1;
        match ch {
            '(' => {
                stack.push(pos);
                pairs[i] = (pos, 0);
            }
            ')' => {
                let opening = stack.pop().ok_or("Unbalanced: too many ')'")?;
                pairs[opening - 1] = (opening, pos);
                pairs[i] = (pos, opening);
            }
            '.' => {
                pairs[i] = (pos, 0);
            }
            _ => return Err("Unexpected character in dot-bracket."),
        }
    }
    if !stack.is_empty() {
        return Err("Unbalanced: some '(' left unmatched.");
    }
    Ok(pairs)
}

/// Check if any coordinate in `stem` is +/-1 from the bulge coords
fn any_difference_of_one(
    stem: ((usize, usize), (usize, usize)),
    bulge: (usize, usize),
) -> bool {
    let ((a1, a2), (b1, b2)) = stem;
    let (x, y) = bulge;
    for s in &[a1, a2, b1, b2] {
        for b in &[x, y] {
            if (*s as isize - *b as isize).abs() == 1 {
                return true;
            }
        }
    }
    false
}

/// Sort each half of the 2×2 stem, then sort the two halves
fn sort_stem_tuple(stem: ((usize, usize), (usize, usize))) -> ((usize, usize), (usize, usize)) {
    let part1 = if stem.0.0 <= stem.0.1 {
        (stem.0.0, stem.0.1)
    } else {
        (stem.0.1, stem.0.0)
    };
    let part2 = if stem.1.0 <= stem.1.1 {
        (stem.1.0, stem.1.1)
    } else {
        (stem.1.1, stem.1.0)
    };
    if part1 <= part2 {
        (part1, part2)
    } else {
        (part2, part1)
    }
}
