use crate::model::{Component, ComponentKind, Node, Connection, SimplificationResult, Rotation};
use iced::Point;
use std::collections::{HashMap, HashSet, VecDeque};

const NODE_TOLERANCE: f32 = 5.0;
const R_WIDTH: f32 = 20.0; // Standard resistor width

// --- Local Utility ---
fn is_nearly_zero(value: f64) -> bool {
    value.abs() < 1e-9
}

// --- Helper Functions for Prefixes ---
fn get_multiplier(prefix: &str) -> f64 {
    match prefix {
        "p" => 1e-12,
        "n" => 1e-9,
        "µ" | "u" => 1e-6,
        "m" => 1e-3,
        "k" => 1e3,
        "M" => 1e6,
        "G" => 1e9,
        _ => 1.0,
    }
}

fn get_base_ohms(comp: &Component) -> f64 {
    comp.value * get_multiplier(&comp.prefix)
}

fn format_to_best_prefix(ohms: f64) -> (f64, String) {
    if is_nearly_zero(ohms) { return (0.0, "".to_string()); }
    let abs_ohms = ohms.abs();

    if abs_ohms >= 1e9 { (ohms / 1e9, "G".to_string()) }
    else if abs_ohms >= 1e6 { (ohms / 1e6, "M".to_string()) }
    else if abs_ohms >= 1e3 { (ohms / 1e3, "k".to_string()) }
    else if abs_ohms >= 1.0 { (ohms, "".to_string()) }
    else if abs_ohms >= 1e-3 { (ohms * 1e3, "m".to_string()) }
    else if abs_ohms >= 1e-6 { (ohms * 1e6, "µ".to_string()) }
    else if abs_ohms >= 1e-9 { (ohms * 1e9, "n".to_string()) }
    else { (ohms * 1e12, "p".to_string()) }
}

// --- CORE LOGIC: Netlist Extraction via Flood Fill ---

// Helper to snap float coordinates to integer grid
fn snap(p: Point) -> (i32, i32) {
    (
        (p.x / NODE_TOLERANCE).round() as i32,
        (p.y / NODE_TOLERANCE).round() as i32
    )
}

// Checks if point p lies on the segment a-b (within tolerance)
fn is_on_segment(p: Point, a: Point, b: Point) -> bool {
    let tolerance = NODE_TOLERANCE; 

    // Quick Bounding box check
    let min_x = a.x.min(b.x) - tolerance;
    let max_x = a.x.max(b.x) + tolerance;
    let min_y = a.y.min(b.y) - tolerance;
    let max_y = a.y.max(b.y) + tolerance;

    if p.x < min_x || p.x > max_x || p.y < min_y || p.y > max_y {
        return false;
    }

    // Distance from point to line segment
    let ab_len_sq = (b.x - a.x).powi(2) + (b.y - a.y).powi(2);
    if is_nearly_zero(ab_len_sq as f64) {
        return (p.x - a.x).hypot(p.y - a.y) < tolerance;
    }

    // Project p onto line ab, clamp to segment
    let t = ((p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y)) / ab_len_sq;
    let t_clamped = t.max(0.0).min(1.0);
    
    let closest_x = a.x + t_clamped * (b.x - a.x);
    let closest_y = a.y + t_clamped * (b.y - a.y);

    let dist = (p.x - closest_x).hypot(p.y - closest_y);
    dist < tolerance
}

// Splits wires if another component's endpoint lies on them
fn split_wires_at_intersections(components: &[Component]) -> (Vec<Component>, bool) {
    let mut split_components = Vec::new();
    let mut points: HashSet<(i32, i32)> = HashSet::new();
    let mut next_id = components.iter().map(|c| c.id).max().unwrap_or(0) + 1;
    let mut changed = false;

    // 1. Collect all endpoints from ALL components
    for comp in components {
        let (p1, p2) = comp.endpoints();
        points.insert(snap(p1));
        points.insert(snap(p2));
    }

    // 2. Iterate wires and split if they cross any of these points
    for comp in components {
        if comp.kind != ComponentKind::Wire {
            split_components.push(comp.clone());
            continue;
        }

        let (p1, p2) = comp.endpoints();
        let mut split_points = Vec::new();

        for &pt_grid in &points {
            // Check if this grid point lies on the wire segment (excluding own endpoints)
            if pt_grid == snap(p1) || pt_grid == snap(p2) { continue; }
            
            let pt_real = Point::new(
                pt_grid.0 as f32 * NODE_TOLERANCE,
                pt_grid.1 as f32 * NODE_TOLERANCE
            );
            
            if is_on_segment(pt_real, p1, p2) {
                split_points.push(pt_real);
            }
        }

        if split_points.is_empty() {
            split_components.push(comp.clone());
        } else {
            changed = true;
            // Sort split points by distance from p1 to ensure orderly splitting
            split_points.sort_by(|a, b| {
                let da = (a.x - p1.x).powi(2) + (a.y - p1.y).powi(2);
                let db = (b.x - p1.x).powi(2) + (b.y - p1.y).powi(2);
                da.partial_cmp(&db).unwrap()
            });

            let mut current_start = p1;
            for pt in split_points {
                let mut new_wire = comp.clone();
                new_wire.id = next_id;
                next_id += 1;
                new_wire.name = format!("W_split_{}", new_wire.id);
                update_wire_geometry(&mut new_wire, current_start, pt);
                split_components.push(new_wire);
                current_start = pt;
            }
            // Final segment
            let mut last_wire = comp.clone();
            last_wire.id = next_id;
            next_id += 1;
            last_wire.name = format!("W_split_{}", last_wire.id);
            update_wire_geometry(&mut last_wire, current_start, p2);
            split_components.push(last_wire);
        }
    }

    (split_components, changed)
}

fn update_wire_geometry(comp: &mut Component, p1: Point, p2: Point) {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let len = dx.hypot(dy);
    
    // Set wire center to midpoint
    comp.center = Point::new((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
    comp.value = len as f64;
    comp.width = len;
    
    // Determine rotation based on dominant axis
    if dy.abs() > dx.abs() {
        comp.rotation = Rotation::Deg90;
    } else {
        comp.rotation = Rotation::Deg0;
    }
    
    // Clear overrides as this is a new Manhattan wire segment
    comp.p1_override = None;
    comp.p2_override = None;
}

// Returns a Map: Component ID -> (NetID 1, NetID 2)
// Also returns a list of Nodes for visualization
pub fn identify_nodes(components: &[Component]) -> (Vec<Node>, Vec<Connection>) {
    // 1. Build Adjacency Graph for WIRES only
    let mut wire_graph: HashMap<(i32, i32), Vec<(i32, i32)>> = HashMap::new();
    let mut all_points: HashSet<(i32, i32)> = HashSet::new();

    let mut grid_to_real: HashMap<(i32, i32), Point> = HashMap::new();

    for comp in components {
        let (p1, p2) = comp.endpoints();
        let gp1 = snap(p1);
        let gp2 = snap(p2);
        all_points.insert(gp1);
        all_points.insert(gp2);
        
        grid_to_real.insert(gp1, p1);
        grid_to_real.insert(gp2, p2);

        if comp.kind == ComponentKind::Wire {
            wire_graph.entry(gp1).or_default().push(gp2);
            wire_graph.entry(gp2).or_default().push(gp1);
        }
    }

    // 2. Flood Fill to assign Net IDs
    let mut point_to_net: HashMap<(i32, i32), u32> = HashMap::new();
    let mut next_net_id = 1;

    for point in all_points {
        if point_to_net.contains_key(&point) { continue; }

        // Start a new Net
        let current_net = next_net_id;
        next_net_id += 1;

        // BFS Flood Fill along wires
        let mut queue = VecDeque::new();
        queue.push_back(point);
        point_to_net.insert(point, current_net);

        while let Some(p) = queue.pop_front() {
            if let Some(neighbors) = wire_graph.get(&p) {
                for &neighbor in neighbors {
                    if !point_to_net.contains_key(&neighbor) {
                        point_to_net.insert(neighbor, current_net);
                        queue.push_back(neighbor);
                    }
                }
            }
        }
    }

    // 3. Build Result Structures
    let mut nodes = Vec::new();
    let mut connections = Vec::new();
    
    // Create visualization nodes (one dot per Net)
    let mut net_visual_pos: HashMap<u32, Point> = HashMap::new();
    
    for (pt, net_id) in &point_to_net {
        net_visual_pos.entry(*net_id).or_insert(Point::new(
            pt.0 as f32 * NODE_TOLERANCE, 
            pt.1 as f32 * NODE_TOLERANCE
        ));
    }

    for (net_id, pos) in net_visual_pos {
        nodes.push(Node { id: net_id, position: pos });
    }

    // Create connections for Non-Wire components
    for comp in components {
        if comp.kind == ComponentKind::Wire || comp.kind == ComponentKind::SimplifiedEquivalent {
            continue;
        }

        let (p1, p2) = comp.endpoints();
        let net1 = *point_to_net.get(&snap(p1)).unwrap_or(&0);
        let net2 = *point_to_net.get(&snap(p2)).unwrap_or(&0);

        if net1 != 0 && net2 != 0 && net1 != net2 {
            connections.push(Connection {
                comp_id: comp.id,
                node1_id: net1,
                node2_id: net2,
            });
        }
    }

    (nodes, connections)
}

fn cleanup_dangling_wires(mut components: Vec<Component>) -> Vec<Component> {
    loop {
        // Simple geometric check for loose ends
        let mut counts: HashMap<(i32, i32), usize> = HashMap::new();
        for comp in &components {
            let (p1, p2) = comp.endpoints();
            *counts.entry(snap(p1)).or_default() += 1;
            *counts.entry(snap(p2)).or_default() += 1;
        }

        let mut to_remove = HashSet::new();
        for comp in &components {
            if comp.kind == ComponentKind::Wire {
                let (p1, p2) = comp.endpoints();
                // A wire is dangling if either end has <= 1 connection (itself)
                if counts[&snap(p1)] <= 1 || counts[&snap(p2)] <= 1 {
                    to_remove.insert(comp.id);
                }
            }
        }

        if to_remove.is_empty() { break; }
        components.retain(|c| !to_remove.contains(&c.id));
    }
    components
}

// NEW: Consolidates co-linear wire segments to clean up post-simplification geometry.
fn consolidate_co_linear_wires(mut components: Vec<Component>) -> Vec<Component> {
    let mut components_map: HashMap<u32, Component> = components.into_iter().map(|c| (c.id, c)).collect();
    let mut current_wire_ids: HashSet<u32> = components_map.keys().cloned().filter(|&id| components_map[&id].kind == ComponentKind::Wire).collect();
    let mut merged_wires_to_remove: HashSet<u32> = HashSet::new();

    let tolerance = NODE_TOLERANCE;

    loop {
        let mut did_merge_in_pass = false;
        let wires_to_check: Vec<u32> = current_wire_ids.iter().cloned().collect();

        'merge_check: for i_id in &wires_to_check {
            if merged_wires_to_remove.contains(i_id) { continue; }
            let i_comp = components_map.get(i_id).unwrap();
            let (i_p1, i_p2) = i_comp.endpoints();
            
            for j_id in &wires_to_check {
                if i_id == j_id || merged_wires_to_remove.contains(j_id) { continue; }
                let j_comp = components_map.get(j_id).unwrap();
                let (j_p1, j_p2) = j_comp.endpoints();

                let mut p_join = None; // The internal joining point (shared node)
                let mut p_ext1 = None; // The external end of wire i
                let mut p_ext2 = None; // The external end of wire j
                
                // Check if the snapped endpoints align
                let snap_eq = |p_a: Point, p_b: Point| snap(p_a) == snap(p_b);

                if snap_eq(i_p2, j_p1) { p_join = Some(i_p2); p_ext1 = Some(i_p1); p_ext2 = Some(j_p2); }
                else if snap_eq(i_p2, j_p2) { p_join = Some(i_p2); p_ext1 = Some(i_p1); p_ext2 = Some(j_p1); }
                else if snap_eq(i_p1, j_p1) { p_join = Some(i_p1); p_ext1 = Some(i_p2); p_ext2 = Some(j_p2); }
                else if snap_eq(i_p1, j_p2) { p_join = Some(i_p1); p_ext1 = Some(i_p2); p_ext2 = Some(j_p1); }
                
                if let (Some(_join), Some(ext1), Some(ext2)) = (p_join, p_ext1, p_ext2) {
                    
                    // Collinearity check using cross product approximation: (x_j - x_1) * (y_2 - y_1) - (x_2 - x_1) * (y_j - y_1) = 0
                    // The joining point is one of the four endpoints (e.g., i_p2)
                    let p_a = ext1; // External 1
                    let p_b = ext2; // External 2
                    
                    // We need the *real* joining point to check if the total span is straight through it.
                    // The shared point is implicitly captured by the fact that the two wires are adjacent.
                    // Since the two wires (A-B and B-C) are adjacent, we just check if A, B, and C are collinear.
                    
                    let dx_total = p_b.x - p_a.x;
                    let dy_total = p_b.y - p_a.y;
                    
                    // If the wires were A->B and B->C, we check if A, B, and C are collinear.
                    // We only have the endpoints of the total segment (A and C, which are ext1 and ext2)
                    // and the two constituent components (i and j).
                    
                    // Since we rely on snapped endpoints aligning, we just need to confirm 
                    // that the two segments are co-linear. We can check the angle between the two vectors.

                    // Check if the total span is nearly Manhattan (0 or 90 degrees)
                    let is_manhattan = dx_total.abs() < tolerance || dy_total.abs() < tolerance;

                    // If the individual components are already marked with p1/p2 override, they are diagonal.
                    // We trust that if their snapped endpoints align, they form a straight path.
                    // Since the main problem is eliminating *unnecessary* angle wires,
                    // we simplify if the total span is Manhattan (straight line).
                    
                    if is_manhattan {
                        // Merge! The merged wire will be a Manhattan wire, so we use update_wire_geometry
                        let merged_comp = components_map.get_mut(i_id).unwrap();
                        
                        // Update the geometry of the remaining wire (i_comp)
                        // If the total span is Manhattan, update_wire_geometry handles it correctly.
                        update_wire_geometry(merged_comp, ext1, ext2); 
                        merged_comp.p1_override = None; // Enforce Manhattan structure
                        merged_comp.p2_override = None;
                        
                        // Mark the second wire (j_comp) for removal
                        merged_wires_to_remove.insert(*j_id);
                        current_wire_ids.remove(j_id);
                        did_merge_in_pass = true;
                        break 'merge_check; // Restart outer loop for the modified 'i' wire
                    }
                }
            }
        }
        
        if !did_merge_in_pass { break; }
    }

    // Final cleanup
    components_map.retain(|id, _| !merged_wires_to_remove.contains(id));
    components_map.into_values().collect()
}


// NEW: Combines both cleanup steps
fn cleanup_circuit_geometry(components: Vec<Component>) -> Vec<Component> {
    let mut cleaned = consolidate_co_linear_wires(components);
    cleaned = cleanup_dangling_wires(cleaned);
    cleaned
}


// --- SOLVER LOGIC ---

pub fn find_next_simplification(components: &[Component]) -> Option<SimplificationResult> {
    // 0. Pre-process: Consolidate Wires and Split at T-junctions
    
    // 0.1 First, clean up circuit geometry (merge co-linear wires, remove dangling wires)
    let cleaned_components = cleanup_circuit_geometry(components.to_vec());

    // 0.2 Then, ensure wires are split at all component endpoints for net identification
    // We only need to run this if the geometry changed significantly (i.e., we modified the input)
    let (normalized_components, _changed) = split_wires_at_intersections(&cleaned_components);
    
    // Map: NetID -> List of Connected Components
    let (nodes, connections) = identify_nodes(&normalized_components);
    
    let mut net_connections: HashMap<u32, Vec<(u32, u32, ComponentKind)>> = HashMap::new();
    for conn in &connections {
        let comp = normalized_components.iter().find(|c| c.id == conn.comp_id)?;
        net_connections.entry(conn.node1_id).or_default().push((conn.comp_id, conn.node2_id, comp.kind));
        net_connections.entry(conn.node2_id).or_default().push((conn.comp_id, conn.node1_id, comp.kind));
    }
    
    // 1. Series Check
    for node in &nodes {
        if let Some(neighbors) = net_connections.get(&node.id) {
            // If a Net connects exactly TWO components, they are in Series.
            if neighbors.len() == 2 {
                let (comp1_id, _n1, kind1) = neighbors[0];
                let (comp2_id, _n2, kind2) = neighbors[1];

                if kind1 == ComponentKind::Resistor && kind2 == ComponentKind::Resistor {
                    let r1 = normalized_components.iter().find(|c| c.id == comp1_id)?;
                    let r2 = normalized_components.iter().find(|c| c.id == comp2_id)?;
                    
                    return Some(create_series_simplification(r1, r2, &normalized_components));
                }
            }
        }
    }
    
    // 2. Parallel Check
    if let Some(res) = find_parallel_simplification(&normalized_components, &connections) {
        return Some(res);
    }

    // 3. Delta-Wye (Pi-T) Check
    find_delta_simplification(&normalized_components, &net_connections, &nodes)
}

/// Helper: build a diagonal wire component spanning two points
fn make_diagonal_wire(id: u32, p1: Point, p2: Point) -> Component {
    let mut wire = Component::new(id, ComponentKind::Wire, Point::new(0.0, 0.0));
    wire.p1_override = Some(p1);
    wire.p2_override = Some(p2);
    wire.center = Point::new((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    wire.width = dx.hypot(dy);
    wire.name = format!("W_diag_{}", id);
    wire
}

fn create_series_simplification(
    r1: &Component,
    r2: &Component,
    components: &[Component],
) -> SimplificationResult {
    let val1 = get_base_ohms(r1);
    let val2 = get_base_ohms(r2);
    let total_ohms = val1 + val2;
    let (new_val, new_prefix) = format_to_best_prefix(total_ohms);
    
    let explanation = format!(
        "Series: {} ({}{}Ω) + {} ({}{}Ω) = {:.2}{}Ω",
        r1.name, r1.value, r1.prefix,
        r2.name, r2.value, r2.prefix,
        new_val, new_prefix
    );
    
    let mut new_circuit = Vec::new();

    // Determine which (if any) is a Wye leg
    let r1_is_wye_leg = r1.p1_override.is_some();
    let r2_is_wye_leg = r2.p1_override.is_some();

    // --- Determine the two external endpoints of the series combination ---
    let r1_endpoints = r1.endpoints();
    let r2_endpoints = r2.endpoints();
    
    let mut all_points = vec![
        (snap(r1_endpoints.0), r1_endpoints.0), 
        (snap(r1_endpoints.1), r1_endpoints.1), 
        (snap(r2_endpoints.0), r2_endpoints.0), 
        (snap(r2_endpoints.1), r2_endpoints.1)
    ];

    let mut snapped_counts: HashMap<(i32, i32), usize> = HashMap::new();
    for &(snapped_p, _) in &all_points {
        *snapped_counts.entry(snapped_p).or_default() += 1;
    }

    let external_points: Vec<Point> = all_points.iter()
        .filter(|(snapped_p, _)| snapped_counts[snapped_p] == 1)
        .map(|(_, real_p)| *real_p)
        .collect();

    // external span of the series combo (fallback if weird)
    let (series_p1, series_p2) = if external_points.len() == 2 {
        (external_points[0], external_points[1])
    } else {
        (r1_endpoints.0, r2_endpoints.0) 
    };

    // Decide which resistor becomes the equivalent, and what span it should use
    // If one is a Wye leg, the equivalent should use the non-Wye's endpoints.
    let (eq_resistor_source, eq_other_id, eq_p1, eq_p2) = if r1_is_wye_leg && !r2_is_wye_leg {
        let (p1, p2) = r2_endpoints;
        (r2, r1.id, p1, p2)
    } else if r2_is_wye_leg && !r1_is_wye_leg {
        let (p1, p2) = r1_endpoints;
        (r1, r2.id, p1, p2)
    } else {
        // neither or both are Wye legs -> fall back to external span
        (r1, r2.id, series_p1, series_p2)
    };

    // Next free ID for any new wire(s) we add
    let mut next_id = components.iter().map(|c| c.id).max().unwrap_or(0) + 1;

    for comp in components {
        if comp.id == r1.id || comp.id == r2.id {
            // Handle the equivalent resistor
            if comp.id == eq_resistor_source.id {
                let mut new_r = comp.clone();
                new_r.value = new_val;
                new_r.prefix = new_prefix.clone();
                new_r.name = format!("Req_{}", new_r.id);

                let dx = eq_p2.x - eq_p1.x;
                let dy = eq_p2.y - eq_p1.y;

                if !is_nearly_zero(dx.abs() as f64) && !is_nearly_zero(dy.abs() as f64) {
                    // diagonal span
                    new_r.p1_override = Some(eq_p1);
                    new_r.p2_override = Some(eq_p2);
                    new_r.center = Point::new((eq_p1.x + eq_p2.x) / 2.0, (eq_p1.y + eq_p2.y) / 2.0);
                    new_r.width = R_WIDTH;
                    new_r.value = new_val;

                    // FIX: Calculate Rotation for the diagonal resistor
                    let angle = dy.atan2(dx).to_degrees();
                    let snapped = (angle / 45.0).round() * 45.0;

                    // We snap to nearest 0, 45, 90, 135 degrees (modulo 180)
                    new_r.rotation = match (snapped.rem_euclid(180.0)).round() as i32 {
                        0 => Rotation::Deg0,
                        90 => Rotation::Deg90,
                        _ => Rotation::Deg0,
                    };
                } else {
                    // manhattan span
                    update_wire_geometry(&mut new_r, eq_p1, eq_p2);
                    new_r.value = new_val;
                    // Reset overrides if it was previously a diagonal wye leg
                    new_r.p1_override = None; 
                    new_r.p2_override = None;
                }

                new_circuit.push(new_r);
                continue;
            }

            // If the other series component was a Wye leg, replace it with a diagonal wire
            if (comp.id == r1.id && r1_is_wye_leg) || (comp.id == r2.id && r2_is_wye_leg) {
                let (wp1, wp2) = comp.endpoints();
                let wire = make_diagonal_wire(next_id, wp1, wp2);
                next_id += 1;
                new_circuit.push(wire);
            }
            // Otherwise, this resistor is fully absorbed and removed (no replacement)
            continue;
        } else {
            new_circuit.push(comp.clone());
        }
    }
    
    // The simplified circuit will be cleaned up in the next iteration of find_next_simplification
    SimplificationResult {
        simplified_circuit: cleanup_dangling_wires(new_circuit),
        explanation,
    }
}

fn find_parallel_simplification(
    components: &[Component],
    connections: &[Connection],
) -> Option<SimplificationResult> {
    let mut net_pair_to_resistors: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
    
    for conn in connections {
        let comp = components.iter().find(|c| c.id == conn.comp_id)?;
        if comp.kind == ComponentKind::Resistor {
            // Only allow simplification of Manhattan resistors (no p1_override)
            if comp.p1_override.is_some() { continue; } 
            
            let key = if conn.node1_id < conn.node2_id {
                (conn.node1_id, conn.node2_id)
            } else {
                (conn.node2_id, conn.node1_id)
            };
            net_pair_to_resistors.entry(key).or_default().push(conn.comp_id);
        }
    }
    
    for (_, resistor_ids) in net_pair_to_resistors.iter() {
        if resistor_ids.len() >= 2 {
            let r1 = components.iter().find(|c| c.id == resistor_ids[0])?;
            let r2 = components.iter().find(|c| c.id == resistor_ids[1])?;
            
            let val1 = get_base_ohms(r1);
            let val2 = get_base_ohms(r2);
            if is_nearly_zero(val1 + val2) { continue; }
            let total_ohms = (val1 * val2) / (val1 + val2);
            let (new_val, new_prefix) = format_to_best_prefix(total_ohms);
            
            let explanation = format!(
                "Parallel: {} ({}{}Ω) || {} ({}{}Ω) = {:.2}{}Ω",
                r1.name, r1.value, r1.prefix,
                r2.name, r2.value, r2.prefix,
                new_val, new_prefix
            );
            
            let mut new_circuit = Vec::new();
            for comp in components {
                if comp.id == r1.id {
                    let mut new_r1 = comp.clone();
                    new_r1.value = new_val;
                    new_r1.prefix = new_prefix.clone();
                    new_r1.name = format!("Req_{}", r1.id);
                    new_circuit.push(new_r1);
                } else if comp.id == r2.id {
                    continue; // Delete R2
                } else {
                    new_circuit.push(comp.clone());
                }
            }

            // The simplified circuit will be cleaned up in the next iteration of find_next_simplification
            let cleaned_circuit = cleanup_dangling_wires(new_circuit);
            return Some(SimplificationResult {
                simplified_circuit: cleaned_circuit,
                explanation,
            });
        }
    }
    
    None
}

// --- DELTA-WYE TRANSFORMATION LOGIC ---

fn find_delta_simplification(
    components: &[Component],
    net_connections: &HashMap<u32, Vec<(u32, u32, ComponentKind)>>,
    nodes: &[Node],
) -> Option<SimplificationResult> {
    // Build node position map
    let mut node_map: HashMap<u32, Point> = HashMap::new();
    for n in nodes {
        node_map.insert(n.id, n.position);
    }

    // Build quick lookup: unordered node-pair -> resistor component id
    let mut resistor_map: HashMap<(u32, u32), u32> = HashMap::new();
    for (&node, neighbors) in net_connections {
        for &(comp_id, other_node, kind) in neighbors {
            if kind == ComponentKind::Resistor {
                let comp = components.iter().find(|c| c.id == comp_id)?;
                // Exclude components that are already part of a Wye (diagonal)
                if comp.p1_override.is_some() { continue; }
                
                let key = if node < other_node { (node, other_node) } else { (other_node, node) };
                resistor_map.entry(key).or_insert(comp_id);
            }
        }
    }
    
    let mut node_ids: Vec<u32> = nodes.iter().map(|n| n.id).collect();
    node_ids.sort_unstable();
    
    for a in 0..node_ids.len() {
        for b in (a + 1)..node_ids.len() {
            for c in (b + 1)..node_ids.len() {
                let n1 = node_ids[a];
                let n2 = node_ids[b];
                let n3 = node_ids[c];

                let k12 = (n1.min(n2), n1.max(n2));
                let k13 = (n1.min(n3), n1.max(n3));
                let k23 = (n2.min(n3), n2.max(n3));

                if let (Some(&comp12_id), Some(&comp13_id), Some(&comp23_id)) = (
                    resistor_map.get(&k12),
                    resistor_map.get(&k13),
                    resistor_map.get(&k23),
                ) {
                    if comp12_id != comp13_id && comp12_id != comp23_id && comp13_id != comp23_id {
                        let r12 = components.iter().find(|c| c.id == comp12_id)?;
                        let r13 = components.iter().find(|c| c.id == comp13_id)?;
                        let r23 = components.iter().find(|c| c.id == comp23_id)?;

                        let result = perform_delta_wye_transform(
                            components,
                            n1, n2, n3,
                            r12, r13, r23,
                            &node_map,
                        );
                        return Some(result);
                    }
                }
            }
        }
    }

    None
}

fn perform_delta_wye_transform(
    components: &[Component],
    n1: u32, n2: u32, n3: u32, // Node IDs
    r12: &Component, r13: &Component, r23: &Component, // Components
    node_map: &HashMap<u32, Point>,
) -> SimplificationResult {
    // Get base values in ohms
    let val12 = get_base_ohms(r12);
    let val13 = get_base_ohms(r13);
    let val23 = get_base_ohms(r23);
    
    // Calculate Wye values
    let denom = val12 + val13 + val23;
    
    if is_nearly_zero(denom) || denom < 0.0 {
        panic!("Sum of resistances must be positive for delta-wye transformation");
    }
    
    let val_a = (val12 * val13) / denom; // Ra
    let val_b = (val12 * val23) / denom; // Rb
    let val_c = (val13 * val23) / denom; // Rc
    
    let (s_a, p_a) = format_to_best_prefix(val_a);
    let (s_b, p_b) = format_to_best_prefix(val_b);
    let (s_c, p_c) = format_to_best_prefix(val_c);

    let explanation = format!(
        "Delta-Wye Transform on Nodes {}, {}, {}.\nOriginal: R{}-{} ({:.1}{}), R{}-{} ({:.1}{}), R{}-{} ({:.1}{})\nNew Wye: Ra ({:.1}{}), Rb ({:.1}{}), Rc ({:.1}{}). Using diagonal geometry.",
        n1, n2, n3,
        n1, n2, r12.value, r12.prefix, n1, n3, r13.value, r13.prefix, n2, n3, r23.value, r23.prefix,
        s_a, p_a, s_b, p_b, s_c, p_c
    );

    let mut new_circuit = Vec::new();
    let removed_ids = [r12.id, r13.id, r23.id];

    for comp in components {
        if !removed_ids.contains(&comp.id) {
            new_circuit.push(comp.clone());
        }
    }
    
    let p1 = *node_map.get(&n1).expect("Node N1 not found for delta-wye transform");
    let p2 = *node_map.get(&n2).expect("Node N2 not found for delta-wye transform");
    let p3 = *node_map.get(&n3).expect("Node N3 not found for delta-wye transform");
    
    let center_x = ((p1.x + p2.x + p3.x) / 3.0 / NODE_TOLERANCE).round() * NODE_TOLERANCE;
    let center_y = ((p1.y + p2.y + p3.y) / 3.0 / NODE_TOLERANCE).round() * NODE_TOLERANCE;
    let center = Point::new(center_x, center_y);

    let mut id_gen = components.iter().map(|c| c.id).max().unwrap_or(0);

    let mut add_leg = |node_pos: Point, val: f64, prefix: String, name_suffix: &str| {
        let dx_r = center.x - node_pos.x;
        let dy_r = center.y - node_pos.y;
        let len_r_total = dx_r.hypot(dy_r);
        
        if len_r_total < 0.1 { return; }

        id_gen += 1;
        let mut r = Component::new(id_gen, ComponentKind::Resistor, Point::new(0.0,0.0));
        r.value = val;
        r.prefix = prefix;
        r.name = format!("R{}_{}", name_suffix, id_gen);
        r.width = R_WIDTH; 
        
        r.p1_override = Some(node_pos);
        r.p2_override = Some(center);
        r.center = Point::new((node_pos.x + center.x) / 2.0, (node_pos.y + center.y) / 2.0);
        
        // Set rotation for visualization
        let angle = dy_r.atan2(dx_r).to_degrees();
        let snapped = (angle / 45.0).round() * 45.0;

        r.rotation = match (snapped.rem_euclid(180.0)).round() as i32 {
            0 => Rotation::Deg0,
            90 => Rotation::Deg90,
            _ => Rotation::Deg0,
        };
        
        new_circuit.push(r); 
    };

    add_leg(p1, s_a, p_a, &format!("R{}_center", n1));
    add_leg(p2, s_b, p_b, &format!("R{}_center", n2));
    add_leg(p3, s_c, p_c, &format!("R{}_center", n3));
    
    // The simplified circuit will be cleaned up in the next iteration of find_next_simplification
    SimplificationResult {
        simplified_circuit: cleanup_dangling_wires(new_circuit),
        explanation,
    }
}