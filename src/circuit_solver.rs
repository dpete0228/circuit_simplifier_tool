use crate::model::{Component, ComponentKind, Node, Connection, SimplificationResult, Rotation};
use iced::Point;
use std::collections::{HashMap, HashSet, VecDeque};

const NODE_TOLERANCE: f32 = 5.0;

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
    if ohms == 0.0 { return (0.0, "".to_string()); }
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
    if ab_len_sq == 0.0 {
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
    
    comp.center = Point::new((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
    comp.value = len as f64;
    comp.width = len;
    
    // Determine rotation based on dominant axis
    if dy.abs() > dx.abs() {
        comp.rotation = Rotation::Deg90;
    } else {
        comp.rotation = Rotation::Deg0;
    }
}

// Returns a Map: Component ID -> (NetID 1, NetID 2)
// Also returns a list of Nodes for visualization
pub fn identify_nodes(components: &[Component]) -> (Vec<Node>, Vec<Connection>) {
    // 1. Build Adjacency Graph for WIRES only
    // Map: GridPoint -> List of connected GridPoints via Wires
    let mut wire_graph: HashMap<(i32, i32), Vec<(i32, i32)>> = HashMap::new();
    let mut all_points: HashSet<(i32, i32)> = HashSet::new();

    // Map snapped points back to real points for geometric checks
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
    // We pick the "first" point we found for that net as the visual position
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
        let net1 = *point_to_net.get(&snap(p1)).unwrap_or(&0); // Should always exist
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

// --- SOLVER LOGIC ---

pub fn find_next_simplification(components: &[Component]) -> Option<SimplificationResult> {
    // 0. Pre-process: Normalize Wires (Split at T-junctions)
    let (normalized_components, _changed) = split_wires_at_intersections(components);

    let (nodes, connections) = identify_nodes(&normalized_components);
    
    // Map: NetID -> List of Connected Components
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

    for comp in components {
        if comp.id == r1.id {
            let mut new_r1 = comp.clone();
            new_r1.value = new_val;
            new_r1.prefix = new_prefix.clone();
            new_r1.name = format!("Req_{}", r1.id); 
            new_circuit.push(new_r1);
        } else if comp.id == r2.id {
            // Replace R2 with a wire
            let mut new_wire = comp.clone();
            new_wire.kind = ComponentKind::Wire;
            new_wire.name = format!("W_replaced_{}", r2.id);
            let (p1, p2) = comp.endpoints();
            update_wire_geometry(&mut new_wire, p1, p2);
            new_circuit.push(new_wire);
        } else {
            new_circuit.push(comp.clone());
        }
    }
    
    SimplificationResult {
        simplified_circuit: new_circuit,
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
            if val1 + val2 == 0.0 { continue; }
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
    let mut node_map: HashMap<u32, Point> = HashMap::new();
    for n in nodes {
        node_map.insert(n.id, n.position);
    }

    // Iterate over all nodes to find a potential 'Node A' of a Delta triangle
    for (&node_a, neighbors_a) in net_connections {
        // Filter neighbors of A that are Resistors
        let r_neighbors_a: Vec<_> = neighbors_a.iter()
            .filter(|(_, _, kind)| *kind == ComponentKind::Resistor)
            .collect();

        // We need at least 2 resistor neighbors to form two sides of a triangle
        if r_neighbors_a.len() < 2 { continue; }

        // Check pairs of neighbors (Node B and Node C)
        for i in 0..r_neighbors_a.len() {
            for j in (i + 1)..r_neighbors_a.len() {
                let (comp_ab_id, node_b, _) = r_neighbors_a[i];
                let (comp_ac_id, node_c, _) = r_neighbors_a[j];

                // Optimization: Enforce ID order to avoid checking same triangle 3 times
                // Triangle: A-B-C. Only process if A < B < C
                if !(node_a < *node_b && *node_b < *node_c) {
                   continue;
                }

                // Check if Node B and Node C are connected by a Resistor (The 3rd side)
                if let Some(neighbors_b) = net_connections.get(node_b) {
                    if let Some((comp_bc_id, _, _)) = neighbors_b.iter().find(|(_, n, k)| *n == *node_c && *k == ComponentKind::Resistor) {
                        
                        // FOUND DELTA (A-B-C) with resistors Rab, Rac, Rbc
                        let rab = components.iter().find(|c| c.id == *comp_ab_id)?;
                        let rac = components.iter().find(|c| c.id == *comp_ac_id)?;
                        let rbc = components.iter().find(|c| c.id == *comp_bc_id)?;

                        return Some(perform_delta_wye_transform(
                            components, 
                            node_a, *node_b, *node_c,
                            rab, rac, rbc,
                            &node_map
                        ));
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
    let val12 = get_base_ohms(r12);
    let val13 = get_base_ohms(r13);
    let val23 = get_base_ohms(r23);
    
    let denom = val12 + val13 + val23;
    
    // Wye Resistor values
    let val_a = (val12 * val13) / denom;
    let val_b = (val12 * val23) / denom;
    let val_c = (val13 * val23) / denom;

    let (s_a, p_a) = format_to_best_prefix(val_a);
    let (s_b, p_b) = format_to_best_prefix(val_b);
    let (s_c, p_c) = format_to_best_prefix(val_c);

    let explanation = format!(
        "Delta-Wye Transform on Nodes {}, {}, {}.\nOriginal: R{}-{} ({:.1}Ω), R{}-{} ({:.1}Ω), R{}-{} ({:.1}Ω)\nNew Wye: Ra ({:.1}{}), Rb ({:.1}{}), Rc ({:.1}{})",
        n1, n2, n3,
        n1, n2, val12, n1, n3, val13, n2, n3, val23,
        s_a, p_a, s_b, p_b, s_c, p_c
    );

    let mut new_circuit = Vec::new();
    let removed_ids = [r12.id, r13.id, r23.id];

    // 1. Keep all components except the delta resistors
    for comp in components {
        if !removed_ids.contains(&comp.id) {
            new_circuit.push(comp.clone());
        }
    }

    // 2. Determine Geometry using Bounding Box Center for cleaner Manhattan layout
    let p1 = node_map.get(&n1).cloned().unwrap_or(Point::new(0.0,0.0));
    let p2 = node_map.get(&n2).cloned().unwrap_or(Point::new(0.0,0.0));
    let p3 = node_map.get(&n3).cloned().unwrap_or(Point::new(0.0,0.0));
    
    let min_x = p1.x.min(p2.x).min(p3.x);
    let max_x = p1.x.max(p2.x).max(p3.x);
    let min_y = p1.y.min(p2.y).min(p3.y);
    let max_y = p1.y.max(p2.y).max(p3.y);

    // Bounding Box Center often aligns better with grid/Manhattan lines 
    // than the geometric centroid (average of points)
    let center = Point::new(
        (min_x + max_x) / 2.0,
        (min_y + max_y) / 2.0
    );

    let mut id_gen = components.iter().map(|c| c.id).max().unwrap_or(0);

    // Helper to add a "Resistor with leads" to ensure connectivity
    // using Manhattan routing (L-shapes) to avoid diagonal wires.
    let mut add_leg = |node_pos: Point, val: f64, prefix: String, name_suffix: &str| {
        id_gen += 1;
        let r_id = id_gen;
        
        let dx = center.x - node_pos.x;
        let dy = center.y - node_pos.y;
        
        let mut r = Component::new(r_id, ComponentKind::Resistor, Point::new(0.0,0.0));
        r.value = val;
        r.prefix = prefix;
        // Truncate name value to reduce clutter (e.g. "R_3.33" instead of "R_3.333333333")
        r.name = format!("R{}_{:.2}", name_suffix, val);

        let corner: Point;
        let r_start: Point;
        let r_end: Point;

        // Decide dominant axis for the resistor placement
        if dx.abs() >= dy.abs() {
            // Horizontal dominant: Resistor is Horizontal.
            // Path: Node -> (Horizontal Wire) -> Resistor -> (Horizontal Wire) -> Corner -> (Vertical Wire) -> Center
            r.rotation = Rotation::Deg0;
            
            // Place resistor at midpoint of the horizontal span
            let mid_x = (node_pos.x + center.x) / 2.0;
            // Y matches node_pos (first segment is horizontal)
            r.center = Point::new(mid_x, node_pos.y);
            
            // Check if we need to shrink the resistor (if segment is too small)
            let segment_len = (center.x - node_pos.x).abs();
            if segment_len < r.width {
                r.width = segment_len.max(10.0) - 2.0; 
            }

            let (p1, p2) = r.endpoints(); // p1 is left, p2 is right
            
            // Orient based on direction
            if node_pos.x < center.x {
                    r_start = p1; r_end = p2;
            } else {
                    r_start = p2; r_end = p1;
            }
            
            corner = Point::new(center.x, node_pos.y);
        } else {
            // Vertical dominant: Resistor is Vertical.
            // Path: Node -> (Vertical Wire) -> Resistor -> (Vertical Wire) -> Corner -> (Horizontal Wire) -> Center
            r.rotation = Rotation::Deg90;
            
            // Place resistor at midpoint of the vertical span
            let mid_y = (node_pos.y + center.y) / 2.0;
            // X matches node_pos (first segment is vertical)
            r.center = Point::new(node_pos.x, mid_y);
            
            // Check if we need to shrink
            let segment_len = (center.y - node_pos.y).abs();
            if segment_len < r.width {
                r.width = segment_len.max(10.0) - 2.0;
            }

            let (p1, p2) = r.endpoints(); // p1 is top, p2 is bottom
            
            if node_pos.y < center.y {
                r_start = p1; r_end = p2;
            } else {
                r_start = p2; r_end = p1;
            }
            
            corner = Point::new(node_pos.x, center.y);
        }

        // Add the resistor
        new_circuit.push(r);
        
        // Add wires connecting the sequence
        // Wire 1: Node -> Resistor Start
        if (node_pos.x - r_start.x).abs() > 0.1 || (node_pos.y - r_start.y).abs() > 0.1 {
            id_gen += 1;
            let mut w1 = Component::new(id_gen, ComponentKind::Wire, Point::new(0.0,0.0));
            w1.name = format!("W_leg_{}_1", name_suffix);
            update_wire_geometry(&mut w1, node_pos, r_start);
            new_circuit.push(w1);
        }

        // Wire 2: Resistor End -> Corner
        if (r_end.x - corner.x).abs() > 0.1 || (r_end.y - corner.y).abs() > 0.1 {
            id_gen += 1;
            let mut w2 = Component::new(id_gen, ComponentKind::Wire, Point::new(0.0,0.0));
            w2.name = format!("W_leg_{}_2", name_suffix);
            update_wire_geometry(&mut w2, r_end, corner);
            new_circuit.push(w2);
        }

        // Wire 3: Corner -> Center (Only needed if Corner != Center)
        if (corner.x - center.x).abs() > 0.1 || (corner.y - center.y).abs() > 0.1 {
            id_gen += 1;
            let mut w3 = Component::new(id_gen, ComponentKind::Wire, Point::new(0.0,0.0));
            w3.name = format!("W_leg_{}_3", name_suffix);
            update_wire_geometry(&mut w3, corner, center);
            new_circuit.push(w3);
        }
    };

    add_leg(p1, val_a, p_a, "a");
    add_leg(p2, val_b, p_b, "b");
    add_leg(p3, val_c, p_c, "c");

    SimplificationResult {
        simplified_circuit: new_circuit,
        explanation,
    }
}