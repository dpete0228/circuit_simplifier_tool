use iced::widget::canvas::{self, Frame, Geometry, Path, Stroke, Text, Cache, Cursor, Event};
use iced::{Color, Point, Size, Rectangle, mouse};
use crate::model::{Component, ComponentKind, Rotation, SYMBOL_FONT};
use crate::circuit_solver::identify_nodes;

// Helper to draw the translucent ghost preview component
fn draw_ghost_component(frame: &mut Frame, kind: ComponentKind, center: Point, rotation: Rotation) {
    let mut temp_comp = Component::new(0, kind, center);
    temp_comp.rotation = rotation;
    
    // NOTE: Ghost components are always Manhattan, so no p1_override is set.
    let (p1, p2) = temp_comp.endpoints();
    let color = Color::from_rgba(0.3, 0.5, 1.0, 0.6); // Blue, semi-transparent
    let stroke = Stroke::default().with_width(2.0).with_color(color);
    
    match kind {
        ComponentKind::Resistor => {
            let dx = p2.x - p1.x;
            let dy = p2.y - p1.y;
            let length = (dx * dx + dy * dy).sqrt();
            let ux = dx / length;
            let uy = dy / length;
            
            let lead_len = 10.0;
            let zigzag_len = length - 2.0 * lead_len;
            let zz_start = Point::new(p1.x + ux * lead_len, p1.y + uy * lead_len);
            let zz_end = Point::new(p2.x - ux * lead_len, p2.y - uy * lead_len);
            
            frame.stroke(&Path::line(p1, zz_start), stroke.clone());
            frame.stroke(&Path::line(zz_end, p2), stroke.clone());
            
            let segments = 8; 
            let seg_len = zigzag_len / segments as f32;
            let perp_x = -uy * 5.0; 
            let perp_y = ux * 5.0;
            
            let mut path = canvas::path::Builder::new();
            path.move_to(zz_start);
            
            for i in 1..segments {
                let mid_point = Point::new(
                    zz_start.x + ux * seg_len * i as f32,
                    zz_start.y + uy * seg_len * i as f32,
                );
                let offset = if i % 2 == 1 { 1.0 } else { -1.0 };
                path.line_to(Point::new(mid_point.x + perp_x * offset, mid_point.y + perp_y * offset));
            }
            path.line_to(zz_end);
            
            frame.stroke(&path.build(), stroke.clone());
        }
        ComponentKind::Capacitor => {
            let dx = p2.x - p1.x;
            let dy = p2.y - p1.y;
            let length = (dx * dx + dy * dy).sqrt();
            let ux = dx / length;
            let uy = dy / length;
            let perp_x = -uy;
            let perp_y = ux;
            
            let plate_height = 12.0;
            let gap = 4.0;
            let plate1_center = Point::new(center.x - ux * gap / 2.0, center.y - uy * gap / 2.0);
            let plate2_center = Point::new(center.x + ux * gap / 2.0, center.y + uy * gap / 2.0);
            
            frame.stroke(&Path::line(p1, plate1_center), stroke.clone());
            frame.stroke(&Path::line(plate2_center, p2), stroke.clone());
            
            let plate1_p1 = Point::new(plate1_center.x - perp_x * plate_height / 2.0, plate1_center.y - perp_y * plate_height / 2.0);
            let plate1_p2 = Point::new(plate1_center.x + perp_x * plate_height / 2.0, plate1_center.y + perp_y * plate_height / 2.0);
            frame.stroke(&Path::line(plate1_p1, plate1_p2), Stroke::default().with_width(3.0).with_color(color));
            
            let plate2_p1 = Point::new(plate2_center.x - perp_x * plate_height / 2.0, plate2_center.y - perp_y * plate_height / 2.0);
            let plate2_p2 = Point::new(plate2_center.x + perp_x * plate_height / 2.0, plate2_center.y + perp_y * plate_height / 2.0);
            frame.stroke(&Path::line(plate2_p1, plate2_p2), Stroke::default().with_width(3.0).with_color(color));
        }
        ComponentKind::Inductor => {
            let dx = p2.x - p1.x;
            let dy = p2.y - p1.y;
            let length = (dx * dx + dy * dy).sqrt();
            let ux = dx / length;
            let uy = dy / length;
            
            let lead_len = 10.0;
            let coil_len = length - 2.0 * lead_len;
            let coil_start = Point::new(p1.x + ux * lead_len, p1.y + uy * lead_len);
            let coil_end = Point::new(p2.x - ux * lead_len, p2.y - uy * lead_len);
            
            frame.stroke(&Path::line(p1, coil_start), stroke.clone());
            frame.stroke(&Path::line(coil_end, p2), stroke.clone());
            
            let num_coils = 4;
            let coil_width = coil_len / num_coils as f32;
            let radius = coil_width / 2.0;
            let perp_x = -uy;
            let perp_y = ux;
            
            for i in 0..num_coils {
                let coil_center = Point::new(
                    coil_start.x + ux * (i as f32 + 0.5) * coil_width,
                    coil_start.y + uy * (i as f32 + 0.5) * coil_width,
                );
                
                let mut path = canvas::path::Builder::new();
                let arc_start = Point::new(coil_center.x - ux * radius, coil_center.y - uy * radius);
                path.move_to(arc_start);
                
                for j in 0..=8 {
                    let angle = std::f32::consts::PI * j as f32 / 8.0;
                    let offset_along = -angle.cos() * radius;
                    let offset_perp = angle.sin() * radius;
                    let pt = Point::new(
                        coil_center.x + ux * offset_along + perp_x * offset_perp,
                        coil_center.y + uy * offset_along + perp_y * offset_perp,
                    );
                    path.line_to(pt);
                }
                frame.stroke(&path.build(), stroke.clone());
            }
        }
        ComponentKind::VoltageSource => {
            let dx = p2.x - p1.x;
            let dy = p2.y - p1.y;
            let length = (dx * dx + dy * dy).sqrt();
            let ux = dx / length;
            let uy = dy / length;
            
            let radius = 10.0;
            let lead_len = (length - 2.0 * radius) / 2.0;
            let circle_start = Point::new(p1.x + ux * lead_len, p1.y + uy * lead_len);
            let circle_end = Point::new(p2.x - ux * lead_len, p2.y - uy * lead_len);
            
            frame.stroke(&Path::line(p1, circle_start), stroke.clone());
            frame.stroke(&Path::line(circle_end, p2), stroke.clone());
            frame.stroke(&Path::circle(center, radius), Stroke::default().with_width(2.0).with_color(color));
            
            let plus_size = 4.0;
            let perp_x = -uy;
            let perp_y = ux;

            let plus_center = Point::new(center.x - ux * radius * 0.3, center.y - uy * radius * 0.3);
            let p_start_par = Point::new(plus_center.x - ux * plus_size, plus_center.y - uy * plus_size);
            let p_end_par = Point::new(plus_center.x + ux * plus_size, plus_center.y + uy * plus_size);
            frame.stroke(&Path::line(p_start_par, p_end_par), Stroke::default().with_width(1.5).with_color(color));

            let p_start_perp = Point::new(plus_center.x - perp_x * plus_size, plus_center.y - perp_y * plus_size);
            let p_end_perp = Point::new(plus_center.x + perp_x * plus_size, plus_center.y + perp_y * plus_size);
            frame.stroke(&Path::line(p_start_perp, p_end_perp), Stroke::default().with_width(1.5).with_color(color));
            
            let minus_center = Point::new(center.x + ux * radius * 0.3, center.y + uy * radius * 0.3);
            let m_start_perp = Point::new(minus_center.x - perp_x * plus_size, minus_center.y - perp_y * plus_size);
            let m_end_perp = Point::new(minus_center.x + perp_x * plus_size, minus_center.y + perp_y * plus_size);
            frame.stroke(&Path::line(m_start_perp, m_end_perp), Stroke::default().with_width(1.5).with_color(color));
        }
        _ => {}
    }
}

pub struct CircuitCanvas<'a> {
    pub components: &'a [Component],
    pub ghost: Option<(ComponentKind, Point, Rotation)>,
    pub interactive_wire: Option<(Point, Point)>,
    pub selected_id: Option<u32>,
    pub cache: &'a Cache,
}

impl<'a> canvas::Program<crate::Message> for CircuitCanvas<'a> {
    type State = ();

    fn draw(
        &self,
        _state: &Self::State,
        _theme: &iced::Theme,
        bounds: Rectangle,
        cursor: Cursor,
    ) -> Vec<Geometry> {
        let cached_geometry = self.cache.draw(bounds.size(), |frame| {
            // Draw light gray background
            frame.fill_rectangle(
                Point::ORIGIN,
                bounds.size(),
                Color::from_rgb(0.9, 0.9, 0.9),
            );
            
            // --- 0. Draw Grid dots every 20 units ---
            let grid_spacing = 20.0;
            let bounds_size = bounds.size();
            let grid_color = Color::from_rgb(0.75, 0.75, 0.75);

            for x in (0..=(bounds_size.width / grid_spacing).ceil() as u32).map(|i| i as f32 * grid_spacing) {
                for y in (0..=(bounds_size.height / grid_spacing).ceil() as u32).map(|i| i as f32 * grid_spacing) {
                    frame.fill(
                        &Path::circle(Point::new(x, y), 0.5),
                        grid_color,
                    );
                }
            }

            // --- 1. Draw Components and Wires ---
            // Draw wires first so components overlap them
            for comp in self.components.iter().filter(|c| c.kind == ComponentKind::Wire) {
                draw_component(frame, comp, self.selected_id);
            }
            // Draw other components, including potential diagonal ones
            for comp in self.components.iter().filter(|c| c.kind != ComponentKind::Wire) {
                draw_component(frame, comp, self.selected_id);
            }

            // --- 2. Draw Interactive Wire ---
            if let Some((p1, p2)) = self.interactive_wire {
                frame.stroke(
                    &Path::line(p1, p2),
                    Stroke::default()
                        .with_width(2.0)
                        .with_color(Color::from_rgb(0.5, 0.5, 0.5)),
                );
            }

            // --- 3. Draw Circuit Nodes ---
            let (nodes, _connections) = identify_nodes(self.components);
            for node in nodes {
                frame.fill(
                    &Path::circle(node.position, 3.0),
                    Color::from_rgb(0.9, 0.2, 0.2),
                );
            }
        });

        // Draw the dynamic ghost component
        let mut ghost_frame = Frame::new(bounds.size());
        if let Some((kind, _, rotation)) = &self.ghost {
            if let Some(cursor_pos) = cursor.position_in(&bounds) {
                let snapped = Point::new(
                    (cursor_pos.x / 20.0).round() * 20.0,
                    (cursor_pos.y / 20.0).round() * 20.0,
                );
                draw_ghost_component(&mut ghost_frame, *kind, snapped, *rotation);
            }
        }
        let ghost_geometry = ghost_frame.into_geometry();

        vec![cached_geometry, ghost_geometry]
    }
    
    fn update(
        &self,
        _state: &mut Self::State,
        event: Event,
        bounds: Rectangle,
        cursor: Cursor,
    ) -> (canvas::event::Status, Option<crate::Message>) {
        match event {
            Event::Mouse(mouse::Event::ButtonPressed(mouse::Button::Left)) => {
                if let Some(position) = cursor.position_in(&bounds) {
                    return (
                        canvas::event::Status::Captured,
                        Some(crate::Message::CanvasClicked(position)),
                    );
                }
            }
            _ => {}
        }
        (canvas::event::Status::Ignored, None)
    }
}

// Helper to draw the actual components
fn draw_component(frame: &mut Frame, comp: &Component, selected_id: Option<u32>) {
    // Endpoints are guaranteed to be the start/end of the physical component path, 
    // whether calculated (Manhattan) or overridden (Diagonal).
    let (p1_end, p2_end) = comp.endpoints();
    let is_selected = selected_id == Some(comp.id);
    
    let color = if is_selected {
        Color::from_rgb(1.0, 0.5, 0.0)
    } else {
        Color::from_rgb(0.1, 0.1, 0.1)
    };
    
    let stroke = Stroke::default().with_width(2.0).with_color(color);
    
    match comp.kind {
        ComponentKind::Wire => {
            frame.stroke(&Path::line(p1_end, p2_end), stroke);
        }
        ComponentKind::Resistor | ComponentKind::Capacitor | 
        ComponentKind::Inductor | ComponentKind::VoltageSource => {
            
            // Calculate the component's unit vector based on its true endpoints
            let dx = p2_end.x - p1_end.x;
            let dy = p2_end.y - p1_end.y;
            let length = (dx * dx + dy * dy).sqrt();
            let (ux, uy) = if length > 0.001 {
                (dx / length, dy / length)
            } else {
                (1.0, 0.0)
            };
            
            // Component geometry is drawn along the vector defined by (p1_end, p2_end)
            
            let mut draw_leads = |p1: Point, p2: Point, comp_len: f32| -> (Point, Point) {
                let lead_len = 10.0;
                let center_len = comp_len - 2.0 * lead_len;
                
                if center_len <= 0.0 {
                    // Draw full wire if too short for symbol
                    frame.stroke(&Path::line(p1, p2), stroke.clone());
                    return (p1, p2); // Return full length leads
                }
                
                let zz_start = Point::new(p1.x + ux * lead_len, p1.y + uy * lead_len);
                let zz_end = Point::new(p2.x - ux * lead_len, p2.y - uy * lead_len);
                
                // Draw leads
                frame.stroke(&Path::line(p1, zz_start), stroke.clone());
                frame.stroke(&Path::line(zz_end, p2), stroke.clone());
                
                (zz_start, zz_end)
            };
            
            match comp.kind {
                ComponentKind::Resistor => {
                    let (zz_start, zz_end) = draw_leads(p1_end, p2_end, length);
                    let zigzag_len = (zz_end.x - zz_start.x).hypot(zz_end.y - zz_start.y);
                    
                    let segments = 8;
                    let seg_len = zigzag_len / segments as f32;
                    let perp_x = -uy * 5.0; 
                    let perp_y = ux * 5.0;
                    
                    let mut path = canvas::path::Builder::new();
                    path.move_to(zz_start);
                    
                    for i in 1..segments {
                        let mid_point = Point::new(
                            zz_start.x + ux * seg_len * i as f32,
                            zz_start.y + uy * seg_len * i as f32,
                        );
                        let offset = if i % 2 == 1 { 1.0 } else { -1.0 };
                        path.line_to(Point::new(mid_point.x + perp_x * offset, mid_point.y + perp_y * offset));
                    }
                    path.line_to(zz_end);
                    frame.stroke(&path.build(), stroke.clone());
                }
                ComponentKind::Capacitor => {
                    let lead_len = (length - 4.0) / 2.0; // Shorten leads slightly more
                    if lead_len > 0.0 {
                        let plate1_center = Point::new(p1_end.x + ux * lead_len, p1_end.y + uy * lead_len);
                        let plate2_center = Point::new(p2_end.x - ux * lead_len, p2_end.y - uy * lead_len);
                        
                        frame.stroke(&Path::line(p1_end, plate1_center), stroke.clone());
                        frame.stroke(&Path::line(plate2_center, p2_end), stroke.clone());
                        
                        let perp_x = -uy;
                        let perp_y = ux;
                        let plate_height = 12.0;
                        
                        let plate1_p1 = Point::new(plate1_center.x - perp_x * plate_height / 2.0, plate1_center.y - perp_y * plate_height / 2.0);
                        let plate1_p2 = Point::new(plate1_center.x + perp_x * plate_height / 2.0, plate1_center.y + perp_y * plate_height / 2.0);
                        frame.stroke(&Path::line(plate1_p1, plate1_p2), Stroke::default().with_width(3.0).with_color(color));
                        
                        let plate2_p1 = Point::new(plate2_center.x - perp_x * plate_height / 2.0, plate2_center.y - perp_y * plate_height / 2.0);
                        let plate2_p2 = Point::new(plate2_center.x + perp_x * plate_height / 2.0, plate2_center.y + perp_y * plate_height / 2.0);
                        frame.stroke(&Path::line(plate2_p1, plate2_p2), Stroke::default().with_width(3.0).with_color(color));
                    } else {
                        frame.stroke(&Path::line(p1_end, p2_end), stroke.clone());
                    }
                }
                ComponentKind::Inductor => {
                    let (coil_start, coil_end) = draw_leads(p1_end, p2_end, length);
                    let coil_len = (coil_end.x - coil_start.x).hypot(coil_end.y - coil_start.y);

                    let num_coils = 4;
                    let coil_width = coil_len / num_coils as f32;
                    let radius = coil_width / 2.0;
                    let perp_x = -uy;
                    let perp_y = ux;
                    
                    for i in 0..num_coils {
                        let coil_center = Point::new(
                            coil_start.x + ux * (i as f32 + 0.5) * coil_width,
                            coil_start.y + uy * (i as f32 + 0.5) * coil_width,
                        );
                        
                        let mut path = canvas::path::Builder::new();
                        let arc_start = Point::new(coil_center.x - ux * radius, coil_center.y - uy * radius);
                        path.move_to(arc_start);
                        
                        for j in 0..=8 {
                            let angle = std::f32::consts::PI * j as f32 / 8.0;
                            let offset_along = -angle.cos() * radius;
                            let offset_perp = angle.sin() * radius;
                            let pt = Point::new(
                                coil_center.x + ux * offset_along + perp_x * offset_perp,
                                coil_center.y + uy * offset_along + perp_y * offset_perp,
                            );
                            path.line_to(pt);
                        }
                        
                        frame.stroke(&path.build(), stroke.clone());
                    }
                }
                ComponentKind::VoltageSource => {
                    let radius = 10.0;
                    let lead_len = (length - 2.0 * radius) / 2.0;
                    if lead_len > 0.0 {
                        let circle_start = Point::new(p1_end.x + ux * lead_len, p1_end.y + uy * lead_len);
                        let circle_end = Point::new(p2_end.x - ux * lead_len, p2_end.y - uy * lead_len);
                        
                        frame.stroke(&Path::line(p1_end, circle_start), stroke.clone());
                        frame.stroke(&Path::line(circle_end, p2_end), stroke.clone());
                    } else {
                        frame.stroke(&Path::line(p1_end, p2_end), stroke.clone());
                    }
                    
                    frame.stroke(&Path::circle(comp.center, radius), Stroke::default().with_width(2.0).with_color(color));
                    
                    let plus_size = 4.0;
                    let perp_x = -uy;
                    let perp_y = ux;
                    
                    let plus_center = Point::new(comp.center.x - ux * radius * 0.3, comp.center.y - uy * radius * 0.3);
                    let p_start_par = Point::new(plus_center.x - ux * plus_size, plus_center.y - uy * plus_size);
                    let p_end_par = Point::new(plus_center.x + ux * plus_size, plus_center.y + uy * plus_size);
                    frame.stroke(&Path::line(p_start_par, p_end_par), Stroke::default().with_width(1.5).with_color(color));

                    let p_start_perp = Point::new(plus_center.x - perp_x * plus_size, plus_center.y - perp_y * plus_size);
                    let p_end_perp = Point::new(plus_center.x + perp_x * plus_size, plus_center.y + perp_y * plus_size);
                    frame.stroke(&Path::line(p_start_perp, p_end_perp), Stroke::default().with_width(1.5).with_color(color));
                    
                    let minus_center = Point::new(comp.center.x + ux * radius * 0.3, comp.center.y + uy * radius * 0.3);
                    let m_start_perp = Point::new(minus_center.x - perp_x * plus_size, minus_center.y - perp_y * plus_size);
                    let m_end_perp = Point::new(minus_center.x + perp_x * plus_size, minus_center.y + perp_y * plus_size);
                    frame.stroke(&Path::line(m_start_perp, m_end_perp), Stroke::default().with_width(1.5).with_color(color));
                }
                _ => {} // Already handled outside this match
            }
            
            draw_component_labels(frame, comp, color);
        }
        _ => {
            // Should be covered by the match block above or Wire.
        }
    }
}

fn draw_component_labels(frame: &mut Frame, comp: &Component, _color: Color) {
    let (p1, p2) = comp.endpoints();

    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let length = (dx * dx + dy * dy).hypot(1.0); // Use hypot with 1.0 floor to avoid div by zero if length is 0

    let (ux, uy) = if length > 0.001 {
        (dx / length, dy / length)
    } else {
        (1.0, 0.0)
    };

    // Calculate perpendicular offset vector for labels
    let perp_x = -uy;
    let perp_y = ux;

    let offset = 17.0; 

    // Calculate position for combined Name and Value text, offset perpendicular to component axis
    let combined_pos = Point::new(
        comp.center.x + perp_x * offset,
        comp.center.y + perp_y * offset,
    );

    // Format the value concisely (Name and Value combined)
    let value_str = format!("{:.2}{}{}", comp.value, comp.prefix, comp.unit_symbol());
    let combined_text_content = format!("{} {}", comp.name, value_str);

    let combined_text = Text {
        content: combined_text_content,
        position: combined_pos,
        size: 10.0.into(),
        color: Color::from_rgb(0.2, 0.2, 0.2),
        font: SYMBOL_FONT,
        horizontal_alignment: iced::alignment::Horizontal::Center,
        vertical_alignment: iced::alignment::Vertical::Center,
        ..Text::default()
    };
    frame.fill_text(combined_text);
    
    // Remove the separate value display to reduce clutter (it's now combined)
}