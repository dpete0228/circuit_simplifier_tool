use iced::{Font, Point};

// Define the custom font for the application
pub const SYMBOL_FONT: Font = Font::Default;

// --- Component Definitions ---

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rotation {
    Deg0,
    Deg90,
    Deg180,
    Deg270,
    
}

impl Rotation {
    pub fn rotate_cw(self) -> Self {
        match self {
            Rotation::Deg0 => Rotation::Deg90,
            Rotation::Deg90 => Rotation::Deg180,
            Rotation::Deg180 => Rotation::Deg270,
            Rotation::Deg270 => Rotation::Deg0,
            
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComponentKind {
    Resistor,
    Capacitor,
    Inductor,
    VoltageSource,
    Wire,
    SimplifiedEquivalent, 
}

impl ComponentKind {
    pub fn symbol(&self) -> char {
        match self {
            ComponentKind::Resistor => 'R',
            ComponentKind::Capacitor => 'C',
            ComponentKind::Inductor => 'L',
            ComponentKind::VoltageSource => 'V',
            ComponentKind::Wire => '/',
            ComponentKind::SimplifiedEquivalent => 'E',
        }
    }
}

#[derive(Debug, Clone)]
pub struct Component {
    pub id: u32,
    pub kind: ComponentKind,
    pub center: Point,
    pub rotation: Rotation,
    pub name: String,
    pub value: f64,
    pub prefix: String,
    pub width: f32,
    pub height: f32,
    /// Overrides calculated endpoints for solver-generated diagonal components.
    pub p1_override: Option<Point>,
    /// Overrides calculated endpoints for solver-generated diagonal components.
    pub p2_override: Option<Point>,
}

impl Component {
    pub fn new(id: u32, kind: ComponentKind, center: Point) -> Self {
        let default_name = format!("{}{}", kind.symbol(), id);
        let default_value = match kind {
            ComponentKind::Resistor => 10.0,
            ComponentKind::Capacitor => 1.0,
            ComponentKind::Inductor => 1.0,
            ComponentKind::VoltageSource => 5.0,
            ComponentKind::Wire => 40.0,
            ComponentKind::SimplifiedEquivalent => 0.0,
        };

        let (width, height) = match kind {
            ComponentKind::Wire => (default_value as f32, 0.0),
            _ => (40.0, 40.0),
        };

        Component {
            id,
            kind,
            center,
            rotation: Rotation::Deg0,
            name: default_name,
            value: default_value,
            prefix: String::from(""),
            width,
            height,
            p1_override: None,
            p2_override: None,
        }
    }

    pub fn endpoints(&self) -> (Point, Point) {
        // Use explicit overrides if set (for solver-generated diagonal components)
        if let (Some(p1), Some(p2)) = (self.p1_override, self.p2_override) {
            return (p1, p2);
        }

        let half_len = if self.kind == ComponentKind::Wire {
            self.value as f32 / 2.0
        } else {
            self.width / 2.0 
        };
        
        let (dx, dy) = match self.rotation {
            Rotation::Deg0 | Rotation::Deg180 => (half_len, 0.0),
            Rotation::Deg90 | Rotation::Deg270 => (0.0, half_len),
        };

        let start = Point::new(self.center.x - dx, self.center.y - dy);
        let end = Point::new(self.center.x + dx, self.center.y + dy);

        (start, end)
    }

    pub fn contains(&self, p: Point) -> bool {
        let (p1, p2) = self.endpoints();
        let min_x = p1.x.min(p2.x);
        let max_x = p1.x.max(p2.x);
        let min_y = p1.y.min(p2.y);
        let max_y = p1.y.max(p2.y);
        
        let tolerance = 10.0;
        
        // This geometric check works correctly even for diagonal components
        p.x >= min_x - tolerance && p.x <= max_x + tolerance && 
        p.y >= min_y - tolerance && p.y <= max_y + tolerance
    }
    
    pub fn unit_symbol(&self) -> &'static str {
        match self.kind {
            ComponentKind::Resistor => "Ω",
            ComponentKind::Capacitor => "F",
            ComponentKind::Inductor => "H",
            ComponentKind::VoltageSource => "V",
            ComponentKind::Wire => "px",
            ComponentKind::SimplifiedEquivalent => "Ω/F/H/V",
        }
    }
}

// --- Solver Definitions ---

#[derive(Debug, Clone)]
pub struct Node {
    pub position: Point,
    pub id: u32,
}

#[derive(Debug, Clone)]
pub struct SimplificationResult {
    pub simplified_circuit: Vec<Component>,
    pub explanation: String,
}

#[derive(Debug, Clone)]
pub struct Connection {
    pub comp_id: u32,
    pub node1_id: u32,
    pub node2_id: u32,
}