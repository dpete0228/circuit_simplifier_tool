use crate::model::{Component, ComponentKind, Rotation, SYMBOL_FONT};
use iced::widget::{column, container, pick_list, row, text, text_input, button, canvas as canvas_widget};
use iced::widget::canvas::Cache;
use iced::{
    executor, keyboard, mouse, Application, Command, Element,
    Event, Length, Point, Settings, Subscription, Theme,
};
use tools::SelectedTool;

mod canvas;
mod model;
mod tools;
mod circuit_solver;

use canvas::CircuitCanvas;

pub fn main() -> iced::Result {
    let mut settings = iced::Settings::default();
    settings.default_font = Some(include_bytes!("../fonts/NotoSans-Regular.ttf"));
    CircuitEditor::run(settings)
}

struct CircuitEditor {
    components: Vec<Component>, 
    next_id: u32,
    selected_tool: SelectedTool,
    temp_rotation: Rotation,
    mouse_position: Option<Point>,
    wire_start: Option<Point>,
    selected_id: Option<u32>,
    temp_value_input: String,
    
    // History Fields
    circuit_history: Vec<(Vec<Component>, String)>,
    current_step: usize,
    
    // Canvas cache
    canvas_cache: Cache,
}

#[derive(Debug, Clone)]
pub enum Message {
    KeyPressed(keyboard::KeyCode),
    MouseMoved(Option<Point>),
    CanvasClicked(Point),
    Ignore,
    UpdateName(String),
    UpdateValue(String),
    UpdatePrefix(String),
    Simplify,
    NextStep,
    PreviousStep,
}

const GRID_SIZE: f32 = 20.0;

fn snap(point: Point) -> Point {
    let x = (point.x / GRID_SIZE).round() * GRID_SIZE;
    let y = (point.y / GRID_SIZE).round() * GRID_SIZE;
    Point::new(x, y)
}

impl CircuitEditor {
    fn sync_circuit_state(&mut self) {
        if self.current_step < self.circuit_history.len() {
            self.components = self.circuit_history[self.current_step].0.clone();
            
            if let Some(id) = self.selected_id {
                if !self.components.iter().any(|c| c.id == id) {
                    self.selected_id = None;
                    self.temp_value_input.clear();
                } else {
                    if let Some(comp) = self.components.iter().find(|c| c.id == id) {
                        self.temp_value_input = comp.value.to_string();
                    }
                }
            }
            
            // Clear cache when circuit changes
            self.canvas_cache.clear();
        }
    }
    
    fn commit_user_change(&mut self) {
        self.circuit_history.truncate(self.current_step + 1);
        
        let new_state = self.components.clone();
        
        let explanation = if new_state.is_empty() {
            "Circuit cleared.".to_string()
        } else {
            "User modification (Add/Edit/Delete component).".to_string()
        };
        
        self.circuit_history.push((new_state, explanation));
        self.current_step += 1;
        
        self.sync_circuit_state();
    }
}

impl Application for CircuitEditor {
    type Executor = executor::Default;
    type Message = Message;
    type Flags = ();
    type Theme = Theme;

    fn new(_flags: Self::Flags) -> (Self, Command<Self::Message>) {
        let initial_components = Vec::new();

        (
            Self {
                components: initial_components.clone(), 
                next_id: 1,
                selected_tool: SelectedTool::None,
                temp_rotation: Rotation::Deg0,
                mouse_position: None,
                wire_start: None,
                selected_id: None,
                temp_value_input: String::new(),
                
                circuit_history: vec![(initial_components, "Initial Circuit".to_string())],
                current_step: 0,
                
                canvas_cache: Cache::default(),
            },
            Command::none(),
        )
    }

    fn title(&self) -> String {
        "Circuit Simplifier".into()
    }

    fn theme(&self) -> Self::Theme {
        Theme::Dark
    }

    fn update(&mut self, message: Self::Message) -> Command<Self::Message> {
        match message {
            Message::UpdateName(new_name) => {
                if let Some(id) = self.selected_id {
                    if let Some(comp) = self.components.iter_mut().find(|c| c.id == id) {
                        comp.name = new_name;
                        self.commit_user_change();
                    }
                }
            }
            Message::UpdateValue(new_str) => {
                self.temp_value_input = new_str.clone();
                if let Some(id) = self.selected_id {
                    if let Some(comp) = self.components.iter_mut().find(|c| c.id == id) {
                        if let Ok(val) = new_str.parse::<f64>() {
                            comp.value = val;
                            self.commit_user_change();
                        }
                    }
                }
            }
            Message::UpdatePrefix(new_prefix) => {
                if let Some(id) = self.selected_id {
                    if let Some(comp) = self.components.iter_mut().find(|c| c.id == id) {
                        comp.prefix = new_prefix;
                        self.commit_user_change();
                    }
                }
            }
            Message::KeyPressed(key) => {
                match key {
                    keyboard::KeyCode::R => {
                        self.selected_tool = SelectedTool::Resistor;
                    }
                    keyboard::KeyCode::C => {
                        self.selected_tool = SelectedTool::Capacitor;
                    }
                    keyboard::KeyCode::L => {
                        self.selected_tool = SelectedTool::Inductor;
                    }
                    keyboard::KeyCode::W => {
                        self.selected_tool = SelectedTool::Wire;
                    }
                    keyboard::KeyCode::V => {
                        self.selected_tool = SelectedTool::VoltageSource;
                    }
                    keyboard::KeyCode::Space => {
                        self.temp_rotation = self.temp_rotation.rotate_cw();
                    }
                    keyboard::KeyCode::Delete | keyboard::KeyCode::Backspace => {
                        if let Some(id) = self.selected_id {
                            self.components.retain(|c| c.id != id);
                            self.selected_id = None;
                            self.temp_value_input.clear();
                            // commit_user_change calls sync_circuit_state, which clears the cache
                            self.commit_user_change();
                        }
                    }
                    keyboard::KeyCode::Escape => {
                        // Clear the cache if we are clearing an active selection
                        if self.selected_id.is_some() {
                            self.canvas_cache.clear();
                        }

                        self.selected_tool = SelectedTool::None;
                        self.wire_start = None;
                        self.selected_id = None;
                        self.temp_value_input.clear();
                    }
                    _ => {}
                }
            },
            Message::MouseMoved(pos_option) => {
                self.mouse_position = pos_option;
            }
            Message::CanvasClicked(raw_pos) => {
                let pos = snap(raw_pos);
                match self.selected_tool {
                    SelectedTool::None => {
                        let old_selected_id = self.selected_id; // Store previous selection
                        let hit = self.components.iter().rev().find(|c| c.contains(raw_pos));
                        if let Some(comp) = hit {
                            self.selected_id = Some(comp.id);
                            self.temp_value_input = comp.value.to_string();
                        } else {
                            self.selected_id = None;
                            self.temp_value_input.clear();
                        }
                        
                        // FIX: If the selection status changed, clear the cache to force redraw
                        if old_selected_id != self.selected_id {
                            self.canvas_cache.clear();
                        }
                    }
                    SelectedTool::Wire => {
                        if let Some(start) = self.wire_start {
                            let dx = pos.x - start.x;
                            let dy = pos.y - start.y;
                            let (rotation, length, end_pos) = if dx.abs() > dy.abs() {
                                (Rotation::Deg0, dx.abs(), Point::new(pos.x, start.y))
                            } else {
                                (Rotation::Deg90, dy.abs(), Point::new(start.x, pos.y))
                            };
                            if length > 0.1 {
                                let mid_x = (start.x + end_pos.x) / 2.0;
                                let mid_y = (start.y + end_pos.y) / 2.0;
                                let center = Point::new(mid_x, mid_y);
                                let mut comp = Component::new(self.next_id, ComponentKind::Wire, center);
                                comp.rotation = rotation;
                                comp.value = length as f64;
                                self.components.push(comp);
                                self.next_id += 1;
                            }
                            self.wire_start = None;
                            self.selected_id = None;
                            self.commit_user_change();
                        } else {
                            self.wire_start = Some(pos);
                        }
                    }
                    _ => {
                        let kind = match self.selected_tool {
                            SelectedTool::Resistor => ComponentKind::Resistor,
                            SelectedTool::Capacitor => ComponentKind::Capacitor,
                            SelectedTool::Inductor => ComponentKind::Inductor,
                            SelectedTool::VoltageSource => ComponentKind::VoltageSource,
                            _ => ComponentKind::Resistor,
                        };
                        let mut comp = Component::new(self.next_id, kind, pos);
                        comp.rotation = self.temp_rotation;
                        self.components.push(comp);
                        self.next_id += 1;
                        self.selected_id = None;
                        self.commit_user_change();
                    }
                }
            }
            Message::Simplify => {
                let (current_circuit, _) = self.circuit_history[self.current_step].clone(); 
                
                if let Some(result) = circuit_solver::find_next_simplification(&current_circuit) {
                    self.circuit_history.truncate(self.current_step + 1);
                    self.circuit_history.push((result.simplified_circuit, result.explanation));
                    self.current_step += 1;
                    self.sync_circuit_state();
                }
            }
            
            Message::NextStep => {
                if self.current_step < self.circuit_history.len() - 1 {
                    self.current_step += 1;
                    self.sync_circuit_state();
                }
            }

            Message::PreviousStep => {
                if self.current_step > 0 {
                    self.current_step -= 1;
                    self.sync_circuit_state();
                }
            }
            
            Message::Ignore => {}
        }
        Command::none()
    }

    fn view(&self) -> Element<Message> {
        let (current_components, current_explanation) = &self.circuit_history[self.current_step];
        let snapped_mouse = self.mouse_position.map(snap);

        let ghost = if let SelectedTool::Wire = self.selected_tool {
            None
        } else {
            match self.selected_tool {
                SelectedTool::None => None,
                _ => snapped_mouse.map(|p| {
                    let kind = match self.selected_tool {
                        SelectedTool::Resistor => ComponentKind::Resistor,
                        SelectedTool::Capacitor => ComponentKind::Capacitor,
                        SelectedTool::Inductor => ComponentKind::Inductor,
                        SelectedTool::VoltageSource => ComponentKind::VoltageSource,
                        _ => ComponentKind::Resistor,
                    };
                    (kind, p, self.temp_rotation)
                }),
            }
        };

       

        let interactive_wire = if let (SelectedTool::Wire, Some(start), Some(current)) =
            (&self.selected_tool, self.wire_start, snapped_mouse)
        {
            let dx = current.x - start.x;
            let dy = current.y - start.y;
            let end = if dx.abs() > dy.abs() {
                Point::new(current.x, start.y)
            } else {
                Point::new(start.x, current.y)
            };
            Some((start, end))
        } else {
            None
        };

        let canvas_layer = canvas_widget::Canvas::new(CircuitCanvas {
            components: current_components,
            ghost,
            interactive_wire,
            selected_id: self.selected_id,
            cache: &self.canvas_cache,
        })
        .width(Length::Fill)
        .height(Length::Fill);

        let prev_button = if self.current_step > 0 {
            button(text("Previous")).on_press(Message::PreviousStep)
        } else {
            button(text("Previous"))
        };
        
        let simplify_button = if self.current_step == self.circuit_history.len() - 1 && matches!(self.selected_tool, SelectedTool::None) {
            button(text("Simplify")).on_press(Message::Simplify)
        } else {
            button(text("Simplify"))
        };
            
        let next_button = if self.current_step < self.circuit_history.len() - 1 {
            button(text("Next Step")).on_press(Message::NextStep)
        } else {
            button(text("Next Step"))
        };
            
        let nav_controls = row![
            prev_button.width(Length::Fill),
            simplify_button.width(Length::Fill),
            next_button.width(Length::Fill)
        ]
        .spacing(10);
        
        let main_content = if let Some(id) = self.selected_id {
            if let Some(comp) = current_components.iter().find(|c| c.id == id) {
                let prefixes = vec!["M", "k", "", "m", "µ", "n", "p"];
                let selected_prefix = prefixes.iter().find(|&&p| p == comp.prefix).copied();

                column![
                    text("Properties").size(20),
                    text("Name:"),
                    text_input("Name", &comp.name).on_input(Message::UpdateName),
                    text("Value:"),
                    row![
                        text_input("Value", &self.temp_value_input)
                            .on_input(Message::UpdateValue)
                            
                            .width(Length::Fill),
                        pick_list(
                            prefixes,
                            selected_prefix,
                            |s| Message::UpdatePrefix(s.to_string())
                        )
                        
                        .width(Length::Fixed(50.0)),
                        text(comp.unit_symbol())
                            .size(18)
                            ,
                    ]
                    .spacing(5),
                    text(format!("Kind: {:?}", comp.kind)).size(12),
                ]
                .spacing(10)
            } else {
                 column![
                    text("Component not visible in current step.").size(14),
                ]
                .spacing(10)
            }
        } else {
            column![
                text("Circuit Simplifier").size(20),
                text("Select a component to edit.").size(14),
                text("R: Resistor"),
                text("C: Capacitor"),
                text("L: Inductor"),
                text("V: Voltage"),
                text("W: Wire"),
                text("Space: Rotate"),
                text("Del: Delete"),
                text("Esc: Cancel"),
            ]
            .spacing(10)
        };

        let sidebar_content = column![
            text("--- History ---").size(12),
            text(format!("Step {}/{}", self.current_step + 1, self.circuit_history.len())).size(16),
            nav_controls,
            text("--- Explanation ---").size(12),
            text(current_explanation).size(14),
            text("---").size(12),
            main_content,
        ]
        .spacing(10);
        
        let sidebar = container(sidebar_content)
            .padding(20)
            .width(Length::Fixed(300.0));

        row![canvas_layer, sidebar].into()
    }

    fn subscription(&self) -> Subscription<Message> {
        iced::subscription::events().map(|event| match event {
            Event::Keyboard(key_event) => {
                if let keyboard::Event::KeyPressed { key_code, .. } = key_event {
                    Message::KeyPressed(key_code)
                } else {
                    Message::Ignore
                }
            }
            Event::Mouse(mouse::Event::CursorMoved { position }) => {
                Message::MouseMoved(Some(position))
            }
            _ => Message::Ignore,
        })
    }
}