use bevy::prelude::*;

#[derive(Resource)]
pub struct Tissue {
    pub segments_per_ring: Vec<usize>,
    pub cell_centers: Vec<Vec<Vec2>>, // [ring][seg]
    pub auxin: Vec<Vec<f32>>,         // [ring][seg]
    pub next_auxin: Vec<Vec<f32>>,    // [ring][seg]
    pub pin_dir: Vec<Vec<Vec2>>,      // [ring][seg]
    pub cell_entities: Vec<Vec<Entity>>, // entity per cell for color updates
}

#[derive(Component)]
pub struct CellVisual;

#[derive(Resource, Default)]
pub struct PrimordiumList {
    pub items: Vec<Primordium>,
    pub steps_since_last: u32,
}

#[derive(Clone)]
pub struct Primordium {
    pub ring_index: usize,
    pub angle: f32, // radians
    pub sink_strength: f32,
    pub radial_velocity: f32, // rings per second (mapped to steps)
}

#[derive(Resource, Clone, Copy)]
pub struct SimParams {
    pub production_rate: f32,      // auxin production per step
    pub decay_rate: f32,           // auxin decay per step
    pub diffusion_rate: f32,       // symmetric diffusion weight
    pub transport_rate: f32,       // PIN-mediated transport weight
    pub pin_alignment_rate: f32,   // how fast pins align to gradient
    pub primordium_threshold: f32, // auxin threshold to initiate new primordium
    pub primordium_min_separation: f32, // radians separation from recent primordia
    pub plastochron_steps: u32,    // minimal steps between initiations
    pub sink_radius: f32,          // influence radius for sink (world units)
}

impl Default for SimParams {
    fn default() -> Self {
        Self {
            production_rate: 0.03,
            decay_rate: 0.01,
            diffusion_rate: 0.05,
            transport_rate: 0.20,
            pin_alignment_rate: 0.25,
            primordium_threshold: 1.2,
            primordium_min_separation: std::f32::consts::PI / 3.0, // 60 deg
            plastochron_steps: 30,
            sink_radius: 28.0,
        }
    }
}

#[derive(Resource, Default)]
pub struct SimulationControl {
    pub paused: bool,
}

pub fn simulation_running(control: Res<SimulationControl>) -> bool {
    !control.paused
}


