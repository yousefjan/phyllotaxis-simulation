use bevy::prelude::*;

pub mod config;
pub mod dynamics;
pub mod helpers;
pub mod input;
pub mod setup;
pub mod types;
pub mod visuals;

use config::GridConfig;
use dynamics::{initiate_primordia, update_auxin, update_pin_polarity};
use input::handle_input;
use setup::{setup_camera, setup_simulation};
use types::{simulation_running, SimParams, SimulationControl};
use visuals::{draw_gizmos, update_colors};

pub struct SimulationPlugin;

impl Plugin for SimulationPlugin {
    fn build(&self, app: &mut App) {
        app.insert_resource(SimulationControl::default())
            .insert_resource(GridConfig::default())
            .insert_resource(SimParams::default())
            .add_systems(Startup, (setup_camera, setup_simulation))
            .add_systems(
                Update,
                (
                    handle_input,
                    update_pin_polarity.run_if(simulation_running),
                    update_auxin.run_if(simulation_running),
                    initiate_primordia.run_if(simulation_running),
                    update_colors,
                    draw_gizmos,
                ),
            );
    }
}


