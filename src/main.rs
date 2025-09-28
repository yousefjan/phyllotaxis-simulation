use bevy::prelude::*;
use phyllotaxis_simulation::simulation::SimulationPlugin;

fn main() {
    App::new()
        .insert_resource(ClearColor(Color::BLACK))
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: String::from("Phyllotaxis Simulation (Smith et al., PNAS 2006)"),
                resolution: (1280.0, 960.0).into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(SimulationPlugin)
        .run();
}
