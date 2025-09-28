use bevy::prelude::*;
use rand::Rng;

use crate::simulation::types::{PrimordiumList, SimulationControl, SimParams, Tissue};

pub fn handle_input(
    keys: Res<ButtonInput<KeyCode>>,
    mut control: ResMut<SimulationControl>,
    mut tissue: ResMut<Tissue>,
    mut primordia: ResMut<PrimordiumList>,
    mut params: ResMut<SimParams>,
) {
    if keys.just_pressed(KeyCode::Space) {
        control.paused = !control.paused;
    }
    if keys.just_pressed(KeyCode::KeyR) {
        let mut rng = rand::thread_rng();
        for r in 0..tissue.auxin.len() {
            for s in 0..tissue.auxin[r].len() {
                tissue.auxin[r][s] = rng.gen_range(0.6..0.8);
                tissue.next_auxin[r][s] = 0.0;
                tissue.pin_dir[r][s] = Vec2::X;
            }
        }
        primordia.items.clear();
        primordia.steps_since_last = 0;
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        params.transport_rate = (params.transport_rate * 0.9).max(0.01);
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        params.transport_rate = (params.transport_rate * 1.1).min(1.0);
    }
}


