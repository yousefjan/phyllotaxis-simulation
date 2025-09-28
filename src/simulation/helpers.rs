use bevy::prelude::*;

use crate::simulation::config::GridConfig;
use crate::simulation::types::PrimordiumList;

pub fn ring_angle_to_pos(config: &GridConfig, ring: usize, angle: f32) -> Vec2 {
    let radius = (ring as f32 + 0.5) * config.ring_radius_step;
    Vec2::from_angle(angle) * radius + config.center
}

pub fn nearest_primordium_direction(
    config: &GridConfig,
    primordia: &PrimordiumList,
    pos: Vec2,
) -> Option<Vec2> {
    let mut best: Option<(f32, Vec2)> = None;
    for p in &primordia.items {
        let pp = ring_angle_to_pos(config, p.ring_index, p.angle);
        let delta = pp - pos;
        let d2 = delta.length_squared();
        if d2 <= 0.0001 {
            continue;
        }
        match best {
            None => best = Some((d2, delta.normalize())),
            Some((bd2, _)) if d2 < bd2 => best = Some((d2, delta.normalize())),
            _ => {}
        }
    }
    best.map(|(_, dir)| dir)
}

pub fn angular_distance(a: f32, b: f32) -> f32 {
    let mut d = (a - b).abs();
    while d > std::f32::consts::PI {
        d -= 2.0 * std::f32::consts::PI;
    }
    d.abs()
}

pub fn auxin_color(a: f32) -> Color {
    let x = (a / 2.0).clamp(0.0, 1.0);
    let r = x;
    let g = (0.2 + x).clamp(0.0, 1.0);
    let b = (1.0 - x).clamp(0.0, 1.0);
    Color::rgba(r, g, b, 1.0)
}


