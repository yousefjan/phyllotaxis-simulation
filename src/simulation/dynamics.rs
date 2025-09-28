use bevy::prelude::*;

use crate::simulation::config::GridConfig;
use crate::simulation::helpers::{angular_distance, nearest_primordium_direction, ring_angle_to_pos};
use crate::simulation::types::{Primordium, PrimordiumList, SimParams, Tissue};

pub fn compute_auxin_gradient(tissue: &Tissue, config: &GridConfig, r: usize, s: usize) -> Vec2 {
    let a_center = tissue.auxin[r][s];
    let pos_center = tissue.cell_centers[r][s];
    let mut grad = Vec2::ZERO;

    let segs = tissue.segments_per_ring[r];
    let s_prev = (s + segs - 1) % segs;
    let s_next = (s + 1) % segs;
    for (rr, ss) in [(r, s_prev), (r, s_next)] {
        let delta = tissue.cell_centers[rr][ss] - pos_center;
        let a_delta = tissue.auxin[rr][ss] - a_center;
        grad += delta.normalize_or_zero() * a_delta;
    }

    if r > 0 {
        let segs_in = tissue.segments_per_ring[r - 1];
        let ss = ((s as f32) * (segs_in as f32) / (segs as f32)).round() as usize % segs_in;
        let delta = tissue.cell_centers[r - 1][ss] - pos_center;
        grad += delta.normalize_or_zero() * (tissue.auxin[r - 1][ss] - a_center);
    }
    if r + 1 < config.ring_count {
        let segs_out = tissue.segments_per_ring[r + 1];
        let ss = ((s as f32) * (segs_out as f32) / (segs as f32)).round() as usize % segs_out;
        let delta = tissue.cell_centers[r + 1][ss] - pos_center;
        grad += delta.normalize_or_zero() * (tissue.auxin[r + 1][ss] - a_center);
    }

    grad.normalize_or_zero()
}

pub fn update_pin_polarity(
    config: Res<GridConfig>,
    primordia: Res<PrimordiumList>,
    params: Res<SimParams>,
    mut tissue: ResMut<Tissue>,
) {
    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            let pos = tissue.cell_centers[r][s];
            let gradient_dir = compute_auxin_gradient(&tissue, &config, r, s);
            let target = if let Some(dir_to_sink) = nearest_primordium_direction(&config, &primordia, pos) {
                (dir_to_sink * 0.7 + gradient_dir * 0.3).normalize_or_zero()
            } else {
                gradient_dir
            };
            let current = tissue.pin_dir[r][s];
            let blended = (current * (1.0 - params.pin_alignment_rate) + target * params.pin_alignment_rate)
                .normalize_or_zero();
            tissue.pin_dir[r][s] = if blended.length_squared() < 1e-6 { Vec2::X } else { blended };
        }
    }
}

pub fn update_auxin(config: Res<GridConfig>, params: Res<SimParams>, mut tissue: ResMut<Tissue>) {
    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            tissue.next_auxin[r][s] = tissue.auxin[r][s];
        }
    }

    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            let a_here = tissue.auxin[r][s];
            let pos_here = tissue.cell_centers[r][s];
            let pin_dir = tissue.pin_dir[r][s];

            for &ss in [((s + segs - 1) % segs), ((s + 1) % segs)].iter() {
                let delta = tissue.cell_centers[r][ss] - pos_here;
                let edge_dir = delta.normalize_or_zero();
                let a_there = tissue.auxin[r][ss];

                let diff = params.diffusion_rate * (a_here - a_there);
                tissue.next_auxin[r][s] -= diff;
                tissue.next_auxin[r][ss] += diff;

                let conductance = (pin_dir.dot(edge_dir)).max(0.0);
                let flux = params.transport_rate * conductance * a_here;
                tissue.next_auxin[r][s] -= flux;
                tissue.next_auxin[r][ss] += flux;
            }

            if r > 0 {
                let segs_in = tissue.segments_per_ring[r - 1];
                let ss = ((s as f32) * (segs_in as f32) / (segs as f32)).round() as usize % segs_in;
                let delta = tissue.cell_centers[r - 1][ss] - pos_here;
                let edge_dir = delta.normalize_or_zero();
                let a_there = tissue.auxin[r - 1][ss];
                let diff = params.diffusion_rate * (a_here - a_there);
                tissue.next_auxin[r][s] -= diff;
                tissue.next_auxin[r - 1][ss] += diff;
                let conductance = (pin_dir.dot(edge_dir)).max(0.0);
                let flux = params.transport_rate * conductance * a_here;
                tissue.next_auxin[r][s] -= flux;
                tissue.next_auxin[r - 1][ss] += flux;
            }
            if r + 1 < config.ring_count {
                let segs_out = tissue.segments_per_ring[r + 1];
                let ss = ((s as f32) * (segs_out as f32) / (segs as f32)).round() as usize % segs_out;
                let delta = tissue.cell_centers[r + 1][ss] - pos_here;
                let edge_dir = delta.normalize_or_zero();
                let a_there = tissue.auxin[r + 1][ss];
                let diff = params.diffusion_rate * (a_here - a_there);
                tissue.next_auxin[r][s] -= diff;
                tissue.next_auxin[r + 1][ss] += diff;
                let conductance = (pin_dir.dot(edge_dir)).max(0.0);
                let flux = params.transport_rate * conductance * a_here;
                tissue.next_auxin[r][s] -= flux;
                tissue.next_auxin[r + 1][ss] += flux;
            }
        }
    }

    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            let mut a = tissue.next_auxin[r][s];
            a += params.production_rate;
            a -= params.decay_rate * a;
            tissue.next_auxin[r][s] = a.max(0.0);
        }
    }

    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            tissue.auxin[r][s] = tissue.next_auxin[r][s];
        }
    }
}

pub fn initiate_primordia(
    time: Res<Time>,
    config: Res<GridConfig>,
    mut primordia: ResMut<PrimordiumList>,
    params: Res<SimParams>,
    mut tissue: ResMut<Tissue>,
) {
    for p in primordia.items.iter_mut() {
        let delta_rings = p.radial_velocity * time.delta_seconds();
        p.ring_index = (p.ring_index as f32 + delta_rings).floor() as usize;
    }
    primordia.items.retain(|p| p.ring_index < config.ring_count);

    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            let pos = tissue.cell_centers[r][s];
            let mut sink = 0.0;
            for p in primordia.items.iter() {
                let pp = ring_angle_to_pos(&config, p.ring_index, p.angle);
                let d = pos.distance(pp);
                if d < params.sink_radius {
                    let w = 1.0 - (d / params.sink_radius);
                    sink += p.sink_strength * w;
                }
            }
            tissue.auxin[r][s] = (tissue.auxin[r][s] - sink).max(0.0);
        }
    }

    primordia.steps_since_last = primordia.steps_since_last.saturating_add(1);
    if primordia.steps_since_last < params.plastochron_steps {
        return;
    }

    let r = config.ring_count - 2;
    let segs = tissue.segments_per_ring[r];
    let mut best_s = 0usize;
    let mut best_a = -1.0f32;
    for s in 0..segs {
        let a = tissue.auxin[r][s];
        if a > best_a {
            best_a = a;
            best_s = s;
        }
    }

    if best_a >= params.primordium_threshold {
        let theta = 2.0 * std::f32::consts::PI * (best_s as f32 + 0.5) / segs as f32;
        let ok = primordia
            .items
            .iter()
            .rev()
            .take(2)
            .all(|p| angular_distance(theta, p.angle) >= params.primordium_min_separation);
        if ok {
            primordia.items.push(Primordium {
                ring_index: r,
                angle: theta,
                sink_strength: 0.20,
                radial_velocity: 0.8,
            });
            primordia.steps_since_last = 0;
        }
    }
}


