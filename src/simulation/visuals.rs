use bevy::prelude::*;

use crate::simulation::config::GridConfig;
use crate::simulation::helpers::auxin_color;
use crate::simulation::types::{CellVisual, Tissue};

pub fn update_colors(
    mut materials: ResMut<Assets<ColorMaterial>>,
    tissue: Res<Tissue>,
) {
    let _ = materials; // material recoloring is done in gizmos for now
    let _ = tissue;
}

pub fn draw_gizmos(
    mut gizmos: Gizmos,
    mut q_visuals: Query<(&mut Handle<ColorMaterial>, Entity), With<CellVisual>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    tissue: Res<Tissue>,
    config: Res<GridConfig>,
) {
    for (mat_handle, entity) in q_visuals.iter_mut() {
        if let Some((r, s)) = find_cell_by_entity(&tissue, entity) {
            let a = tissue.auxin[r][s];
            let color = auxin_color(a);
            if let Some(mat) = materials.get_mut(mat_handle.id()) {
                mat.color = color;
            }
        }
    }

    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            let pos = tissue.cell_centers[r][s];
            let dir = tissue.pin_dir[r][s];
            let p0 = Vec3::new(pos.x, pos.y, 1.0);
            let p1 = Vec3::new(pos.x + dir.x * 12.0, pos.y + dir.y * 12.0, 1.0);
            gizmos.line(p0, p1, Color::WHITE);
        }
    }
}

fn find_cell_by_entity(tissue: &Tissue, entity: Entity) -> Option<(usize, usize)> {
    for (r, row) in tissue.cell_entities.iter().enumerate() {
        for (s, &e) in row.iter().enumerate() {
            if e == entity {
                return Some((r, s));
            }
        }
    }
    None
}


