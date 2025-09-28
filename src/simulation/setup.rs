use bevy::prelude::*;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy::render::render_asset::RenderAssetUsages;
use bevy::sprite::MaterialMesh2dBundle;
use rand::Rng;

use crate::simulation::config::GridConfig;
use crate::simulation::helpers::auxin_color;
use crate::simulation::types::{CellVisual, PrimordiumList, Tissue};

pub fn setup_camera(mut commands: Commands) {
    commands.spawn(Camera2dBundle::default());
}

pub fn setup_simulation(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    config: Res<GridConfig>,
) {
    let mut segments_per_ring = Vec::with_capacity(config.ring_count);
    let mut segments = config.base_segments as f32;
    for _ in 0..config.ring_count {
        segments_per_ring.push(segments.round().max(6.0) as usize);
        segments *= config.segment_growth;
    }

    let mut cell_centers: Vec<Vec<Vec2>> = Vec::with_capacity(config.ring_count);
    for r in 0..config.ring_count {
        let radius = (r as f32 + 0.5) * config.ring_radius_step;
        let segs = segments_per_ring[r];
        let mut centers = Vec::with_capacity(segs);
        for s in 0..segs {
            let theta = 2.0 * std::f32::consts::PI * (s as f32 + 0.5) / segs as f32;
            let pos = config.center + Vec2::from_angle(theta) * radius;
            centers.push(pos);
        }
        cell_centers.push(centers);
    }

    let mut rng = rand::thread_rng();
    let mut auxin = Vec::with_capacity(config.ring_count);
    let mut next_auxin = Vec::with_capacity(config.ring_count);
    let mut pin_dir = Vec::with_capacity(config.ring_count);
    for r in 0..config.ring_count {
        let segs = segments_per_ring[r];
        let mut a_row = Vec::with_capacity(segs);
        let mut n_row = Vec::with_capacity(segs);
        let mut p_row = Vec::with_capacity(segs);
        for _ in 0..segs {
            let a0 = rng.gen_range(0.6..0.8);
            a_row.push(a0);
            n_row.push(0.0);
            p_row.push(Vec2::X);
        }
        auxin.push(a_row);
        next_auxin.push(n_row);
        pin_dir.push(p_row);
    }

    let mut cell_entities: Vec<Vec<Entity>> = Vec::with_capacity(config.ring_count);
    for r in 0..config.ring_count {
        let inner_r = r as f32 * config.ring_radius_step;
        let outer_r = (r as f32 + 1.0) * config.ring_radius_step;
        let segs = segments_per_ring[r];
        let mut row_entities = Vec::with_capacity(segs);
        for s in 0..segs {
            let theta0 = 2.0 * std::f32::consts::PI * (s as f32) / segs as f32;
            let theta1 = 2.0 * std::f32::consts::PI * (s as f32 + 1.0) / segs as f32;

            let vertices = vec![
                Vec3::new((inner_r * theta0.cos()), (inner_r * theta0.sin()), 0.0),
                Vec3::new((outer_r * theta0.cos()), (outer_r * theta0.sin()), 0.0),
                Vec3::new((outer_r * theta1.cos()), (outer_r * theta1.sin()), 0.0),
                Vec3::new((inner_r * theta1.cos()), (inner_r * theta1.sin()), 0.0),
            ];
            let indices: Vec<u32> = vec![0, 1, 2, 0, 2, 3];
            let normals = vec![Vec3::Z; 4];
            let uvs = vec![
                Vec2::new(0.0, 0.0),
                Vec2::new(0.0, 1.0),
                Vec2::new(1.0, 1.0),
                Vec2::new(1.0, 0.0),
            ];
            let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::default());
            mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, vertices);
            mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
            mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uvs);
            mesh.insert_indices(Indices::U32(indices));

            let color = auxin_color(auxin[r][s]);
            let entity = commands
                .spawn(MaterialMesh2dBundle {
                    mesh: meshes.add(mesh).into(),
                    material: materials.add(ColorMaterial::from(color)),
                    transform: Transform::from_translation(Vec3::new(config.center.x, config.center.y, 0.0)),
                    ..default()
                })
                .insert(CellVisual)
                .id();
            row_entities.push(entity);
        }
        cell_entities.push(row_entities);
    }

    commands.insert_resource(Tissue {
        segments_per_ring,
        cell_centers,
        auxin,
        next_auxin,
        pin_dir,
        cell_entities,
    });

    commands.insert_resource(PrimordiumList::default());
}


