use bevy::prelude::*;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy::render::render_asset::RenderAssetUsages;
use bevy::sprite::MaterialMesh2dBundle;
use rand::Rng;

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
        .insert_resource(SimulationControl::default())
        .insert_resource(GridConfig::default())
        .insert_resource(SimParams::default())
        .add_systems(Startup, (setup_camera, setup_simulation))
        .add_systems(Update, (
            handle_input,
            update_pin_polarity.run_if(simulation_running),
            update_auxin.run_if(simulation_running),
            initiate_primordia.run_if(simulation_running),
            update_colors,
            draw_gizmos,
        ))
        .run();
}

fn setup_camera(mut commands: Commands) {
    commands.spawn(Camera2dBundle::default());
}

#[derive(Resource, Clone, Copy)]
struct GridConfig {
    ring_count: usize,
    ring_radius_step: f32,
    base_segments: usize,
    segment_growth: f32,
    center: Vec2,
}

impl Default for GridConfig {
    fn default() -> Self {
        Self {
            ring_count: 12,
            ring_radius_step: 22.0,
            base_segments: 12,
            segment_growth: 1.22,
            center: Vec2::ZERO,
        }
    }
}

#[derive(Resource)]
struct Tissue {
    segments_per_ring: Vec<usize>,
    cell_centers: Vec<Vec<Vec2>>, // [ring][seg]
    auxin: Vec<Vec<f32>>,         // [ring][seg]
    next_auxin: Vec<Vec<f32>>,    // [ring][seg]
    pin_dir: Vec<Vec<Vec2>>,      // [ring][seg]
    cell_entities: Vec<Vec<Entity>>, // entity per cell for color updates
}

#[derive(Component)]
struct CellVisual;

#[derive(Resource, Default)]
struct PrimordiumList {
    items: Vec<Primordium>,
    steps_since_last: u32,
}

#[derive(Clone)]
struct Primordium {
    ring_index: usize,
    angle: f32, // radians
    sink_strength: f32,
    radial_velocity: f32, // rings per second (mapped to steps)
}

#[derive(Resource, Clone, Copy)]
struct SimParams {
    production_rate: f32,      // auxin production per step
    decay_rate: f32,           // auxin decay per step
    diffusion_rate: f32,       // symmetric diffusion weight
    transport_rate: f32,       // PIN-mediated transport weight
    pin_alignment_rate: f32,   // how fast pins align to gradient
    primordium_threshold: f32, // auxin threshold to initiate new primordium
    primordium_min_separation: f32, // radians separation from recent primordia
    plastochron_steps: u32,    // minimal steps between initiations
    sink_radius: f32,          // influence radius for sink (world units)
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
struct SimulationControl {
    paused: bool,
}

fn simulation_running(control: Res<SimulationControl>) -> bool {
    !control.paused
}

fn setup_simulation(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    config: Res<GridConfig>,
) {
    // Derive segments per ring
    let mut segments_per_ring = Vec::with_capacity(config.ring_count);
    let mut segments = config.base_segments as f32;
    for _ in 0..config.ring_count {
        segments_per_ring.push(segments.round().max(6.0) as usize);
        segments *= config.segment_growth;
    }

    // Precompute centers
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

    // Initialize fields
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

    // Create cell meshes and entities
    let mut cell_entities: Vec<Vec<Entity>> = Vec::with_capacity(config.ring_count);
    for r in 0..config.ring_count {
        let inner_r = r as f32 * config.ring_radius_step;
        let outer_r = (r as f32 + 1.0) * config.ring_radius_step;
        let segs = segments_per_ring[r];
        let mut row_entities = Vec::with_capacity(segs);
        for s in 0..segs {
            let theta0 = 2.0 * std::f32::consts::PI * (s as f32) / segs as f32;
            let theta1 = 2.0 * std::f32::consts::PI * (s as f32 + 1.0) / segs as f32;

            // Build a wedge quad as two triangles
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

fn ring_angle_to_pos(config: &GridConfig, ring: usize, angle: f32) -> Vec2 {
    let radius = (ring as f32 + 0.5) * config.ring_radius_step;
    Vec2::from_angle(angle) * radius + config.center
}

fn nearest_primordium_direction(
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

fn update_pin_polarity(
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
                // Bias towards sinks but also consider local auxin gradient
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

fn compute_auxin_gradient(tissue: &Tissue, config: &GridConfig, r: usize, s: usize) -> Vec2 {
    let a_center = tissue.auxin[r][s];
    let pos_center = tissue.cell_centers[r][s];
    let mut grad = Vec2::ZERO;

    // Angular neighbors within ring
    let segs = tissue.segments_per_ring[r];
    let s_prev = (s + segs - 1) % segs;
    let s_next = (s + 1) % segs;
    for (rr, ss) in [(r, s_prev), (r, s_next)] {
        let delta = tissue.cell_centers[rr][ss] - pos_center;
        let a_delta = tissue.auxin[rr][ss] - a_center;
        grad += delta.normalize_or_zero() * a_delta;
    }

    // Radial neighbors (nearest angular index mapping)
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

fn update_auxin(config: Res<GridConfig>, params: Res<SimParams>, mut tissue: ResMut<Tissue>) {
    // Compute fluxes and next auxin
    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            tissue.next_auxin[r][s] = tissue.auxin[r][s];
        }
    }

    // Transport and diffusion
    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            let a_here = tissue.auxin[r][s];
            let pos_here = tissue.cell_centers[r][s];
            let pin_dir = tissue.pin_dir[r][s];

            // Angular neighbors
            for &ss in [((s + segs - 1) % segs), ((s + 1) % segs)].iter() {
                let delta = tissue.cell_centers[r][ss] - pos_here;
                let edge_dir = delta.normalize_or_zero();
                let a_there = tissue.auxin[r][ss];

                // Symmetric diffusion
                let diff = params.diffusion_rate * (a_here - a_there);
                tissue.next_auxin[r][s] -= diff;
                tissue.next_auxin[r][ss] += diff;

                // PIN-mediated transport: project pin onto edge direction
                let conductance = (pin_dir.dot(edge_dir)).max(0.0);
                let flux = params.transport_rate * conductance * a_here;
                tissue.next_auxin[r][s] -= flux;
                tissue.next_auxin[r][ss] += flux;
            }

            // Radial neighbors: inward
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
            // outward neighbor
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

    // Sources and sinks, decay
    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            let mut a = tissue.next_auxin[r][s];
            a += params.production_rate;
            a -= params.decay_rate * a;
            tissue.next_auxin[r][s] = a.max(0.0);
        }
    }

    // Swap
    for r in 0..config.ring_count {
        let segs = tissue.segments_per_ring[r];
        for s in 0..segs {
            tissue.auxin[r][s] = tissue.next_auxin[r][s];
        }
    }
}

fn initiate_primordia(
    time: Res<Time>,
    config: Res<GridConfig>,
    mut primordia: ResMut<PrimordiumList>,
    params: Res<SimParams>,
    mut tissue: ResMut<Tissue>,
) {
    // Move existing primordia outward over time
    for p in primordia.items.iter_mut() {
        let delta_rings = p.radial_velocity * time.delta_seconds();
        p.ring_index = (p.ring_index as f32 + delta_rings).floor() as usize;
    }
    // Cull primordia outside grid
    primordia.items.retain(|p| p.ring_index < config.ring_count);

    // Apply sink effect to nearby cells
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

    // Initiation rule: find max auxin in outer ring
    primordia.steps_since_last = primordia.steps_since_last.saturating_add(1);
    if primordia.steps_since_last < params.plastochron_steps {
        return;
    }

    let r = config.ring_count - 2; // peripheral zone (avoid border)
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
        // Enforce minimal angular separation from the two most recent primordia
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
                radial_velocity: 0.8, // rings per second
            });
            primordia.steps_since_last = 0;
        }
    }
}

fn angular_distance(a: f32, b: f32) -> f32 {
    let mut d = (a - b).abs();
    while d > std::f32::consts::PI {
        d -= 2.0 * std::f32::consts::PI;
    }
    d.abs()
}

fn update_colors(
    mut materials: ResMut<Assets<ColorMaterial>>,
    tissue: Res<Tissue>,
    params: Res<SimParams>,
) {
    let _ = params; // currently unused in color mapping, reserved for future
    for (r, row) in tissue.cell_entities.iter().enumerate() {
        for (s, &entity) in row.iter().enumerate() {
            // SAFETY: we only set per-entity material at spawn time; to update, fetch handle via world. Simpler: lookup material by entity is cumbersome here.
            // Instead, store ColorMaterial handles alongside entities. For simplicity in this prototype, we do a world query in draw_gizmos for overlays, and re-color by changing Color on material via Commands.
            // As a compromise: skip here; colors updated via a separate system using queries.
            let _ = (r, s, entity);
        }
    }
}

// Recolor by querying materials and remapping by auxin value
fn draw_gizmos(
    mut gizmos: Gizmos,
    mut q_visuals: Query<(&mut Handle<ColorMaterial>, Entity), With<CellVisual>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    tissue: Res<Tissue>,
    config: Res<GridConfig>,
) {
    // Update cell colors
    for (mat_handle, entity) in q_visuals.iter_mut() {
        // Find cell indices by scanning mapping (small grid, acceptable). In production, store a direct map.
        if let Some((r, s)) = find_cell_by_entity(&tissue, entity) {
            let a = tissue.auxin[r][s];
            let color = auxin_color(a);
            if let Some(mat) = materials.get_mut(mat_handle.id()) {
                mat.color = color;
            }
        }
    }

    // Draw PIN directions as arrows
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

    // Draw primordia markers
    let red = Color::rgb(1.0, 0.2, 0.2);
    // SAFETY: access primordia via world query would be cleaner; here we read resource by system param order in schedule.
    // This system signature does not include PrimordiumList; add separate system parameter if needed.
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

fn auxin_color(a: f32) -> Color {
    // Map [0, ~2] to a blue->green->yellow gradient
    let x = (a / 2.0).clamp(0.0, 1.0);
    let r = x;
    let g = (0.2 + x).clamp(0.0, 1.0);
    let b = (1.0 - x).clamp(0.0, 1.0);
    Color::rgba(r, g, b, 1.0)
}

fn handle_input(
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
        // Reset auxin and primordia
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
