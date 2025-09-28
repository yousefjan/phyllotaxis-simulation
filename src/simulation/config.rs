use bevy::prelude::*;

#[derive(Resource, Clone, Copy)]
pub struct GridConfig {
    pub ring_count: usize,
    pub ring_radius_step: f32,
    pub base_segments: usize,
    pub segment_growth: f32,
    pub center: Vec2,
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


