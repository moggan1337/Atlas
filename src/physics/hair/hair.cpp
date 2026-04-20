/**
 * @file hair.cpp
 * @brief Hair Simulation Implementation
 */

#include "hair.hpp"
#include <algorithm>
#include <cmath>
#include <random>

namespace atlas {
namespace physics {

HairSimulation::HairSimulation(const HairConfig& cfg) 
    : config(cfg)
    , rng(42)
    , dist(0.0f, 1.0f) 
{}

int HairSimulation::add_strand(const Vec3& root, const Vec3& normal, float length) {
    HairStrand strand;
    strand.root_position = root;
    strand.root_normal = normal;
    strand.segment_count = config.segments_per_strand;
    strand.length = length > 0.0f ? length : 
                    config.length_min + dist(rng) * (config.length_max - config.length_min);
    strand.thickness = config.thickness;
    strand.curl = config.curl;
    strand.color = config.base_color + Vec3(
        (dist(rng) - 0.5f) * config.color_variation,
        (dist(rng) - 0.5f) * config.color_variation,
        (dist(rng) - 0.5f) * config.color_variation
    );
    
    // Create segments
    float segment_length = strand.length / strand.segment_count;
    strand.root_index = static_cast<int>(all_segments.size());
    
    for (int i = 0; i < strand.segment_count; ++i) {
        HairSegment seg;
        float t = static_cast<float>(i + 1) / strand.segment_count;
        
        // Curled position along strand
        float curl_offset = std::sin(t * M_PI) * strand.curl * segment_length * i;
        Vec3 tangent = normal;
        Vec3 bitangent = std::abs(tangent.y) < 0.9f ? 
                        Vec3(0.0f, 1.0f, 0.0f).cross(tangent).normalized() :
                        Vec3(1.0f, 0.0f, 0.0f).cross(tangent).normalized();
        
        seg.local_offset = tangent * (t * strand.length) + 
                          bitangent * curl_offset * (dist(rng) - 0.5f);
        seg.position = root + seg.local_offset;
        seg.previous_position = seg.position;
        seg.mass = config.segment_mass;
        seg.inv_mass = 1.0f / config.segment_mass;
        
        all_segments.push_back(seg);
    }
    
    strands.push_back(strand);
    return static_cast<int>(strands.size()) - 1;
}

void HairSimulation::add_strands_on_surface(const std::vector<Vec3>& positions,
                                           const std::vector<Vec3>& normals,
                                           float density) {
    float d = density > 0.0f ? density : config.density;
    
    for (size_t i = 0; i < positions.size(); ++i) {
        int strands_to_add = static_cast<int>(d);
        for (int j = 0; j < strands_to_add; ++j) {
            Vec3 pos = positions[i];
            Vec3 normal = normals[i];
            
            // Add small random offset
            pos += normal * 0.001f;  // Slightly above surface
            add_strand(pos, normal);
        }
    }
}

void HairSimulation::generate_follicles(const Vec3& center, const Vec3& size, float density) {
    std::uniform_real_distribution<float> dist_x(-size.x/2, size.x/2);
    std::uniform_real_distribution<float> dist_y(-size.y/2, size.y/2);
    std::uniform_real_distribution<float> dist_z(-size.z/2, size.z/2);
    
    int count = static_cast<int>(density * size.x * size.y);
    
    for (int i = 0; i < count; ++i) {
        HairFollicle follicle;
        follicle.position = center + Vec3(dist_x(rng), dist_y(rng), dist_z(rng));
        follicle.normal = Vec3::unit_y();
        follicle.density = config.density;
        follicle.length = config.length_min + dist(rng) * (config.length_max - config.length_min);
        follicle.thickness = config.thickness;
        follicle.active = true;
        
        int strand_idx = add_strand(follicle.position, follicle.normal, follicle.length);
        follicle.strand_index = strand_idx;
        
        follicles.push_back(follicle);
    }
}

void HairSimulation::update(float dt, const Vec3& gravity) {
    // Apply forces
    apply_forces(dt, gravity);
    
    // Integration
    integrate_verlet(dt);
    
    // Constraints
    solve_constraints();
    
    // Wind
    if (config.wind_strength > 0.0f) {
        apply_wind(0.0f, dt);
    }
    
    // Update bounds
    update_bounds();
}

void HairSimulation::integrate_verlet(float dt) {
    float dt_sq = dt * dt;
    
    for (auto& seg : all_segments) {
        Vec3 velocity = seg.position - seg.previous_position;
        seg.velocity = velocity / dt;
        
        Vec3 temp = seg.position;
        seg.position = seg.position * 2.0f - seg.previous_position 
                     + seg.velocity * config.damping * dt
                     + seg.velocity * dt;  // Damping on velocity
        seg.previous_position = temp;
    }
}

void HairSimulation::solve_constraints() {
    float stiffness = config.stiffness;
    
    // Iterate through all strands
    for (auto& strand : strands) {
        Vec3 root_pos = strand.root_position;
        
        for (int i = 0; i < strand.segment_count; ++i) {
            int seg_idx = strand.root_index + i;
            HairSegment& seg = all_segments[seg_idx];
            
            // Constraint to parent (root or previous segment)
            Vec3 parent_pos = (i == 0) ? root_pos : all_segments[seg_idx - 1].position;
            float segment_length = strand.length / strand.segment_count;
            
            Vec3 delta = seg.position - parent_pos;
            float dist = delta.magnitude();
            
            if (dist < 1e-6f) continue;
            
            float diff = (dist - segment_length) / dist;
            Vec3 correction = delta * diff * 0.5f;
            
            // Stiffness factor
            float k = 1.0f - std::pow(1.0f - stiffness / 10000.0f, 0.25f);
            
            seg.position -= correction * k;
        }
    }
}

void HairSimulation::apply_forces(float dt, const Vec3& gravity) {
    for (auto& seg : all_segments) {
        seg.velocity += gravity * config.gravity_scale * dt;
    }
}

void HairSimulation::apply_wind(float time, float dt) {
    for (auto& strand : strands) {
        for (int i = 0; i < strand.segment_count; ++i) {
            int seg_idx = strand.root_index + i;
            HairSegment& seg = all_segments[seg_idx];
            
            // Wind force with turbulence
            float t = static_cast<float>(i) / strand.segment_count;
            float turbulence = std::sin(time * 10.0f + seg.position.x * 5.0f) * config.wind_turbulence;
            
            Vec3 wind = config.wind_direction * (config.wind_strength * (1.0f + turbulence));
            
            // Bend hair in wind direction based on height
            wind *= t * t;
            
            seg.velocity += wind * dt;
        }
    }
}

void HairSimulation::apply_gravity(const Vec3& gravity) {
    for (auto& seg : all_segments) {
        seg.velocity += gravity * config.gravity_scale;
    }
}

void HairSimulation::collide_with_sphere(const Vec3& center, float radius) {
    float collision_radius = radius + config.collision_radius;
    float radius_sq = collision_radius * collision_radius;
    
    for (auto& seg : all_segments) {
        Vec3 delta = seg.position - center;
        float dist_sq = delta.magnitude_squared();
        
        if (dist_sq < radius_sq && dist_sq > 1e-6f) {
            float dist = std::sqrt(dist_sq);
            Vec3 normal = delta / dist;
            seg.position = center + normal * collision_radius;
        }
    }
}

void HairSimulation::collide_with_body(const std::vector<Vec3>& vertices,
                                      const std::vector<uint32_t>& indices) {
    // Simple sphere collision with bounding sphere
    if (vertices.empty()) return;
    
    Vec3 center = Vec3::zero();
    for (const auto& v : vertices) {
        center += v;
    }
    center /= static_cast<float>(vertices.size());
    
    float max_dist = 0.0f;
    for (const auto& v : vertices) {
        max_dist = std::max(max_dist, (v - center).magnitude());
    }
    
    collide_with_sphere(center, max_dist);
}

void HairSimulation::update_bounds() {
    if (all_segments.empty()) return;
    
    bounds_min = all_segments[0].position;
    bounds_max = all_segments[0].position;
    
    for (const auto& seg : all_segments) {
        bounds_min = bounds_min.min(seg.position);
        bounds_max = bounds_max.max(seg.position);
    }
}

void HairSimulation::cut_hair(float max_length) {
    for (auto& strand : strands) {
        if (strand.length > max_length) {
            strand.length = max_length;
            
            // Update segment positions
            float segment_length = strand.length / strand.segment_count;
            for (int i = 0; i < strand.segment_count; ++i) {
                int seg_idx = strand.root_index + i;
                float t = static_cast<float>(i + 1) / strand.segment_count;
                all_segments[seg_idx].local_offset = strand.root_normal * (t * strand.length);
                all_segments[seg_idx].position = strand.root_position + all_segments[seg_idx].local_offset;
                all_segments[seg_idx].previous_position = all_segments[seg_idx].position;
            }
        }
    }
}

void HairSimulation::style_hair(const Vec3& direction) {
    Vec3 dir = direction.normalized();
    
    for (auto& strand : strands) {
        strand.root_normal = dir;
        
        float segment_length = strand.length / strand.segment_count;
        for (int i = 0; i < strand.segment_count; ++i) {
            int seg_idx = strand.root_index + i;
            float t = static_cast<float>(i + 1) / strand.segment_count;
            all_segments[seg_idx].local_offset = dir * (t * strand.length);
            all_segments[seg_idx].position = strand.root_position + all_segments[seg_idx].local_offset;
            all_segments[seg_idx].previous_position = all_segments[seg_idx].position;
        }
    }
}

void HairSimulation::apply_curl(float amount) {
    for (auto& strand : strands) {
        strand.curl = std::clamp(amount, 0.0f, 1.0f);
    }
}

std::vector<Vec3> HairSimulation::get_render_positions() const {
    return std::vector<Vec3>();
}

std::vector<Vec3> HairSimulation::get_render_normals() const {
    return std::vector<Vec3>();
}

std::vector<float> HairSimulation::get_render_thicknesses() const {
    std::vector<float> thicknesses;
    for (const auto& strand : strands) {
        for (int i = 0; i < strand.segment_count; ++i) {
            float t = static_cast<float>(i) / strand.segment_count;
            thicknesses.push_back(strand.thickness * (1.0f - t * 0.5f));
        }
    }
    return thicknesses;
}

std::vector<Vec3> HairSimulation::get_render_colors() const {
    std::vector<Vec3> colors;
    for (const auto& strand : strands) {
        for (int i = 0; i < strand.segment_count; ++i) {
            colors.push_back(strand.color);
        }
    }
    return colors;
}

std::vector<uint32_t> HairSimulation::get_render_indices() const {
    std::vector<uint32_t> indices;
    for (const auto& strand : strands) {
        for (int i = 0; i < strand.segment_count - 1; ++i) {
            int idx = strand.root_index + i;
            indices.push_back(idx);
            indices.push_back(idx + 1);
        }
    }
    return indices;
}

// ============== FurRenderer Implementation ==============

std::vector<Vec3> FurRenderer::generate_shell(const HairSimulation& hair, int layer_index) {
    const auto& config = hair.config;
    std::vector<Vec3> positions;
    
    ShellLayer& layer = layers[layer_index];
    
    for (size_t s = 0; s < hair.get_strand_count(); ++s) {
        const auto& strand = hair.get_strand(s);
        
        for (int seg = 0; seg < strand.segment_count; ++seg) {
            const auto& segment = hair.get_segment(strand.root_index + seg);
            
            float t = static_cast<float>(seg) / strand.segment_count;
            Vec3 shell_pos = strand.root_position + segment.local_offset * layer.scale;
            
            positions.push_back(shell_pos);
        }
    }
    
    return positions;
}

} // namespace physics
} // namespace atlas
