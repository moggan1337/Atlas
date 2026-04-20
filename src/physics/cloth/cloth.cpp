/**
 * @file cloth.cpp
 * @brief Cloth Simulation Implementation
 */

#include "cloth.hpp"
#include <algorithm>
#include <cmath>

namespace atlas {
namespace physics {

int Cloth::add_particle(const Vec3& position, float mass) {
    ClothParticle p;
    p.position = position;
    p.previous_position = position;
    p.mass = mass > 0.0f ? mass : config.mass;
    p.inv_mass = 1.0f / p.mass;
    particles.push_back(p);
    return static_cast<int>(particles.size()) - 1;
}

int Cloth::add_particle(const Vec3& position, const Vec2& uv, float mass) {
    int idx = add_particle(position, mass);
    particles[idx].uv = uv;
    return idx;
}

void Cloth::add_constraint(int a, int b, ClothConstraintType type, float stiffness) {
    float rest_length = (particles[a].position - particles[b].position).magnitude();
    float k = stiffness > 0.0f ? stiffness : 
              (type == ClothConstraintType::Structural ? config.structural_stiffness :
               type == ClothConstraintType::Shear ? config.shear_stiffness :
               config.bend_stiffness);
    constraints.emplace_back(a, b, rest_length, k, type);
}

void Cloth::add_structural_constraint(int a, int b) {
    add_constraint(a, b, ClothConstraintType::Structural);
}

void Cloth::add_shear_constraint(int a, int b) {
    add_constraint(a, b, ClothConstraintType::Shear);
}

void Cloth::add_bend_constraint(int a, int b) {
    add_constraint(a, b, ClothConstraintType::Bend);
}

int Cloth::add_triangle(int i0, int i1, int i2) {
    triangles.emplace_back(i0, i1, i2);
    return static_cast<int>(triangles.size()) - 1;
}

void Cloth::pin_particle(int index) {
    if (index >= 0 && index < static_cast<int>(particles.size())) {
        particles[index].pinned = true;
        particles[index].inv_mass = 0.0f;
    }
}

void Cloth::unpin_particle(int index) {
    if (index >= 0 && index < static_cast<int>(particles.size())) {
        particles[index].pinned = false;
        particles[index].inv_mass = 1.0f / particles[index].mass;
    }
}

std::shared_ptr<Cloth> Cloth::create_grid(const Vec3& origin, float width, float height,
                                          int segments_x, int segments_y,
                                          bool pin_top) {
    auto cloth = std::make_shared<Cloth>();
    cloth->name = "cloth_grid";
    
    float dx = width / segments_x;
    float dy = height / segments_y;
    
    // Create particles
    for (int j = 0; j <= segments_y; ++j) {
        for (int i = 0; i <= segments_x; ++i) {
            Vec3 pos = origin + Vec3(i * dx, 0.0f, j * dy);
            Vec2 uv(static_cast<float>(i) / segments_x, 
                   static_cast<float>(j) / segments_y);
            int idx = cloth->add_particle(pos, uv);
            
            if (pin_top && j == 0) {
                cloth->pin_particle(idx);
            }
        }
    }
    
    // Create constraints
    int verts_per_row = segments_x + 1;
    for (int j = 0; j <= segments_y; ++j) {
        for (int i = 0; i <= segments_x; ++i) {
            int idx = j * verts_per_row + i;
            
            // Structural
            if (i < segments_x) cloth->add_structural_constraint(idx, idx + 1);
            if (j < segments_y) cloth->add_structural_constraint(idx, idx + verts_per_row);
            
            // Shear
            if (i < segments_x && j < segments_y) {
                cloth->add_shear_constraint(idx, idx + verts_per_row + 1);
                cloth->add_shear_constraint(idx + 1, idx + verts_per_row);
            }
            
            // Bend
            if (i < segments_x - 1) cloth->add_bend_constraint(idx, idx + 2);
            if (j < segments_y - 1) cloth->add_bend_constraint(idx, idx + 2 * verts_per_row);
        }
    }
    
    // Create triangles
    for (int j = 0; j < segments_y; ++j) {
        for (int i = 0; i < segments_x; ++i) {
            int idx = j * verts_per_row + i;
            cloth->add_triangle(idx, idx + 1, idx + verts_per_row + 1);
            cloth->add_triangle(idx, idx + verts_per_row + 1, idx + verts_per_row);
        }
    }
    
    return cloth;
}

std::shared_ptr<Cloth> Cloth::create_curtain(const Vec3& origin, float width, float height,
                                             int segments_x, int segments_y,
                                             int pin_spacing) {
    auto cloth = create_grid(origin, width, height, segments_x, segments_y, false);
    cloth->name = "curtain";
    
    // Pin at intervals
    for (int i = 0; i <= segments_x; i += pin_spacing) {
        cloth->pin_particle(i);
    }
    
    return cloth;
}

void Cloth::update(float dt, const Vec3& gravity) {
    // Apply forces
    for (auto& p : particles) {
        if (!p.pinned) {
            p.acceleration = gravity * config.gravity_scale;
        }
    }
    
    // Wind
    if (config.wind_strength > 0.0f) {
        apply_wind(dt);
    }
    
    // Integration
    integrate_verlet(dt);
    
    // Constraint solving
    for (int iter = 0; iter < config.constraint_iterations; ++iter) {
        solve_constraints(1);
    }
    
    // Self-collision
    if (config.self_collision) {
        self_collide(dt);
    }
    
    // Normals
    update_normals();
    
    // AABB
    update_aabb();
}

void Cloth::integrate_verlet(float dt) {
    float dt_sq = dt * dt;
    
    for (auto& p : particles) {
        if (p.pinned) continue;
        
        // Velocity estimate
        p.position = p.position * 2.0f - p.previous_position 
                   + p.acceleration * dt_sq;
        p.position *= (1.0f - config.damping * dt);
        
        p.previous_position = p.position;
        p.acceleration = Vec3::zero();
    }
}

void Cloth::solve_constraints(int iterations) {
    for (int iter = 0; iter < iterations; ++iter) {
        for (const auto& c : constraints) {
            ClothParticle& p_a = particles[c.particle_a];
            ClothParticle& p_b = particles[c.particle_b];
            
            Vec3 delta = p_b.position - p_a.position;
            float current_length = delta.magnitude();
            
            if (current_length < 1e-6f) continue;
            
            float diff = (current_length - c.rest_length) / current_length;
            Vec3 correction = delta * diff * 0.5f;
            
            // Stiffness factor
            float k = 1.0f - std::pow(1.0f - c.stiffness / 1000.0f, 0.25f);
            
            if (!p_a.pinned) {
                p_a.position += correction * k * (p_a.inv_mass / (p_a.inv_mass + p_b.inv_mass));
            }
            if (!p_b.pinned) {
                p_b.position -= correction * k * (p_b.inv_mass / (p_a.inv_mass + p_b.inv_mass));
            }
        }
    }
}

void Cloth::apply_wind(float dt) {
    // Apply wind force based on triangle normals
    for (const auto& tri : triangles) {
        ClothParticle& p0 = particles[tri.indices[0]];
        ClothParticle& p1 = particles[tri.indices[1]];
        ClothParticle& p2 = particles[tri.indices[2]];
        
        if (p0.pinned && p1.pinned && p2.pinned) continue;
        
        Vec3 edge_a = p1.position - p0.position;
        Vec3 edge_b = p2.position - p0.position;
        Vec3 normal = edge_a.cross(edge_b);
        
        float area = normal.magnitude() * 0.5f;
        normal = normal.normalized();
        
        // Wind force
        float wind_dot = normal.dot(config.wind_direction);
        Vec3 wind_force = config.wind_direction * (wind_dot * area * config.wind_strength);
        
        float inv_mass_sum = p0.inv_mass + p1.inv_mass + p2.inv_mass;
        if (inv_mass_sum < 1e-6f) continue;
        
        Vec3 force_per_mass = wind_force / inv_mass_sum;
        
        if (!p0.pinned) p0.acceleration += force_per_mass;
        if (!p1.pinned) p1.acceleration += force_per_mass;
        if (!p2.pinned) p2.acceleration += force_per_mass;
    }
}

void Cloth::apply_damping(float dt) {
    for (auto& p : particles) {
        if (p.pinned) continue;
        // Damping is applied during integration
    }
}

void Cloth::collide_with_sphere(const Vec3& center, float radius) {
    float collision_radius = radius + config.thickness;
    
    for (auto& p : particles) {
        if (p.pinned) continue;
        
        Vec3 delta = p.position - center;
        float dist = delta.magnitude();
        
        if (dist < collision_radius && dist > 1e-6f) {
            Vec3 normal = delta / dist;
            p.position = center + normal * collision_radius;
        }
    }
}

void Cloth::collide_with_box(const Vec3& min, const Vec3& max) {
    float margin = config.thickness;
    
    for (auto& p : particles) {
        if (p.pinned) continue;
        
        if (p.position.x < min.x + margin) p.position.x = min.x + margin;
        if (p.position.x > max.x - margin) p.position.x = max.x - margin;
        if (p.position.y < min.y + margin) p.position.y = min.y + margin;
        if (p.position.y > max.y - margin) p.position.y = max.y - margin;
        if (p.position.z < min.z + margin) p.position.z = min.z + margin;
        if (p.position.z > max.z - margin) p.position.z = max.z - margin;
    }
}

void Cloth::collide_with_plane(const Vec3& normal, float offset) {
    for (auto& p : particles) {
        if (p.pinned) continue;
        
        float dist = p.position.dot(normal) - offset;
        if (dist < config.thickness) {
            p.position = p.position + normal * (config.thickness - dist);
        }
    }
}

void Cloth::self_collide(float dt) {
    float min_dist = config.thickness * 2.0f;
    float min_dist_sq = min_dist * min_dist;
    
    // Simple O(n²) self-collision
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            // Skip connected particles
            bool connected = false;
            for (const auto& c : constraints) {
                if ((c.particle_a == static_cast<int>(i) && c.particle_b == static_cast<int>(j)) ||
                    (c.particle_b == static_cast<int>(i) && c.particle_a == static_cast<int>(j))) {
                    connected = true;
                    break;
                }
            }
            if (connected) continue;
            
            Vec3 delta = particles[j].position - particles[i].position;
            float dist_sq = delta.magnitude_squared();
            
            if (dist_sq < min_dist_sq && dist_sq > 1e-6f) {
                float dist = std::sqrt(dist_sq);
                Vec3 normal = delta / dist;
                
                float overlap = min_dist - dist;
                
                if (!particles[i].pinned && !particles[j].pinned) {
                    particles[i].position -= normal * (overlap * 0.5f);
                    particles[j].position += normal * (overlap * 0.5f);
                } else if (!particles[i].pinned) {
                    particles[i].position -= normal * overlap;
                } else if (!particles[j].pinned) {
                    particles[j].position += normal * overlap;
                }
            }
        }
    }
}

void Cloth::update_normals() {
    // Reset normals
    for (auto& p : particles) {
        p.normal = Vec3::zero();
    }
    
    // Accumulate face normals
    for (const auto& tri : triangles) {
        Vec3& p0 = particles[tri.indices[0]].position;
        Vec3& p1 = particles[tri.indices[1]].position;
        Vec3& p2 = particles[tri.indices[2]].position;
        
        Vec3 edge_a = p1 - p0;
        Vec3 edge_b = p2 - p0;
        Vec3 normal = edge_a.cross(edge_b);
        
        for (int i = 0; i < 3; ++i) {
            particles[tri.indices[i]].normal += normal;
        }
    }
    
    // Normalize
    for (auto& p : particles) {
        p.normal.normalize();
    }
}

void Cloth::update_aabb() {
    if (particles.empty()) return;
    
    aabb_min = particles[0].position;
    aabb_max = particles[0].position;
    
    for (const auto& p : particles) {
        aabb_min = aabb_min.min(p.position);
        aabb_max = aabb_max.max(p.position);
    }
}

Vec3 Cloth::get_center() const {
    Vec3 center = Vec3::zero();
    for (const auto& p : particles) {
        center += p.position;
    }
    return particles.empty() ? Vec3::zero() : center / static_cast<float>(particles.size());
}

std::vector<Vec3> Cloth::get_positions() const {
    std::vector<Vec3> positions;
    positions.reserve(particles.size());
    for (const auto& p : particles) {
        positions.push_back(p.position);
    }
    return positions;
}

std::vector<Vec3> Cloth::get_normals() const {
    std::vector<Vec3> normals;
    normals.reserve(particles.size());
    for (const auto& p : particles) {
        normals.push_back(p.normal);
    }
    return normals;
}

std::vector<Vec2> Cloth::get_uvs() const {
    std::vector<Vec2> uvs;
    uvs.reserve(particles.size());
    for (const auto& p : particles) {
        uvs.push_back(p.uv);
    }
    return uvs;
}

std::vector<uint32_t> Cloth::get_indices() const {
    std::vector<uint32_t> indices;
    indices.reserve(triangles.size() * 3);
    for (const auto& tri : triangles) {
        indices.push_back(tri.indices[0]);
        indices.push_back(tri.indices[1]);
        indices.push_back(tri.indices[2]);
    }
    return indices;
}

void Cloth::check_tearing() {
    if (config.tear_threshold >= INFINITY) return;
    
    for (auto it = constraints.begin(); it != constraints.end(); ) {
        const Vec3& p_a = particles[it->particle_a].position;
        const Vec3& p_b = particles[it->particle_b].position;
        
        float stretch = (p_a - p_b).magnitude() / it->rest_length;
        
        if (stretch > config.tear_threshold) {
            it = constraints.erase(it);
        } else {
            ++it;
        }
    }
}

// ============== ClothWorld Implementation ==============

void ClothWorld::step(float dt) {
    dt = std::min(dt, config.max_dt);
    float sub_dt = dt / static_cast<float>(config.substeps);
    
    for (int i = 0; i < config.substeps; ++i) {
        step_fixed(sub_dt);
    }
}

void ClothWorld::step_fixed(float dt) {
    for (auto& cloth : cloths) {
        cloth->update(dt, config.gravity);
    }
}

} // namespace physics
} // namespace atlas
