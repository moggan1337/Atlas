/**
 * @file soft_body.cpp
 * @brief Soft Body Implementation
 */

#include "soft_body.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace atlas {
namespace physics {

// ============== SoftBody Implementation ==============

int SoftBody::add_particle(const Vec3& position, float mass) {
    SoftBodyParticle p(position, mass);
    p.previous_position = position;
    particles.push_back(p);
    rest_positions.push_back(position);
    return static_cast<int>(particles.size()) - 1;
}

void SoftBody::pin_particle(int index) {
    if (index >= 0 && index < static_cast<int>(particles.size())) {
        particles[index].pinned = true;
        particles[index].inv_mass = 0.0f;
    }
}

void SoftBody::unpin_particle(int index) {
    if (index >= 0 && index < static_cast<int>(particles.size())) {
        particles[index].pinned = false;
        particles[index].inv_mass = 1.0f / particles[index].mass;
    }
}

int SoftBody::add_spring(int particle_a, int particle_b, float stiffness, float damping) {
    if (particle_a < 0 || particle_a >= static_cast<int>(particles.size()) ||
        particle_b < 0 || particle_b >= static_cast<int>(particles.size())) {
        return -1;
    }
    
    float k = stiffness > 0.0f ? stiffness : config.stiffness;
    float d = damping > 0.0f ? damping : config.damping;
    float rest = (particles[particle_a].position - particles[particle_b].position).magnitude();
    
    springs.emplace_back(particle_a, particle_b, rest, k, d);
    return static_cast<int>(springs.size()) - 1;
}

void SoftBody::remove_spring(int index) {
    if (index >= 0 && index < static_cast<int>(springs.size())) {
        springs.erase(springs.begin() + index);
    }
}

int SoftBody::add_face(int i0, int i1, int i2) {
    faces.emplace_back(i0, i1, i2);
    return static_cast<int>(faces.size()) - 1;
}

std::shared_ptr<SoftBody> SoftBody::create_ball(const Vec3& center, float radius, int segments) {
    auto ball = std::make_shared<SoftBody>();
    ball->name = "ball";
    
    // Create particles in a sphere pattern
    int stacks = segments;
    int slices = segments * 2;
    
    for (int i = 0; i <= stacks; ++i) {
        float phi = M_PI * i / stacks;
        for (int j = 0; j <= slices; ++j) {
            float theta = 2.0f * M_PI * j / slices;
            
            float x = center.x + radius * std::sin(phi) * std::cos(theta);
            float y = center.y + radius * std::cos(phi);
            float z = center.z + radius * std::sin(phi) * std::sin(theta);
            
            ball->add_particle(Vec3(x, y, z));
        }
    }
    
    // Create faces
    for (int i = 0; i < stacks; ++i) {
        for (int j = 0; j < slices; ++j) {
            int first = i * (slices + 1) + j;
            int second = first + slices + 1;
            
            ball->add_face(first, second, first + 1);
            ball->add_face(second, second + 1, first + 1);
        }
    }
    
    // Create springs (structural)
    int verts_per_row = slices + 1;
    for (int i = 0; i <= stacks; ++i) {
        for (int j = 0; j <= slices; ++j) {
            int idx = i * verts_per_row + j;
            
            // Horizontal spring
            if (j < slices) {
                ball->add_spring(idx, idx + 1);
            }
            
            // Vertical spring
            if (i < stacks) {
                ball->add_spring(idx, idx + verts_per_row);
            }
            
            // Diagonal springs (shear)
            if (i < stacks && j < slices) {
                ball->add_spring(idx, idx + verts_per_row + 1);
                ball->add_spring(idx + 1, idx + verts_per_row);
            }
        }
    }
    
    return ball;
}

std::shared_ptr<SoftBody> SoftBody::create_box(const Vec3& center, const Vec3& size) {
    auto box = std::make_shared<SoftBody>();
    box->name = "box";
    
    float hx = size.x * 0.5f;
    float hy = size.y * 0.5f;
    float hz = size.z * 0.5f;
    
    // 8 corners
    std::vector<int> corner_indices;
    Vec3 corners[] = {
        Vec3(center.x - hx, center.y - hy, center.z - hz),
        Vec3(center.x + hx, center.y - hy, center.z - hz),
        Vec3(center.x - hx, center.y + hy, center.z - hz),
        Vec3(center.x + hx, center.y + hy, center.z - hz),
        Vec3(center.x - hx, center.y - hy, center.z + hz),
        Vec3(center.x + hx, center.y - hy, center.z + hz),
        Vec3(center.x - hx, center.y + hy, center.z + hz),
        Vec3(center.x + hx, center.y + hy, center.z + hz)
    };
    
    for (const auto& corner : corners) {
        corner_indices.push_back(box->add_particle(corner));
    }
    
    // 12 edges (structural springs)
    box->add_spring(0, 1); box->add_spring(2, 3);
    box->add_spring(4, 5); box->add_spring(6, 7);
    box->add_spring(0, 2); box->add_spring(1, 3);
    box->add_spring(4, 6); box->add_spring(5, 7);
    box->add_spring(0, 4); box->add_spring(1, 5);
    box->add_spring(2, 6); box->add_spring(3, 7);
    
    // Face diagonals (shear springs)
    box->add_spring(0, 3); box->add_spring(1, 2);
    box->add_spring(4, 7); box->add_spring(5, 6);
    box->add_spring(0, 5); box->add_spring(1, 4);
    box->add_spring(2, 7); box->add_spring(3, 6);
    box->add_spring(0, 6); box->add_spring(2, 4);
    box->add_spring(1, 7); box->add_spring(3, 5);
    
    // 6 faces
    box->add_face(0, 2, 1); box->add_face(1, 2, 3);  // Front
    box->add_face(4, 5, 6); box->add_face(5, 7, 6);  // Back
    box->add_face(0, 4, 2); box->add_face(2, 4, 6);  // Left
    box->add_face(1, 3, 5); box->add_face(3, 7, 5);  // Right
    box->add_face(2, 6, 3); box->add_face(3, 6, 7);  // Top
    box->add_face(0, 1, 4); box->add_face(1, 5, 4);  // Bottom
    
    return box;
}

std::shared_ptr<SoftBody> SoftBody::create_cloth(const Vec3& origin, float width, float height,
                                                  int segments_x, int segments_y) {
    auto cloth = std::make_shared<SoftBody>();
    cloth->name = "cloth";
    
    float dx = width / segments_x;
    float dy = height / segments_y;
    
    // Create particles
    for (int j = 0; j <= segments_y; ++j) {
        for (int i = 0; i <= segments_x; ++i) {
            Vec3 pos = origin + Vec3(i * dx, -j * dy, 0.0f);
            cloth->add_particle(pos);
            
            // Pin top row
            if (j == 0) {
                cloth->pin_particle(cloth->get_particle_count() - 1);
            }
        }
    }
    
    // Create springs
    int verts_per_row = segments_x + 1;
    for (int j = 0; j <= segments_y; ++j) {
        for (int i = 0; i <= segments_x; ++i) {
            int idx = j * verts_per_row + i;
            
            // Structural springs
            if (i < segments_x) cloth->add_spring(idx, idx + 1);
            if (j < segments_y) cloth->add_spring(idx, idx + verts_per_row);
            
            // Shear springs
            if (i < segments_x && j < segments_y) {
                cloth->add_spring(idx, idx + verts_per_row + 1);
                cloth->add_spring(idx + 1, idx + verts_per_row);
            }
            
            // Bend springs (skip one)
            if (i < segments_x - 1) cloth->add_spring(idx, idx + 2);
            if (j < segments_y - 1) cloth->add_spring(idx, idx + 2 * verts_per_row);
        }
    }
    
    // Create faces
    for (int j = 0; j < segments_y; ++j) {
        for (int i = 0; i < segments_x; ++i) {
            int idx = j * verts_per_row + i;
            cloth->add_face(idx, idx + 1, idx + verts_per_row + 1);
            cloth->add_face(idx, idx + verts_per_row + 1, idx + verts_per_row);
        }
    }
    
    return cloth;
}

std::shared_ptr<SoftBody> SoftBody::create_chain(const Vec3& start, const Vec3& end, int segments) {
    auto chain = std::make_shared<SoftBody>();
    chain->name = "chain";
    
    Vec3 delta = end - start;
    float segment_length = delta.magnitude() / segments;
    Vec3 direction = delta.normalized();
    
    for (int i = 0; i <= segments; ++i) {
        Vec3 pos = start + direction * (i * segment_length);
        chain->add_particle(pos);
        
        if (i > 0) {
            chain->add_spring(i - 1, i);
        }
    }
    
    return chain;
}

std::shared_ptr<SoftBody> SoftBody::create_rope(const Vec3& start, float length, int segments) {
    auto rope = std::make_shared<SoftBody>();
    rope->name = "rope";
    
    float segment_length = length / segments;
    
    for (int i = 0; i <= segments; ++i) {
        Vec3 pos = start + Vec3(0.0f, -i * segment_length, 0.0f);
        rope->add_particle(pos);
        
        if (i > 0) {
            rope->add_spring(i - 1, i);
        }
    }
    
    return rope;
}

void SoftBody::update(float dt, const Vec3& gravity) {
    // Apply gravity
    for (auto& particle : particles) {
        if (!particle.pinned) {
            particle.acceleration = gravity * config.gravity_scale;
        }
    }
    
    // Integration
    integrate_verlet(dt);
    
    // Constraint solving iterations
    for (int i = 0; i < 3; ++i) {
        solve_springs(dt);
        solve_constraints(dt);
    }
    
    // Pressure
    if (config.pressure > 0.0f && !faces.empty()) {
        apply_pressure();
    }
    
    // Shape matching
    if (config.shape_matching && config.shape_stiffness > 0.0f) {
        enforce_shape_matching(dt);
    }
    
    // Update normals
    update_normals();
    
    // Update AABB
    update_aabb();
}

void SoftBody::integrate_verlet(float dt) {
    float dt_sq = dt * dt;
    
    for (auto& particle : particles) {
        if (particle.pinned) continue;
        
        // Calculate velocity from positions
        particle.velocity = (particle.position - particle.previous_position) / dt;
        
        // Store current position
        Vec3 temp = particle.position;
        
        // Verlet integration: x_new = 2*x - x_old + a*dt²
        particle.position = particle.position * 2.0f - particle.previous_position 
                          + particle.acceleration * dt_sq;
        
        // Apply damping
        particle.position = particle.position * (1.0f - config.damping * dt);
        
        particle.previous_position = temp;
        particle.acceleration = Vec3::zero();
    }
}

void SoftBody::solve_springs(float dt) {
    for (auto& spring : springs) {
        SoftBodyParticle& p_a = particles[spring.particle_a];
        SoftBodyParticle& p_b = particles[spring.particle_b];
        
        Vec3 delta = p_b.position - p_a.position;
        float current_length = delta.magnitude();
        
        if (current_length < 1e-6f) continue;
        
        // Check for tearing
        if (config.tear_threshold < INFINITY) {
            if (current_length > spring.rest_length * config.tear_threshold) {
                // Spring should break - mark for removal (handled externally)
                continue;
            }
        }
        
        float diff = (current_length - spring.rest_length) / current_length;
        Vec3 correction = delta * diff * 0.5f;
        
        // Stiffness factor
        float k = spring.stiffness / (p_a.inv_mass + p_b.inv_mass);
        k = 1.0f - std::pow(1.0f - k, 0.25f);
        
        if (p_a.inv_mass > 0.0f) {
            p_a.position += correction * k * p_a.inv_mass * 2.0f;
        }
        if (p_b.inv_mass > 0.0f) {
            p_b.position -= correction * k * p_b.inv_mass * 2.0f;
        }
    }
}

void SoftBody::solve_constraints(float dt) {
    // Simple position constraints can be added here
}

void SoftBody::apply_pressure() {
    if (faces.empty()) return;
    
    // Calculate current volume
    float volume = calculate_volume();
    if (volume < 1e-6f) return;
    
    // Calculate pressure
    float pressure = config.pressure * (volume / 1.0f); // Target volume = 1
    
    // Apply pressure force to each face
    for (auto& face : faces) {
        Vec3& p0 = particles[face.indices[0]].position;
        Vec3& p1 = particles[face.indices[1]].position;
        Vec3& p2 = particles[face.indices[2]].position;
        
        Vec3 edge_a = p1 - p0;
        Vec3 edge_b = p2 - p0;
        Vec3 normal = edge_a.cross(edge_b);
        
        float area = normal.magnitude() * 0.5f;
        normal = normal.normalized();
        
        // Force = pressure * area
        Vec3 force = normal * (pressure * area / 3.0f);
        
        for (int i = 0; i < 3; ++i) {
            SoftBodyParticle& p = particles[face.indices[i]];
            if (!p.pinned) {
                p.position += force * p.inv_mass * dt * dt;
            }
        }
    }
}

void SoftBody::apply_damping(float dt) {
    for (auto& particle : particles) {
        if (particle.pinned) continue;
        particle.velocity *= (1.0f - config.damping * dt);
    }
}

void SoftBody::enforce_shape_matching(float dt) {
    if (rest_positions.empty() || particles.size() != rest_positions.size()) return;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        SoftBodyParticle& p = particles[i];
        if (p.pinned) continue;
        
        Vec3 rest_pos = rest_positions[i];
        Vec3 target_pos = rest_transform.transform_point(rest_pos);
        
        Vec3 delta = target_pos - p.position;
        p.position += delta * config.shape_stiffness;
    }
}

void SoftBody::collide_with_plane(const Vec3& normal, float offset, float restitution) {
    for (auto& particle : particles) {
        if (particle.pinned) continue;
        
        float dist = particle.position.dot(normal) - offset;
        
        if (dist < 0.0f) {
            // Push particle out of plane
            particle.position = particle.position + normal * (-dist + 0.01f);
            
            // Reflect velocity
            Vec3 velocity = particle.position - particle.previous_position;
            float vel_normal = velocity.dot(normal);
            
            if (vel_normal < 0.0f) {
                particle.previous_position = particle.position + 
                    (velocity - normal * vel_normal * (1.0f + restitution));
            }
        }
    }
}

void SoftBody::collide_with_sphere(const Vec3& center, float radius, float restitution) {
    for (auto& particle : particles) {
        if (particle.pinned) continue;
        
        Vec3 delta = particle.position - center;
        float dist = delta.magnitude();
        
        if (dist < radius) {
            Vec3 normal = delta / dist;
            particle.position = center + normal * radius;
            
            // Reflect velocity
            Vec3 velocity = particle.position - particle.previous_position;
            float vel_normal = velocity.dot(normal);
            
            if (vel_normal < 0.0f) {
                particle.previous_position = particle.position + 
                    (velocity - normal * vel_normal * (1.0f + restitution));
            }
        }
    }
}

void SoftBody::collide_with_box(const Vec3& min, const Vec3& max, float restitution) {
    for (auto& particle : particles) {
        if (particle.pinned) continue;
        
        // Check if particle is inside box
        if (particle.position.x < min.x) {
            particle.position.x = min.x;
            particle.previous_position.x = particle.position.x + 
                (particle.position.x - particle.previous_position.x) * (1.0f + restitution);
        }
        if (particle.position.x > max.x) {
            particle.position.x = max.x;
            particle.previous_position.x = particle.position.x + 
                (particle.position.x - particle.previous_position.x) * (1.0f + restitution);
        }
        
        if (particle.position.y < min.y) {
            particle.position.y = min.y;
            particle.previous_position.y = particle.position.y + 
                (particle.position.y - particle.previous_position.y) * (1.0f + restitution);
        }
        if (particle.position.y > max.y) {
            particle.position.y = max.y;
            particle.previous_position.y = particle.position.y + 
                (particle.position.y - particle.previous_position.y) * (1.0f + restitution);
        }
        
        if (particle.position.z < min.z) {
            particle.position.z = min.z;
            particle.previous_position.z = particle.position.z + 
                (particle.position.z - particle.previous_position.z) * (1.0f + restitution);
        }
        if (particle.position.z > max.z) {
            particle.position.z = max.z;
            particle.previous_position.z = particle.position.z + 
                (particle.position.z - particle.previous_position.z) * (1.0f + restitution);
        }
    }
}

void SoftBody::update_aabb() {
    if (particles.empty()) return;
    
    aabb_min = particles[0].position;
    aabb_max = particles[0].position;
    
    for (const auto& particle : particles) {
        aabb_min = aabb_min.min(particle.position);
        aabb_max = aabb_max.max(particle.position);
    }
}

Vec3 SoftBody::get_center_of_mass() const {
    if (particles.empty()) return Vec3::zero();
    
    Vec3 center = Vec3::zero();
    float total_mass = 0.0f;
    
    for (const auto& particle : particles) {
        center += particle.position * particle.mass;
        total_mass += particle.mass;
    }
    
    return total_mass > 0.0f ? center / total_mass : Vec3::zero();
}

Vec3 SoftBody::get_velocity_at_point(const Vec3& point) const {
    Vec3 center = get_center_of_mass();
    Vec3 r = point - center;
    
    // For a simple soft body, assume rigid-like rotation
    // In a full implementation, would interpolate from nearby particles
    return Vec3::zero();
}

void SoftBody::apply_force(const Vec3& force) {
    for (auto& particle : particles) {
        if (!particle.pinned) {
            particle.acceleration += force * particle.inv_mass;
        }
    }
}

void SoftBody::apply_force_at_point(const Vec3& force, const Vec3& point) {
    // Find nearest particle and apply force there
    float min_dist = INFINITY;
    int nearest_idx = -1;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        float dist = (particles[i].position - point).magnitude_squared();
        if (dist < min_dist) {
            min_dist = dist;
            nearest_idx = static_cast<int>(i);
        }
    }
    
    if (nearest_idx >= 0) {
        particles[nearest_idx].acceleration += force * particles[nearest_idx].inv_mass;
    }
}

void SoftBody::apply_impulse(const Vec3& impulse) {
    for (auto& particle : particles) {
        if (!particle.pinned) {
            particle.previous_position -= impulse * particle.inv_mass;
        }
    }
}

void SoftBody::apply_impulse_at_point(const Vec3& impulse, const Vec3& point) {
    // Distribute impulse to nearby particles
    float min_dist = INFINITY;
    int nearest_idx = -1;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        float dist = (particles[i].position - point).magnitude_squared();
        if (dist < min_dist) {
            min_dist = dist;
            nearest_idx = static_cast<int>(i);
        }
    }
    
    if (nearest_idx >= 0) {
        particles[nearest_idx].previous_position -= impulse;
    }
}

void SoftBody::update_normals() {
    // Reset normals
    for (auto& particle : particles) {
        particle.normal = Vec3::zero();
    }
    
    // Accumulate face normals
    for (const auto& face : faces) {
        Vec3& p0 = particles[face.indices[0]].position;
        Vec3& p1 = particles[face.indices[1]].position;
        Vec3& p2 = particles[face.indices[2]].position;
        
        Vec3 normal = (p1 - p0).cross(p2 - p0);
        
        for (int i = 0; i < 3; ++i) {
            particles[face.indices[i]].normal += normal;
        }
    }
    
    // Normalize
    for (auto& particle : particles) {
        particle.normal.normalize();
    }
}

float SoftBody::calculate_volume() const {
    if (faces.empty()) return 0.0f;
    
    float volume = 0.0f;
    
    for (const auto& face : faces) {
        const Vec3& p0 = particles[face.indices[0]].position;
        const Vec3& p1 = particles[face.indices[1]].position;
        const Vec3& p2 = particles[face.indices[2]].position;
        
        // Signed volume of tetrahedron with origin
        volume += p0.dot(p1.cross(p2)) / 6.0f;
    }
    
    return std::abs(volume);
}

void SoftBody::check_tearing() {
    for (auto it = springs.begin(); it != springs.end(); ) {
        const Vec3& p_a = particles[it->particle_a].position;
        const Vec3& p_b = particles[it->particle_b].position;
        
        float stretch = (p_a - p_b).magnitude() / it->rest_length;
        
        if (stretch > config.tear_threshold) {
            it = springs.erase(it);
        } else {
            ++it;
        }
    }
}

// ============== SoftBodyWorld Implementation ==============

void SoftBodyWorld::step(float dt) {
    dt = std::min(dt, config.max_dt);
    float sub_dt = dt / static_cast<float>(config.substeps);
    
    for (int i = 0; i < config.substeps; ++i) {
        step_fixed(sub_dt);
    }
}

void SoftBodyWorld::step_fixed(float dt) {
    for (auto& body : bodies) {
        body->update(dt, config.gravity);
    }
}

// ============== ParticleSystem Implementation ==============

void ParticleSystem::add_particle(const Vec3& position, float mass) {
    SoftBodyParticle p(position, mass);
    p.previous_position = position;
    particles.push_back(p);
}

void ParticleSystem::step(float dt, const Vec3& gravity) {
    // Simplified particle system step
    for (auto& particle : particles) {
        if (particle.inv_mass > 0.0f) {
            particle.acceleration = gravity * config.gravity_scale;
            particle.velocity += particle.acceleration * dt;
            particle.position += particle.velocity * dt;
        }
    }
}

} // namespace physics
} // namespace atlas
