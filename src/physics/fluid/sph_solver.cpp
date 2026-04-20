/**
 * @file sph_solver.cpp
 * @brief SPH Solver Implementation
 */

#include "sph_solver.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace atlas {
namespace physics {

// ============== SPHSolver Implementation ==============

SPHSolver::SPHSolver(const SPHConfig& cfg) : config(cfg) {
    gravity = Vec3(0.0f, -9.81f, 0.0f);
    kernel_radius_sq = config.kernel_radius * config.kernel_radius;
    avg_particle_density = config.rest_density;
}

int SPHSolver::add_particle(const Vec3& position) {
    int index = static_cast<int>(particles.size());
    particles.emplace_back(position, config.particle_mass);
    neighbor_lists.emplace_back();
    return index;
}

int SPHSolver::add_particle(const Vec3& position, const Vec3& velocity) {
    int index = add_particle(position);
    particles[index].velocity = velocity;
    return index;
}

void SPHSolver::remove_particle(int index) {
    if (index >= 0 && index < static_cast<int>(particles.size())) {
        particles.erase(particles.begin() + index);
        neighbor_lists.erase(neighbor_lists.begin() + index);
    }
}

void SPHSolver::clear() {
    particles.clear();
    neighbor_lists.clear();
}

void SPHSolver::update(float dt) {
    if (particles.empty()) return;
    
    // Find neighbors
    find_neighbors();
    
    // Compute density and pressure
    compute_density();
    compute_pressure();
    
    // Compute forces
    compute_forces(dt);
    
    // Integrate
    integrate(dt);
    
    // Handle boundaries
    handle_boundaries();
}

void SPHSolver::step_fixed(float dt) {
    update(dt);
}

void SPHSolver::find_neighbors() {
    float h_sq = kernel_radius_sq;
    
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        auto& neighbors = neighbor_lists[i];
        neighbors.clear();
        
        const Vec3& pos_i = particles[i].position;
        
        for (int j = 0; j < static_cast<int>(particles.size()); ++j) {
            if (i == j) continue;
            
            Vec3 delta = particles[j].position - pos_i;
            if (delta.magnitude_squared() < h_sq) {
                neighbors.push_back(j);
                
                // Limit neighbors for performance
                if (neighbors.size() >= config.max_neighbors) break;
            }
        }
    }
}

void SPHSolver::compute_density() {
    float h_sq = kernel_radius_sq;
    float particle_mass = config.particle_mass;
    
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        float density = 0.0f;
        const Vec3& pos_i = particles[i].position;
        
        // Self-contribution
        density += particle_mass * kernel_poly6(0.0f, config.kernel_radius);
        
        // Neighbor contributions
        for (int j : neighbor_lists[i]) {
            Vec3 delta = particles[j].position - pos_i;
            float r_sq = delta.magnitude_squared();
            density += particle_mass * kernel_poly6(r_sq, config.kernel_radius);
        }
        
        particles[i].density = std::max(density, 1.0f);
    }
    
    // Calculate average density
    if (!particles.empty()) {
        float total_density = 0.0f;
        for (const auto& p : particles) {
            total_density += p.density;
        }
        avg_particle_density = total_density / particles.size();
    }
}

void SPHSolver::compute_pressure() {
    float rest_density = config.rest_density;
    float gas_constant = config.gas_constant;
    
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        // Equation of state: P = k * (rho - rho_0)
        float density_error = particles[i].density - rest_density;
        particles[i].pressure = gas_constant * density_error;
        
        // Clamp pressure
        particles[i].pressure = std::max(particles[i].pressure, 0.0f);
        
        // Color based on pressure for rendering
        particles[i].pressure_color = std::min(1.0f, particles[i].pressure / 5000.0f);
    }
}

void SPHSolver::compute_forces(float dt) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        SPHParticle& particle = particles[i];
        
        // Gravity
        particle.force = gravity * particle.density * config.gravity_scale;
        
        // Pressure force
        particle.force += compute_pressure_force(particle, neighbor_lists[i]);
        
        // Viscosity force
        particle.force += compute_viscosity_force(particle, neighbor_lists[i]);
        
        // Surface tension
        particle.force += compute_surface_tension(particle, neighbor_lists[i]);
        
        // Compute acceleration
        if (particle.density > 0.0f) {
            particle.acceleration = particle.force / particle.density;
        } else {
            particle.acceleration = Vec3::zero();
        }
    }
}

Vec3 SPHSolver::compute_pressure_force(const SPHParticle& particle, 
                                        const std::vector<int>& neighbors) {
    float h = config.kernel_radius;
    float particle_mass = config.particle_mass;
    Vec3 pressure_force = Vec3::zero();
    
    for (int j : neighbors) {
        const SPHParticle& neighbor = particles[j];
        
        if (neighbor.density < 1e-6f) continue;
        
        Vec3 delta = particle.position - neighbor.position;
        float r = delta.magnitude();
        
        if (r < h && r > 1e-6f) {
            // Pressure contribution
            float pressure = (particle.pressure + neighbor.pressure) / (2.0f * neighbor.density);
            
            // Spiky kernel gradient
            float kernel_grad = -45.0f / (M_PI * std::pow(h, 6.0f)) * 
                               std::pow(h - r, 2.0f);
            
            pressure_force += delta * (particle_mass * pressure * kernel_grad / r);
        }
    }
    
    return pressure_force;
}

Vec3 SPHSolver::compute_viscosity_force(const SPHParticle& particle,
                                         const std::vector<int>& neighbors) {
    float h = config.kernel_radius;
    float mu = config.viscosity;
    float particle_mass = config.particle_mass;
    Vec3 viscosity_force = Vec3::zero();
    
    for (int j : neighbors) {
        const SPHParticle& neighbor = particles[j];
        
        if (neighbor.density < 1e-6f) continue;
        
        Vec3 delta = particle.position - neighbor.position;
        float r = delta.magnitude();
        
        if (r < h && r > 1e-6f) {
            // Viscosity kernel Laplacian
            float kernel_lap = 45.0f / (M_PI * std::pow(h, 6.0f)) * (h - r);
            
            Vec3 velocity_diff = neighbor.velocity - particle.velocity;
            viscosity_force += velocity_diff * (particle_mass * mu * kernel_lap / (neighbor.density * r));
        }
    }
    
    return viscosity_force;
}

Vec3 SPHSolver::compute_surface_tension(const SPHParticle& particle,
                                         const std::vector<int>& neighbors) {
    // Simplified surface tension - move toward center of mass of neighbors
    float h = config.kernel_radius;
    float gamma = config.surface_tension;
    float particle_mass = config.particle_mass;
    
    Vec3 tension_force = Vec3::zero();
    Vec3 color_field_gradient = Vec3::zero();
    float color_field_laplacian = 0.0f;
    
    for (int j : neighbors) {
        const SPHParticle& neighbor = particles[j];
        
        if (neighbor.density < 1e-6f) continue;
        
        Vec3 delta = particle.position - neighbor.position;
        float r = delta.magnitude();
        
        if (r < h && r > 1e-6f) {
            float r_sq = r * r;
            
            // Color field (normalized density)
            float density_ratio = particle_mass / (neighbor.density * r_sq);
            
            // Gradient
            float kernel_deriv = -945.0f / (32.0f * M_PI * std::pow(h, 9.0f)) * 
                                std::pow(h * h - r_sq, 2.0f) * r;
            color_field_gradient += delta * density_ratio * kernel_deriv;
            
            // Laplacian
            float kernel_lap = 945.0f / (32.0f * M_PI * std::pow(h, 9.0f)) * 
                              (h * h - r_sq) * (3.0f * h * h - 7.0f * r_sq);
            color_field_laplacian += density_ratio * kernel_lap;
        }
    }
    
    // Surface tension force
    if (color_field_gradient.magnitude() > 1e-6f) {
        tension_force = -gamma * color_field_laplacian * color_field_gradient.normalized();
    }
    
    return tension_force;
}

void SPHSolver::integrate(float dt) {
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        SPHParticle& p = particles[i];
        
        // Semi-implicit Euler
        p.velocity += p.acceleration * dt;
        
        // Speed limit
        float speed = p.velocity.magnitude();
        if (speed > config.max_speed) {
            p.velocity *= config.max_speed / speed;
        }
        
        // Update position
        p.position += p.velocity * dt;
    }
}

void SPHSolver::handle_boundaries() {
    float margin = config.bounds_margin;
    
    for (auto& p : particles) {
        // Simple boundary handling
        if (p.position.x < bounds_min.x + margin) {
            p.position.x = bounds_min.x + margin;
            p.velocity.x *= -0.3f;  // Damping on bounce
        }
        if (p.position.x > bounds_max.x - margin) {
            p.position.x = bounds_max.x - margin;
            p.velocity.x *= -0.3f;
        }
        
        if (p.position.y < bounds_min.y + margin) {
            p.position.y = bounds_min.y + margin;
            p.velocity.y *= -0.3f;
        }
        if (p.position.y > bounds_max.y - margin) {
            p.position.y = bounds_max.y - margin;
            p.velocity.y *= -0.3f;
        }
        
        if (p.position.z < bounds_min.z + margin) {
            p.position.z = bounds_min.z + margin;
            p.velocity.z *= -0.3f;
        }
        if (p.position.z > bounds_max.z - margin) {
            p.position.z = bounds_max.z - margin;
            p.velocity.z *= -0.3f;
        }
    }
}

float SPHSolver::kernel_poly6(float r_squared, float h) {
    if (r_squared >= h * h) return 0.0f;
    
    float h_sq = h * h;
    float diff = h_sq - r_squared;
    
    return 315.0f / (64.0f * M_PI * std::pow(h, 9.0f)) * diff * diff * diff;
}

float SPHSolver::kernel_spiky_gradient(float r, float h) {
    if (r >= h || r < 1e-6f) return 0.0f;
    
    float diff = h - r;
    return -45.0f / (M_PI * std::pow(h, 6.0f)) * diff * diff;
}

float SPHSolver::kernel_viscosity_laplacian(float r, float h) {
    if (r >= h) return 0.0f;
    
    return 45.0f / (M_PI * std::pow(h, 6.0f)) * (h - r);
}

Vec3 SPHSolver::kernel_poly6_gradient(const Vec3& r, float r_length, float h) {
    if (r_length >= h || r_length < 1e-6f) return Vec3::zero();
    
    float h_sq = h * h;
    float r_sq = r_length * r_length;
    float diff = h_sq - r_sq;
    float coeff = -945.0f / (32.0f * M_PI * std::pow(h, 9.0f)) * diff * diff;
    
    return r * coeff;
}

void SPHSolver::emit_particles(const Vec3& position, int count, const Vec3& velocity) {
    for (int i = 0; i < count; ++i) {
        Vec3 offset(
            (static_cast<float>(rand()) / RAND_MAX - 0.5f) * 0.1f,
            (static_cast<float>(rand()) / RAND_MAX - 0.5f) * 0.1f,
            (static_cast<float>(rand()) / RAND_MAX - 0.5f) * 0.1f
        );
        add_particle(position + offset, velocity);
    }
}

void SPHSolver::emit_particle_grid(const Vec3& origin, const Vec3& size, float spacing,
                                   const Vec3& initial_velocity) {
    int count_x = static_cast<int>(size.x / spacing) + 1;
    int count_y = static_cast<int>(size.y / spacing) + 1;
    int count_z = static_cast<int>(size.z / spacing) + 1;
    
    for (int z = 0; z < count_z; ++z) {
        for (int y = 0; y < count_y; ++y) {
            for (int x = 0; x < count_x; ++x) {
                Vec3 pos = origin + Vec3(x * spacing, y * spacing, z * spacing);
                add_particle(pos, initial_velocity);
            }
        }
    }
}

std::vector<int> SPHSolver::query_neighbors(const Vec3& point, float radius) {
    std::vector<int> neighbors;
    float radius_sq = radius * radius;
    
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        Vec3 delta = particles[i].position - point;
        if (delta.magnitude_squared() < radius_sq) {
            neighbors.push_back(i);
        }
    }
    
    return neighbors;
}

float SPHSolver::compute_density_at(const Vec3& point) {
    float density = 0.0f;
    float particle_mass = config.particle_mass;
    float h = config.kernel_radius;
    float h_sq = h * h;
    
    for (const auto& p : particles) {
        Vec3 delta = p.position - point;
        float r_sq = delta.magnitude_squared();
        
        if (r_sq < h_sq) {
            density += particle_mass * kernel_poly6(r_sq, h);
        }
    }
    
    return density;
}

void SPHSolver::print_stats() {
    printf("SPH Particles: %zu\n", particles.size());
    printf("Average Density: %.2f (target: %.2f)\n", 
           avg_particle_density, config.rest_density);
    
    if (!particles.empty()) {
        float total_ke = get_total_kinetic_energy();
        printf("Total Kinetic Energy: %.2f\n", total_ke);
    }
}

float SPHSolver::get_total_kinetic_energy() const {
    float ke = 0.0f;
    for (const auto& p : particles) {
        ke += 0.5f * p.mass * p.velocity.magnitude_squared();
    }
    return ke;
}

std::vector<Vec3> SPHSolver::get_particle_positions() const {
    std::vector<Vec3> positions;
    positions.reserve(particles.size());
    for (const auto& p : particles) {
        positions.push_back(p.position);
    }
    return positions;
}

std::vector<float> SPHSolver::get_particle_colors() const {
    std::vector<float> colors;
    colors.reserve(particles.size());
    for (const auto& p : particles) {
        colors.push_back(p.pressure_color);
    }
    return colors;
}

// ============== SPHWorld Implementation ==============

SPHWorld::SPHWorld(const Config& cfg) 
    : config(cfg)
    , solver(cfg.sph)
    , emission_accumulator(0.0f) 
{
    solver.set_gravity(cfg.gravity);
    solver.set_bounds(cfg.bounds_min, cfg.bounds_max);
}

void SPHWorld::update(float dt) {
    // Handle sources
    for (auto& source : sources) {
        if (!source.active) continue;
        
        emission_accumulator += source.emission_rate * dt;
        int to_emit = static_cast<int>(emission_accumulator);
        
        if (to_emit > 0) {
            solver.emit_particle(source.position, to_emit, source.velocity);
            emission_accumulator -= to_emit;
        }
    }
    
    // Handle sinks
    for (auto& sink : sinks) {
        if (!sink.active) continue;
        
        // Remove particles in sink
        auto& particles = const_cast<std::vector<SPHParticle>&>(solver.get_particle(0));  // Hack
    }
    
    // Update fluid
    solver.update(dt);
}

} // namespace physics
} // namespace atlas
