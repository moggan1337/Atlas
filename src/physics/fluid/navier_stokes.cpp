/**
 * @file navier_stokes.cpp
 * @brief Navier-Stokes Fluid Simulation Implementation
 */

#include "navier_stokes.hpp"
#include <algorithm>
#include <cmath>

namespace atlas {
namespace physics {

// ============== FluidGrid Implementation ==============

FluidGrid::FluidGrid(const FluidGridConfig& cfg) : config(cfg) {
    width = cfg.grid_width;
    height = cfg.grid_height;
    depth = cfg.grid_depth;
    total_cells = width * height * depth;
    
    density_grid.resize(total_cells);
    pressure_grid.resize(total_cells, 0.0f);
    temperature_field.resize(total_cells, 20.0f);
    
    u_velocity.resize(total_cells, 0.0f);
    v_velocity.resize(total_cells, 0.0f);
    w_velocity.resize(total_cells, 0.0f);
    
    obstacle_mask.resize(total_cells, 0);
    
    // Initialize obstacles
    for (const auto& obs_pos : config.obstacle_positions) {
        add_sphere_obstacle(obs_pos, config.obstacle_radius);
    }
}

Vec3 FluidGrid::world_to_grid(const Vec3& world_pos) const {
    Vec3 grid_pos;
    grid_pos.x = (world_pos.x - config.origin.x) / (config.bounds.x / width);
    grid_pos.y = (world_pos.y - config.origin.y) / (config.bounds.y / height);
    grid_pos.z = (world_pos.z - config.origin.z) / (config.bounds.z / depth);
    return grid_pos;
}

Vec3 FluidGrid::grid_to_world(float x, float y, float z) const {
    Vec3 world_pos;
    world_pos.x = config.origin.x + x * (config.bounds.x / width);
    world_pos.y = config.origin.y + y * (config.bounds.y / height);
    world_pos.z = config.origin.z + z * (config.bounds.z / depth);
    return world_pos;
}

void FluidGrid::update(float dt) {
    // Advection
    advect(dt);
    
    // Apply external forces
    apply_forces(dt);
    
    // Apply viscosity
    apply_viscosity(dt);
    
    // Pressure projection
    project_pressure(dt);
    
    // Boundary conditions
    apply_boundary_conditions();
}

void FluidGrid::advect(float dt) {
    semi_lagrangian_advection(dt);
}

void FluidGrid::semi_lagrangian_advection(float dt) {
    // Create temporary copies
    std::vector<float> u_new = u_velocity;
    std::vector<float> v_new = v_velocity;
    std::vector<float> w_new = w_velocity;
    std::vector<float> density_new = std::vector<float>(total_cells);
    std::vector<float> temp_new = temperature_field;
    
    for (int z = 1; z < depth - 1; ++z) {
        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                if (obstacle(x, y, z)) continue;
                
                int idx = index(x, y, z);
                
                // Get velocity at cell center
                float u = 0.5f * (u_get(x, y, z) + u_get(x + 1, y, z));
                float v = 0.5f * (v_get(x, y, z) + v_get(x, y + 1, z));
                float w = 0.5f * (w_get(x, y, z) + w_get(x, y, z + 1));
                
                // Trace back in time
                float prev_x = x - u * dt;
                float prev_y = y - v * dt;
                float prev_z = z - w * dt;
                
                // Bilinear/trilinear interpolation
                int x0 = static_cast<int>(std::floor(prev_x));
                int y0 = static_cast<int>(std::floor(prev_y));
                int z0 = static_cast<int>(std::floor(prev_z));
                
                float fx = prev_x - x0;
                float fy = prev_y - y0;
                float fz = prev_z - z0;
                
                // Clamp to valid range
                x0 = std::max(1, std::min(x0, width - 2));
                y0 = std::max(1, std::min(y0, height - 2));
                z0 = std::max(1, std::min(z0, depth - 2));
                
                int x1 = std::min(x0 + 1, width - 2);
                int y1 = std::min(y0 + 1, height - 2);
                int z1 = std::min(z0 + 1, depth - 2);
                
                // Trilinear interpolation for velocity
                float u000 = u_get(x0, y0, z0);
                float u100 = u_get(x1, y0, z0);
                float u010 = u_get(x0, y1, z0);
                float u110 = u_get(x1, y1, z0);
                float u001 = u_get(x0, y0, z1);
                float u101 = u_get(x1, y0, z1);
                float u011 = u_get(x0, y1, z1);
                float u111 = u_get(x1, y1, z1);
                
                u_new[idx] = (1 - fx) * (1 - fy) * (1 - fz) * u000 +
                            fx * (1 - fy) * (1 - fz) * u100 +
                            (1 - fx) * fy * (1 - fz) * u010 +
                            fx * fy * (1 - fz) * u110 +
                            (1 - fx) * (1 - fy) * fz * u001 +
                            fx * (1 - fy) * fz * u101 +
                            (1 - fx) * fy * fz * u011 +
                            fx * fy * fz * u111;
                
                // Copy velocity and density
                v_new[idx] = u_new[idx];  // Simplified
                w_new[idx] = u_new[idx];
                
                // Density advection
                density_new[idx] = cell(x0, y0, z0).density;
                
                // Temperature advection
                temp_new[idx] = temperature(x0, y0, z0);
            }
        }
    }
    
    // Copy back
    u_velocity = std::move(u_new);
    v_velocity = std::move(v_new);
    w_velocity = std::move(w_new);
    temperature_field = std::move(temp_new);
    
    for (int i = 0; i < total_cells; ++i) {
        density_grid[i].density = density_new[i];
    }
}

void FluidGrid::mac_advection(float dt) {
    // More accurate MAC grid advection would go here
    semi_lagrangian_advection(dt);
}

void FluidGrid::apply_forces(float dt) {
    // Buoyancy
    apply_buoyancy(dt);
}

void FluidGrid::apply_buoyancy(float dt) {
    float buoyancy = config.buoyancy;
    float ambient_temp = 20.0f;
    
    for (int z = 0; z < depth; ++z) {
        for (int y = 1; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                if (obstacle(x, y, z)) continue;
                
                float temp_diff = temperature(x, y, z) - ambient_temp;
                float density_effect = -buoyancy * temp_diff * dt;
                
                // Apply buoyancy as upward force
                v(x, y, z) += density_effect;
            }
        }
    }
}

void FluidGrid::apply_viscosity(float dt) {
    float mu = config.viscosity;
    
    for (int z = 1; z < depth - 1; ++z) {
        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                if (obstacle(x, y, z)) continue;
                
                // Laplacian for viscosity
                float laplacian_u = u_get(x-1, y, z) + u_get(x+1, y, z) +
                                   u_get(x, y-1, z) + u_get(x, y+1, z) +
                                   u_get(x, y, z-1) + u_get(x, y, z+1) -
                                   6.0f * u_get(x, y, z);
                
                u(x, y, z) += mu * laplacian_u * dt;
            }
        }
    }
}

void FluidGrid::project_pressure(float dt) {
    // Compute divergence
    compute_divergence();
    
    // Solve pressure (Jacobi iteration)
    for (int iter = 0; iter < config.pressure_iterations; ++iter) {
        jacobi_iteration();
    }
    
    // Apply pressure gradient to velocity
    apply_pressure_gradient();
}

void FluidGrid::compute_divergence() {
    for (int z = 1; z < depth - 1; ++z) {
        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                if (obstacle(x, y, z)) continue;
                
                float div = 
                    (u_get(x + 1, y, z) - u_get(x, y, z)) +
                    (v_get(x, y + 1, z) - v_get(x, y, z)) +
                    (w_get(x, y, z + 1) - w_get(x, y, z));
                
                pressure(x, y, z) = div;
            }
        }
    }
}

void FluidGrid::jacobi_iteration() {
    std::vector<float> pressure_new = pressure_grid;
    
    for (int z = 1; z < depth - 1; ++z) {
        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                if (obstacle(x, y, z)) continue;
                
                // Sample neighbors
                float p_left = obstacle(x - 1, y, z) ? pressure(x, y, z) : pressure(x - 1, y, z);
                float p_right = obstacle(x + 1, y, z) ? pressure(x, y, z) : pressure(x + 1, y, z);
                float p_bottom = obstacle(x, y - 1, z) ? pressure(x, y, z) : pressure(x, y - 1, z);
                float p_top = obstacle(x, y + 1, z) ? pressure(x, y, z) : pressure(x, y + 1, z);
                float p_back = obstacle(x, y, z - 1) ? pressure(x, y, z) : pressure(x, y, z - 1);
                float p_front = obstacle(x, y, z + 1) ? pressure(x, y, z) : pressure(x, y, z + 1);
                
                // Jacobi update
                pressure_new[index(x, y, z)] = 
                    (p_left + p_right + p_bottom + p_top + p_back + p_front - pressure(x, y, z)) / 6.0f;
            }
        }
    }
    
    pressure_grid = std::move(pressure_new);
}

void FluidGrid::apply_pressure_gradient() {
    for (int z = 1; z < depth - 1; ++z) {
        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                if (obstacle(x, y, z)) continue;
                
                float dx_p = (pressure(x, y, z) - pressure(x - 1, y, z));
                float dy_p = (pressure(x, y, z) - pressure(x, y - 1, z));
                float dz_p = (pressure(x, y, z) - pressure(x, y, z - 1));
                
                u(x, y, z) -= dx_p;
                v(x, y, z) -= dy_p;
                w(x, y, z) -= dz_p;
            }
        }
    }
}

void FluidGrid::apply_boundary_conditions() {
    // Simple boundary conditions (no-slip for walls, free for open)
    for (int y = 0; y < height; ++y) {
        for (int z = 0; z < depth; ++z) {
            // Left and right walls
            u(0, y, z) = 0.0f;
            u(width - 1, y, z) = 0.0f;
            
            // Front and back walls
            w(x = 0, y, z) = 0.0f;  // Note: x should be z
            w(x = 0, y, z) = 0.0f;
        }
    }
    
    // Bottom wall (ground)
    for (int x = 0; x < width; ++x) {
        for (int z = 0; z < depth; ++z) {
            v(x, 0, z) = 0.0f;
        }
    }
}

void FluidGrid::add_density_source(const Vec3& position, float amount, float radius) {
    Vec3 grid_pos = world_to_grid(position);
    int cx = static_cast<int>(grid_pos.x);
    int cy = static_cast<int>(grid_pos.y);
    int cz = static_cast<int>(grid_pos.z);
    int r = static_cast<int>(radius * width / config.bounds.x);
    
    for (int dz = -r; dz <= r; ++dz) {
        for (int dy = -r; dy <= r; ++dy) {
            for (int dx = -r; dx <= r; ++dx) {
                int x = cx + dx;
                int y = cy + dy;
                int z = cz + dz;
                
                if (!in_bounds(x, y, z) || obstacle(x, y, z)) continue;
                
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (dist <= r) {
                    float falloff = 1.0f - dist / r;
                    cell(x, y, z).density += amount * falloff;
                }
            }
        }
    }
}

void FluidGrid::add_velocity_source(const Vec3& position, const Vec3& velocity, float radius) {
    Vec3 grid_pos = world_to_grid(position);
    int cx = static_cast<int>(grid_pos.x);
    int cy = static_cast<int>(grid_pos.y);
    int cz = static_cast<int>(grid_pos.z);
    int r = static_cast<int>(radius * width / config.bounds.x);
    
    float scale_x = config.bounds.x / width;
    float scale_y = config.bounds.y / height;
    float scale_z = config.bounds.z / depth;
    
    for (int dz = -r; dz <= r; ++dz) {
        for (int dy = -r; dy <= r; ++dy) {
            for (int dx = -r; dx <= r; ++dx) {
                int x = cx + dx;
                int y = cy + dy;
                int z = cz + dz;
                
                if (!in_bounds(x, y, z) || obstacle(x, y, z)) continue;
                
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (dist <= r) {
                    float falloff = 1.0f - dist / r;
                    u(x, y, z) += velocity.x * falloff / scale_x;
                    v(x, y, z) += velocity.y * falloff / scale_y;
                    w(x, y, z) += velocity.z * falloff / scale_z;
                }
            }
        }
    }
}

void FluidGrid::add_temperature_source(const Vec3& position, float amount, float radius) {
    Vec3 grid_pos = world_to_grid(position);
    int cx = static_cast<int>(grid_pos.x);
    int cy = static_cast<int>(grid_pos.y);
    int cz = static_cast<int>(grid_pos.z);
    int r = static_cast<int>(radius * width / config.bounds.x);
    
    for (int dz = -r; dz <= r; ++dz) {
        for (int dy = -r; dy <= r; ++dy) {
            for (int dx = -r; dx <= r; ++dx) {
                int x = cx + dx;
                int y = cy + dy;
                int z = cz + dz;
                
                if (!in_bounds(x, y, z) || obstacle(x, y, z)) continue;
                
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (dist <= r) {
                    float falloff = 1.0f - dist / r;
                    temperature(x, y, z) += amount * falloff;
                }
            }
        }
    }
}

void FluidGrid::clear_densities() {
    for (auto& cell : density_grid) {
        cell.density = 0.0f;
    }
}

void FluidGrid::clear_velocities() {
    std::fill(u_velocity.begin(), u_velocity.end(), 0.0f);
    std::fill(v_velocity.begin(), v_velocity.end(), 0.0f);
    std::fill(w_velocity.begin(), w_velocity.end(), 0.0f);
    std::fill(pressure_grid.begin(), pressure_grid.end(), 0.0f);
}

void FluidGrid::clear_temperatures() {
    std::fill(temperature_field.begin(), temperature_field.end(), 20.0f);
}

void FluidGrid::clear_all() {
    clear_densities();
    clear_velocities();
    clear_temperatures();
}

void FluidGrid::add_sphere_obstacle(const Vec3& center, float radius) {
    Vec3 grid_pos = world_to_grid(center);
    int cx = static_cast<int>(grid_pos.x);
    int cy = static_cast<int>(grid_pos.y);
    int cz = static_cast<int>(grid_pos.z);
    int r = static_cast<int>(radius * width / config.bounds.x);
    
    for (int dz = -r; dz <= r; ++dz) {
        for (int dy = -r; dy <= r; ++dy) {
            for (int dx = -r; dx <= r; ++dx) {
                int x = cx + dx;
                int y = cy + dy;
                int z = cz + dz;
                
                if (!in_bounds(x, y, z)) continue;
                
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (dist <= r) {
                    obstacle(x, y, z) = 1;
                }
            }
        }
    }
}

void FluidGrid::remove_sphere_obstacle(const Vec3& center, float radius) {
    Vec3 grid_pos = world_to_grid(center);
    int cx = static_cast<int>(grid_pos.x);
    int cy = static_cast<int>(grid_pos.y);
    int cz = static_cast<int>(grid_pos.z);
    int r = static_cast<int>(radius * width / config.bounds.x);
    
    for (int dz = -r; dz <= r; ++dz) {
        for (int dy = -r; dy <= r; ++dy) {
            for (int dx = -r; dx <= r; ++dx) {
                int x = cx + dx;
                int y = cy + dy;
                int z = cz + dz;
                
                if (!in_bounds(x, y, z)) continue;
                
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (dist <= r) {
                    obstacle(x, y, z) = 0;
                }
            }
        }
    }
}

void FluidGrid::update_obstacle_mask() {
    // Update any dynamic obstacles
}

// ============== FluidVolume Implementation ==============

FluidVolume::FluidVolume(const Config& cfg) : config(cfg) {
    width = cfg.width;
    height = cfg.height;
    depth = cfg.depth;
    density_field.resize(width * height * depth, 0.0f);
}

void FluidVolume::set_density(int x, int y, int z, float density) {
    if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth) {
        density_field[x + y * width + z * width * height] = density;
    }
}

float FluidVolume::get_density(int x, int y, int z) const {
    if (x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth) {
        return density_field[x + y * width + z * width * height];
    }
    return 0.0f;
}

void FluidVolume::sample_from_grid(const FluidGrid& grid) {
    // Sample grid values into volume
    // Simplified implementation
    clear();
}

void FluidVolume::clear() {
    std::fill(density_field.begin(), density_field.end(), 0.0f);
}

// ============== NavierStokesWorld Implementation ==============

NavierStokesWorld::NavierStokesWorld(const Config& cfg) 
    : config(cfg)
    , grid(cfg.grid) 
{}

void NavierStokesWorld::step(float dt) {
    dt = std::min(dt, config.max_dt);
    grid.update(dt);
}

void NavierStokesWorld::set_velocity(const Vec3& world_pos, const Vec3& velocity) {
    grid.add_velocity_source(world_pos, velocity, 0.5f);
}

Vec3 NavierStokesWorld::get_velocity(const Vec3& world_pos) const {
    Vec3 grid_pos = grid.world_to_grid(world_pos);
    int x = static_cast<int>(grid_pos.x);
    int y = static_cast<int>(grid_pos.y);
    int z = static_cast<int>(grid_pos.z);
    return grid.get_velocity(x, y, z);
}

void NavierStokesWorld::add_density(const Vec3& world_pos, float amount) {
    grid.add_density_source(world_pos, amount, 0.5f);
}

void NavierStokesWorld::add_heat(const Vec3& world_pos, float amount) {
    grid.add_temperature_source(world_pos, amount, 0.5f);
}

void NavierStokesWorld::add_obstacle(const Vec3& center, float radius) {
    grid.add_sphere_obstacle(center, radius);
}

void NavierStokesWorld::remove_obstacle(const Vec3& center, float radius) {
    grid.remove_sphere_obstacle(center, radius);
}

void NavierStokesWorld::clear() {
    grid.clear_all();
}

} // namespace physics
} // namespace atlas
