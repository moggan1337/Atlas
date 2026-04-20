/**
 * @file navier_stokes.hpp
 * @brief Grid-based Eulerian Fluid Simulation using Navier-Stokes Equations
 * 
 * Implements:
 * - MAC (Marker-and-Cell) grid velocity field
 * - Advection (semi-Lagrangian)
 * - Pressure projection (Jacobi iteration)
 * - Buoyancy and external forces
 */

#pragma once

#include "../../core/vec3.hpp"
#include <vector>
#include <memory>

namespace atlas {
namespace physics {

using core::Vec3;

/**
 * @brief MAC Grid Face
 * Velocities are stored on cell faces for staggered grid
 */
struct MACCell {
    // Velocities on faces
    float u_left;     // x- face
    float u_right;    // x+ face
    float v_bottom;   // y- face
    float v_top;      // y+ face
    float w_back;     // z- face
    float w_front;    // z+ face
    
    // Cell-centered quantities
    float pressure;
    float divergence;
    float density;
    float temperature;
    
    MACCell() 
        : u_left(0), u_right(0)
        , v_bottom(0), v_top(0)
        , w_back(0), w_front(0)
        , pressure(0), divergence(0)
        , density(1.0f), temperature(20.0f) {}
};

/**
 * @brief Grid Cell for density/temperature
 */
struct GridCell {
    float density;
    float temperature;
    float pressure;
    float divergence;
    
    GridCell() : density(0.0f), temperature(20.0f), pressure(0.0f), divergence(0.0f) {}
};

/**
 * @brief Eulerian Fluid Grid Configuration
 */
struct FluidGridConfig {
    Vec3 origin = Vec3(0.0f, 0.0f, 0.0f);
    Vec3 bounds = Vec3(10.0f, 10.0f, 10.0f);
    int grid_width = 64;
    int grid_height = 64;
    int grid_depth = 64;
    
    float viscosity = 0.001f;
    float dissipation = 0.001f;
    float buoyancy = 1.0f;
    float temperature_diffusion = 0.001f;
    
    int pressure_iterations = 50;
    float pressure_tolerance = 1e-5f;
    
    // Obstacles
    std::vector<Vec3> obstacle_positions;
    float obstacle_radius = 0.5f;
};

/**
 * @brief Eulerian Fluid Grid
 * 
 * Uses MAC (staggered) grid for velocity storage.
 * Implements semi-Lagrangian advection and pressure projection.
 */
class FluidGrid {
public:
    FluidGridConfig config;
    
private:
    int width, height, depth;
    int total_cells;
    std::vector<GridCell> density_grid;
    std::vector<float> pressure_grid;
    std::vector<Vec3> velocity_field;  // Stored as Vec3 for simplicity
    std::vector<float> temperature_field;
    
    // MAC grid velocities (staggered)
    std::vector<float> u_velocity;  // x-direction
    std::vector<float> v_velocity;  // y-direction
    std::vector<float> w_velocity;  // z-direction
    
    std::vector<char> obstacle_mask;
    
public:
    FluidGrid(const FluidGridConfig& cfg = FluidGridConfig());
    ~FluidGrid() = default;
    
    // Grid indexing
    int index(int x, int y, int z) const {
        return x + y * width + z * width * height;
    }
    
    bool in_bounds(int x, int y, int z) const {
        return x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth;
    }
    
    // Cell access
    GridCell& cell(int x, int y, int z) { return density_grid[index(x, y, z)]; }
    const GridCell& cell(int x, int y, int z) const { return density_grid[index(x, y, z)]; }
    
    float& pressure(int x, int y, int z) { return pressure_grid[index(x, y, z)]; }
    float pressure(int x, int y, int z) const { return pressure_grid[index(x, y, z)]; }
    
    float& u(int x, int y, int z) { return u_velocity[index(x, y, z)]; }
    float& v(int x, int y, int z) { return v_velocity[index(x, y, z)]; }
    float& w(int x, int y, int z) { return w_velocity[index(x, y, z)]; }
    
    float u_get(int x, int y, int z) const { 
        return in_bounds(x, y, z) ? u_velocity[index(x, y, z)] : 0.0f; 
    }
    float v_get(int x, int y, int z) const { 
        return in_bounds(x, y, z) ? v_velocity[index(x, y, z)] : 0.0f; 
    }
    float w_get(int x, int y, int z) const { 
        return in_bounds(x, y, z) ? w_velocity[index(x, y, z)] : 0.0f; 
    }
    
    Vec3 get_velocity(int x, int y, int z) const {
        return Vec3(u_get(x, y, z), v_get(x, y, z), w_get(x, y, z));
    }
    
    float& temperature(int x, int y, int z) { return temperature_field[index(x, y, z)]; }
    float temperature(int x, int y, int z) const { 
        return in_bounds(x, y, z) ? temperature_field[index(x, y, z)] : 0.0f; 
    }
    
    bool& obstacle(int x, int y, int z) { return obstacle_mask[index(x, y, z)]; }
    bool obstacle(int x, int y, int z) const { 
        return in_bounds(x, y, z) ? obstacle_mask[index(x, y, z)] : true; 
    }
    
    // Size queries
    int get_width() const { return width; }
    int get_height() const { return height; }
    int get_depth() const { return depth; }
    int get_total_cells() const { return total_cells; }
    
    // World position to grid
    Vec3 world_to_grid(const Vec3& world_pos) const;
    Vec3 grid_to_world(float x, float y, float z) const;
    
    // Simulation
    void update(float dt);
    void advect(float dt);
    void apply_forces(float dt);
    void project_pressure(float dt);
    void apply_boundary_conditions();
    
    // Advection schemes
    void semi_lagrangian_advection(float dt);
    void mac_advection(float dt);
    
    // Pressure solver
    void jacobi_iteration();
    void compute_divergence();
    void apply_pressure_gradient();
    
    // Forces
    void apply_buoyancy(float dt);
    void apply_viscosity(float dt);
    
    // Sources
    void add_density_source(const Vec3& position, float amount, float radius);
    void add_velocity_source(const Vec3& position, const Vec3& velocity, float radius);
    void add_temperature_source(const Vec3& position, float amount, float radius);
    
    // Clearing
    void clear_densities();
    void clear_velocities();
    void clear_temperatures();
    void clear_all();
    
    // Rendering data
    std::vector<float> get_density_volume() const { return std::vector<float>(); }
    
    // Obstacles
    void add_sphere_obstacle(const Vec3& center, float radius);
    void remove_sphere_obstacle(const Vec3& center, float radius);
    void update_obstacle_mask();
};

/**
 * @brief Fluid Volume for rendering (marching cubes ready)
 */
class FluidVolume {
public:
    struct Config {
        int width = 64;
        int height = 64;
        int depth = 64;
        float isosurface_level = 0.5f;
    } config;
    
private:
    int width, height, depth;
    std::vector<float> density_field;
    
public:
    FluidVolume(const Config& cfg = Config());
    ~FluidVolume() = default;
    
    void set_density(int x, int y, int z, float density);
    float get_density(int x, int y, int z) const;
    
    void sample_from_grid(const FluidGrid& grid);
    
    void clear();
    
    int get_width() const { return width; }
    int get_height() const { return height; }
    int get_depth() const { return depth; }
    
    const float* get_data() const { return density_field.data(); }
};

/**
 * @brief Navier-Stokes World
 */
class NavierStokesWorld {
public:
    struct Config {
        FluidGridConfig grid;
        Vec3 gravity = Vec3(0.0f, -9.81f, 0.0f);
        float max_dt = 1.0f / 60.0f;
    } config;
    
private:
    FluidGrid grid;
    
public:
    NavierStokesWorld(const Config& cfg = Config());
    ~NavierStokesWorld() = default;
    
    FluidGrid& get_grid() { return grid; }
    const FluidGrid& get_grid() const { return grid; }
    
    void step(float dt);
    
    void set_velocity(const Vec3& world_pos, const Vec3& velocity);
    Vec3 get_velocity(const Vec3& world_pos) const;
    
    void add_density(const Vec3& world_pos, float amount);
    void add_heat(const Vec3& world_pos, float amount);
    
    void add_obstacle(const Vec3& center, float radius);
    void remove_obstacle(const Vec3& center, float radius);
    
    void clear();
};

} // namespace physics
} // namespace atlas
