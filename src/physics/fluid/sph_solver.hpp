/**
 * @file sph_solver.hpp
 * @brief SPH (Smoothed Particle Hydrodynamics) Fluid Simulation
 * 
 * Implements particle-based fluid dynamics:
 * - Kernel functions (Poly6, Spiky, Viscosity)
 * - Pressure and viscosity forces
 * - Surface tension
 * - Fluid-rigid body interaction
 */

#pragma once

#include "../../core/vec3.hpp"
#include <vector>
#include <memory>

namespace atlas {
namespace physics {

using core::Vec3;

/**
 * @brief SPH Particle
 */
struct SPHParticle {
    Vec3 position;
    Vec3 velocity;
    Vec3 acceleration;
    Vec3 force;
    float density;
    float pressure;
    float mass;
    float density_error;  // For pressure solver
    
    // Color for rendering
    float pressure_color;  // 0-1 based on pressure
    
    SPHParticle() 
        : position(Vec3::zero())
        , velocity(Vec3::zero())
        , acceleration(Vec3::zero())
        , force(Vec3::zero())
        , density(0.0f)
        , pressure(0.0f)
        , mass(0.01f)
        , density_error(0.0f)
        , pressure_color(0.0f) {}
    
    SPHParticle(const Vec3& pos, float m = 0.01f) 
        : position(pos)
        , velocity(Vec3::zero())
        , acceleration(Vec3::zero())
        , force(Vec3::zero())
        , density(0.0f)
        , pressure(0.0f)
        , mass(m)
        , density_error(0.0f)
        , pressure_color(0.0f) {}
};

/**
 * @brief SPH Configuration
 */
struct SPHConfig {
    float particle_radius = 0.1f;
    float particle_mass = 0.01f;
    float rest_density = 1000.0f;
    float gas_constant = 2000.0f;      // Pressure constant
    float viscosity = 0.018f;          // Viscosity coefficient
    float surface_tension = 0.0728f;   // Surface tension coefficient
    float kernel_radius = 0.3f;         // Smoothing radius
    float gravity_scale = 1.0f;
    
    // Solver params
    int max_neighbors = 100;
    float max_speed = 10.0f;
    float bounds_margin = 0.5f;
};

/**
 * @brief SPH Solver
 */
class SPHSolver {
public:
    SPHConfig config;
    
private:
    std::vector<SPHParticle> particles;
    std::vector<std::vector<int>> neighbor_lists;
    Vec3 gravity;
    Vec3 bounds_min;
    Vec3 bounds_max;
    
    // Performance metrics
    float avg_particle_density;
    
public:
    SPHSolver(const SPHConfig& cfg = SPHConfig());
    ~SPHSolver() = default;
    
    // Particle management
    int add_particle(const Vec3& position);
    int add_particle(const Vec3& position, const Vec3& velocity);
    void remove_particle(int index);
    void clear();
    
    SPHParticle& get_particle(int index) { return particles[index]; }
    const SPHParticle& get_particle(int index) const { return particles[index]; }
    size_t get_particle_count() const { return particles.size(); }
    
    // Setup
    void set_gravity(const Vec3& g) { gravity = g; }
    Vec3 get_gravity() const { return gravity; }
    void set_bounds(const Vec3& min, const Vec3& max) { bounds_min = min; bounds_max = max; }
    
    // Simulation
    void update(float dt);
    void step_fixed(float dt);
    
    // Simulation steps
    void find_neighbors();
    void compute_density();
    void compute_pressure();
    void compute_forces(float dt);
    void integrate(float dt);
    void handle_boundaries();
    
    // Force computations
    Vec3 compute_pressure_force(const SPHParticle& particle, const std::vector<int>& neighbors);
    Vec3 compute_viscosity_force(const SPHParticle& particle, const std::vector<int>& neighbors);
    Vec3 compute_surface_tension(const SPHParticle& particle, const std::vector<int>& neighbors);
    
    // Kernel functions
    float kernel_poly6(float r_squared, float h);
    float kernel_spiky_gradient(float r, float h);
    float kernel_viscosity_laplacian(float r, float h);
    Vec3 kernel_poly6_gradient(const Vec3& r, float r_length, float h);
    
    // Particle emission
    void emit_particles(const Vec3& position, int count, const Vec3& velocity = Vec3::zero());
    void emit_particle_grid(const Vec3& origin, const Vec3& size, float spacing, 
                           const Vec3& initial_velocity = Vec3::zero());
    
    // Queries
    std::vector<int> query_neighbors(const Vec3& point, float radius);
    float compute_density_at(const Vec3& point);
    
    // Utilities
    void print_stats();
    float get_average_density() const { return avg_particle_density; }
    float get_total_kinetic_energy() const;
    
    // Rendering data
    std::vector<Vec3> get_particle_positions() const;
    std::vector<float> get_particle_colors() const;
    
private:
    float kernel_radius_sq;
};

/**
 * @brief Fluid Source (emitter)
 */
class FluidSource {
public:
    Vec3 position;
    Vec3 velocity;
    float emission_rate;  // particles per second
    float particle_mass;
    bool active;
    
    FluidSource() 
        : position(Vec3::zero())
        , velocity(0.0f, -1.0f, 0.0f)
        , emission_rate(100.0f)
        , particle_mass(0.01f)
        , active(true) {}
    
    FluidSource(const Vec3& pos, const Vec3& vel, float rate = 100.0f)
        : position(pos)
        , velocity(vel)
        , emission_rate(rate)
        , particle_mass(0.01f)
        , active(true) {}
};

/**
 * @brief Fluid Sink (drain)
 */
class FluidSink {
public:
    Vec3 position;
    float radius;
    bool active;
    
    FluidSink() 
        : position(Vec3::zero())
        , radius(0.5f)
        , active(true) {}
    
    FluidSink(const Vec3& pos, float r = 0.5f)
        : position(pos)
        , radius(r)
        , active(true) {}
};

/**
 * @brief SPH World with sources and sinks
 */
class SPHWorld {
public:
    struct Config {
        SPHConfig sph;
        Vec3 gravity = Vec3(0.0f, -9.81f, 0.0f);
        Vec3 bounds_min = Vec3(-5.0f, -5.0f, -5.0f);
        Vec3 bounds_max = Vec3(5.0f, 10.0f, 5.0f);
        float max_dt = 1.0f / 60.0f;
    } config;
    
private:
    SPHSolver solver;
    std::vector<FluidSource> sources;
    std::vector<FluidSink> sinks;
    float emission_accumulator;
    
public:
    SPHWorld(const Config& cfg = Config());
    ~SPHWorld() = default;
    
    SPHSolver& get_solver() { return solver; }
    const SPHSolver& get_solver() const { return solver; }
    
    void add_source(const FluidSource& source) { sources.push_back(source); }
    void add_sink(const FluidSink& sink) { sinks.push_back(sink); }
    void remove_source(size_t index);
    void remove_sink(size_t index);
    
    void update(float dt);
    
    size_t get_particle_count() const { return solver.get_particle_count(); }
    
    void clear() { solver.clear(); }
};

} // namespace physics
} // namespace atlas
