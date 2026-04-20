/**
 * @file benchmark_main.cpp
 * @brief Performance benchmarks for Atlas physics engine
 */

#include "atlas.hpp"
#include <chrono>
#include <iostream>
#include <vector>
#include <random>

using namespace atlas;
using namespace atlas::physics;

// High-resolution timer
class Timer {
public:
    void start() { start_time = std::chrono::high_resolution_clock::now(); }
    void stop() { end_time = std::chrono::high_resolution_clock::now(); }
    double elapsed_ms() const {
        return std::chrono::duration<double, std::milli>(end_time - start_time).count();
    }
private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point end_time;
};

void benchmark_rigid_bodies(int num_bodies) {
    std::cout << "Benchmarking " << num_bodies << " rigid bodies..." << std::endl;
    
    WorldConfig config;
    config.gravity = {0, -9.81f, 0};
    config.substeps = 4;
    
    RigidBodyWorld world(config);
    
    // Create bodies
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-10, 10);
    
    for (int i = 0; i < num_bodies; ++i) {
        auto sphere = RigidBody::create_sphere(0.5f, 1.0f);
        sphere->set_position({dist(rng), 5 + dist(rng), dist(rng)});
        world.add_body(sphere);
    }
    
    // Add ground
    auto ground = RigidBody::create_static_plane({0, 1, 0}, 0);
    world.add_body(ground);
    
    Timer timer;
    timer.start();
    
    // Simulate 1 second
    for (int i = 0; i < 60; ++i) {
        world.step(1.0f / 60.0f);
    }
    
    timer.stop();
    
    std::cout << "  Time: " << timer.elapsed_ms() << " ms" << std::endl;
    std::cout << "  Bodies/sec: " << num_bodies / (timer.elapsed_ms() / 1000.0) << std::endl;
}

void benchmark_soft_body(int particles) {
    std::cout << "Benchmarking soft body with " << particles << " particles..." << std::endl;
    
    SoftBodyConfig config;
    config.stiffness = 500;
    config.damping = 2;
    
    // Create a large cloth
    int segments = static_cast<int>(std::sqrt(particles));
    auto cloth = SoftBody::create_cloth({0, 5, 0}, 5, 5, segments, segments);
    
    SoftBodyWorld world;
    world.config.gravity = {0, -9.81f, 0};
    world.add_body(cloth);
    
    Timer timer;
    timer.start();
    
    // Simulate
    for (int i = 0; i < 120; ++i) {
        world.step(1.0f / 60.0f);
    }
    
    timer.stop();
    
    std::cout << "  Time: " << timer.elapsed_ms() << " ms" << std::endl;
}

void benchmark_sph(int num_particles) {
    std::cout << "Benchmarking SPH with " << num_particles << " particles..." << std::endl;
    
    SPHWorld::Config config;
    config.sph.particle_radius = 0.1f;
    config.sph.rest_density = 1000.0f;
    
    SPHWorld world(config);
    
    // Emit particles
    int side = static_cast<int>(std::cbrt(num_particles));
    world.get_solver().emit_particle_grid(
        {-1, 0, -1}, 
        {2.0f, 2.0f, 2.0f}, 
        0.2f
    );
    
    Timer timer;
    timer.start();
    
    // Simulate
    for (int i = 0; i < 60; ++i) {
        world.update(1.0f / 60.0f);
    }
    
    timer.stop();
    
    std::cout << "  Particles: " << world.get_particle_count() << std::endl;
    std::cout << "  Time: " << timer.elapsed_ms() << " ms" << std::endl;
}

void benchmark_raycast(int num_primitives, int num_queries) {
    std::cout << "Benchmarking raycasting with " << num_primitives << " primitives, " 
              << num_queries << " queries..." << std::endl;
    
    BVHRaycast accelerator;
    
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-100, 100);
    
    std::vector<BoundingVolume> bounds;
    for (int i = 0; i < num_primitives; ++i) {
        bounds.emplace_back(Vec3(dist(rng), dist(rng), dist(rng)), Vec3(1, 1, 1));
    }
    accelerator.set_primitives(bounds);
    accelerator.build();
    
    Timer timer;
    timer.start();
    
    int hits = 0;
    for (int i = 0; i < num_queries; ++i) {
        Ray ray(
            Vec3(dist(rng), dist(rng), dist(rng)),
            Vec3(1, 0, 0).normalized(),
            0, 1000
        );
        RayHit hit;
        if (accelerator.raycast(ray, hit)) hits++;
    }
    
    timer.stop();
    
    std::cout << "  Hits: " << hits << std::endl;
    std::cout << "  Time: " << timer.elapsed_ms() << " ms" << std::endl;
    std::cout << "  Queries/sec: " << num_queries / (timer.elapsed_ms() / 1000.0) << std::endl;
}

void benchmark_terrain(int resolution) {
    std::cout << "Benchmarking terrain generation at " << resolution << "x" << resolution << "..." << std::endl;
    
    TerrainConfig config;
    config.width = resolution;
    config.depth = resolution;
    config.max_height = 50.0f;
    config.octaves = 6;
    
    Timer timer;
    timer.start();
    
    TerrainGenerator terrain(config);
    terrain.generate();
    
    timer.stop();
    
    std::cout << "  Time: " << timer.elapsed_ms() << " ms" << std::endl;
    std::cout << "  Vertices: " << terrain.get_vertices().size() << std::endl;
}

int main() {
    std::cout << "╔══════════════════════════════════════════╗" << std::endl;
    std::cout << "║     Atlas Physics Benchmarks v1.0        ║" << std::endl;
    std::cout << "╚══════════════════════════════════════════╝" << std::endl << std::endl;
    
    // Rigid Body Benchmarks
    std::cout << "═══ Rigid Body Benchmarks ═══" << std::endl;
    benchmark_rigid_bodies(100);
    benchmark_rigid_bodies(500);
    benchmark_rigid_bodies(1000);
    std::cout << std::endl;
    
    // Soft Body Benchmarks
    std::cout << "═══ Soft Body Benchmarks ═══" << std::endl;
    benchmark_soft_body(100);
    benchmark_soft_body(500);
    std::cout << std::endl;
    
    // SPH Benchmarks
    std::cout << "═══ SPH Fluid Benchmarks ═══" << std::endl;
    benchmark_sph(1000);
    benchmark_sph(5000);
    std::cout << std::endl;
    
    // Raycast Benchmarks
    std::cout << "═══ Raycasting Benchmarks ═══" << std::endl;
    benchmark_raycast(10000, 100);
    benchmark_raycast(100000, 100);
    std::cout << std::endl;
    
    // Terrain Benchmarks
    std::cout << "═══ Terrain Generation Benchmarks ═══" << std::endl;
    benchmark_terrain(128);
    benchmark_terrain(256);
    benchmark_terrain(512);
    std::cout << std::endl;
    
    std::cout << "All benchmarks complete!" << std::endl;
    
    return 0;
}
