/**
 * @file example_main.cpp
 * @brief Atlas Physics Engine - Example Usage
 */

#include "atlas.hpp"
#include <iostream>

using namespace atlas;

void example_rigid_bodies() {
    std::cout << "=== Rigid Body Example ===" << std::endl;
    
    physics::WorldConfig config;
    config.gravity = {0, -9.81f, 0};
    config.substeps = 4;
    
    physics::RigidBodyWorld world(config);
    
    // Create ground plane
    auto ground = physics::RigidBody::create_static_plane({0, 1, 0}, 0);
    ground->name = "Ground";
    world.add_body(ground);
    
    // Create some spheres
    for (int i = 0; i < 5; ++i) {
        auto sphere = physics::RigidBody::create_sphere(0.5f, 1.0f);
        sphere->set_position({i * 1.2f, 5 + i, 0});
        sphere->name = "Sphere_" + std::to_string(i);
        world.add_body(sphere);
    }
    
    // Create a box
    auto box = physics::RigidBody::create_box({1, 1, 1}, 2.0f);
    box->set_position({0, 8, 0});
    box->name = "Box";
    world.add_body(box);
    
    // Simulation loop
    std::cout << "Simulating 3 seconds..." << std::endl;
    for (int i = 0; i < 180; ++i) {
        world.step(1.0f / 60.0f);
        
        if (i % 60 == 0) {
            std::cout << "Time: " << (i / 60) << "s" << std::endl;
            for (auto& body : world.get_bodies()) {
                if (body->name != "Ground") {
                    auto pos = body->get_position();
                    std::cout << "  " << body->name << ": (" 
                              << pos.x << ", " << pos.y << ", " << pos.z << ")" << std::endl;
                }
            }
        }
    }
    std::cout << std::endl;
}

void example_soft_body() {
    std::cout << "=== Soft Body Example ===" << std::endl;
    
    physics::SoftBodyWorld world;
    world.config.gravity = {0, -9.81f, 0};
    
    // Create a jelly ball
    auto jelly = physics::SoftBody::create_ball({0, 5, 0}, 1.0f, 12);
    jelly->name = "Jelly";
    world.add_body(jelly);
    
    std::cout << "Soft body created with " << jelly->get_particle_count() 
              << " particles and " << jelly->get_spring_count() << " springs" << std::endl;
    
    // Simulate
    for (int i = 0; i < 120; ++i) {
        world.step(1.0f / 60.0f);
    }
    
    std::cout << "Soft body simulation complete" << std::endl << std::endl;
}

void example_cloth() {
    std::cout << "=== Cloth Example ===" << std::endl;
    
    physics::ClothWorld world;
    world.config.gravity = {0, -9.81f, 0};
    
    // Create cloth
    auto cloth = physics::Cloth::create_grid({0, 5, 0}, 3.0f, 3.0f, 10, 10, true);
    cloth->name = "Curtain";
    cloth->config.wind_strength = 2.0f;
    cloth->config.wind_direction = {1, 0, 0.5f};
    world.add_cloth(cloth);
    
    std::cout << "Cloth created with " << cloth->get_particle_count() 
              << " particles" << std::endl;
    
    // Simulate
    for (int i = 0; i < 180; ++i) {
        world.step(1.0f / 60.0f);
        
        if (i % 60 == 0) {
            auto center = cloth->get_center();
            std::cout << "Time: " << (i / 60) << "s, Cloth center: (" 
                      << center.x << ", " << center.y << ", " << center.z << ")" << std::endl;
        }
    }
    std::cout << std::endl;
}

void example_fluid() {
    std::cout << "=== SPH Fluid Example ===" << std::endl;
    
    physics::SPHWorld::Config config;
    config.sph.particle_radius = 0.1f;
    config.sph.rest_density = 1000.0f;
    config.gravity = {0, -9.81f, 0};
    
    physics::SPHWorld fluid_world(config);
    
    // Add source
    physics::FluidSource source({0, 5, 0}, {0, -1, 0}, 50);
    fluid_world.add_source(source);
    
    // Emit initial particles
    fluid_world.get_solver().emit_particle_grid({-1, 4, -1}, {2, 0.5f, 2}, 0.2f, {0, 0, 0});
    
    std::cout << "Initial particles: " << fluid_world.get_particle_count() << std::endl;
    
    // Simulate
    for (int i = 0; i < 120; ++i) {
        fluid_world.update(1.0f / 60.0f);
        
        if (i % 30 == 0) {
            std::cout << "Time: " << (i / 30) << "0.5s, Particles: " 
                      << fluid_world.get_particle_count() << std::endl;
        }
    }
    std::cout << std::endl;
}

void example_terrain() {
    std::cout << "=== Terrain Generation Example ===" << std::endl;
    
    physics::TerrainConfig config;
    config.width = 128;
    config.depth = 128;
    config.max_height = 30.0f;
    config.octaves = 6;
    config.frequency = 0.02f;
    config.enable_thermal_erosion = true;
    
    physics::TerrainGenerator terrain(config);
    
    std::cout << "Terrain generated with " << terrain.get_vertices().size() 
              << " vertices" << std::endl;
    
    // Get height at a point
    float h = terrain.get_height_at(64, 64);
    std::cout << "Height at center: " << h << std::endl;
    
    // Get biome
    auto biome = terrain.get_biome_at(64, 64);
    std::cout << "Biome at center: " << static_cast<int>(biome) << std::endl;
    std::cout << std::endl;
}

void example_constraints() {
    std::cout << "=== Constraints Example ===" << std::endl;
    
    physics::WorldConfig rb_config;
    physics::RigidBodyWorld rb_world(rb_config);
    
    physics::ConstraintWorld constraint_world;
    
    // Create two bodies
    auto body_a = physics::RigidBody::create_box({0.5f, 0.5f, 0.5f}, 1.0f);
    body_a->set_position({0, 5, 0});
    rb_world.add_body(body_a);
    
    auto body_b = physics::RigidBody::create_box({0.5f, 0.5f, 0.5f}, 1.0f);
    body_b->set_position({2, 5, 0});
    rb_world.add_body(body_b);
    
    // Create distance constraint
    auto constraint = constraint_world.add_distance_constraint(
        body_a.get(), body_b.get(),
        {0, 0, 0}, {0, 0, 0},
        2.0f
    );
    
    std::cout << "Constraint added between bodies" << std::endl;
    std::cout << "Target distance: " << constraint->target_distance << std::endl;
    std::cout << std::endl;
}

void example_raycasting() {
    std::cout << "=== Raycasting Example ===" << std::endl;
    
    // Create BVH accelerator
    physics::BVHRaycast accelerator;
    
    // Add some bounding volumes
    accelerator.set_primitives({
        physics::BoundingVolume({0, 0, 0}, {1, 1, 1}),
        physics::BoundingVolume({3, 0, 0}, {0.5f, 0.5f, 0.5f}),
        physics::BoundingVolume({0, 3, 0}, {0.5f, 0.5f, 0.5f})
    });
    
    accelerator.build();
    
    // Cast a ray
    physics::Ray ray({-5, 0, 0}, {1, 0, 0}, 0, 20);
    physics::RayHit hit;
    
    if (accelerator.raycast(ray, hit)) {
        std::cout << "Hit at t=" << hit.t << std::endl;
        std::cout << "Point: (" << hit.point.x << ", " 
                  << hit.point.y << ", " << hit.point.z << ")" << std::endl;
    } else {
        std::cout << "No hit" << std::endl;
    }
    
    // Test ray-sphere
    physics::Ray sphere_ray({0, 5, 0}, {0, -1, 0});
    physics::RayHit sphere_hit;
    
    if (physics::Raycast::ray_sphere(sphere_ray, {0, 0, 0}, 1.0f, sphere_hit)) {
        std::cout << "Sphere hit at t=" << sphere_hit.t << std::endl;
    }
    std::cout << std::endl;
}

void example_combined_world() {
    std::cout << "=== Combined Physics World Example ===" << std::endl;
    
    PhysicsWorld world;
    
    // Add rigid bodies
    auto ground = physics::RigidBody::create_static_plane({0, 1, 0}, 0);
    world.rigid_world.add_body(ground);
    
    auto ball = physics::RigidBody::create_sphere(0.5f, 1.0f);
    ball->set_position({0, 10, 0});
    world.rigid_world.add_body(ball);
    
    // Add soft body
    auto jelly = physics::SoftBody::create_ball({0, 7, 0}, 0.5f, 8);
    world.soft_world.add_body(jelly);
    
    // Add cloth
    auto cloth = physics::Cloth::create_grid({3, 5, 0}, 2, 2, 8, 8, true);
    world.cloth_world.add_body(cloth);
    
    std::cout << "Combined world created:" << std::endl;
    std::cout << "  Rigid bodies: " << world.rigid_world.get_body_count() << std::endl;
    std::cout << "  Soft bodies: " << world.soft_world.get_body_count() << std::endl;
    std::cout << "  Cloth objects: " << world.cloth_world.get_cloth_count() << std::endl;
    
    // Simulate
    for (int i = 0; i < 60; ++i) {
        world.update(1.0f / 60.0f);
    }
    
    std::cout << "Simulation complete!" << std::endl;
}

int main() {
    std::cout << "╔══════════════════════════════════════════╗" << std::endl;
    std::cout << "║     Atlas Physics Simulator v1.0        ║" << std::endl;
    std::cout << "╚══════════════════════════════════════════╝" << std::endl << std::endl;
    
    example_rigid_bodies();
    example_soft_body();
    example_cloth();
    example_fluid();
    example_terrain();
    example_constraints();
    example_raycasting();
    example_combined_world();
    
    std::cout << "All examples completed successfully!" << std::endl;
    
    return 0;
}
