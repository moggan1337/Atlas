/**
 * @file test_rigid_body.cpp
 * @brief Unit tests for rigid body physics
 */

#include "atlas.hpp"
#include <iostream>
#include <cassert>

using namespace atlas;
using namespace atlas::physics;

void test_rigid_body_creation() {
    std::cout << "Testing rigid body creation..." << std::endl;
    
    // Sphere
    auto sphere = RigidBody::create_sphere(1.0f, 5.0f);
    assert(sphere != nullptr);
    assert(std::abs(sphere->mass - 5.0f) < 1e-6f);
    assert(std::abs(sphere->inv_mass - 0.2f) < 1e-6f);
    
    // Box
    auto box = RigidBody::create_box({1, 1, 1}, 2.0f);
    assert(box != nullptr);
    assert(std::abs(box->mass - 2.0f) < 1e-6f);
    
    // Static plane
    auto plane = RigidBody::create_static_plane({0, 1, 0}, 0);
    assert(plane->mass == 0);
    assert(plane->inv_mass == 0);
    assert(plane->state == BodyState::Static);
    
    std::cout << "  Rigid body creation tests passed!" << std::endl;
}

void test_rigid_body_forces() {
    std::cout << "Testing rigid body forces..." << std::endl;
    
    auto sphere = RigidBody::create_sphere(0.5f, 1.0f);
    Vec3 initial_pos = sphere->get_position();
    
    // Apply force
    sphere->apply_force({10, 0, 0});
    
    // Integrate
    sphere->integrate(1.0f / 60.0f);
    
    // Check velocity changed
    assert(sphere->linear_velocity.x > 0);
    
    std::cout << "  Rigid body forces tests passed!" << std::endl;
}

void test_rigid_body_world() {
    std::cout << "Testing rigid body world..." << std::endl;
    
    WorldConfig config;
    config.gravity = {0, -9.81f, 0};
    
    RigidBodyWorld world(config);
    
    // Add ground
    auto ground = RigidBody::create_static_plane({0, 1, 0}, 0);
    world.add_body(ground);
    
    // Add sphere
    auto sphere = RigidBody::create_sphere(0.5f, 1.0f);
    sphere->set_position({0, 5, 0});
    world.add_body(sphere);
    
    // Simulate
    for (int i = 0; i < 60; ++i) {
        world.step(1.0f / 60.0f);
    }
    
    // Sphere should have fallen
    assert(sphere->get_position().y < 5.0f);
    
    // Sphere should have bounced (velocity should be changing)
    assert(sphere->linear_velocity.y != 0);
    
    std::cout << "  Rigid body world tests passed!" << std::endl;
}

void test_collision_shapes() {
    std::cout << "Testing collision shapes..." << std::endl;
    
    // Sphere shape
    SphereShape sphere(1.0f);
    assert(std::abs(sphere.get_volume() - 4.18879f) < 0.01f);
    
    // Box shape
    BoxShape box({1, 1, 1});
    assert(std::abs(box.get_volume() - 8.0f) < 0.01f);
    
    // AABB
    Transform t;
    t.position = {0, 0, 0};
    Vec3 aabb_min = box.get_aabb_min(t);
    Vec3 aabb_max = box.get_aabb_max(t);
    assert((aabb_min - Vec3(-1, -1, -1)).magnitude() < 0.01f);
    assert((aabb_max - Vec3(1, 1, 1)).magnitude() < 0.01f);
    
    std::cout << "  Collision shapes tests passed!" << std::endl;
}

int main() {
    std::cout << "═══ Atlas Rigid Body Tests ═══" << std::endl << std::endl;
    
    test_rigid_body_creation();
    test_rigid_body_forces();
    test_rigid_body_world();
    test_collision_shapes();
    
    std::cout << std::endl << "All rigid body tests passed!" << std::endl;
    return 0;
}
