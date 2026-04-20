/**
 * @file test_soft_body.cpp
 * @brief Unit tests for soft body physics
 */

#include "atlas.hpp"
#include <iostream>
#include <cassert>

using namespace atlas;
using namespace atlas::physics;

void test_soft_body_creation() {
    std::cout << "Testing soft body creation..." << std::endl;
    
    // Create ball
    auto ball = SoftBody::create_ball({0, 0, 0}, 1.0f, 8);
    assert(ball != nullptr);
    assert(ball->get_particle_count() > 0);
    assert(ball->get_spring_count() > 0);
    
    // Create box
    auto box = SoftBody::create_box({0, 0, 0}, {1, 1, 1});
    assert(box != nullptr);
    assert(box->get_particle_count() == 8);
    assert(box->get_spring_count() > 0);
    
    // Create cloth
    auto cloth = SoftBody::create_cloth({0, 0, 0}, 2, 2, 5, 5);
    assert(cloth != nullptr);
    assert(cloth->get_particle_count() == 36);
    
    std::cout << "  Soft body creation tests passed!" << std::endl;
}

void test_soft_body_simulation() {
    std::cout << "Testing soft body simulation..." << std::endl;
    
    auto ball = SoftBody::create_ball({0, 5, 0}, 0.5f, 6);
    
    Vec3 initial_pos = ball->get_center();
    
    // Simulate
    for (int i = 0; i < 60; ++i) {
        ball->update(1.0f / 60.0f, {0, -9.81f, 0});
    }
    
    Vec3 final_pos = ball->get_center();
    
    // Ball should have moved down
    assert(final_pos.y < initial_pos.y);
    
    std::cout << "  Soft body simulation tests passed!" << std::endl;
}

void test_soft_body_collisions() {
    std::cout << "Testing soft body collisions..." << std::endl;
    
    auto ball = SoftBody::create_ball({0, 3, 0}, 0.5f, 8);
    
    // Simulate with ground collision
    for (int i = 0; i < 120; ++i) {
        ball->update(1.0f / 60.0f, {0, -9.81f, 0});
        ball->collide_with_plane({0, 1, 0}, 0);
    }
    
    // Ball should have settled near ground
    assert(ball->get_center().y < 1.0f);
    
    std::cout << "  Soft body collision tests passed!" << std::endl;
}

int main() {
    std::cout << "═══ Atlas Soft Body Tests ═══" << std::endl << std::endl;
    
    test_soft_body_creation();
    test_soft_body_simulation();
    test_soft_body_collisions();
    
    std::cout << std::endl << "All soft body tests passed!" << std::endl;
    return 0;
}
