/**
 * @file test_collision.cpp
 * @brief Unit tests for collision detection
 */

#include "atlas.hpp"
#include <iostream>
#include <cassert>

using namespace atlas;
using namespace atlas::physics;

void test_ray_sphere() {
    std::cout << "Testing ray-sphere intersection..." << std::endl;
    
    Ray ray({0, 0, 0}, {1, 0, 0}, 0, 100);
    RayHit hit;
    
    // Miss
    bool result1 = Raycast::ray_sphere(ray, {5, 0, 0}, 1.0f, hit);
    assert(!result1);
    
    // Hit
    bool result2 = Raycast::ray_sphere(ray, {5, 0, 0}, 6.0f, hit);
    assert(result2);
    assert(hit.t > 0);
    
    // Origin inside
    Ray ray2({0, 0, 0}, {1, 0, 0});
    bool result3 = Raycast::ray_sphere(ray2, {0, 0, 0}, 1.0f, hit);
    assert(result3);
    
    std::cout << "  Ray-sphere tests passed!" << std::endl;
}

void test_ray_box() {
    std::cout << "Testing ray-box intersection..." << std::endl;
    
    Ray ray({-5, 0, 0}, {1, 0, 0});
    RayHit hit;
    
    // Hit
    bool result = Raycast::ray_box(ray, {-1, -1, -1}, {1, 1, 1}, hit);
    assert(result);
    assert(hit.t > 0);
    
    // Miss
    Ray ray2({0, 0, 5}, {1, 0, 0});
    bool result2 = Raycast::ray_box(ray2, {-1, -1, -1}, {1, 1, 1}, hit);
    assert(!result2);
    
    std::cout << "  Ray-box tests passed!" << std::endl;
}

void test_ray_plane() {
    std::cout << "Testing ray-plane intersection..." << std::endl;
    
    Ray ray({0, 5, 0}, {0, -1, 0});
    RayHit hit;
    
    bool result = Raycast::ray_plane(ray, {0, 0, 0}, {0, 1, 0}, hit);
    assert(result);
    assert(std::abs(hit.t - 5.0f) < 0.01f);
    
    // Parallel ray
    Ray ray2({0, 5, 0}, {1, 0, 0});
    bool result2 = Raycast::ray_plane(ray2, {0, 0, 0}, {0, 1, 0}, hit);
    assert(!result2);
    
    std::cout << "  Ray-plane tests passed!" << std::endl;
}

void test_ray_triangle() {
    std::cout << "Testing ray-triangle intersection..." << std::endl;
    
    Vec3 v0(0, 0, 0), v1(1, 0, 0), v2(0, 1, 0);
    
    Ray ray({0.2f, 0.2f, 1}, {0, 0, -1});
    RayHit hit;
    
    bool result = Raycast::ray_triangle(ray, v0, v1, v2, hit);
    assert(result);
    assert(hit.t > 0);
    
    // Outside triangle
    Ray ray2({0.6f, 0.6f, 1}, {0, 0, -1});
    bool result2 = Raycast::ray_triangle(ray2, v0, v1, v2, hit);
    assert(!result2);
    
    std::cout << "  Ray-triangle tests passed!" << std::endl;
}

void test_bvh() {
    std::cout << "Testing BVH..." << std::endl;
    
    BVHRaycast accelerator;
    accelerator.set_primitives({
        BoundingVolume({0, 0, 0}, {1, 1, 1}),
        BoundingVolume({5, 0, 0}, {0.5f, 0.5f, 0.5f}),
        BoundingVolume({0, 5, 0}, {0.5f, 0.5f, 0.5f})
    });
    accelerator.build();
    
    Ray ray({-5, 0, 0}, {1, 0, 0}, 0, 20);
    RayHit hit;
    
    bool result = accelerator.raycast(ray, hit);
    assert(result);
    assert(hit.t < 10);
    
    std::cout << "  BVH tests passed!" << std::endl;
}

int main() {
    std::cout << "═══ Atlas Collision Detection Tests ═══" << std::endl << std::endl;
    
    test_ray_sphere();
    test_ray_box();
    test_ray_plane();
    test_ray_triangle();
    test_bvh();
    
    std::cout << std::endl << "All collision tests passed!" << std::endl;
    return 0;
}
