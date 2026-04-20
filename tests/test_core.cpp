/**
 * @file test_core.cpp
 * @brief Unit tests for core math library
 */

#include "atlas.hpp"
#include <cmath>
#include <cassert>
#include <iostream>

using namespace atlas;

void test_vec3() {
    std::cout << "Testing Vec3..." << std::endl;
    
    // Construction
    math::Vec3 v1(1, 2, 3);
    assert(v1.x == 1 && v1.y == 2 && v1.z == 3);
    
    // Addition
    math::Vec3 v2 = v1 + math::Vec3(1, 1, 1);
    assert(v2.x == 2 && v2.y == 3 && v2.z == 4);
    
    // Subtraction
    math::Vec3 v3 = v1 - math::Vec3(1, 0, 0);
    assert(v3.x == 0 && v3.y == 2 && v3.z == 3);
    
    // Multiplication
    math::Vec3 v4 = v1 * 2.0f;
    assert(v4.x == 2 && v4.y == 4 && v4.z == 6);
    
    // Dot product
    float dot = v1.dot(math::Vec3(1, 0, 0));
    assert(std::abs(dot - 1.0f) < 1e-6f);
    
    // Cross product
    math::Vec3 cross = math::Vec3(1, 0, 0).cross(math::Vec3(0, 1, 0));
    assert(std::abs(cross.x - 0) < 1e-6f);
    assert(std::abs(cross.y - 0) < 1e-6f);
    assert(std::abs(cross.z - 1) < 1e-6f);
    
    // Magnitude
    float mag = math::Vec3(3, 4, 0).magnitude();
    assert(std::abs(mag - 5.0f) < 1e-6f);
    
    // Normalize
    math::Vec3 v5 = math::Vec3(3, 4, 0).normalized();
    assert(std::abs(v5.x - 0.6f) < 1e-5f);
    assert(std::abs(v5.y - 0.8f) < 1e-5f);
    
    // Lerp
    math::Vec3 v6 = math::Vec3::lerp(math::Vec3(0), math::Vec3(10), 0.5f);
    assert(std::abs(v6.x - 5) < 1e-6f);
    
    std::cout << "  Vec3 tests passed!" << std::endl;
}

void test_quaternion() {
    std::cout << "Testing Quaternion..." << std::endl;
    
    // Identity
    math::Quaternion q1 = math::Quaternion::identity();
    assert(q1.w == 1 && q1.x == 0 && q1.y == 0 && q1.z == 0);
    
    // From axis-angle
    math::Quaternion q2 = math::Quaternion::from_axis_angle(math::Vec3::unit_y(), M_PI_2);
    assert(std::abs(q2.w - std::cos(M_PI_4)) < 1e-5f);
    
    // To matrix and back
    math::Quaternion q3 = math::Quaternion::from_euler(0.5f, 0.3f, 0.1f);
    math::Mat4 m = q3.to_matrix();
    math::Quaternion q4 = math::Quaternion::from_matrix(m);
    
    // Normalize
    q4.normalize();
    assert(std::abs(q3.dot(q4) - 1.0f) < 0.01f);
    
    // Rotate vector
    math::Vec3 v1(1, 0, 0);
    math::Quaternion q5 = math::Quaternion::from_axis_angle(math::Vec3::unit_y(), M_PI_2);
    math::Vec3 v2 = q5.rotate(v1);
    assert(std::abs(v2.x) < 1e-5f);
    assert(std::abs(v2.z - 1) < 1e-5f);
    
    // SLERP
    math::Quaternion q6 = math::Quaternion::identity();
    math::Quaternion q7 = math::Quaternion::from_axis_angle(math::Vec3::unit_y(), M_PI);
    math::Quaternion q8 = math::Quaternion::slerp(q6, q7, 0.5f);
    assert(std::abs(q8.w - std::cos(M_PI_4)) < 1e-5f);
    
    std::cout << "  Quaternion tests passed!" << std::endl;
}

void test_matrix() {
    std::cout << "Testing Mat4..." << std::endl;
    
    // Identity
    math::Mat4 m1 = math::Mat4::identity();
    assert(m1(0, 0) == 1 && m1(1, 1) == 1 && m1(2, 2) == 1 && m1(3, 3) == 1);
    
    // Translation
    math::Mat4 m2 = math::Mat4::translation(1, 2, 3);
    assert(m2(0, 3) == 1 && m2(1, 3) == 2 && m2(2, 3) == 3);
    
    // Scale
    math::Mat4 m3 = math::Mat4::scale(2, 3, 4);
    assert(m3(0, 0) == 2 && m3(1, 1) == 3 && m3(2, 2) == 4);
    
    // Rotation
    math::Mat4 m4 = math::Mat4::rotation_x(M_PI_2);
    math::Vec3 v1(0, 1, 0);
    math::Vec3 v2 = m4.transform_direction(v1);
    assert(std::abs(v2.y) < 1e-5f);
    assert(std::abs(v2.z - 1) < 1e-5f);
    
    // Inverse
    math::Mat4 m5 = math::Mat4::translation(1, 2, 3);
    math::Mat4 m6 = m5.inverse();
    math::Mat4 product = m5 * m6;
    assert(product.is_identity());
    
    std::cout << "  Mat4 tests passed!" << std::endl;
}

void test_transform() {
    std::cout << "Testing Transform..." << std::endl;
    
    // Create transform
    scene::Transform t1;
    t1.position = {1, 2, 3};
    t1.rotation = math::Quaternion::from_euler(0, M_PI_2, 0);
    
    // To matrix
    math::Mat4 m = t1.to_matrix();
    
    // Transform point
    math::Vec3 p1(1, 0, 0);
    math::Vec3 p2 = t1.transform_point(p1);
    
    // Point should be rotated 90 degrees around Y
    assert(std::abs(p2.x) < 1e-3f);
    assert(std::abs(p2.z - 1) < 1e-3f);
    
    // Inverse transform
    math::Vec3 p3 = t1.inverse_transform_point(p2);
    assert((p3 - p1).magnitude() < 1e-5f);
    
    std::cout << "  Transform tests passed!" << std::endl;
}

int main() {
    std::cout << "═══ Atlas Core Math Tests ═══" << std::endl << std::endl;
    
    test_vec3();
    test_quaternion();
    test_matrix();
    test_transform();
    
    std::cout << std::endl << "All core tests passed!" << std::endl;
    return 0;
}
