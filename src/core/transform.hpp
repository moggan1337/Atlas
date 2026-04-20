/**
 * @file transform.hpp
 * @brief Transform Component for 3D Objects
 * 
 * Combines position, rotation (quaternion), and scale:
 * - World and local space transformations
 * - Parent-child hierarchy support
 * - Efficient transform operations
 */

#pragma once

#include "vec3.hpp"
#include "quaternion.hpp"
#include "matrix4x4.hpp"

namespace atlas {
namespace scene {

using math::Vec3;
using math::Quaternion;
using math::Mat4;

/**
 * @brief Transform Component
 * 
 * Represents the transformation of an object in 3D space.
 * Supports hierarchical transformations (parent-child).
 */
class Transform {
public:
    Vec3 position;
    Quaternion rotation;
    Vec3 scale;

    // Constructors
    Transform() 
        : position(Vec3::zero())
        , rotation(Quaternion::identity())
        , scale(Vec3::one()) {}

    Transform(const Vec3& position) 
        : position(position)
        , rotation(Quaternion::identity())
        , scale(Vec3::one()) {}

    Transform(const Vec3& position, const Quaternion& rotation) 
        : position(position)
        , rotation(rotation)
        , scale(Vec3::one()) {}

    Transform(const Vec3& position, const Quaternion& rotation, const Vec3& scale) 
        : position(position)
        , rotation(rotation)
        , scale(scale) {}

    Transform(const Vec3& position, const Vec3& euler) 
        : position(position)
        , rotation(Quaternion::from_euler(euler))
        , scale(Vec3::one()) {}

    Transform(const Vec3& position, const Vec3& euler, const Vec3& scale) 
        : position(position)
        , rotation(Quaternion::from_euler(euler))
        , scale(scale) {}

    // Static factories
    static Transform identity() { return Transform(); }

    static Transform from_position(const Vec3& position) {
        return Transform(position);
    }

    static Transform from_rotation(const Quaternion& rotation) {
        return Transform(Vec3::zero(), rotation);
    }

    static Transform from_position_rotation(const Vec3& position, const Quaternion& rotation) {
        return Transform(position, rotation);
    }

    static Transform from_position_euler(const Vec3& position, const Vec3& euler) {
        return Transform(position, euler);
    }

    static Transform look_at(const Vec3& eye, const Vec3& target, const Vec3& up = Vec3::up()) {
        Transform t;
        t.position = eye;
        Vec3 forward = (target - eye).normalized();
        t.rotation = Quaternion::from_axis_angle(Vec3::up(), 
            std::atan2(forward.x, forward.z));
        return t;
    }

    // Matrix conversion
    Mat4 to_matrix() const {
        Mat4 result = Mat4::translation(position);
        result = result * rotation.to_matrix();
        result = result * Mat4::scale(scale);
        return result;
    }

    Mat4 to_inverse_matrix() const {
        return inverse().to_matrix();
    }

    // Extract from matrix
    static Transform from_matrix(const Mat4& m) {
        Transform t;
        
        // Extract translation
        t.position = Vec3(m(0, 3), m(1, 3), m(2, 3));
        
        // Extract scale
        t.scale.x = Vec3(m(0, 0), m(1, 0), m(2, 0)).magnitude();
        t.scale.y = Vec3(m(0, 1), m(1, 1), m(2, 1)).magnitude();
        t.scale.z = Vec3(m(0, 2), m(1, 2), m(2, 2)).magnitude();
        
        // Extract rotation
        Mat4 rotation_matrix = m;
        if (t.scale.x > 1e-6f) {
            rotation_matrix(0, 0) /= t.scale.x;
            rotation_matrix(1, 0) /= t.scale.x;
            rotation_matrix(2, 0) /= t.scale.x;
        }
        if (t.scale.y > 1e-6f) {
            rotation_matrix(0, 1) /= t.scale.y;
            rotation_matrix(1, 1) /= t.scale.y;
            rotation_matrix(2, 1) /= t.scale.y;
        }
        if (t.scale.z > 1e-6f) {
            rotation_matrix(0, 2) /= t.scale.z;
            rotation_matrix(1, 2) /= t.scale.z;
            rotation_matrix(2, 2) /= t.scale.z;
        }
        
        t.rotation = Quaternion::from_matrix(rotation_matrix);
        
        return t;
    }

    // Decomposition
    Vec3 get_euler() const {
        return rotation.to_euler();
    }

    Vec3 get_forward() const {
        return rotation.rotate(Vec3::forward());
    }

    Vec3 get_backward() const {
        return rotation.rotate(Vec3::backward());
    }

    Vec3 get_right() const {
        return rotation.rotate(Vec3::right());
    }

    Vec3 get_left() const {
        return rotation.rotate(Vec3::left());
    }

    Vec3 get_up() const {
        return rotation.rotate(Vec3::up());
    }

    Vec3 get_down() const {
        return rotation.rotate(Vec3::down());
    }

    // Operations
    Transform operator*(const Transform& other) const {
        Mat4 combined = to_matrix() * other.to_matrix();
        return from_matrix(combined);
    }

    Transform inverse() const {
        return from_matrix(to_inverse_matrix());
    }

    // Transform points and directions
    Vec3 transform_point(const Vec3& point) const {
        return to_matrix().transform_point(point);
    }

    Vec3 transform_vector(const Vec3& vector) const {
        return rotation.rotate(vector) * scale;
    }

    Vec3 transform_direction(const Vec3& direction) const {
        return rotation.rotate(direction);
    }

    Vec3 inverse_transform_point(const Vec3& point) const {
        return inverse().transform_point(point);
    }

    Vec3 inverse_transform_vector(const Vec3& vector) const {
        return rotation.rotate_inverse(vector) / scale;
    }

    Vec3 inverse_transform_direction(const Vec3& direction) const {
        return rotation.rotate_inverse(direction);
    }

    // Local to world / World to local
    Vec3 world_to_local(const Vec3& world_point) const {
        return inverse_transform_point(world_point);
    }

    Vec3 local_to_world(const Vec3& local_point) const {
        return transform_point(local_point);
    }

    Vec3 world_to_local_vector(const Vec3& world_vector) const {
        return inverse_transform_vector(world_vector);
    }

    Vec3 local_to_world_vector(const Vec3& local_vector) const {
        return transform_vector(local_vector);
    }

    Vec3 world_to_local_direction(const Vec3& world_direction) const {
        return inverse_transform_direction(world_direction);
    }

    Vec3 local_to_world_direction(const Vec3& local_direction) const {
        return transform_direction(local_direction);
    }

    // Interpolation
    static Transform lerp(const Transform& a, const Transform& b, float t) {
        Transform result;
        result.position = Vec3::lerp(a.position, b.position, t);
        result.rotation = Quaternion::lerp(a.rotation, b.rotation, t);
        result.scale = Vec3::lerp(a.scale, b.scale, t);
        return result;
    }

    Transform lerp(const Transform& target, float t) const {
        return lerp(*this, target, t);
    }

    static Transform slerp(const Transform& a, const Transform& b, float t) {
        Transform result;
        result.position = Vec3::lerp(a.position, b.position, t);
        result.rotation = Quaternion::slerp(a.rotation, b.rotation, t);
        result.scale = Vec3::lerp(a.scale, b.scale, t);
        return result;
    }

    Transform slerp(const Transform& target, float t) const {
        return slerp(*this, target, t);
    }

    // Smoothing
    void smooth_follow(const Transform& target, float smooth_factor) {
        position = Vec3::lerp(position, target.position, smooth_factor);
        rotation = Quaternion::slerp(rotation, target.rotation, smooth_factor);
        scale = Vec3::lerp(scale, target.scale, smooth_factor);
    }

    // Utility
    bool is_identity(float epsilon = 1e-6f) const {
        return position.magnitude_squared() < epsilon &&
               rotation.is_identity(epsilon) &&
               (scale - Vec3::one()).magnitude_squared() < epsilon;
    }

    void set_identity() {
        position = Vec3::zero();
        rotation = Quaternion::identity();
        scale = Vec3::one();
    }

    void print() const {
        printf("Transform:\n");
        printf("  Position: %.3f, %.3f, %.3f\n", position.x, position.y, position.z);
        printf("  Euler:    %.3f, %.3f, %.3f\n", 
               get_euler().x, get_euler().y, get_euler().z);
        printf("  Scale:    %.3f, %.3f, %.3f\n", scale.x, scale.y, scale.z);
    }

    // Operators
    bool operator==(const Transform& other) const {
        return position == other.position &&
               rotation == other.rotation &&
               scale == other.scale;
    }

    bool operator!=(const Transform& other) const {
        return !(*this == other);
    }

    Transform operator+(const Vec3& translation) const {
        Transform result = *this;
        result.position += translation;
        return result;
    }

    Transform operator-(const Vec3& translation) const {
        Transform result = *this;
        result.position -= translation;
        return result;
    }
};

// Stream output
inline std::ostream& operator<<(std::ostream& os, const Transform& t) {
    os << "Transform(pos=" << t.position << ", euler=" << t.get_euler() 
       << ", scale=" << t.scale << ")";
    return os;
}

} // namespace scene
} // namespace atlas
