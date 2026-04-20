/**
 * @file quaternion.hpp
 * @brief Quaternion Mathematics for 3D Rotation
 * 
 * Comprehensive quaternion implementation for efficient 3D rotation:
 * - All standard quaternion operations
 * - Conversion to/from matrices and Euler angles
 * - Spherical linear interpolation (SLERP)
 * - Angular velocity integration
 */

#pragma once

#include "vec3.hpp"
#include <cmath>
#include <ostream>

namespace atlas {
namespace math {

/**
 * @brief Quaternion Class for 3D Rotation
 * 
 * Quaternions are represented as (w, x, y, z) where:
 * - w is the scalar (cos(angle/2))
 * - (x, y, z) is the vector (axis * sin(angle/2))
 */
class Quaternion {
public:
    float w, x, y, z;

    // Constructors
    Quaternion() : w(1.0f), x(0.0f), y(0.0f), z(0.0f) {}

    Quaternion(float w, float x, float y, float z) : w(w), x(x), y(y), z(z) {}

    explicit Quaternion(float scalar) : w(scalar), x(scalar), y(scalar), z(scalar) {}

    explicit Quaternion(const Vec3& axis, float angle_rad) {
        float half_angle = angle_rad * 0.5f;
        float sin_half = std::sin(half_angle);
        float cos_half = std::cos(half_angle);
        
        Vec3 axis_normalized = axis.normalized();
        w = cos_half;
        x = axis_normalized.x * sin_half;
        y = axis_normalized.y * sin_half;
        z = axis_normalized.z * sin_half;
    }

    // Create from Euler angles (XYZ convention)
    static Quaternion from_euler(float pitch, float yaw, float roll) {
        float cy = std::cos(yaw * 0.5f);
        float sy = std::sin(yaw * 0.5f);
        float cp = std::cos(pitch * 0.5f);
        float sp = std::sin(pitch * 0.5f);
        float cr = std::cos(roll * 0.5f);
        float sr = std::sin(roll * 0.5f);

        Quaternion q;
        q.w = cr * cp * cy + sr * sp * sy;
        q.x = sr * cp * cy - cr * sp * sy;
        q.y = cr * sp * cy + sr * cp * sy;
        q.z = cr * cp * sy - sr * sp * cy;
        return q;
    }

    static Quaternion from_euler(const Vec3& euler) {
        return from_euler(euler.x, euler.y, euler.z);
    }

    // Create from axis-angle
    static Quaternion from_axis_angle(const Vec3& axis, float angle_rad) {
        return Quaternion(axis, angle_rad);
    }

    // Create from rotation matrix
    static Quaternion from_matrix(const Mat4& m) {
        float trace = m(0, 0) + m(1, 1) + m(2, 2);
        Quaternion q;

        if (trace > 0.0f) {
            float s = std::sqrt(trace + 1.0f) * 2.0f; // 4 * q.w
            q.w = 0.25f * s;
            q.x = (m(2, 1) - m(1, 2)) / s;
            q.y = (m(0, 2) - m(2, 0)) / s;
            q.z = (m(1, 0) - m(0, 1)) / s;
        } else if (m(0, 0) > m(1, 1) && m(0, 0) > m(2, 2)) {
            float s = std::sqrt(1.0f + m(0, 0) - m(1, 1) - m(2, 2)) * 2.0f;
            q.w = (m(2, 1) - m(1, 2)) / s;
            q.x = 0.25f * s;
            q.y = (m(0, 1) + m(1, 0)) / s;
            q.z = (m(0, 2) + m(2, 0)) / s;
        } else if (m(1, 1) > m(2, 2)) {
            float s = std::sqrt(1.0f + m(1, 1) - m(0, 0) - m(2, 2)) * 2.0f;
            q.w = (m(0, 2) - m(2, 0)) / s;
            q.x = (m(0, 1) + m(1, 0)) / s;
            q.y = 0.25f * s;
            q.z = (m(1, 2) + m(2, 1)) / s;
        } else {
            float s = std::sqrt(1.0f + m(2, 2) - m(0, 0) - m(1, 1)) * 2.0f;
            q.w = (m(1, 0) - m(0, 1)) / s;
            q.x = (m(0, 2) + m(2, 0)) / s;
            q.y = (m(1, 2) + m(2, 1)) / s;
            q.z = 0.25f * s;
        }

        return q.normalized();
    }

    // Static factories
    static Quaternion identity() { return Quaternion(); }
    static Quaternion zero() { return Quaternion(0.0f, 0.0f, 0.0f, 0.0f); }

    // Access
    float& operator[](int index) { return (&w)[index]; }
    const float& operator[](int index) const { return (&w)[index]; }

    Vec3 get_axis() const {
        float sin_sq = 1.0f - w * w;
        if (sin_sq < 1e-6f) return Vec3(1.0f, 0.0f, 0.0f);
        float inv_sin = 1.0f / std::sqrt(sin_sq);
        return Vec3(x * inv_sin, y * inv_sin, z * inv_sin);
    }

    float get_angle() const {
        return 2.0f * std::acos(std::clamp(w, -1.0f, 1.0f));
    }

    // Basic operations
    Quaternion operator*(const Quaternion& q) const {
        return Quaternion(
            w * q.w - x * q.x - y * q.y - z * q.z,
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x,
            w * q.z + x * q.y - y * q.x + z * q.w
        );
    }

    Quaternion operator*(float scalar) const {
        return Quaternion(w * scalar, x * scalar, y * scalar, z * scalar);
    }

    Quaternion operator+(const Quaternion& q) const {
        return Quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
    }

    Quaternion operator-(const Quaternion& q) const {
        return Quaternion(w - q.w, x - q.x, y - q.y, z - q.z);
    }

    Quaternion operator-() const {
        return Quaternion(-w, -x, -y, -z);
    }

    Quaternion& operator*=(const Quaternion& q) {
        *this = *this * q;
        return *this;
    }

    Quaternion& operator*=(float scalar) {
        w *= scalar; x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }

    // Comparison
    bool operator==(const Quaternion& q) const {
        return w == q.w && x == q.x && y == q.y && z == q.z;
    }

    bool operator!=(const Quaternion& q) const {
        return !(*this == q);
    }

    // Quaternion operations
    float dot(const Quaternion& q) const {
        return w * q.w + x * q.x + y * q.y + z * q.z;
    }

    float magnitude() const {
        return std::sqrt(w * w + x * x + y * y + z * z);
    }

    float magnitude_squared() const {
        return w * w + x * x + y * y + z * z;
    }

    Quaternion normalized() const {
        float mag = magnitude();
        if (mag > 1e-6f) {
            float inv = 1.0f / mag;
            return Quaternion(w * inv, x * inv, y * inv, z * inv);
        }
        return identity();
    }

    Quaternion& normalize() {
        *this = normalized();
        return *this;
    }

    Quaternion conjugate() const {
        return Quaternion(w, -x, -y, -z);
    }

    Quaternion inverse() const {
        float mag_sq = magnitude_squared();
        if (mag_sq < 1e-6f) return identity();
        float inv_mag_sq = 1.0f / mag_sq;
        return conjugate() * inv_mag_sq;
    }

    // Vector rotation
    Vec3 rotate(const Vec3& v) const {
        // nVidia SDK implementation
        Vec3 qvec(x, y, z);
        Vec3 uv = qvec.cross(v) * 2.0f;
        Vec3 uuv = qvec.cross(uv);
        return v + uv * w + uuv;
    }

    Vec3 rotate_inverse(const Vec3& v) const {
        Quaternion q_inv = inverse();
        return q_inv.rotate(v);
    }

    // Rotation using the quaternion as a rotation matrix
    Vec3 transform(const Vec3& v) const {
        return rotate(v);
    }

    // Convert to rotation matrix
    Mat4 to_matrix() const {
        Mat4 m;
        
        float xx = x * x;
        float xy = x * y;
        float xz = x * z;
        float xw = x * w;
        float yy = y * y;
        float yz = y * z;
        float yw = y * w;
        float zz = z * z;
        float zw = z * w;

        m(0, 0) = 1.0f - 2.0f * (yy + zz);
        m(0, 1) = 2.0f * (xy - zw);
        m(0, 2) = 2.0f * (xz + yw);
        
        m(1, 0) = 2.0f * (xy + zw);
        m(1, 1) = 1.0f - 2.0f * (xx + zz);
        m(1, 2) = 2.0f * (yz - xw);
        
        m(2, 0) = 2.0f * (xz - yw);
        m(2, 1) = 2.0f * (yz + xw);
        m(2, 2) = 1.0f - 2.0f * (xx + yy);
        
        m(3, 3) = 1.0f;

        return m;
    }

    // Convert to Euler angles (XYZ convention)
    Vec3 to_euler() const {
        Vec3 euler;
        
        // Roll (x-axis rotation)
        float sinr_cosp = 2.0f * (w * x + y * z);
        float cosr_cosp = 1.0f - 2.0f * (x * x + y * y);
        euler.x = std::atan2(sinr_cosp, cosr_cosp);

        // Pitch (y-axis rotation)
        float sinp = 2.0f * (w * y - z * x);
        if (std::abs(sinp) >= 1.0f) {
            euler.y = std::copysign(M_PI / 2.0f, sinp); // use 90 degrees if out of range
        } else {
            euler.y = std::asin(sinp);
        }

        // Yaw (z-axis rotation)
        float siny_cosp = 2.0f * (w * z + x * y);
        float cosy_cosp = 1.0f - 2.0f * (y * y + z * z);
        euler.z = std::atan2(siny_cosp, cosy_cosp);

        return euler;
    }

    // Interpolation
    static Quaternion lerp(const Quaternion& a, const Quaternion& b, float t) {
        float cos_angle = a.dot(b);
        float sign = (cos_angle < 0.0f) ? -1.0f : 1.0f;
        return (a + (b * sign - a) * t).normalized();
    }

    Quaternion lerp(const Quaternion& target, float t) const {
        return lerp(*this, target, t);
    }

    // Spherical linear interpolation
    static Quaternion slerp(const Quaternion& a, const Quaternion& b, float t) {
        float cos_angle = a.dot(b);
        
        // Ensure we take the shortest path
        Quaternion b_adjusted = b;
        if (cos_angle < 0.0f) {
            b_adjusted = b * -1.0f;
            cos_angle = -cos_angle;
        }

        // If very close, use linear interpolation
        if (cos_angle > 0.9995f) {
            return lerp(a, b_adjusted, t);
        }

        // Slerp
        float sin_angle = std::sqrt(1.0f - cos_angle * cos_angle);
        float angle = std::acos(cos_angle);
        float sin_t = std::sin(angle * t);
        float cos_t = std::cos(angle * t);

        float s0 = std::cos(angle * t) - cos_angle * sin_t / sin_angle;
        float s1 = sin_t / sin_angle;

        return a * s0 + b_adjusted * s1;
    }

    Quaternion slerp(const Quaternion& target, float t) const {
        return slerp(*this, target, t);
    }

    // Squad (smooth quaternion interpolation)
    static Quaternion squad(
        const Quaternion& q0, const Quaternion& q1,
        const Quaternion& q2, const Quaternion& q3,
        float t
    ) {
        Quaternion a = slerp(q0, q1, t);
        Quaternion b = slerp(q2, q3, t);
        return slerp(a, b, 2.0f * t * (1.0f - t));
    }

    // Normalize angle to [-PI, PI]
    static float normalize_angle(float angle) {
        while (angle > M_PI) angle -= 2.0f * M_PI;
        while (angle < -M_PI) angle += 2.0f * M_PI;
        return angle;
    }

    // Angular velocity to quaternion derivative
    static Quaternion angular_velocity_to_derivative(
        const Quaternion& q, const Vec3& angular_velocity
    ) {
        float half_w = -0.5f * (q.x * angular_velocity.x + 
                               q.y * angular_velocity.y + 
                               q.z * angular_velocity.z);
        float half_x = 0.5f * (q.w * angular_velocity.x + 
                               q.y * angular_velocity.z - 
                               q.z * angular_velocity.y);
        float half_y = 0.5f * (q.w * angular_velocity.y + 
                               q.z * angular_velocity.x - 
                               q.x * angular_velocity.z);
        float half_z = 0.5f * (q.w * angular_velocity.z + 
                               q.x * angular_velocity.y - 
                               q.y * angular_velocity.x);
        return Quaternion(half_w, half_x, half_y, half_z);
    }

    // Utility
    bool is_normalized(float epsilon = 1e-6f) const {
        return std::abs(magnitude_squared() - 1.0f) < epsilon;
    }

    bool is_identity(float epsilon = 1e-6f) const {
        return std::abs(w - 1.0f) < epsilon && 
               std::abs(x) < epsilon && 
               std::abs(y) < epsilon && 
               std::abs(z) < epsilon;
    }

    float* data() { return &w; }
    const float* data() const { return &w; }

    void print() const {
        printf("Quaternion(w=%.4f, x=%.4f, y=%.4f, z=%.4f)\n", w, x, y, z);
    }
};

// Stream output
inline std::ostream& operator<<(std::ostream& os, const Quaternion& q) {
    os << "Quaternion(" << q.w << ", " << q.x << ", " << q.y << ", " << q.z << ")";
    return os;
}

// Scalar * Quaternion
inline Quaternion operator*(float scalar, const Quaternion& q) {
    return q * scalar;
}

// Vec3 * Quaternion (rotate vector)
inline Vec3 operator*(const Vec3& v, const Quaternion& q) {
    return q.rotate(v);
}

} // namespace math
} // namespace atlas

// Hash function
namespace std {
    template<>
    struct hash<atlas::math::Quaternion> {
        size_t operator()(const atlas::math::Quaternion& q) const {
            size_t h1 = std::hash<float>{}(q.w);
            size_t h2 = std::hash<float>{}(q.x);
            size_t h3 = std::hash<float>{}(q.y);
            size_t h4 = std::hash<float>{}(q.z);
            return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
        }
    };
}
