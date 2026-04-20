/**
 * @file vec2.hpp
 * @brief 2D Vector Class
 */

#pragma once

#include <cmath>
#include <ostream>

namespace atlas {
namespace math {

class Vec2 {
public:
    float x, y;

    Vec2() : x(0.0f), y(0.0f) {}
    Vec2(float x, float y) : x(x), y(y) {}
    explicit Vec2(float scalar) : x(scalar), y(scalar) {}

    float& operator[](int index) { return (&x)[index]; }
    const float& operator[](int index) const { return (&x)[index]; }

    Vec2 operator+(const Vec2& v) const { return Vec2(x + v.x, y + v.y); }
    Vec2 operator-(const Vec2& v) const { return Vec2(x - v.x, y - v.y); }
    Vec2 operator*(float scalar) const { return Vec2(x * scalar, y * scalar); }
    Vec2 operator/(float scalar) const { float inv = 1.0f / scalar; return Vec2(x * inv, y * inv); }

    Vec2& operator+=(const Vec2& v) { x += v.x; y += v.y; return *this; }
    Vec2& operator-=(const Vec2& v) { x -= v.x; y -= v.y; return *this; }
    Vec2& operator*=(float scalar) { x *= scalar; y *= scalar; return *this; }
    Vec2& operator/=(float scalar) { float inv = 1.0f / scalar; x *= inv; y *= inv; return *this; }

    float dot(const Vec2& v) const { return x * v.x + y * v.y; }
    float cross(const Vec2& v) const { return x * v.y - y * v.x; }
    float magnitude() const { return std::sqrt(x * x + y * y); }
    float magnitude_squared() const { return x * x + y * y; }

    Vec2 normalized() const {
        float mag = magnitude();
        if (mag > 1e-6f) return *this * (1.0f / mag);
        return Vec2(0.0f, 0.0f);
    }

    Vec2& normalize() { *this = normalized(); return *this; }

    static Vec2 lerp(const Vec2& a, const Vec2& b, float t) {
        return a + (b - a) * t;
    }

    static Vec2 zero() { return Vec2(0.0f, 0.0f); }
    static Vec2 one() { return Vec2(1.0f, 1.0f); }
    static Vec2 unit_x() { return Vec2(1.0f, 0.0f); }
    static Vec2 unit_y() { return Vec2(0.0f, 1.0f); }

    float* data() { return &x; }
    const float* data() const { return &x; }
};

inline Vec2 operator*(float scalar, const Vec2& v) { return v * scalar; }

inline std::ostream& operator<<(std::ostream& os, const Vec2& v) {
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}

} // namespace math
} // namespace atlas
