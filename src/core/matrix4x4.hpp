/**
 * @file matrix4x4.hpp
 * @brief 4x4 Matrix Mathematics Library
 * 
 * Comprehensive 4x4 matrix operations for 3D graphics and physics:
 * - Transformation matrices (translation, rotation, scale)
 * - Matrix multiplication and inversion
 * - Vector transformation (points and directions)
 * - Projection matrices
 * - Decomposition and extraction
 */

#pragma once

#include "vec3.hpp"
#include <cstring>
#include <ostream>
#include <initializer_list>

namespace atlas {
namespace math {

/**
 * @brief Column-major 4x4 Matrix Class
 * 
 * Uses column-major storage for OpenGL compatibility.
 * Matrix layout:
 *   m[0] m[4] m[8]  m[12]
 *   m[1] m[5] m[9]  m[13]
 *   m[2] m[6] m[10] m[14]
 *   m[3] m[7] m[11] m[15]
 */
class Mat4 {
public:
    float m[16];

    // Constructors
    Mat4() { set_identity(); }

    Mat4(float diagonal) { set_diagonal(diagonal); }

    Mat4(std::initializer_list<float> values) {
        int i = 0;
        for (float v : values) {
            if (i < 16) m[i++] = v;
        }
        for (; i < 16; ++i) m[i] = 0.0f;
    }

    Mat4(const float* values) {
        std::memcpy(m, values, 16 * sizeof(float));
    }

    // Element access
    float& operator()(int row, int col) {
        return m[col * 4 + row];
    }

    const float& operator()(int row, int col) const {
        return m[col * 4 + row];
    }

    float& at(int row, int col) { return (*this)(row, col); }
    const float& at(int row, int col) const { return (*this)(row, col); }

    // Column access
    Vec3 get_column(int col) const {
        return Vec3(m[col * 4], m[col * 4 + 1], m[col * 4 + 2]);
    }

    void set_column(int col, const Vec3& v) {
        m[col * 4] = v.x;
        m[col * 4 + 1] = v.y;
        m[col * 4 + 2] = v.z;
    }

    // Row access
    Vec3 get_row(int row) const {
        return Vec3(m[row + 4], m[row + 8], m[row + 12]);
    }

    void set_row(int row, const Vec3& v) {
        m[row + 4] = v.x;
        m[row + 8] = v.y;
        m[row + 12] = v.z;
    }

    // Array access
    float* data() { return m; }
    const float* data() const { return m; }

    // Basic setters
    void set_identity() {
        std::memset(m, 0, 16 * sizeof(float));
        m[0] = m[5] = m[10] = m[15] = 1.0f;
    }

    void set_diagonal(float value) {
        std::memset(m, 0, 16 * sizeof(float));
        m[0] = m[5] = m[10] = m[15] = value;
    }

    void set_diagonal(const Vec3& v) {
        std::memset(m, 0, 16 * sizeof(float));
        m[0] = v.x;
        m[5] = v.y;
        m[10] = v.z;
        m[15] = 1.0f;
    }

    void set_zero() {
        std::memset(m, 0, 16 * sizeof(float));
    }

    // Static factory methods - Identity and Zero
    static Mat4 identity() {
        Mat4 mat;
        mat.set_identity();
        return mat;
    }

    static Mat4 zero() {
        Mat4 mat;
        mat.set_zero();
        return mat;
    }

    // Translation matrix
    static Mat4 translation(const Vec3& t) {
        Mat4 mat;
        mat.set_identity();
        mat(0, 3) = t.x;
        mat(1, 3) = t.y;
        mat(2, 3) = t.z;
        return mat;
    }

    static Mat4 translation(float x, float y, float z) {
        return translation(Vec3(x, y, z));
    }

    // Scale matrix
    static Mat4 scale(const Vec3& s) {
        Mat4 mat;
        mat.set_zero();
        mat(0, 0) = s.x;
        mat(1, 1) = s.y;
        mat(2, 2) = s.z;
        mat(3, 3) = 1.0f;
        return mat;
    }

    static Mat4 scale(float uniform) {
        return scale(Vec3(uniform, uniform, uniform));
    }

    static Mat4 scale(float x, float y, float z) {
        return scale(Vec3(x, y, z));
    }

    // Rotation matrices
    static Mat4 rotation_x(float angle_rad) {
        Mat4 mat;
        mat.set_identity();
        float c = std::cos(angle_rad);
        float s = std::sin(angle_rad);
        mat(1, 1) = c;  mat(1, 2) = -s;
        mat(2, 1) = s;  mat(2, 2) = c;
        return mat;
    }

    static Mat4 rotation_y(float angle_rad) {
        Mat4 mat;
        mat.set_identity();
        float c = std::cos(angle_rad);
        float s = std::sin(angle_rad);
        mat(0, 0) = c;  mat(0, 2) = s;
        mat(2, 0) = -s; mat(2, 2) = c;
        return mat;
    }

    static Mat4 rotation_z(float angle_rad) {
        Mat4 mat;
        mat.set_identity();
        float c = std::cos(angle_rad);
        float s = std::sin(angle_rad);
        mat(0, 0) = c;  mat(0, 1) = -s;
        mat(1, 0) = s;  mat(1, 1) = c;
        return mat;
    }

    static Mat4 rotation(float angle_rad, const Vec3& axis) {
        Mat4 mat;
        float c = std::cos(angle_rad);
        float s = std::sin(angle_rad);
        float t = 1.0f - c;
        
        Vec3 axis_normalized = axis.normalized();
        float x = axis_normalized.x;
        float y = axis_normalized.y;
        float z = axis_normalized.z;

        mat(0, 0) = t*x*x + c;     mat(0, 1) = t*x*y - s*z; mat(0, 2) = t*x*z + s*y;
        mat(1, 0) = t*x*y + s*z;   mat(1, 1) = t*y*y + c;   mat(1, 2) = t*y*z - s*x;
        mat(2, 0) = t*x*z - s*y;   mat(2, 1) = t*y*z + s*x; mat(2, 2) = t*z*z + c;
        mat(3, 3) = 1.0f;

        return mat;
    }

    // Axis-aligned rotation (no normalization)
    static Mat4 rotation(float angle_rad, float ax, float ay, float az) {
        return rotation(angle_rad, Vec3(ax, ay, az));
    }

    // Rotation from direction vectors
    static Mat4 rotation(const Vec3& forward, const Vec3& up) {
        Vec3 f = forward.normalized();
        Vec3 r = up.normalized().cross(f).normalized();
        Vec3 u = f.cross(r);

        Mat4 mat;
        mat(0, 0) = r.x; mat(0, 1) = r.y; mat(0, 2) = r.z;
        mat(1, 0) = u.x; mat(1, 1) = u.y; mat(1, 2) = u.z;
        mat(2, 0) = f.x; mat(2, 1) = f.y; mat(2, 2) = f.z;
        mat(3, 3) = 1.0f;
        return mat;
    }

    // Transform matrix from components
    static Mat4 transform(const Vec3& position, const Vec3& rotation_euler, const Vec3& scale) {
        Mat4 result = translation(position);
        result = result * rotation_x(rotation_euler.x);
        result = result * rotation_y(rotation_euler.y);
        result = result * rotation_z(rotation_euler.z);
        result = result * scale(scale);
        return result;
    }

    // Perspective projection matrices
    static Mat4 perspective(float fov_rad, float aspect, float near, float far) {
        Mat4 mat;
        float tan_half_fov = std::tan(fov_rad * 0.5f);
        float range = far - near;

        mat.set_zero();
        mat(0, 0) = 1.0f / (aspect * tan_half_fov);
        mat(1, 1) = 1.0f / tan_half_fov;
        mat(2, 2) = -(far + near) / range;
        mat(2, 3) = -2.0f * far * near / range;
        mat(3, 2) = -1.0f;

        return mat;
    }

    static Mat4 perspective_gl(float left, float right, float bottom, float top, float near, float far) {
        Mat4 mat;
        mat.set_zero();

        mat(0, 0) = 2.0f * near / (right - left);
        mat(1, 1) = 2.0f * near / (top - bottom);
        mat(2, 0) = (right + left) / (right - left);
        mat(2, 1) = (top + bottom) / (top - bottom);
        mat(2, 2) = -(far + near) / (far - near);
        mat(2, 3) = -1.0f;
        mat(3, 2) = -2.0f * far * near / (far - near);

        return mat;
    }

    static Mat4 perspective_vulkan(float fov_rad, float aspect, float near, float far) {
        Mat4 mat;
        float tan_half_fov = std::tan(fov_rad * 0.5f);

        mat.set_zero();
        mat(0, 0) = 1.0f / (aspect * tan_half_fov);
        mat(1, 1) = 1.0f / tan_half_fov;
        mat(2, 2) = far / (near - far);
        mat(2, 3) = -1.0f;
        mat(3, 2) = far * near / (near - far);

        return mat;
    }

    // Orthographic projection
    static Mat4 orthographic(float left, float right, float bottom, float top, float near, float far) {
        Mat4 mat;
        mat.set_identity();

        mat(0, 0) = 2.0f / (right - left);
        mat(1, 1) = 2.0f / (top - bottom);
        mat(2, 2) = -2.0f / (far - near);
        mat(3, 0) = -(right + left) / (right - left);
        mat(3, 1) = -(top + bottom) / (top - bottom);
        mat(3, 2) = -(far + near) / (far - near);

        return mat;
    }

    static Mat4 orthographic(float fov_rad, float aspect, float near, float far) {
        float half_height = near * std::tan(fov_rad * 0.5f);
        float half_width = half_height * aspect;
        return orthographic(-half_width, half_width, -half_height, half_height, near, far);
    }

    // LookAt matrices
    static Mat4 look_at(const Vec3& eye, const Vec3& target, const Vec3& up) {
        Vec3 forward = (target - eye).normalized();
        Vec3 right = forward.cross(up).normalized();
        Vec3 new_up = right.cross(forward);

        Mat4 mat;
        mat.set_identity();
        mat(0, 0) = right.x;  mat(0, 1) = new_up.x;  mat(0, 2) = -forward.x;
        mat(1, 0) = right.y;  mat(1, 1) = new_up.y;  mat(1, 2) = -forward.y;
        mat(2, 0) = right.z;  mat(2, 1) = new_up.z;  mat(2, 2) = -forward.z;
        mat(3, 0) = -right.dot(eye);
        mat(3, 1) = -new_up.dot(eye);
        mat(3, 2) = forward.dot(eye);

        return mat;
    }

    // Matrix multiplication
    Mat4 operator*(const Mat4& other) const {
        Mat4 result;
        for (int col = 0; col < 4; ++col) {
            for (int row = 0; row < 4; ++row) {
                float sum = 0.0f;
                for (int k = 0; k < 4; ++k) {
                    sum += (*this)(row, k) * other(k, col);
                }
                result(row, col) = sum;
            }
        }
        return result;
    }

    Mat4& operator*=(const Mat4& other) {
        *this = *this * other;
        return *this;
    }

    // Vector transformations
    Vec3 transform_point(const Vec3& p) const {
        return Vec3(
            m[0] * p.x + m[4] * p.y + m[8] * p.z + m[12],
            m[1] * p.x + m[5] * p.y + m[9] * p.z + m[13],
            m[2] * p.x + m[6] * p.y + m[10] * p.z + m[14]
        );
    }

    Vec3 transform_vector(const Vec3& v) const {
        return Vec3(
            m[0] * v.x + m[4] * v.y + m[8] * v.z,
            m[1] * v.x + m[5] * v.y + m[9] * v.z,
            m[2] * v.x + m[6] * v.y + m[8] * v.z
        );
    }

    Vec3 transform_direction(const Vec3& d) const {
        return Vec3(
            m[0] * d.x + m[4] * d.y + m[8] * d.z,
            m[1] * d.x + m[5] * d.y + m[9] * d.z,
            m[2] * d.x + m[6] * d.y + m[10] * d.z
        ).normalized();
    }

    Vec3 transform_normal(const Vec3& n) const {
        // For normals, use inverse transpose
        Mat4 inv = inverse();
        return Vec3(
            inv(0, 0) * n.x + inv(1, 0) * n.y + inv(2, 0) * n.z,
            inv(0, 1) * n.x + inv(1, 1) * n.y + inv(2, 1) * n.z,
            inv(0, 2) * n.x + inv(1, 2) * n.y + inv(2, 2) * n.z
        ).normalized();
    }

    // Homogeneous transformation
    Vec3 transform(const Vec3& v, float w) const {
        return Vec3(
            m[0] * v.x + m[4] * v.y + m[8] * v.z + m[12] * w,
            m[1] * v.x + m[5] * v.y + m[9] * v.z + m[13] * w,
            m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14] * w
        );
    }

    Vec4 transform(const Vec4& v) const {
        return Vec4(
            m[0] * v.x + m[4] * v.y + m[8] * v.z + m[12] * v.w,
            m[1] * v.x + m[5] * v.y + m[9] * v.z + m[13] * v.w,
            m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14] * v.w,
            m[3] * v.x + m[7] * v.y + m[11] * v.z + m[15] * v.w
        );
    }

    // Transpose
    Mat4 transposed() const {
        Mat4 result;
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                result(row, col) = (*this)(col, row);
            }
        }
        return result;
    }

    Mat4& transpose() {
        *this = transposed();
        return *this;
    }

    // Determinant
    float determinant() const {
        float det = 0.0f;
        for (int col = 0; col < 4; ++col) {
            det += (col % 2 == 0 ? 1.0f : -1.0f) * (*this)(0, col) * minor(0, col);
        }
        return det;
    }

    float minor(int row, int col) const {
        float minors[9];
        int idx = 0;
        for (int r = 0; r < 4; ++r) {
            if (r == row) continue;
            for (int c = 0; c < 4; ++c) {
                if (c == col) continue;
                minors[idx++] = (*this)(r, c);
            }
        }
        return minors[0] * (minors[4] * minors[8] - minors[5] * minors[7])
             - minors[1] * (minors[3] * minors[8] - minors[5] * minors[6])
             + minors[2] * (minors[3] * minors[7] - minors[4] * minors[6]);
    }

    // Inverse
    Mat4 inverse() const {
        Mat4 inv;
        float inv_det = 1.0f / determinant();
        
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                inv(col, row) = ((row + col) % 2 == 0 ? 1.0f : -1.0f) 
                                * minor(row, col) * inv_det;
            }
        }
        return inv;
    }

    Mat4& invert() {
        *this = inverse();
        return *this;
    }

    // Decomposition
    void decompose(Vec3& translation, Vec3& rotation, Vec3& scale) const {
        translation = Vec3(m[12], m[13], m[14]);
        
        scale.x = Vec3(m[0], m[1], m[2]).magnitude();
        scale.y = Vec3(m[4], m[5], m[6]).magnitude();
        scale.z = Vec3(m[8], m[9], m[10]).magnitude();

        // Extract rotation (simplified - assumes orthogonal matrix)
        rotation.x = std::atan2(m[6], m[10]);
        rotation.y = std::atan2(-m[2], std::sqrt(m[0]*m[0] + m[1]*m[1]));
        rotation.z = 0.0f;
    }

    // Interpolation
    static Mat4 lerp(const Mat4& a, const Mat4& b, float t) {
        Mat4 result;
        for (int i = 0; i < 16; ++i) {
            result.m[i] = a.m[i] * (1.0f - t) + b.m[i] * t;
        }
        return result;
    }

    Mat4 lerp(const Mat4& target, float t) const {
        return lerp(*this, target, t);
    }

    // Comparison
    bool operator==(const Mat4& other) const {
        for (int i = 0; i < 16; ++i) {
            if (m[i] != other.m[i]) return false;
        }
        return true;
    }

    bool operator!=(const Mat4& other) const {
        return !(*this == other);
    }

    bool is_identity(float epsilon = 1e-6f) const {
        for (int i = 0; i < 16; ++i) {
            float expected = (i % 5 == 0) ? 1.0f : 0.0f;
            if (std::abs(m[i] - expected) > epsilon) return false;
        }
        return true;
    }

    bool is_symmetric(float epsilon = 1e-6f) const {
        for (int row = 0; row < 4; ++row) {
            for (int col = row + 1; col < 4; ++col) {
                if (std::abs((*this)(row, col) - (*this)(col, row)) > epsilon) {
                    return false;
                }
            }
        }
        return true;
    }

    // Utility
    void print() const {
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col) {
                printf("%8.3f ", (*this)(row, col));
            }
            printf("\n");
        }
    }
};

// Stream output
inline std::ostream& operator<<(std::ostream& os, const Mat4& mat) {
    for (int row = 0; row < 4; ++row) {
        os << "[";
        for (int col = 0; col < 4; ++col) {
            os << mat(row, col);
            if (col < 3) os << ", ";
        }
        os << "]\n";
    }
    return os;
}

// Additional Vec4 class for homogeneous coordinates
class Vec4 {
public:
    float x, y, z, w;

    Vec4() : x(0), y(0), z(0), w(0) {}
    Vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
    explicit Vec4(const Vec3& v, float w = 1.0f) : x(v.x), y(v.y), z(v.z), w(w) {}

    Vec4 operator*(float scalar) const {
        return Vec4(x * scalar, y * scalar, z * scalar, w * scalar);
    }

    Vec4 operator+(const Vec4& v) const {
        return Vec4(x + v.x, y + v.y, z + v.z, w + v.w);
    }

    Vec4 operator-(const Vec4& v) const {
        return Vec4(x - v.x, y - v.y, z - v.z, w - v.w);
    }

    Vec3 to_vec3() const {
        if (std::abs(w) < 1e-6f) return Vec3(x, y, z);
        return Vec3(x / w, y / w, z / w);
    }

    float dot(const Vec4& v) const {
        return x * v.x + y * v.y + z * v.z + w * v.w;
    }

    float magnitude() const {
        return std::sqrt(x*x + y*y + z*z + w*w);
    }

    Vec4 normalized() const {
        float mag = magnitude();
        if (mag > 1e-6f) {
            return *this * (1.0f / mag);
        }
        return Vec4(0, 0, 0, 0);
    }

    float* data() { return &x; }
    const float* data() const { return &x; }
};

} // namespace math
} // namespace atlas
