/**
 * @file raycasting.hpp
 * @brief Real-time Raycasting System
 * 
 * Implements:
 * - Ray-sphere, ray-box, ray-plane intersection
 * - Ray-mesh intersection (triangles)
 * - BVH-accelerated raycasting
 * - Spatial partitioning
 * - Scene queries (point, ray, volume)
 */

#pragma once

#include "../../core/vec3.hpp"
#include "../../core/transform.hpp"
#include <vector>
#include <memory>
#include <functional>

namespace atlas {
namespace physics {

using core::Vec3;
using core::Transform;

/**
 * @brief Ray
 */
struct Ray {
    Vec3 origin;
    Vec3 direction;
    float min_t;  // Minimum distance
    float max_t;  // Maximum distance
    
    Ray() : origin(Vec3::zero()), direction(Vec3::unit_z()), 
            min_t(0.0f), max_t(INFINITY) {}
    
    Ray(const Vec3& o, const Vec3& d, float min = 0.0f, float max = INFINITY)
        : origin(o), direction(d.normalized()), min_t(min), max_t(max) {}
    
    Vec3 get_point(float t) const { return origin + direction * t; }
};

/**
 * @brief Ray Hit Result
 */
struct RayHit {
    bool hit;
    float t;               // Distance along ray
    Vec3 point;            // Hit point in world space
    Vec3 normal;           // Surface normal at hit
    Vec3 barycentric;      // Barycentric coords for triangles
    int triangle_index;    // Triangle that was hit
    void* user_data;      // User-defined data
    float u, v;           // UV coordinates on triangle
    
    RayHit() : hit(false), t(INFINITY), point(Vec3::zero()), 
               normal(Vec3::unit_y()), barycentric(Vec3::zero()),
               triangle_index(-1), user_data(nullptr), u(0), v(0) {}
};

/**
 * @brief Bounding Volume
 */
struct BoundingVolume {
    Vec3 center;
    Vec3 half_extents;
    
    BoundingVolume() : center(Vec3::zero()), half_extents(Vec3::zero()) {}
    BoundingVolume(const Vec3& c, const Vec3& h) : center(c), half_extents(h) {}
    
    Vec3 min() const { return center - half_extents; }
    Vec3 max() const { return center + half_extents; }
    
    bool contains(const Vec3& point) const;
    bool intersects(const BoundingVolume& other) const;
    bool intersects_ray(const Ray& ray, float& tmin, float& tmax) const;
};

/**
 * @brief BVH Node
 */
struct BVHNode {
    BoundingVolume bounds;
    int left_child;   // -1 if leaf
    int right_child;  // -1 if leaf
    int start_index;  // For leaf nodes
    int count;        // Number of primitives in leaf
    
    bool is_leaf() const { return left_child == -1; }
};

/**
 * @brief BVH Builder
 */
class BVHBuilder {
public:
    struct Config {
        int max_primitives_per_leaf = 4;
        int max_depth = 64;
    } config;
    
private:
    std::vector<BVHNode> nodes;
    std::vector<BoundingVolume> primitive_bounds;
    
public:
    BVHBuilder() {}
    
    void build(const std::vector<BoundingVolume>& primitives);
    
    const std::vector<BVHNode>& get_nodes() const { return nodes; }
    int get_root() const { return nodes.empty() ? -1 : 0; }
    
private:
    int build_recursive(int start, int end, int depth);
    int partition_median(int start, int end, int axis);
    BoundingVolume compute_bounds(int start, int end);
};

/**
 * @brief Raycast Accelerator
 */
class RaycastAccelerator {
public:
    virtual ~RaycastAccelerator() = default;
    virtual void build() = 0;
    virtual bool raycast(const Ray& ray, RayHit& hit) = 0;
    virtual void raycast_all(const Ray& ray, std::vector<RayHit>& hits) = 0;
};

/**
 * @brief Brute Force Raycast Accelerator
 */
class BruteForceRaycast : public RaycastAccelerator {
public:
    std::vector<BoundingVolume> bounds;
    std::vector<void*> user_data;
    
    BruteForceRaycast() = default;
    
    void add(const BoundingVolume& bound, void* data = nullptr) {
        bounds.push_back(bound);
        user_data.push_back(data);
    }
    
    void build() override {}
    
    bool raycast(const Ray& ray, RayHit& hit) override;
    void raycast_all(const Ray& ray, std::vector<RayHit>& hits) override;
};

/**
 * @brief BVH Raycast Accelerator
 */
class BVHRaycast : public RaycastAccelerator {
private:
    std::vector<BVHNode> nodes;
    std::vector<BoundingVolume> primitive_bounds;
    std::vector<void*> primitive_data;
    int root;
    
public:
    BVHRaycast() : root(-1) {}
    
    void set_primitives(const std::vector<BoundingVolume>& bounds,
                        const std::vector<void*>& data = {});
    void build() override;
    
    bool raycast(const Ray& ray, RayHit& hit) override;
    void raycast_all(const Ray& ray, std::vector<RayHit>& hits) override;
    
private:
    bool raycast_node(int node_index, const Ray& ray, RayHit& hit);
    void raycast_node_all(int node_index, const Ray& ray, std::vector<RayHit>& hits);
};

/**
 * @brief Raycast Functions
 */
class Raycast {
public:
    // Ray-Geometry intersections
    static bool ray_sphere(const Ray& ray, const Vec3& center, float radius, 
                          RayHit& hit);
    static bool ray_box(const Ray& ray, const Vec3& min, const Vec3& max,
                        RayHit& hit);
    static bool ray_box(const Ray& ray, const Vec3& center, const Vec3& half_extents,
                        const Transform& transform, RayHit& hit);
    static bool ray_plane(const Ray& ray, const Vec3& plane_point, 
                         const Vec3& plane_normal, RayHit& hit);
    static bool ray_triangle(const Ray& ray, const Vec3& v0, const Vec3& v1, 
                             const Vec3& v2, RayHit& hit);
    static bool ray_disk(const Ray& ray, const Vec3& center, const Vec3& normal,
                        float radius, RayHit& hit);
    static bool ray_capsule(const Ray& ray, const Vec3& p1, const Vec3& p2,
                           float radius, RayHit& hit);
    
    // Multiple triangles
    static bool ray_mesh(const Ray& ray, const std::vector<Vec3>& vertices,
                        const std::vector<uint32_t>& indices,
                        std::vector<RayHit>& hits, bool closest_only = true);
    
    // Scene queries
    static std::vector<int> query_sphere(const Vec3& center, float radius,
                                        const std::vector<Vec3>& points);
    static std::vector<int> query_box(const Vec3& min, const Vec3& max,
                                      const std::vector<Vec3>& points);
    
    // Sweep tests
    static bool sweep_sphere(const Vec3& start, const Vec3& end, float radius,
                            const Vec3& obstacle_center, float obstacle_radius);
    static bool sweep_box(const Vec3& start, const Vec3& end,
                         const Vec3& box_half_extents,
                         const Vec3& obstacle_center, const Vec3& obstacle_half_extents);
};

/**
 * @brief Scene Query System
 */
class SceneQuery {
public:
    struct Sphere {
        Vec3 center;
        float radius;
    };
    
    struct Box {
        Vec3 center;
        Vec3 half_extents;
        Transform transform;
    };
    
    struct Frustum {
        Vec3 corners[8];
        Vec3 planes[6];  // Normal pointing inward
    };
    
    // Overlap tests
    static bool sphere_sphere(const Sphere& a, const Sphere& b);
    static bool sphere_box(const Sphere& sphere, const Box& box);
    static bool box_box(const Box& a, const Box& b);
    static bool sphere_frustum(const Sphere& sphere, const Frustum& frustum);
    
    // Containment tests
    static bool point_in_sphere(const Vec3& point, const Sphere& sphere);
    static bool point_in_box(const Vec3& point, const Box& box);
    static bool point_in_frustum(const Vec3& point, const Frustum& frustum);
    
    // Ray intersections
    static bool ray_sphere(const Ray& ray, const Sphere& sphere, RayHit& hit);
    static bool ray_box(const Ray& ray, const Box& box, RayHit& hit);
};

/**
 * @brief Spatial Hash Grid
 */
class SpatialHashGrid {
public:
    struct Config {
        float cell_size = 1.0f;
        int max_objects_per_cell = 16;
    } config;
    
private:
    struct Cell {
        std::vector<int> object_indices;
    };
    
    std::unordered_map<int64_t, Cell> cells;
    std::vector<void*> objects;
    std::vector<Vec3> object_positions;
    
    int64_t hash_position(int x, int y, int z) const;
    void get_cell_coords(const Vec3& pos, int& x, int& y, int& z) const;
    
public:
    SpatialHashGrid(const Config& cfg = Config()) : config(cfg) {}
    
    void clear();
    void insert(int object_id, const Vec3& position);
    void remove(int object_id);
    void update(int object_id, const Vec3& new_position);
    
    std::vector<int> query_point(const Vec3& point);
    std::vector<int> query_radius(const Vec3& center, float radius);
    std::vector<int> query_box(const Vec3& min, const Vec3& max);
};

} // namespace physics
} // namespace atlas
