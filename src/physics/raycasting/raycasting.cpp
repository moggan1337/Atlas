/**
 * @file raycasting.cpp
 * @brief Raycasting Implementation
 */

#include "raycasting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace atlas {
namespace physics {

// ============== BoundingVolume Implementation ==============

bool BoundingVolume::contains(const Vec3& point) const {
    Vec3 p = (point - center).abs();
    return p.x <= half_extents.x && p.y <= half_extents.y && p.z <= half_extents.z;
}

bool BoundingVolume::intersects(const BoundingVolume& other) const {
    Vec3 delta = other.center - center;
    Vec3 combined = half_extents + other.half_extents;
    return std::abs(delta.x) <= combined.x &&
           std::abs(delta.y) <= combined.y &&
           std::abs(delta.z) <= combined.z;
}

bool BoundingVolume::intersects_ray(const Ray& ray, float& tmin, float& tmax) const {
    Vec3 min_b = min();
    Vec3 max_b = max();
    
    tmin = 0.0f;
    tmax = ray.max_t;
    
    for (int i = 0; i < 3; ++i) {
        float origin = (i == 0) ? ray.origin.x : (i == 1) ? ray.origin.y : ray.origin.z;
        float dir = (i == 0) ? ray.direction.x : (i == 1) ? ray.direction.y : ray.direction.z;
        float min_v = (i == 0) ? min_b.x : (i == 1) ? min_b.y : min_b.z;
        float max_v = (i == 0) ? max_b.x : (i == 1) ? max_b.y : max_b.z;
        
        if (std::abs(dir) < 1e-6f) {
            if (origin < min_v || origin > max_v) return false;
        } else {
            float t1 = (min_v - origin) / dir;
            float t2 = (max_v - origin) / dir;
            
            if (t1 > t2) std::swap(t1, t2);
            
            tmin = std::max(tmin, t1);
            tmax = std::min(tmax, t2);
            
            if (tmin > tmax) return false;
        }
    }
    
    return tmin <= ray.max_t && tmax >= ray.min_t;
}

// ============== BVHBuilder Implementation ==============

void BVHBuilder::build(const std::vector<BoundingVolume>& primitives) {
    primitive_bounds = primitives;
    nodes.clear();
    
    if (primitives.empty()) return;
    
    nodes.reserve(primitives.size() * 2);
    build_recursive(0, static_cast<int>(primitives.size()), 0);
}

int BVHBuilder::build_recursive(int start, int end, int depth) {
    if (depth > config.max_depth || end - start <= config.max_primitives_per_leaf) {
        BVHNode leaf;
        leaf.bounds = compute_bounds(start, end);
        leaf.left_child = -1;
        leaf.right_child = -1;
        leaf.start_index = start;
        leaf.count = end - start;
        nodes.push_back(leaf);
        return static_cast<int>(nodes.size()) - 1;
    }
    
    // Find best split axis
    BoundingVolume bounds = compute_bounds(start, end);
    Vec3 extent = bounds.half_extents * 2.0f;
    int best_axis = (extent.x > extent.y && extent.x > extent.z) ? 0 :
                   (extent.y > extent.z) ? 1 : 2;
    
    // Partition
    int mid = partition_median(start, end, best_axis);
    
    // Ensure valid split
    if (mid == start || mid == end) {
        mid = (start + end) / 2;
    }
    
    // Build children
    BVHNode node;
    node.bounds = bounds;
    node.left_child = build_recursive(start, mid, depth + 1);
    node.right_child = build_recursive(mid, end, depth + 1);
    node.start_index = -1;
    node.count = 0;
    
    nodes.push_back(node);
    return static_cast<int>(nodes.size()) - 1;
}

int BVHBuilder::partition_median(int start, int end, int axis) {
    int mid = (start + end) / 2;
    
    std::nth_element(
        primitive_bounds.begin() + start,
        primitive_bounds.begin() + mid,
        primitive_bounds.begin() + end,
        [axis](const BoundingVolume& a, const BoundingVolume& b) {
            float ca = (axis == 0) ? a.center.x : (axis == 1) ? a.center.y : a.center.z;
            float cb = (axis == 0) ? b.center.x : (axis == 1) ? b.center.y : b.center.z;
            return ca < cb;
        }
    );
    
    return mid;
}

BoundingVolume BVHBuilder::compute_bounds(int start, int end) {
    if (start >= end) return BoundingVolume();
    
    Vec3 min_p = primitive_bounds[start].min();
    Vec3 max_p = primitive_bounds[start].max();
    
    for (int i = start + 1; i < end; ++i) {
        Vec3 pmin = primitive_bounds[i].min();
        Vec3 pmax = primitive_bounds[i].max();
        min_p = min_p.min(pmin);
        max_p = max_p.max(pmax);
    }
    
    Vec3 center = (min_p + max_p) * 0.5f;
    Vec3 half_extents = (max_p - min_p) * 0.5f;
    
    return BoundingVolume(center, half_extents);
}

// ============== BVHRaycast Implementation ==============

void BVHRaycast::set_primitives(const std::vector<BoundingVolume>& bounds,
                                const std::vector<void*>& data) {
    primitive_bounds = bounds;
    primitive_data = data;
}

void BVHRaycast::build() {
    if (primitive_bounds.empty()) {
        root = -1;
        return;
    }
    
    BVHBuilder builder;
    builder.build(primitive_bounds);
    nodes = builder.get_nodes();
    root = builder.get_root();
}

bool BVHRaycast::raycast(const Ray& ray, RayHit& hit) {
    hit = RayHit();
    
    if (root < 0 || nodes.empty()) return false;
    
    return raycast_node(root, ray, hit);
}

bool BVHRaycast::raycast_node(int node_index, const Ray& ray, RayHit& hit) {
    if (node_index < 0 || node_index >= static_cast<int>(nodes.size())) {
        return false;
    }
    
    const BVHNode& node = nodes[node_index];
    
    float tmin, tmax;
    if (!node.bounds.intersects_ray(ray, tmin, tmax)) {
        return false;
    }
    
    if (node.is_leaf()) {
        bool any_hit = false;
        
        for (int i = 0; i < node.count; ++i) {
            int prim_idx = node.start_index + i;
            
            if (prim_idx < static_cast<int>(primitive_bounds.size())) {
                // Simple center-radius test
                const auto& bound = primitive_bounds[prim_idx];
                RayHit temp_hit;
                
                if (Raycast::ray_sphere(ray, bound.center, bound.half_extents.x, temp_hit)) {
                    if (temp_hit.t < hit.t) {
                        hit = temp_hit;
                        if (!primitive_data.empty()) {
                            hit.user_data = primitive_data[prim_idx];
                        }
                        any_hit = true;
                    }
                }
            }
        }
        
        return any_hit;
    } else {
        // Traverse children
        RayHit left_hit, right_hit;
        bool left = raycast_node(node.left_child, ray, left_hit);
        bool right = raycast_node(node.right_child, ray, right_hit);
        
        if (left && right) {
            hit = (left_hit.t < right_hit.t) ? left_hit : right_hit;
            return true;
        } else if (left) {
            hit = left_hit;
            return true;
        } else if (right) {
            hit = right_hit;
            return true;
        }
        
        return false;
    }
}

void BVHRaycast::raycast_all(const Ray& ray, std::vector<RayHit>& hits) {
    hits.clear();
    if (root >= 0) {
        raycast_node_all(root, ray, hits);
    }
}

void BVHRaycast::raycast_node_all(int node_index, const Ray& ray, std::vector<RayHit>& hits) {
    if (node_index < 0 || node_index >= static_cast<int>(nodes.size())) return;
    
    const BVHNode& node = nodes[node_index];
    
    float tmin, tmax;
    if (!node.bounds.intersects_ray(ray, tmin, tmax)) return;
    
    if (node.is_leaf()) {
        for (int i = 0; i < node.count; ++i) {
            int prim_idx = node.start_index + i;
            
            if (prim_idx < static_cast<int>(primitive_bounds.size())) {
                const auto& bound = primitive_bounds[prim_idx];
                RayHit hit;
                
                if (Raycast::ray_sphere(ray, bound.center, bound.half_extents.x, hit)) {
                    if (!primitive_data.empty()) {
                        hit.user_data = primitive_data[prim_idx];
                    }
                    hits.push_back(hit);
                }
            }
        }
    } else {
        raycast_node_all(node.left_child, ray, hits);
        raycast_node_all(node.right_child, ray, hits);
    }
}

// ============== Raycast Implementation ==============

bool Raycast::ray_sphere(const Ray& ray, const Vec3& center, float radius, RayHit& hit) {
    Vec3 oc = ray.origin - center;
    float a = ray.direction.dot(ray.direction);
    float b = 2.0f * oc.dot(ray.direction);
    float c = oc.dot(oc) - radius * radius;
    float discriminant = b * b - 4 * a * c;
    
    if (discriminant < 0.0f) {
        hit.hit = false;
        return false;
    }
    
    float t1 = (-b - std::sqrt(discriminant)) / (2.0f * a);
    float t2 = (-b + std::sqrt(discriminant)) / (2.0f * a);
    
    float t = t1;
    if (t < ray.min_t) t = t2;
    
    if (t < ray.min_t || t > ray.max_t) {
        hit.hit = false;
        return false;
    }
    
    hit.hit = true;
    hit.t = t;
    hit.point = ray.get_point(t);
    hit.normal = (hit.point - center).normalized();
    
    return true;
}

bool Raycast::ray_box(const Ray& ray, const Vec3& min, const Vec3& max, RayHit& hit) {
    float tmin = ray.min_t;
    float tmax = ray.max_t;
    
    for (int i = 0; i < 3; ++i) {
        float origin = (i == 0) ? ray.origin.x : (i == 1) ? ray.origin.y : ray.origin.z;
        float dir = (i == 0) ? ray.direction.x : (i == 1) ? ray.direction.y : ray.direction.z;
        float min_v = (i == 0) ? min.x : (i == 1) ? min.y : min.z;
        float max_v = (i == 0) ? max.x : (i == 1) ? max.y : max.z;
        
        if (std::abs(dir) < 1e-6f) {
            if (origin < min_v || origin > max_v) {
                hit.hit = false;
                return false;
            }
        } else {
            float t1 = (min_v - origin) / dir;
            float t2 = (max_v - origin) / dir;
            
            if (t1 > t2) std::swap(t1, t2);
            
            tmin = std::max(tmin, t1);
            tmax = std::min(tmax, t2);
            
            if (tmin > tmax) {
                hit.hit = false;
                return false;
            }
        }
    }
    
    hit.hit = true;
    hit.t = tmin;
    hit.point = ray.get_point(tmin);
    
    // Calculate normal
    Vec3 p = hit.point;
    Vec3 d = (p - min) / (max - min);
    
    float max_comp = std::max(std::abs(d.x - 0.5f), 
                             std::max(std::abs(d.y - 0.5f), std::abs(d.z - 0.5f)));
    
    if (max_comp == std::abs(d.x - 0.5f)) {
        hit.normal = Vec3(p.x < (min.x + max.x) * 0.5f ? -1.0f : 1.0f, 0, 0);
    } else if (max_comp == std::abs(d.y - 0.5f)) {
        hit.normal = Vec3(0, p.y < (min.y + max.y) * 0.5f ? -1.0f : 1.0f, 0);
    } else {
        hit.normal = Vec3(0, 0, p.z < (min.z + max.z) * 0.5f ? -1.0f : 1.0f);
    }
    
    return true;
}

bool Raycast::ray_box(const Ray& ray, const Vec3& center, const Vec3& half_extents,
                      const Transform& transform, RayHit& hit) {
    // Transform ray to local space
    Transform inv = transform.inverse();
    Ray local_ray;
    local_ray.origin = inv.transform_point(ray.origin);
    local_ray.direction = inv.transform_direction(ray.direction);
    local_ray.min_t = ray.min_t;
    local_ray.max_t = ray.max_t;
    
    Vec3 local_min = -half_extents;
    Vec3 local_max = half_extents;
    
    RayHit local_hit;
    if (ray_box(local_ray, local_min, local_max, local_hit)) {
        // Transform hit back to world space
        hit.hit = true;
        hit.t = local_hit.t;
        hit.point = transform.transform_point(local_hit.point);
        hit.normal = transform.transform_direction(local_hit.normal);
        hit.normal.normalize();
        return true;
    }
    
    hit.hit = false;
    return false;
}

bool Raycast::ray_plane(const Ray& ray, const Vec3& plane_point, 
                       const Vec3& plane_normal, RayHit& hit) {
    float denom = ray.direction.dot(plane_normal);
    
    if (std::abs(denom) < 1e-6f) {
        hit.hit = false;
        return false;
    }
    
    float t = (plane_point - ray.origin).dot(plane_normal) / denom;
    
    if (t < ray.min_t || t > ray.max_t) {
        hit.hit = false;
        return false;
    }
    
    hit.hit = true;
    hit.t = t;
    hit.point = ray.get_point(t);
    hit.normal = denom > 0.0f ? -plane_normal : plane_normal;
    
    return true;
}

bool Raycast::ray_triangle(const Ray& ray, const Vec3& v0, const Vec3& v1,
                          const Vec3& v2, RayHit& hit) {
    // Möller–Trumbore algorithm
    Vec3 e1 = v1 - v0;
    Vec3 e2 = v2 - v0;
    Vec3 pvec = ray.direction.cross(e2);
    
    float det = e1.dot(pvec);
    
    if (std::abs(det) < 1e-6f) {
        hit.hit = false;
        return false;
    }
    
    float inv_det = 1.0f / det;
    Vec3 tvec = ray.origin - v0;
    
    float u = tvec.dot(pvec) * inv_det;
    if (u < 0.0f || u > 1.0f) {
        hit.hit = false;
        return false;
    }
    
    Vec3 qvec = tvec.cross(e1);
    float v = ray.direction.dot(qvec) * inv_det;
    if (v < 0.0f || u + v > 1.0f) {
        hit.hit = false;
        return false;
    }
    
    float t = e2.dot(qvec) * inv_det;
    
    if (t < ray.min_t || t > ray.max_t) {
        hit.hit = false;
        return false;
    }
    
    hit.hit = true;
    hit.t = t;
    hit.point = ray.get_point(t);
    hit.normal = e1.cross(e2).normalized();
    hit.barycentric = Vec3(1.0f - u - v, u, v);
    hit.u = u;
    hit.v = v;
    
    return true;
}

bool Raycast::ray_disk(const Ray& ray, const Vec3& center, const Vec3& normal,
                       float radius, RayHit& hit) {
    float denom = ray.direction.dot(normal);
    
    if (std::abs(denom) < 1e-6f) {
        hit.hit = false;
        return false;
    }
    
    float t = (center - ray.origin).dot(normal) / denom;
    
    if (t < ray.min_t || t > ray.max_t) {
        hit.hit = false;
        return false;
    }
    
    Vec3 point = ray.get_point(t);
    if ((point - center).magnitude_squared() > radius * radius) {
        hit.hit = false;
        return false;
    }
    
    hit.hit = true;
    hit.t = t;
    hit.point = point;
    hit.normal = denom > 0.0f ? -normal : normal;
    
    return true;
}

bool Raycast::ray_capsule(const Ray& ray, const Vec3& p1, const Vec3& p2,
                         float radius, RayHit& hit) {
    Vec3 axis = (p2 - p1).normalized();
    Vec3 oc = ray.origin - p1;
    
    float a = ray.direction.dot(ray.direction) - std::pow(ray.direction.dot(axis), 2.0f);
    float b = 2.0f * (ray.direction.dot(oc) - ray.direction.dot(axis) * oc.dot(axis));
    float c = oc.dot(oc) - oc.dot(axis) * oc.dot(axis) - radius * radius;
    
    float discriminant = b * b - 4.0f * a * c;
    
    if (discriminant >= 0.0f) {
        float t1 = (-b - std::sqrt(discriminant)) / (2.0f * a);
        float t2 = (-b + std::sqrt(discriminant)) / (2.0f * a);
        
        float t = t1;
        if (t < ray.min_t) t = t2;
        
        if (t >= ray.min_t && t <= ray.max_t) {
            Vec3 point = ray.get_point(t);
            float h = (point - p1).dot(axis);
            float len = (p2 - p1).magnitude();
            
            if (h >= 0.0f && h <= len) {
                hit.hit = true;
                hit.t = t;
                hit.point = point;
                hit.normal = (point - (p1 + axis * h)).normalized();
                return true;
            }
        }
    }
    
    // Check end caps
    RayHit sphere_hit;
    if (ray_sphere(ray, p1, radius, sphere_hit) && sphere_hit.t < hit.t) {
        hit = sphere_hit;
        return true;
    }
    if (ray_sphere(ray, p2, radius, sphere_hit) && sphere_hit.t < hit.t) {
        hit = sphere_hit;
        return true;
    }
    
    hit.hit = false;
    return false;
}

bool Raycast::ray_mesh(const Ray& ray, const std::vector<Vec3>& vertices,
                      const std::vector<uint32_t>& indices,
                      std::vector<RayHit>& hits, bool closest_only) {
    hits.clear();
    RayHit best_hit;
    best_hit.t = INFINITY;
    
    for (size_t i = 0; i < indices.size(); i += 3) {
        const Vec3& v0 = vertices[indices[i]];
        const Vec3& v1 = vertices[indices[i + 1]];
        const Vec3& v2 = vertices[indices[i + 2]];
        
        RayHit hit;
        if (ray_triangle(ray, v0, v1, v2, hit)) {
            if (closest_only) {
                if (hit.t < best_hit.t) {
                    best_hit = hit;
                    best_hit.triangle_index = static_cast<int>(i / 3);
                }
            } else {
                hit.triangle_index = static_cast<int>(i / 3);
                hits.push_back(hit);
            }
        }
    }
    
    if (closest_only && best_hit.t < INFINITY) {
        hits.push_back(best_hit);
        return true;
    }
    
    return !hits.empty();
}

std::vector<int> Raycast::query_sphere(const Vec3& center, float radius,
                                       const std::vector<Vec3>& points) {
    std::vector<int> result;
    float radius_sq = radius * radius;
    
    for (size_t i = 0; i < points.size(); ++i) {
        if ((points[i] - center).magnitude_squared() <= radius_sq) {
            result.push_back(static_cast<int>(i));
        }
    }
    
    return result;
}

std::vector<int> Raycast::query_box(const Vec3& min, const Vec3& max,
                                   const std::vector<Vec3>& points) {
    std::vector<int> result;
    
    for (size_t i = 0; i < points.size(); ++i) {
        const Vec3& p = points[i];
        if (p.x >= min.x && p.x <= max.x &&
            p.y >= min.y && p.y <= max.y &&
            p.z >= min.z && p.z <= max.z) {
            result.push_back(static_cast<int>(i));
        }
    }
    
    return result;
}

bool Raycast::sweep_sphere(const Vec3& start, const Vec3& end, float radius,
                          const Vec3& obstacle_center, float obstacle_radius) {
    Vec3 movement = end - start;
    float movement_length = movement.magnitude();
    Vec3 ray_dir = movement_length > 0.0f ? movement / movement_length : Vec3::zero();
    
    Ray ray(start, ray_dir, 0.0f, movement_length);
    RayHit hit;
    
    return ray_sphere(ray, obstacle_center, radius + obstacle_radius, hit) &&
           hit.t <= movement_length;
}

bool Raycast::sweep_box(const Vec3& start, const Vec3& end,
                       const Vec3& box_half_extents,
                       const Vec3& obstacle_center, const Vec3& obstacle_half_extents) {
    Vec3 movement = end - start;
    float movement_length = movement.magnitude();
    Vec3 ray_dir = movement_length > 0.0f ? movement / movement_length : Vec3::zero();
    
    // Simple AABB sweep
    Vec3 box_min = start - box_half_extents;
    Vec3 box_max = start + box_half_extents;
    Vec3 obs_min = obstacle_center - obstacle_half_extents;
    Vec3 obs_max = obstacle_center + obstacle_half_extents;
    
    Ray ray(start, ray_dir, 0.0f, movement_length);
    RayHit hit;
    
    return ray_box(ray, obs_min, obs_max, hit) && hit.t <= movement_length;
}

// ============== SceneQuery Implementation ==============

bool SceneQuery::sphere_sphere(const Sphere& a, const Sphere& b) {
    float dist_sq = (a.center - b.center).magnitude_squared();
    float radius_sum = a.radius + b.radius;
    return dist_sq <= radius_sum * radius_sum;
}

bool SceneQuery::sphere_box(const Sphere& sphere, const Box& box) {
    Vec3 closest = sphere.center;
    
    if (box.transform.rotation.x != 0 || box.transform.rotation.y != 0 ||
        box.transform.rotation.z != 0 || box.transform.rotation.w != 1) {
        // For rotated boxes, transform sphere to box space
        Transform inv = box.transform.inverse();
        closest = inv.transform_point(sphere.center);
    } else {
        closest.x = std::max(box.center.x - box.half_extents.x,
                            std::min(closest.x, box.center.x + box.half_extents.x));
        closest.y = std::max(box.center.y - box.half_extents.y,
                            std::min(closest.y, box.center.y + box.half_extents.y));
        closest.z = std::max(box.center.z - box.half_extents.z,
                            std::min(closest.z, box.center.z + box.half_extents.z));
    }
    
    return (closest - sphere.center).magnitude_squared() <= sphere.radius * sphere.radius;
}

bool SceneQuery::box_box(const Box& a, const Box& b) {
    // Simplified - assumes axis-aligned
    Vec3 a_min = a.center - a.half_extents;
    Vec3 a_max = a.center + a.half_extents;
    Vec3 b_min = b.center - b.half_extents;
    Vec3 b_max = b.center + b.half_extents;
    
    return (a_min.x <= b_max.x && a_max.x >= b_min.x) &&
           (a_min.y <= b_max.y && a_max.y >= b_min.y) &&
           (a_min.z <= b_max.z && a_max.z >= b_min.z);
}

bool SceneQuery::sphere_frustum(const Sphere& sphere, const Frustum& frustum) {
    for (int i = 0; i < 6; ++i) {
        float dist = (sphere.center - frustum.planes[i]).dot(frustum.planes[i]);
        if (dist < -sphere.radius) return false;
    }
    return true;
}

bool SceneQuery::point_in_sphere(const Vec3& point, const Sphere& sphere) {
    return (point - sphere.center).magnitude_squared() <= sphere.radius * sphere.radius;
}

bool SceneQuery::point_in_box(const Vec3& point, const Box& box) {
    Vec3 local = point - box.center;
    return std::abs(local.x) <= box.half_extents.x &&
           std::abs(local.y) <= box.half_extents.y &&
           std::abs(local.z) <= box.half_extents.z;
}

bool SceneQuery::point_in_frustum(const Vec3& point, const Frustum& frustum) {
    for (int i = 0; i < 6; ++i) {
        if ((point - frustum.planes[i]).dot(frustum.planes[i]) < 0.0f) {
            return false;
        }
    }
    return true;
}

bool SceneQuery::ray_sphere(const Ray& ray, const Sphere& sphere, RayHit& hit) {
    return Raycast::ray_sphere(ray, sphere.center, sphere.radius, hit);
}

bool SceneQuery::ray_box(const Ray& ray, const Box& box, RayHit& hit) {
    Vec3 min_b = box.center - box.half_extents;
    Vec3 max_b = box.center + box.half_extents;
    return Raycast::ray_box(ray, min_b, max_b, hit);
}

// ============== SpatialHashGrid Implementation ==============

int64_t SpatialHashGrid::hash_position(int x, int y, int z) const {
    return (static_cast<int64_t>(x) * 73856093LL) ^
           (static_cast<int64_t>(y) * 19349663LL) ^
           (static_cast<int64_t>(z) * 83492791LL);
}

void SpatialHashGrid::get_cell_coords(const Vec3& pos, int& x, int& y, int& z) const {
    x = static_cast<int>(std::floor(pos.x / config.cell_size));
    y = static_cast<int>(std::floor(pos.y / config.cell_size));
    z = static_cast<int>(std::floor(pos.z / config.cell_size));
}

void SpatialHashGrid::clear() {
    cells.clear();
    objects.clear();
    object_positions.clear();
}

void SpatialHashGrid::insert(int object_id, const Vec3& position) {
    int x, y, z;
    get_cell_coords(position, x, y, z);
    int64_t hash = hash_position(x, y, z);
    
    if (object_id >= static_cast<int>(objects.size())) {
        objects.resize(object_id + 1, nullptr);
        object_positions.resize(object_id + 1, Vec3::zero());
    }
    
    objects[object_id] = reinterpret_cast<void*>(object_id + 1);
    object_positions[object_id] = position;
    
    cells[hash].object_indices.push_back(object_id);
}

void SpatialHashGrid::remove(int object_id) {
    if (object_id >= static_cast<int>(object_positions.size())) return;
    
    Vec3 pos = object_positions[object_id];
    int x, y, z;
    get_cell_coords(pos, x, y, z);
    int64_t hash = hash_position(x, y, z);
    
    auto& cell = cells[hash];
    cell.object_indices.erase(
        std::remove(cell.object_indices.begin(), cell.object_indices.end(), object_id),
        cell.object_indices.end()
    );
    
    objects[object_id] = nullptr;
}

void SpatialHashGrid::update(int object_id, const Vec3& new_position) {
    if (object_id >= static_cast<int>(object_positions.size())) {
        insert(object_id, new_position);
        return;
    }
    
    Vec3 old_pos = object_positions[object_id];
    int ox, oy, oz, nx, ny, nz;
    get_cell_coords(old_pos, ox, oy, oz);
    get_cell_coords(new_position, nx, ny, nz);
    
    if (ox == nx && oy == ny && oz == nz) {
        object_positions[object_id] = new_position;
        return;
    }
    
    // Remove from old cell
    remove(object_id);
    
    // Insert into new cell
    insert(object_id, new_position);
}

std::vector<int> SpatialHashGrid::query_point(const Vec3& point) {
    std::vector<int> result;
    int x, y, z;
    get_cell_coords(point, x, y, z);
    
    int64_t hash = hash_position(x, y, z);
    auto it = cells.find(hash);
    if (it != cells.end()) {
        for (int id : it->second.object_indices) {
            if (objects[id] != nullptr) {
                result.push_back(id);
            }
        }
    }
    
    return result;
}

std::vector<int> SpatialHashGrid::query_radius(const Vec3& center, float radius) {
    int min_x, min_y, min_z;
    int max_x, max_y, max_z;
    
    get_cell_coords(center - Vec3(radius), min_x, min_y, min_z);
    get_cell_coords(center + Vec3(radius), max_x, max_y, max_z);
    
    std::vector<int> result;
    float radius_sq = radius * radius;
    
    for (int x = min_x; x <= max_x; ++x) {
        for (int y = min_y; y <= max_y; ++y) {
            for (int z = min_z; z <= max_z; ++z) {
                int64_t hash = hash_position(x, y, z);
                auto it = cells.find(hash);
                
                if (it != cells.end()) {
                    for (int id : it->second.object_indices) {
                        if (objects[id] != nullptr) {
                            if ((object_positions[id] - center).magnitude_squared() <= radius_sq) {
                                result.push_back(id);
                            }
                        }
                    }
                }
            }
        }
    }
    
    return result;
}

std::vector<int> SpatialHashGrid::query_box(const Vec3& min, const Vec3& max) {
    int min_x, min_y, min_z;
    int max_x, max_y, max_z;
    
    get_cell_coords(min, min_x, min_y, min_z);
    get_cell_coords(max, max_x, max_y, max_z);
    
    std::vector<int> result;
    
    for (int x = min_x; x <= max_x; ++x) {
        for (int y = min_y; y <= max_y; ++y) {
            for (int z = min_z; z <= max_z; ++z) {
                int64_t hash = hash_position(x, y, z);
                auto it = cells.find(hash);
                
                if (it != cells.end()) {
                    for (int id : it->second.object_indices) {
                        if (objects[id] != nullptr) {
                            const Vec3& p = object_positions[id];
                            if (p.x >= min.x && p.x <= max.x &&
                                p.y >= min.y && p.y <= max.y &&
                                p.z >= min.z && p.z <= max.z) {
                                result.push_back(id);
                            }
                        }
                    }
                }
            }
        }
    }
    
    return result;
}

} // namespace physics
} // namespace atlas
