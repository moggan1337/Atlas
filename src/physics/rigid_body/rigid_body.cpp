/**
 * @file rigid_body.cpp
 * @brief Rigid Body Implementation
 */

#include "rigid_body.hpp"
#include <algorithm>
#include <cmath>

namespace atlas {
namespace physics {

// ============== RigidBody Implementation ==============

RigidBody::RigidBody() 
    : linear_velocity(Vec3::zero())
    , angular_velocity(Vec3::zero())
    , force_accumulator(Vec3::zero())
    , torque_accumulator(Vec3::zero())
{
    static uint32_t next_id = 0;
    id = next_id++;
}

RigidBody::RigidBody(float m, std::shared_ptr<CollisionShape> s)
    : RigidBody()
{
    shape = s;
    set_mass(m);
    if (m > 0) {
        local_inertia = shape->compute_inertia(m);
        inv_inertia = Vec3(1.0f / local_inertia.x,
                          1.0f / local_inertia.y,
                          1.0f / local_inertia.z);
    }
}

// Factory methods
std::shared_ptr<RigidBody> RigidBody::create_sphere(float radius, float mass) {
    auto shape = std::make_shared<SphereShape>(radius);
    return std::make_shared<RigidBody>(mass, shape);
}

std::shared_ptr<RigidBody> RigidBody::create_box(const Vec3& half_extents, float mass) {
    auto shape = std::make_shared<BoxShape>(half_extents);
    return std::make_shared<RigidBody>(mass, shape);
}

std::shared_ptr<RigidBody> RigidBody::create_capsule(float radius, float height, float mass) {
    auto shape = std::make_shared<CapsuleShape>(radius, height);
    return std::make_shared<RigidBody>(mass, shape);
}

std::shared_ptr<RigidBody> RigidBody::create_plane(const Vec3& normal, float offset, float mass) {
    auto shape = std::make_shared<PlaneShape>(normal, offset);
    return std::make_shared<RigidBody>(mass, shape);
}

std::shared_ptr<RigidBody> RigidBody::create_static_plane(const Vec3& normal, float offset) {
    auto body = create_plane(normal, offset, 0.0f);
    body->state = BodyState::Static;
    return body;
}

// Mass properties
void RigidBody::set_mass(float m) {
    mass = m;
    if (m > 0.0f) {
        inv_mass = 1.0f / m;
    } else {
        inv_mass = 0.0f;
        local_inertia = Vec3::zero();
        inv_inertia = Vec3::zero();
    }
}

void RigidBody::set_density(float density) {
    float vol = shape ? shape->get_volume() : 1.0f;
    set_mass(vol * density);
}

void RigidBody::set_mass_properties(float m, const Vec3& inertia) {
    mass = m;
    local_inertia = inertia;
    if (m > 0.0f) {
        inv_mass = 1.0f / m;
        inv_inertia = Vec3(1.0f / std::max(inertia.x, 1e-6f),
                          1.0f / std::max(inertia.y, 1e-6f),
                          1.0f / std::max(inertia.z, 1e-6f));
    } else {
        inv_mass = 0.0f;
        inv_inertia = Vec3::zero();
    }
}

void RigidBody::set_local_inertia(const Vec3& inertia) {
    local_inertia = inertia;
    inv_inertia = Vec3(1.0f / std::max(inertia.x, 1e-6f),
                       1.0f / std::max(inertia.y, 1e-6f),
                       1.0f / std::max(inertia.z, 1e-6f));
}

void RigidBody::set_infinite_mass() {
    mass = 0.0f;
    inv_mass = 0.0f;
    local_inertia = Vec3::zero();
    inv_inertia = Vec3::zero();
    state = BodyState::Static;
}

// State management
void RigidBody::activate() {
    is_sleeping = false;
    sleep_timer = 0.0f;
    state = BodyState::Active;
}

void RigidBody::deactivate() {
    is_sleeping = true;
    sleep_timer = 0.0f;
    state = BodyState::Sleeping;
    linear_velocity = Vec3::zero();
    angular_velocity = Vec3::zero();
}

void RigidBody::set_kinematic() {
    state = BodyState::Kinematic;
    inv_mass = 0.0f;
}

void RigidBody::set_static() {
    state = BodyState::Static;
    inv_mass = 0.0f;
    inv_inertia = Vec3::zero();
}

void RigidBody::update_sleep_state(float dt) {
    if (!allow_sleep || is_static() || is_kinematic()) return;

    float speed = linear_velocity.magnitude();
    float ang_speed = angular_velocity.magnitude();

    if (speed < sleep_threshold && ang_speed < sleep_threshold) {
        sleep_timer += dt;
        if (sleep_timer > sleep_threshold * 2.0f) {
            deactivate();
        }
    } else {
        sleep_timer = 0.0f;
        if (is_sleeping) {
            activate();
        }
    }
}

// Force application
void RigidBody::apply_force(const Vec3& force) {
    if (!has_finite_mass()) return;
    force_accumulator += force;
    if (allow_sleep) activate();
}

void RigidBody::apply_force_at_point(const Vec3& force, const Vec3& point) {
    if (!has_finite_mass()) return;
    force_accumulator += force;
    torque_accumulator += (point - transform.position).cross(force);
    if (allow_sleep) activate();
}

void RigidBody::apply_force_at_body_point(const Vec3& force, const Vec3& local_point) {
    apply_force_at_point(force, transform.position + transform.rotate(local_point));
}

void RigidBody::apply_torque(const Vec3& torque) {
    if (!has_finite_mass()) return;
    torque_accumulator += torque;
    if (allow_sleep) activate();
}

void RigidBody::apply_impulse(const Vec3& impulse) {
    if (!has_finite_mass()) return;
    linear_velocity += impulse * inv_mass;
    if (allow_sleep) activate();
}

void RigidBody::apply_impulse_at_point(const Vec3& impulse, const Vec3& point) {
    if (!has_finite_mass()) return;
    linear_velocity += impulse * inv_mass;
    
    Vec3 r = point - transform.position;
    Vec3 angular_impulse = r.cross(impulse);
    angular_velocity += angular_impulse.component_mul(inv_inertia);
    
    if (allow_sleep) activate();
}

void RigidBody::apply_impulse_at_body_point(const Vec3& impulse, const Vec3& local_point) {
    apply_impulse_at_point(impulse, transform.position + transform.rotate(local_point));
}

void RigidBody::apply_angular_impulse(const Vec3& impulse) {
    if (!has_finite_mass()) return;
    angular_velocity += impulse.component_mul(inv_inertia);
    if (allow_sleep) activate();
}

void RigidBody::clear_forces() {
    force_accumulator = Vec3::zero();
}

void RigidBody::clear_torques() {
    torque_accumulator = Vec3::zero();
}

void RigidBody::clear_accumulators() {
    force_accumulator = Vec3::zero();
    torque_accumulator = Vec3::zero();
}

void RigidBody::apply_gravity(const Vec3& gravity) {
    if (!has_finite_mass()) return;
    force_accumulator += gravity * mass * gravity_scale;
}

Vec3 RigidBody::get_gravity_force(const Vec3& gravity) const {
    if (!has_finite_mass()) return Vec3::zero();
    return gravity * mass * gravity_scale;
}

// Integration
void RigidBody::integrate(float dt) {
    integrate_velocities(dt);
    integrate_positions(dt);
}

void RigidBody::integrate_velocities(float dt) {
    if (!has_finite_mass() || is_sleeping) return;

    // Linear velocity integration: v += (F/m) * dt
    Vec3 acceleration = force_accumulator * inv_mass;
    linear_velocity += acceleration * dt;

    // Angular velocity integration: ω += I⁻¹ * τ * dt
    Vec3 angular_acceleration = torque_accumulator.component_mul(inv_inertia);
    angular_velocity += angular_acceleration * dt;

    // Apply damping
    apply_damping(dt);
}

void RigidBody::integrate_positions(float dt) {
    if (!has_finite_mass() || is_sleeping) return;

    // Position integration: p += v * dt
    transform.position += linear_velocity * dt;

    // Orientation integration using quaternion
    // dq/dt = 0.5 * ω * q
    float half_w = 0.5f * angular_velocity.dot(Vec3(transform.rotation.x,
                                                      transform.rotation.y,
                                                      transform.rotation.z));
    float half_x = 0.5f * (transform.rotation.w * angular_velocity.x +
                          transform.rotation.y * angular_velocity.z -
                          transform.rotation.z * angular_velocity.y);
    float half_y = 0.5f * (transform.rotation.w * angular_velocity.y +
                          transform.rotation.z * angular_velocity.x -
                          transform.rotation.x * angular_velocity.z);
    float half_z = 0.5f * (transform.rotation.w * angular_velocity.z +
                          transform.rotation.x * angular_velocity.y -
                          transform.rotation.y * angular_velocity.x);

    float w_new = transform.rotation.w + half_w * dt;
    float x_new = transform.rotation.x + half_x * dt;
    float y_new = transform.rotation.y + half_y * dt;
    float z_new = transform.rotation.z + half_z * dt;

    transform.rotation = Quaternion(w_new, x_new, y_new, z_new).normalized();
}

void RigidBody::apply_damping(float dt) {
    float ld = std::pow(1.0f - linear_damping, dt);
    float ad = std::pow(1.0f - angular_damping, dt);
    linear_velocity *= ld;
    angular_velocity *= ad;
}

Vec3 RigidBody::get_world_center_of_mass() const {
    return transform.position;
}

Vec3 RigidBody::get_world_point_velocity(const Vec3& local_point) const {
    Vec3 world_point = transform.position + transform.rotate(local_point);
    return linear_velocity + angular_velocity.cross(world_point - transform.position);
}

Vec3 RigidBody::get_local_point_velocity(const Vec3& world_point) const {
    Vec3 r = world_point - transform.position;
    return linear_velocity + angular_velocity.cross(r);
}

void RigidBody::set_linear_velocity(const Vec3& velocity) {
    linear_velocity = velocity;
    if (allow_sleep && velocity.magnitude() > sleep_threshold) {
        activate();
    }
}

void RigidBody::set_angular_velocity(const Vec3& velocity) {
    angular_velocity = velocity;
    if (allow_sleep && velocity.magnitude() > sleep_threshold) {
        activate();
    }
}

void RigidBody::add_linear_velocity(const Vec3& velocity) {
    set_linear_velocity(linear_velocity + velocity);
}

void RigidBody::add_angular_velocity(const Vec3& velocity) {
    set_angular_velocity(angular_velocity + velocity);
}

void RigidBody::set_position(const Vec3& position) {
    transform.position = position;
}

void RigidBody::set_rotation(const Quaternion& rotation) {
    transform.rotation = rotation.normalized();
}

void RigidBody::set_transform(const Transform& t) {
    transform = t;
}

Vec3 RigidBody::get_aabb_min() const {
    if (shape) {
        return shape->get_aabb_min(transform);
    }
    return transform.position;
}

Vec3 RigidBody::get_aabb_max() const {
    if (shape) {
        return shape->get_aabb_max(transform);
    }
    return transform.position;
}

float RigidBody::get_kinetic_energy() const {
    return get_linear_kinetic_energy() + get_angular_kinetic_energy();
}

float RigidBody::get_linear_kinetic_energy() const {
    return 0.5f * mass * linear_velocity.magnitude_squared();
}

float RigidBody::get_angular_kinetic_energy() const {
    Vec3 omega_squared = angular_velocity.component_mul(angular_velocity);
    Vec3 inertia_tensor = local_inertia.component_mul(inv_inertia);
    return 0.5f * (omega_squared.x * local_inertia.x +
                   omega_squared.y * local_inertia.y +
                   omega_squared.z * local_inertia.z);
}

Vec3 RigidBody::get_angular_momentum() const {
    return local_inertia.component_mul(angular_velocity);
}

void RigidBody::set_linear_momentum(const Vec3& momentum) {
    if (has_finite_mass()) {
        linear_velocity = momentum / mass;
    }
}

bool RigidBody::is_moving(float threshold) const {
    return linear_velocity.magnitude() > threshold || 
           angular_velocity.magnitude() > threshold;
}

void RigidBody::print_state() const {
    printf("RigidBody #%u: %s\n", id, name.c_str());
    printf("  Position: %.3f, %.3f, %.3f\n", 
           transform.position.x, transform.position.y, transform.position.z);
    printf("  Velocity: %.3f, %.3f, %.3f\n",
           linear_velocity.x, linear_velocity.y, linear_velocity.z);
    printf("  Angular Vel: %.3f, %.3f, %.3f\n",
           angular_velocity.x, angular_velocity.y, angular_velocity.z);
    printf("  Mass: %.3f (inv: %.3f)\n", mass, inv_mass);
    printf("  Sleeping: %s\n", is_sleeping ? "Yes" : "No");
}

// ============== SphereShape Implementation ==============

bool SphereShape::raycast(const Vec3& origin, const Vec3& direction,
                          float max_distance, Vec3& hit_point,
                          Vec3& hit_normal, float& hit_fraction) const {
    Vec3 oc = origin;
    float a = direction.dot(direction);
    float b = 2.0f * oc.dot(direction);
    float c = oc.dot(oc) - radius * radius;
    float discriminant = b * b - 4.0f * a * c;

    if (discriminant < 0.0f) {
        return false;
    }

    float t1 = (-b - std::sqrt(discriminant)) / (2.0f * a);
    float t2 = (-b + std::sqrt(discriminant)) / (2.0f * a);

    if (t1 > 0.0f && t1 < max_distance) {
        hit_point = origin + direction * t1;
        hit_normal = (hit_point).normalized();
        hit_fraction = t1 / max_distance;
        return true;
    }

    if (t2 > 0.0f && t2 < max_distance) {
        hit_point = origin + direction * t2;
        hit_normal = -(hit_point).normalized();
        hit_fraction = t2 / max_distance;
        return true;
    }

    return false;
}

// ============== BoxShape Implementation ==============

Vec3 BoxShape::get_aabb_min(const Transform& transform) const {
    Vec3 corners[8] = {
        Vec3(-half_extents.x, -half_extents.y, -half_extents.z),
        Vec3(half_extents.x, -half_extents.y, -half_extents.z),
        Vec3(-half_extents.x, half_extents.y, -half_extents.z),
        Vec3(half_extents.x, half_extents.y, -half_extents.z),
        Vec3(-half_extents.x, -half_extents.y, half_extents.z),
        Vec3(half_extents.x, -half_extents.y, half_extents.z),
        Vec3(-half_extents.x, half_extents.y, half_extents.z),
        Vec3(half_extents.x, half_extents.y, half_extents.z)
    };

    Vec3 min_corner = transform.transform_point(corners[0]);
    Vec3 max_corner = min_corner;

    for (int i = 1; i < 8; ++i) {
        Vec3 corner = transform.transform_point(corners[i]);
        min_corner = min_corner.min(corner);
        max_corner = max_corner.max(corner);
    }

    return min_corner;
}

Vec3 BoxShape::get_aabb_max(const Transform& transform) const {
    Vec3 corners[8] = {
        Vec3(-half_extents.x, -half_extents.y, -half_extents.z),
        Vec3(half_extents.x, -half_extents.y, -half_extents.z),
        Vec3(-half_extents.x, half_extents.y, -half_extents.z),
        Vec3(half_extents.x, half_extents.y, -half_extents.z),
        Vec3(-half_extents.x, -half_extents.y, half_extents.z),
        Vec3(half_extents.x, -half_extents.y, half_extents.z),
        Vec3(-half_extents.x, half_extents.y, half_extents.z),
        Vec3(half_extents.x, half_extents.y, half_extents.z)
    };

    Vec3 min_corner = transform.transform_point(corners[0]);
    Vec3 max_corner = min_corner;

    for (int i = 1; i < 8; ++i) {
        Vec3 corner = transform.transform_point(corners[i]);
        min_corner = min_corner.min(corner);
        max_corner = max_corner.max(corner);
    }

    return max_corner;
}

bool BoxShape::raycast(const Vec3& origin, const Vec3& direction,
                        float max_distance, Vec3& hit_point,
                        Vec3& hit_normal, float& hit_fraction) const {
    Vec3 box_min = -half_extents;
    Vec3 box_max = half_extents;

    float tmin = 0.0f;
    float tmax = max_distance;
    hit_normal = Vec3::zero();

    for (int i = 0; i < 3; ++i) {
        if (std::abs(direction[i]) < 1e-6f) {
            if (origin[i] < box_min[i] || origin[i] > box_max[i]) {
                return false;
            }
        } else {
            float inv_d = 1.0f / direction[i];
            float t1 = (box_min[i] - origin[i]) * inv_d;
            float t2 = (box_max[i] - origin[i]) * inv_d;

            if (t1 > t2) std::swap(t1, t2);

            tmin = std::max(tmin, t1);
            tmax = std::min(tmax, t2);

            if (tmin > tmax) return false;
        }
    }

    if (tmin > 0.0f && tmin < max_distance) {
        hit_point = origin + direction * tmin;
        
        // Calculate hit normal
        Vec3 local_hit = hit_point.component_div(half_extents);
        float max_comp = std::max(std::abs(local_hit.x), 
                                  std::max(std::abs(local_hit.y), std::abs(local_hit.z)));
        
        if (max_comp == std::abs(local_hit.x)) {
            hit_normal = Vec3(local_hit.x > 0 ? 1.0f : -1.0f, 0, 0);
        } else if (max_comp == std::abs(local_hit.y)) {
            hit_normal = Vec3(0, local_hit.y > 0 ? 1.0f : -1.0f, 0);
        } else {
            hit_normal = Vec3(0, 0, local_hit.z > 0 ? 1.0f : -1.0f);
        }

        hit_fraction = tmin / max_distance;
        return true;
    }

    return false;
}

// ============== CapsuleShape Implementation ==============

Vec3 CapsuleShape::compute_inertia(float mass) const {
    float r = radius;
    float h = height;

    // Capsule inertia (approximation)
    float cylinder_height = std::max(0.0f, h - 2.0f * r);
    
    // Cylinder part
    float I_cylinder = (1.0f / 12.0f) * mass * (3.0f * r * r + cylinder_height * cylinder_height);
    float I_cylinder_radial = (1.0f / 4.0f) * mass * r * r + (1.0f / 12.0f) * mass * cylinder_height * cylinder_height;
    
    // Hemisphere ends
    float I_hemisphere = (2.0f / 5.0f) * mass * r * r;
    float I_hemisphere_axis = (1.0f / 10.0f) * mass * r * r;
    
    float I_axis = I_cylinder + 2.0f * I_hemisphere_axis;
    float I_radial = I_cylinder_radial + 2.0f * I_hemisphere;

    return Vec3(I_radial, I_radial, I_axis);
}

Vec3 CapsuleShape::get_aabb_min(const Transform& transform) const {
    Vec3 p = transform.position;
    float half_h = height * 0.5f;
    float r = radius;
    
    // Capsule AABB (simplified - assuming up axis)
    return Vec3(p.x - r, p.y - half_h - r, p.z - r);
}

Vec3 CapsuleShape::get_aabb_max(const Transform& transform) const {
    Vec3 p = transform.position;
    float half_h = height * 0.5f;
    float r = radius;
    
    return Vec3(p.x + r, p.y + half_h + r, p.z + r);
}

bool CapsuleShape::raycast(const Vec3& origin, const Vec3& direction,
                            float max_distance, Vec3& hit_point,
                            Vec3& hit_normal, float& hit_fraction) const {
    // Simplified capsule raycast
    // Check cylinder and spheres at ends
    
    float half_h = height * 0.5f - radius;
    
    // Transform to local space
    Vec3 local_origin = origin;
    Vec3 local_dir = direction;

    // Check spheres at ends
    SphereShape top_sphere(radius);
    top_sphere.radius = radius;  // Already set in constructor
    
    Vec3 sphere_center_top(0, half_h, 0);
    Vec3 sphere_center_bottom(0, -half_h, 0);

    Vec3 oc_top = local_origin - sphere_center_top;
    Vec3 oc_bottom = local_origin - sphere_center_bottom;

    for (const Vec3& sphere_center : {sphere_center_top, sphere_center_bottom}) {
        Vec3 oc = local_origin - sphere_center;
        float a = local_dir.dot(local_dir);
        float b = 2.0f * oc.dot(local_dir);
        float c = oc.dot(oc) - radius * radius;
        float discriminant = b * b - 4.0f * a * c;

        if (discriminant >= 0.0f) {
            float t = (-b - std::sqrt(discriminant)) / (2.0f * a);
            if (t > 0.0f && t < max_distance) {
                hit_point = local_origin + local_dir * t;
                hit_normal = (hit_point - sphere_center).normalized();
                hit_fraction = t / max_distance;
                return true;
            }
        }
    }

    // Check cylinder (simplified)
    Vec3 dxy = Vec3(local_dir.x, 0, local_dir.z);
    Vec3 oxy = Vec3(local_origin.x, 0, local_origin.z);
    
    float a = dxy.dot(dxy);
    if (a > 1e-6f) {
        float b = 2.0f * oxy.dot(dxy);
        float c = oxy.dot(oxy) - radius * radius;
        float discriminant = b * b - 4.0f * a * c;

        if (discriminant >= 0.0f) {
            float t = (-b - std::sqrt(discriminant)) / (2.0f * a);
            Vec3 hit = local_origin + local_dir * t;
            if (t > 0.0f && t < max_distance && hit.y >= -half_h && hit.y <= half_h) {
                hit_point = hit;
                hit_normal = Vec3(hit.x, 0, hit.z).normalized();
                hit_fraction = t / max_distance;
                return true;
            }
        }
    }

    return false;
}

// ============== PlaneShape Implementation ==============

Vec3 PlaneShape::get_aabb_min(const Transform& transform) const {
    // Infinite plane - return very large values in negative direction
    return Vec3(-1e6f, -1e6f, -1e6f);
}

Vec3 PlaneShape::get_aabb_max(const Transform& transform) const {
    return Vec3(1e6f, 1e6f, 1e6f);
}

bool PlaneShape::raycast(const Vec3& origin, const Vec3& direction,
                          float max_distance, Vec3& hit_point,
                          Vec3& hit_normal, float& hit_fraction) const {
    float denom = direction.dot(normal);
    if (std::abs(denom) < 1e-6f) return false;

    float t = (offset - origin.dot(normal)) / denom;
    if (t >= 0.0f && t <= max_distance) {
        hit_point = origin + direction * t;
        hit_normal = denom > 0.0f ? -normal : normal;
        hit_fraction = t / max_distance;
        return true;
    }

    return false;
}

} // namespace physics
} // namespace atlas
