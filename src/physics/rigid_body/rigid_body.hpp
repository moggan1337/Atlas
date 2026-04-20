/**
 * @file rigid_body.hpp
 * @brief Rigid Body Dynamics using Newton-Euler Equations
 * 
 * Implements full rigid body physics simulation:
 * - Mass and inertia tensor management
 * - Linear and angular velocity
 * - Force and torque accumulation
 * - Impulse-based dynamics
 * - Sleep optimization for performance
 */

#pragma once

#include "../../core/vec3.hpp"
#include "../../core/quaternion.hpp"
#include "../../core/transform.hpp"
#include <vector>
#include <memory>

namespace atlas {
namespace physics {

using core::Vec3;
using core::Quaternion;
using core::Transform;

// Shape types for collision
enum class ShapeType {
    Sphere,
    Box,
    Capsule,
    ConvexHull,
    Plane,
    Compound
};

// Forward declaration
class CollisionShape;

/**
 * @brief Material properties for rigid bodies
 */
struct RigidBodyMaterial {
    float density = 1000.0f;           // kg/m³
    float restitution = 0.3f;          // Bounciness (0-1)
    float friction = 0.5f;             // Coulomb friction coefficient
    float linear_damping = 0.01f;      // Velocity damping
    float angular_damping = 0.01f;     // Angular velocity damping
    float sleep_threshold = 0.01f;     // Velocity threshold for sleeping

    static RigidBodyMaterial default_material() {
        return RigidBodyMaterial();
    }

    static RigidBodyMaterial metal() {
        RigidBodyMaterial m;
        m.density = 7000.0f;
        m.restitution = 0.2f;
        m.friction = 0.6f;
        return m;
    }

    static RigidBodyMaterial rubber() {
        RigidBodyMaterial m;
        m.density = 1100.0f;
        m.restitution = 0.8f;
        m.friction = 0.9f;
        return m;
    }

    static RigidBodyMaterial ice() {
        RigidBodyMaterial m;
        m.density = 917.0f;
        m.restitution = 0.1f;
        m.friction = 0.02f;
        return m;
    }
};

/**
 * @brief Rigid Body State
 */
enum class BodyState {
    Active,
    Sleeping,
    Kinematic,
    Static
};

/**
 * @brief Rigid Body Class
 * 
 * Represents a rigid body with mass, inertia, position, orientation,
 * linear velocity, and angular velocity. Uses Newton-Euler equations.
 */
class RigidBody {
public:
    uint32_t id;
    std::string name;

    // Shape
    std::shared_ptr<CollisionShape> shape;

    // Material
    RigidBodyMaterial material;

    // State
    BodyState state = BodyState::Active;

    // Mass properties
    float mass = 1.0f;
    float inv_mass = 1.0f;
    Vec3 local_inertia;              // Principal moments of inertia
    Vec3 inv_inertia;                // Inverse principal moments

    // Transform
    Transform transform;

    // Velocity
    Vec3 linear_velocity;
    Vec3 angular_velocity;

    // Accumulated forces and torques (cleared each frame)
    Vec3 force_accumulator;
    Vec3 torque_accumulator;

    // Damping
    float linear_damping = 0.01f;
    float angular_damping = 0.01f;

    // Gravity scale (useful for buoyant objects)
    float gravity_scale = 1.0f;

    // Sleep management
    float sleep_timer = 0.0f;
    float sleep_threshold = 0.5f;
    bool allow_sleep = true;
    bool is_sleeping = false;

    // Collision filtering
    uint32_t collision_group = 0;
    uint32_t collision_mask = 0xFFFFFFFF;

    // Contact/collision callbacks
    std::function<void(const RigidBody& other, const Vec3& point, const Vec3& normal)> on_collision;
    std::function<void(const Vec3& point, const Vec3& normal, float depth)> on_trigger;

public:
    RigidBody();
    RigidBody(float mass, std::shared_ptr<CollisionShape> shape);
    ~RigidBody() = default;

    // Factory methods
    static std::shared_ptr<RigidBody> create_sphere(float radius, float mass = 1.0f);
    static std::shared_ptr<RigidBody> create_box(const Vec3& half_extents, float mass = 1.0f);
    static std::shared_ptr<RigidBody> create_capsule(float radius, float height, float mass = 1.0f);
    static std::shared_ptr<RigidBody> create_plane(const Vec3& normal, float mass = 0.0f);
    static std::shared_ptr<RigidBody> create_static_plane(const Vec3& normal, float offset);

    // Mass properties
    void set_mass(float m);
    void set_density(float density);
    void set_mass_properties(float mass, const Vec3& inertia);
    void set_local_inertia(const Vec3& inertia);
    void set_infinite_mass();
    bool has_finite_mass() const { return mass > 0.0f; }
    bool has_infinite_mass() const { return mass <= 0.0f; }

    // State management
    void activate();
    void deactivate();
    void set_kinematic();
    void set_static();
    void update_sleep_state(float dt);

    // Force and torque application
    void apply_force(const Vec3& force);
    void apply_force_at_point(const Vec3& force, const Vec3& point);
    void apply_force_at_body_point(const Vec3& force, const Vec3& local_point);
    void apply_torque(const Vec3& torque);
    void apply_impulse(const Vec3& impulse);
    void apply_impulse_at_point(const Vec3& impulse, const Vec3& point);
    void apply_impulse_at_body_point(const Vec3& impulse, const Vec3& local_point);
    void apply_angular_impulse(const Vec3& impulse);

    // Clear accumulators
    void clear_forces();
    void clear_torques();
    void clear_accumulators();

    // Gravity
    void apply_gravity(const Vec3& gravity);
    Vec3 get_gravity_force(const Vec3& gravity) const;

    // Integration
    void integrate(float dt);
    void integrate_velocities(float dt);
    void integrate_positions(float dt);

    // World space queries
    Vec3 get_world_center_of_mass() const;
    Vec3 get_world_point_velocity(const Vec3& local_point) const;
    Vec3 get_local_point_velocity(const Vec3& world_point) const;

    // Velocity constraints
    void set_linear_velocity(const Vec3& velocity);
    void set_angular_velocity(const Vec3& velocity);
    void add_linear_velocity(const Vec3& velocity);
    void add_angular_velocity(const Vec3& velocity);

    // Damping
    void apply_damping(float dt);

    // Transform
    void set_position(const Vec3& position);
    Vec3 get_position() const { return transform.position; }
    void set_rotation(const Quaternion& rotation);
    Quaternion get_rotation() const { return transform.rotation; }
    void set_transform(const Transform& t);
    Transform get_transform() const { return transform; }

    // Orientation
    Vec3 get_forward() const { return transform.get_forward(); }
    Vec3 get_up() const { return transform.get_up(); }
    Vec3 get_right() const { return transform.get_right(); }

    // AABB (Axis-Aligned Bounding Box)
    Vec3 get_aabb_min() const;
    Vec3 get_aabb_max() const;

    // Kinetic energy
    float get_kinetic_energy() const;
    float get_linear_kinetic_energy() const;
    float get_angular_kinetic_energy() const;

    // Momentum
    Vec3 get_linear_momentum() const { return linear_velocity * mass; }
    Vec3 get_angular_momentum() const;
    void set_linear_momentum(const Vec3& momentum);

    // Utility
    bool is_moving(float threshold = 0.01f) const;
    void print_state() const;
};

/**
 * @brief Collision Shape Base Class
 */
class CollisionShape {
public:
    virtual ~CollisionShape() = default;

    virtual ShapeType get_type() const = 0;
    virtual Vec3 get_center_of_mass() const = 0;
    virtual float get_volume() const = 0;
    virtual Vec3 compute_inertia(float mass) const = 0;
    virtual Vec3 get_aabb_min(const Transform& transform) const = 0;
    virtual Vec3 get_aabb_max(const Transform& transform) const = 0;
    virtual bool raycast(const Vec3& origin, const Vec3& direction, 
                        float max_distance, Vec3& hit_point, 
                        Vec3& hit_normal, float& hit_fraction) const = 0;
    virtual CollisionShape* clone() const = 0;
};

/**
 * @brief Sphere Collision Shape
 */
class SphereShape : public CollisionShape {
public:
    float radius;

    explicit SphereShape(float r) : radius(r) {}

    ShapeType get_type() const override { return ShapeType::Sphere; }
    
    Vec3 get_center_of_mass() const override {
        return Vec3::zero();
    }

    float get_volume() const override {
        return (4.0f / 3.0f) * M_PI * radius * radius * radius;
    }

    Vec3 compute_inertia(float mass) const override {
        float i = 0.4f * mass * radius * radius;
        return Vec3(i, i, i);
    }

    Vec3 get_aabb_min(const Transform& transform) const override {
        return transform.position - Vec3(radius);
    }

    Vec3 get_aabb_max(const Transform& transform) const override {
        return transform.position + Vec3(radius);
    }

    bool raycast(const Vec3& origin, const Vec3& direction, 
                 float max_distance, Vec3& hit_point,
                 Vec3& hit_normal, float& hit_fraction) const override;

    CollisionShape* clone() const override { return new SphereShape(radius); }
};

/**
 * @brief Box Collision Shape
 */
class BoxShape : public CollisionShape {
public:
    Vec3 half_extents;

    explicit BoxShape(const Vec3& half) : half_extents(half) {}
    BoxShape(float x, float y, float z) : half_extents(x, y, z) {}

    ShapeType get_type() const override { return ShapeType::Box; }

    Vec3 get_center_of_mass() const override { return Vec3::zero(); }

    float get_volume() const override {
        return 8.0f * half_extents.x * half_extents.y * half_extents.z;
    }

    Vec3 compute_inertia(float mass) const override {
        return Vec3(
            (1.0f / 12.0f) * mass * (half_extents.y * half_extents.y * 4 + half_extents.z * half_extents.z * 4),
            (1.0f / 12.0f) * mass * (half_extents.x * half_extents.x * 4 + half_extents.z * half_extents.z * 4),
            (1.0f / 12.0f) * mass * (half_extents.x * half_extents.x * 4 + half_extents.y * half_extents.y * 4)
        );
    }

    Vec3 get_aabb_min(const Transform& transform) const override;
    Vec3 get_aabb_max(const Transform& transform) const override;

    bool raycast(const Vec3& origin, const Vec3& direction,
                 float max_distance, Vec3& hit_point,
                 Vec3& hit_normal, float& hit_fraction) const override;

    CollisionShape* clone() const override { return new BoxShape(half_extents); }
};

/**
 * @brief Capsule Collision Shape
 */
class CapsuleShape : public CollisionShape {
public:
    float radius;
    float height;

    CapsuleShape(float r, float h) : radius(r), height(h) {}

    ShapeType get_type() const override { return ShapeType::Capsule; }

    Vec3 get_center_of_mass() const override { return Vec3::zero(); }

    float get_volume() const override {
        float cylinder_vol = M_PI * radius * radius * (height - 2.0f * radius);
        float sphere_vol = (4.0f / 3.0f) * M_PI * radius * radius * radius;
        return cylinder_vol + sphere_vol;
    }

    Vec3 compute_inertia(float mass) const override;

    Vec3 get_aabb_min(const Transform& transform) const override;
    Vec3 get_aabb_max(const Transform& transform) const override;

    bool raycast(const Vec3& origin, const Vec3& direction,
                 float max_distance, Vec3& hit_point,
                 Vec3& hit_normal, float& hit_fraction) const override;

    CollisionShape* clone() const override { return new CapsuleShape(radius, height); }
};

/**
 * @brief Plane Collision Shape (infinite)
 */
class PlaneShape : public CollisionShape {
public:
    Vec3 normal;
    float offset;  // Distance from origin

    PlaneShape(const Vec3& n, float o) : normal(n.normalized()), offset(o) {}

    ShapeType get_type() const override { return ShapeType::Plane; }

    Vec3 get_center_of_mass() const override { return Vec3::zero(); }

    float get_volume() const override { return INFINITY; }

    Vec3 compute_inertia(float mass) const override { return Vec3(INFINITY); }

    Vec3 get_aabb_min(const Transform& transform) const override;
    Vec3 get_aabb_max(const Transform& transform) const override;

    bool raycast(const Vec3& origin, const Vec3& direction,
                 float max_distance, Vec3& hit_point,
                 Vec3& hit_normal, float& hit_fraction) const override;

    CollisionShape* clone() const override { return new PlaneShape(normal, offset); }
};

} // namespace physics
} // namespace atlas
