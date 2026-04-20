/**
 * @file rigid_body_world.hpp
 * @brief Rigid Body Physics World
 * 
 * Manages all rigid bodies and handles:
 * - Broad phase collision detection (AABB)
 * - Narrow phase collision detection (SAT, GJK)
 * - Contact generation and resolution
 * - Impulse-based constraint solving
 * - Sleep optimization
 */

#pragma once

#include "rigid_body.hpp"
#include <vector>
#include <memory>
#include <unordered_map>
#include <functional>

namespace atlas {
namespace physics {

using core::Vec3;
using core::Quaternion;
using core::Transform;

/**
 * @brief Contact Point
 */
struct ContactPoint {
    Vec3 position;          // World space contact position
    Vec3 normal;           // Contact normal (from body A to B)
    float penetration;     // Penetration depth
    float restitution;     // Combined restitution
    float friction;        // Combined friction
    
    // Cached impulses for warm starting
    Vec3 normal_impulse;
    Vec3 tangent_impulse;
};

/**
 * @brief Contact Manifold
 */
struct ContactManifold {
    RigidBody* body_a;
    RigidBody* body_b;
    std::vector<ContactPoint> contacts;
    Vec3 normal;
    int contact_count;
    
    void add_contact(const ContactPoint& contact) {
        contacts.push_back(contact);
    }
};

/**
 * @brief Collision Pair
 */
struct CollisionPair {
    uint32_t body_a_id;
    uint32_t body_b_id;
    
    CollisionPair(uint32_t a, uint32_t b) : body_a_id(a), body_b_id(b) {}
    
    bool operator==(const CollisionPair& other) const {
        return body_a_id == other.body_a_id && body_b_id == other.body_b_id;
    }
};

/**
 * @brief Broad Phase Result
 */
struct BroadPhasePair {
    RigidBody* body_a;
    RigidBody* body_b;
};

/**
 * @brief Physics Solver Configuration
 */
struct SolverConfig {
    int max_iterations = 10;
    float velocity_iterations = 6.0f;
    float position_iterations = 2.0f;
    float baumgarte_scale = 0.2f;     // Position correction factor
    float slop = 0.005f;              // Allowed penetration tolerance
    float max_velocity = 50.0f;
    float max_angular_velocity = 100.0f;
};

/**
 * @brief Physics World Configuration
 */
struct WorldConfig {
    Vec3 gravity = Vec3(0.0f, -9.81f, 0.0f);
    int substeps = 4;
    float max_dt = 1.0f / 60.0f;
    bool allow_sleep = true;
    SolverConfig solver;
};

/**
 * @brief Rigid Body World
 */
class RigidBodyWorld {
public:
    WorldConfig config;

private:
    std::vector<std::shared_ptr<RigidBody>> bodies;
    std::unordered_map<uint32_t, size_t> body_id_map;
    std::vector<ContactManifold> manifolds;
    std::vector<BroadPhasePair> broadphase_pairs;
    
    uint32_t next_body_id = 0;
    
    // Collision detection
    std::function<void(RigidBody&, RigidBody&, const ContactManifold&)> collision_callback;

public:
    RigidBodyWorld(const WorldConfig& cfg = WorldConfig());
    ~RigidBodyWorld() = default;

    // Body management
    std::shared_ptr<RigidBody> add_body(const std::shared_ptr<RigidBody>& body);
    void remove_body(uint32_t id);
    void remove_body(const std::shared_ptr<RigidBody>& body);
    RigidBody* get_body(uint32_t id);
    const std::vector<std::shared_ptr<RigidBody>>& get_bodies() const { return bodies; }
    size_t get_body_count() const { return bodies.size(); }
    void clear();

    // Simulation
    void step(float dt);
    void step_fixed(float dt);

    // Gravity
    void set_gravity(const Vec3& gravity);
    Vec3 get_gravity() const { return config.gravity; }

    // Collision callback
    void set_collision_callback(std::function<void(RigidBody&, RigidBody&, const ContactManifold&)> callback);

    // Queries
    void query_point(const Vec3& point, std::vector<RigidBody*>& results);
    void query_aabb(const Vec3& min, const Vec3& max, std::vector<RigidBody*>& results);
    void query_ray(const Vec3& origin, const Vec3& direction, float max_distance,
                   std::vector<RaycastHit>& hits);

    // Raycast result
    struct RaycastHit {
        RigidBody* body;
        Vec3 point;
        Vec3 normal;
        float fraction;
    };

    // Debug visualization data
    std::vector<std::pair<Vec3, Vec3>> get_debug_lines();
    std::vector<Vec3> get_debug_points();

private:
    // Integration
    void integrate_bodies(float dt);
    
    // Broad phase
    void broadphase();
    bool aabb_intersect(const RigidBody* a, const RigidBody* b) const;
    
    // Narrow phase
    void narrowphase();
    bool generate_contacts(RigidBody* a, RigidBody* b, ContactManifold& manifold);
    
    // Contact generation helpers
    bool sphere_vs_sphere(RigidBody* a, RigidBody* b, ContactManifold& manifold);
    bool box_vs_box(RigidBody* a, RigidBody* b, ContactManifold& manifold);
    bool sphere_vs_box(RigidBody* sphere, RigidBody* box, ContactManifold& manifold);
    bool sphere_vs_plane(RigidBody* sphere, RigidBody* plane, ContactManifold& manifold);
    bool box_vs_plane(RigidBody* box, RigidBody* plane, ContactManifold& manifold);
    
    // SAT collision detection
    bool sat_test(const std::vector<Vec3>& axes, RigidBody* a, RigidBody* b,
                 Vec3& normal, float& penetration);
    std::vector<Vec3> get_box_axes(RigidBody* box);
    
    // GJK collision detection
    bool gjk_collision(RigidBody* a, RigidBody* b, Vec3& simplex_out);
    Vec3 gjk_support(RigidBody* a, RigidBody* b, const Vec3& direction);
    
    // EPA (Expanding Polytope Algorithm) for penetration depth
    bool epa(CollisionShape* a, CollisionShape* b, const Vec3& init_direction,
            Vec3& normal, float& depth);
    
    // Solving
    void solve_contacts(float dt);
    void solve_velocity_constraints(float dt);
    void solve_position_constraints(float dt);
    
    // Contact resolution
    void resolve_contact(ContactManifold& manifold, float dt);
    void apply_contact_impulse(ContactManifold& manifold, ContactPoint& contact, float dt);
    
    // Warm starting
    void warm_start_contacts();
    
    // Position correction
    void apply_position_correction(ContactManifold& manifold);
    
    // Sleep management
    void update_sleep_states(float dt);
    
    // Utility
    void clamp_velocity(RigidBody* body);
};

} // namespace physics
} // namespace atlas
