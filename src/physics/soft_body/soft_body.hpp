/**
 * @file soft_body.hpp
 * @brief Soft Body Simulation using Mass-Spring System
 * 
 * Implements deformable body physics:
 * - Mass-spring systems with structural, shear, and bend springs
 * - Verlet integration for stable simulation
 * - Pressure constraints for volume preservation
 * - Shape matching for rest shape recovery
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
 * @brief Particle for soft body simulation
 */
struct SoftBodyParticle {
    Vec3 position;           // Current position
    Vec3 previous_position; // Previous position (for Verlet)
    Vec3 velocity;           // Velocity
    Vec3 acceleration;      // Accumulated acceleration
    Vec3 normal;             // Surface normal
    float mass;              // Particle mass
    float inv_mass;         // Inverse mass (0 for pinned)
    float pressure;         // Pressure contribution
    bool pinned;             // Is particle fixed in space
    
    SoftBodyParticle() 
        : position(Vec3::zero())
        , previous_position(Vec3::zero())
        , velocity(Vec3::zero())
        , acceleration(Vec3::zero())
        , normal(Vec3::unit_y())
        , mass(1.0f)
        , inv_mass(1.0f)
        , pressure(0.0f)
        , pinned(false) {}
    
    SoftBodyParticle(const Vec3& pos, float m = 1.0f) 
        : position(pos)
        , previous_position(pos)
        , velocity(Vec3::zero())
        , acceleration(Vec3::zero())
        , normal(Vec3::unit_y())
        , mass(m)
        , inv_mass(m > 0.0f ? 1.0f / m : 0.0f)
        , pressure(0.0f)
        , pinned(false) {}
};

/**
 * @brief Spring connection between particles
 */
struct Spring {
    int particle_a;      // First particle index
    int particle_b;      // Second particle index
    float rest_length;   // Natural length
    float stiffness;     // Spring stiffness (k)
    float damping;        // Damping coefficient
    
    Spring(int a, int b, float rest, float k = 100.0f, float d = 1.0f)
        : particle_a(a), particle_b(b), rest_length(rest), 
          stiffness(k), damping(d) {}
};

/**
 * @brief Face/Triangle for mesh-based soft bodies
 */
struct SoftBodyFace {
    int indices[3];      // Particle indices
    Vec3 normal;         // Face normal
    float area;          // Face area
    
    SoftBodyFace(int i0, int i1, int i2) {
        indices[0] = i0;
        indices[1] = i1;
        indices[2] = i2;
        normal = Vec3::unit_y();
        area = 0.0f;
    }
};

/**
 * @brief Soft Body Configuration
 */
struct SoftBodyConfig {
    float mass = 1.0f;
    float stiffness = 500.0f;          // Structural spring stiffness
    float shear_stiffness = 400.0f;    // Shear spring stiffness
    float bend_stiffness = 100.0f;     // Bend spring stiffness
    float damping = 5.0f;              // Global damping
    float pressure = 100.0f;          // Internal pressure
    float gravity_scale = 1.0f;
    bool shape_matching = true;
    float shape_stiffness = 0.5f;
    bool self_collision = false;
    float tear_threshold = INFINITY;   // Distance at which spring breaks
};

/**
 * @brief Soft Body Class
 */
class SoftBody {
public:
    std::string name;
    SoftBodyConfig config;
    
    std::vector<SoftBodyParticle> particles;
    std::vector<Spring> springs;
    std::vector<SoftBodyFace> faces;
    
    // Rest shape (for shape matching)
    std::vector<Vec3> rest_positions;
    Transform rest_transform;
    
    // Bounding box
    Vec3 aabb_min;
    Vec3 aabb_max;
    
public:
    SoftBody() = default;
    SoftBody(const SoftBodyConfig& cfg) : config(cfg) {}
    ~SoftBody() = default;

    // Particle management
    int add_particle(const Vec3& position, float mass = 1.0f);
    void pin_particle(int index);
    void unpin_particle(int index);
    SoftBodyParticle& get_particle(int index) { return particles[index]; }
    const SoftBodyParticle& get_particle(int index) const { return particles[index]; }
    
    // Spring management
    int add_spring(int particle_a, int particle_b, float stiffness = -1.0f, float damping = -1.0f);
    void remove_spring(int index);
    
    // Face management
    int add_face(int i0, int i1, int i2);
    
    // Shape creation utilities
    static std::shared_ptr<SoftBody> create_ball(const Vec3& center, float radius, int segments = 12);
    static std::shared_ptr<SoftBody> create_box(const Vec3& center, const Vec3& size);
    static std::shared_ptr<SoftBody> create_cloth(const Vec3& origin, float width, float height, 
                                                  int segments_x, int segments_y);
    static std::shared_ptr<SoftBody> create_chain(const Vec3& start, const Vec3& end, int segments);
    static std::shared_ptr<SoftBody> create_rope(const Vec3& start, float length, int segments);
    
    // Simulation
    void update(float dt, const Vec3& gravity);
    void integrate_verlet(float dt);
    void solve_springs(float dt);
    void solve_constraints(float dt);
    void apply_pressure();
    void apply_damping(float dt);
    void enforce_shape_matching(float dt);
    
    // Collision
    void collide_with_plane(const Vec3& normal, float offset, float restitution = 0.3f);
    void collide_with_sphere(const Vec3& center, float radius, float restitution = 0.3f);
    void collide_with_box(const Vec3& min, const Vec3& max, float restitution = 0.3f);
    
    // Utilities
    void update_aabb();
    Vec3 get_center_of_mass() const;
    Vec3 get_velocity_at_point(const Vec3& point) const;
    void apply_force(const Vec3& force);
    void apply_force_at_point(const Vec3& force, const Vec3& point);
    void apply_impulse(const Vec3& impulse);
    void apply_impulse_at_point(const Vec3& impulse, const Vec3& point);
    
    // Normal calculation
    void update_normals();
    
    // Volume calculation (for pressure)
    float calculate_volume() const;
    
    // Tear springs that are stretched too far
    void check_tearing();
    
    // State
    size_t get_particle_count() const { return particles.size(); }
    size_t get_spring_count() const { return springs.size(); }
    size_t get_face_count() const { return faces.size(); }
    bool is_valid() const { return !particles.empty(); }
};

/**
 * @brief Soft Body World
 */
class SoftBodyWorld {
public:
    struct Config {
        Vec3 gravity = Vec3(0.0f, -9.81f, 0.0f);
        int substeps = 4;
        float max_dt = 1.0f / 60.0f;
        int solver_iterations = 3;
        bool enable_pressure = true;
        bool enable_shape_matching = true;
        bool enable_self_collision = false;
    } config;

private:
    std::vector<std::shared_ptr<SoftBody>> bodies;
    
public:
    SoftBodyWorld(const Config& cfg = Config()) : config(cfg) {}
    ~SoftBodyWorld() = default;
    
    std::shared_ptr<SoftBody> add_body(const std::shared_ptr<SoftBody>& body) {
        bodies.push_back(body);
        return body;
    }
    
    void remove_body(const std::shared_ptr<SoftBody>& body) {
        bodies.erase(std::remove(bodies.begin(), bodies.end(), body), bodies.end());
    }
    
    void remove_body(const std::string& name) {
        bodies.erase(
            std::remove_if(bodies.begin(), bodies.end(), 
                [&name](const std::shared_ptr<SoftBody>& b) { return b->name == name; }),
            bodies.end()
        );
    }
    
    SoftBody* get_body(const std::string& name) {
        for (auto& body : bodies) {
            if (body->name == name) return body.get();
        }
        return nullptr;
    }
    
    const std::vector<std::shared_ptr<SoftBody>>& get_bodies() const { return bodies; }
    size_t get_body_count() const { return bodies.size(); }
    
    void clear() { bodies.clear(); }
    
    void step(float dt);
    void step_fixed(float dt);
    
    void set_gravity(const Vec3& gravity) { config.gravity = gravity; }
    Vec3 get_gravity() const { return config.gravity; }
};

/**
 * @brief Soft Body Particle System
 * For large-scale particle-based soft bodies (jelly, blobs)
 */
class ParticleSystem {
public:
    struct Config {
        float rest_density = 1000.0f;
        float stiffness = 50.0f;
        float viscosity = 0.1f;
        float gravity_scale = 1.0f;
    } config;
    
private:
    std::vector<SoftBodyParticle> particles;
    float kernel_radius;
    
public:
    ParticleSystem(float radius = 0.5f) : kernel_radius(radius) {}
    
    void add_particle(const Vec3& position, float mass = 1.0f);
    void step(float dt, const Vec3& gravity);
    
    const std::vector<SoftBodyParticle>& get_particles() const { return particles; }
    size_t get_particle_count() const { return particles.size(); }
};

} // namespace physics
} // namespace atlas
