/**
 * @file cloth.hpp
 * @brief Cloth Simulation using Verlet Integration
 * 
 * Implements realistic cloth simulation:
 * - Particle-based cloth with structural, shear, and bend constraints
 * - Verlet integration for stability
 * - Self-collision detection
 * - Wind and external forces
 */

#pragma once

#include "../../core/vec3.hpp"
#include <vector>
#include <memory>

namespace atlas {
namespace physics {

using core::Vec3;

/**
 * @brief Cloth Constraint Types
 */
enum class ClothConstraintType {
    Structural,  // Direct neighbors (stretch)
    Shear,       // Diagonal neighbors
    Bend         // Skip-one neighbors (bending stiffness)
};

/**
 * @brief Cloth Constraint
 */
struct ClothConstraint {
    int particle_a;
    int particle_b;
    float rest_length;
    float stiffness;
    float compliance;       // 1/stiffness
    ClothConstraintType type;
    
    ClothConstraint(int a, int b, float rest, float stiff, ClothConstraintType t)
        : particle_a(a), particle_b(b), rest_length(rest), stiffness(stiff),
          compliance(1.0f / stiff), type(t) {}
};

/**
 * @brief Cloth Particle
 */
struct ClothParticle {
    Vec3 position;
    Vec3 previous_position;
    Vec3 acceleration;
    Vec3 normal;
    Vec2 uv;              // Texture coordinates
    float mass;
    float inv_mass;
    bool pinned;
    
    ClothParticle() 
        : position(Vec3::zero())
        , previous_position(Vec3::zero())
        , acceleration(Vec3::zero())
        , normal(Vec3::unit_y())
        , uv(0, 0)
        , mass(1.0f)
        , inv_mass(1.0f)
        , pinned(false) {}
};

/**
 * @brief Cloth Triangle
 */
struct ClothTriangle {
    int indices[3];
    Vec3 normal;
    float area;
    
    ClothTriangle(int i0, int i1, int i2) {
        indices[0] = i0;
        indices[1] = i1;
        indices[2] = i2;
        normal = Vec3::unit_y();
        area = 0.0f;
    }
};

/**
 * @brief Cloth Configuration
 */
struct ClothConfig {
    float mass = 0.1f;
    float structural_stiffness = 1000.0f;
    float shear_stiffness = 500.0f;
    float bend_stiffness = 100.0f;
    float damping = 0.01f;
    float gravity_scale = 1.0f;
    float wind_strength = 0.0f;
    Vec3 wind_direction = Vec3(1.0f, 0.0f, 0.0f);
    bool self_collision = true;
    float thickness = 0.01f;
    float tear_threshold = INFINITY;
};

/**
 * @brief Cloth Simulation
 */
class Cloth {
public:
    std::string name;
    ClothConfig config;
    
    std::vector<ClothParticle> particles;
    std::vector<ClothConstraint> constraints;
    std::vector<ClothTriangle> triangles;
    
    // Bounds
    Vec3 aabb_min;
    Vec3 aabb_max;
    
public:
    Cloth() = default;
    Cloth(const ClothConfig& cfg) : config(cfg) {}
    ~Cloth() = default;

    // Creation
    int add_particle(const Vec3& position, float mass = -1.0f);
    int add_particle(const Vec3& position, const Vec2& uv, float mass = -1.0f);
    
    void add_constraint(int a, int b, ClothConstraintType type, float stiffness = -1.0f);
    void add_structural_constraint(int a, int b);
    void add_shear_constraint(int a, int b);
    void add_bend_constraint(int a, int b);
    
    int add_triangle(int i0, int i1, int i2);
    
    void pin_particle(int index);
    void unpin_particle(int index);
    
    // Factory methods
    static std::shared_ptr<Cloth> create_grid(const Vec3& origin, float width, float height,
                                               int segments_x, int segments_y,
                                               bool pin_top = true);
    static std::shared_ptr<Cloth> create_curtain(const Vec3& origin, float width, float height,
                                                  int segments_x, int segments_y,
                                                  int pin_spacing = 3);
    
    // Simulation
    void update(float dt, const Vec3& gravity);
    void integrate_verlet(float dt);
    void solve_constraints(int iterations);
    void apply_wind(float dt);
    void apply_damping(float dt);
    
    // Collisions
    void collide_with_sphere(const Vec3& center, float radius);
    void collide_with_box(const Vec3& min, const Vec3& max);
    void collide_with_plane(const Vec3& normal, float offset);
    void self_collide(float dt);
    
    // Utilities
    void update_normals();
    void update_aabb();
    void check_tearing();
    Vec3 get_center() const;
    
    // Access
    size_t get_particle_count() const { return particles.size(); }
    size_t get_constraint_count() const { return constraints.size(); }
    size_t get_triangle_count() const { return triangles.size(); }
    
    ClothParticle& get_particle(int i) { return particles[i]; }
    const ClothParticle& get_particle(int i) const { return particles[i]; }
    
    // Mesh data for rendering
    std::vector<Vec3> get_positions() const;
    std::vector<Vec3> get_normals() const;
    std::vector<Vec2> get_uvs() const;
    std::vector<uint32_t> get_indices() const;
};

/**
 * @brief Cloth World
 */
class ClothWorld {
public:
    struct Config {
        Vec3 gravity = Vec3(0.0f, -9.81f, 0.0f);
        int substeps = 4;
        int constraint_iterations = 3;
        float max_dt = 1.0f / 60.0f;
        bool enable_self_collision = true;
    } config;
    
private:
    std::vector<std::shared_ptr<Cloth>> cloths;
    
public:
    ClothWorld(const Config& cfg = Config()) : config(cfg) {}
    ~ClothWorld() = default;
    
    std::shared_ptr<Cloth> add_cloth(const std::shared_ptr<Cloth>& cloth) {
        cloths.push_back(cloth);
        return cloth;
    }
    
    void remove_cloth(const std::shared_ptr<Cloth>& cloth);
    
    const std::vector<std::shared_ptr<Cloth>>& get_cloths() const { return cloths; }
    size_t get_cloth_count() const { return cloths.size(); }
    
    void clear() { cloths.clear(); }
    
    void step(float dt);
    void step_fixed(float dt);
};

} // namespace physics
} // namespace atlas
