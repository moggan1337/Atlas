/**
 * @file hair.hpp
 * @brief Hair/Fur Simulation using Super-Segment Method
 * 
 * Implements realistic hair simulation:
 * - Multiple segments per strand (physics chain)
 * - Constraints for hair behavior
 * - Collision with body
 * - Gravity and wind forces
 */

#pragma once

#include "../../core/vec3.hpp"
#include <vector>
#include <memory>
#include <random>

namespace atlas {
namespace physics {

using core::Vec3;

/**
 * @brief Hair Segment
 */
struct HairSegment {
    Vec3 position;
    Vec3 previous_position;
    Vec3 velocity;
    Vec3 local_offset;  // Rest position relative to root
    float mass;
    float inv_mass;
    
    HairSegment() 
        : position(Vec3::zero())
        , previous_position(Vec3::zero())
        , velocity(Vec3::zero())
        , local_offset(Vec3::zero())
        , mass(0.001f)
        , inv_mass(1000.0f) {}
};

/**
 * @brief Hair Strand
 */
struct HairStrand {
    int root_index;           // Index of root particle
    std::vector<HairSegment> segments;
    Vec3 root_position;       // Base position in world space
    Vec3 root_normal;         // Surface normal at root
    float length;             // Total strand length
    float thickness;          // Strand thickness
    int segment_count;
    float curl;               // Curl factor (0-1)
    Vec3 color;              // Hair color
    bool visible;
    
    HairStrand() 
        : root_index(-1)
        , length(0.1f)
        , thickness(0.002f)
        , segment_count(5)
        , curl(0.0f)
        , color(0.1f, 0.05f, 0.02f)
        , visible(true) {}
};

/**
 * @brief Hair Follicle
 */
struct HairFollicle {
    Vec3 position;
    Vec3 normal;
    float density;
    float length;
    float thickness;
    int strand_index;
    bool active;
    
    HairFollicle() 
        : position(Vec3::zero())
        , normal(Vec3::unit_y())
        , density(1.0f)
        , length(0.1f)
        , thickness(0.002f)
        , strand_index(-1)
        , active(true) {}
};

/**
 * @brief Hair Configuration
 */
struct HairConfig {
    int segments_per_strand = 5;
    float segment_mass = 0.001f;
    float stiffness = 1000.0f;
    float damping = 0.1f;
    float gravity_scale = 0.8f;
    float curl = 0.2f;
    float thickness = 0.002f;
    float length_min = 0.05f;
    float length_max = 0.15f;
    float density = 50000.0f;  // strands per square meter
    
    // Forces
    float wind_strength = 0.0f;
    Vec3 wind_direction = Vec3(1.0f, 0.0f, 0.0f);
    float wind_turbulence = 0.0f;
    
    // Collision
    float collision_radius = 0.01f;
    bool self_collision = false;
    
    // Rendering
    Vec3 base_color = Vec3(0.1f, 0.05f, 0.02f);
    float color_variation = 0.1f;
};

/**
 * @brief Hair Simulation
 */
class HairSimulation {
public:
    HairConfig config;
    
private:
    std::vector<HairStrand> strands;
    std::vector<HairFollicle> follicles;
    std::vector<HairSegment> all_segments;
    
    // For fast access
    std::vector<int> strand_offsets;
    
    // Random generator
    std::mt19937 rng;
    std::uniform_real_distribution<float> dist;
    
    // Bounds
    Vec3 bounds_min;
    Vec3 bounds_max;
    
public:
    HairSimulation(const HairConfig& cfg = HairConfig());
    ~HairSimulation() = default;

    // Strand creation
    int add_strand(const Vec3& root, const Vec3& normal, float length = -1.0f);
    void add_strands_on_surface(const std::vector<Vec3>& positions, 
                                const std::vector<Vec3>& normals,
                                float density = -1.0f);
    void generate_follicles(const Vec3& center, const Vec3& size, float density);
    
    // Simulation
    void update(float dt, const Vec3& gravity);
    void integrate_verlet(float dt);
    void solve_constraints();
    void apply_forces(float dt, const Vec3& gravity);
    
    // Forces
    void apply_wind(float time, float dt);
    void apply_gravity(const Vec3& gravity);
    
    // Collisions
    void collide_with_sphere(const Vec3& center, float radius);
    void collide_with_body(const std::vector<Vec3>& vertices,
                           const std::vector<uint32_t>& indices);
    
    // Utilities
    void update_bounds();
    
    // Access
    size_t get_strand_count() const { return strands.size(); }
    size_t get_segment_count() const { return all_segments.size(); }
    size_t get_follicle_count() const { return follicles.size(); }
    
    const HairStrand& get_strand(int i) const { return strands[i]; }
    const HairSegment& get_segment(int i) const { return all_segments[i]; }
    
    // Mesh generation for rendering
    std::vector<Vec3> get_render_positions() const;
    std::vector<Vec3> get_render_normals() const;
    std::vector<float> get_render_thicknesses() const;
    std::vector<Vec3> get_render_colors() const;
    std::vector<uint32_t> get_render_indices() const;
    
    // Grooming
    void cut_hair(float length);
    void style_hair(const Vec3& direction);
    void apply_curl(float amount);
};

/**
 * @brief Fur Rendering Helper
 * Creates shell orfin shell rendering data for fur
 */
class FurRenderer {
public:
    struct ShellLayer {
        float height;  // 0-1, how far along strand
        float scale;   // Expansion factor
    };
    
    std::vector<ShellLayer> layers;
    int layers_count = 16;
    float fur_length = 0.02f;
    
    FurRenderer() {
        for (int i = 0; i < layers_count; ++i) {
            ShellLayer layer;
            layer.height = static_cast<float>(i) / layers_count;
            layer.scale = 1.0f + layer.height * 2.0f;  // Expand outward
            layers.push_back(layer);
        }
    }
    
    // Generate shell mesh positions
    std::vector<Vec3> generate_shell(const HairSimulation& hair, int layer_index);
    
    // Generate fin mesh (cross sections)
    std::vector<Vec3> generate_fins(const HairSimulation& hair);
};

} // namespace physics
} // namespace atlas
