/**
 * @file terrain.hpp
 * @brief Procedural Terrain Generation
 * 
 * Implements various terrain generation algorithms:
 * - Perlin noise
 * - Simplex noise
 * - Diamond-square
 * - Erosion simulation
 * - Biome generation
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
 * @brief Terrain Cell/Vertex
 */
struct TerrainVertex {
    Vec3 position;
    Vec3 normal;
    Vec2 uv;
    float height;
    float moisture;
    float temperature;
    int biome_type;
    
    TerrainVertex() 
        : position(Vec3::zero())
        , normal(Vec3::unit_y())
        , uv(0, 0)
        , height(0.0f)
        , moisture(0.0f)
        , temperature(0.5f)
        , biome_type(0) {}
};

/**
 * @brief Terrain Configuration
 */
struct TerrainConfig {
    int width = 256;
    int depth = 256;
    float world_size = 100.0f;
    float max_height = 30.0f;
    
    // Noise settings
    int octaves = 6;
    float frequency = 0.01f;
    float lacunarity = 2.0f;
    float persistence = 0.5f;
    int seed = 12345;
    
    // Erosion
    bool enable_thermal_erosion = true;
    bool enable_hydraulic_erosion = true;
    int erosion_iterations = 50;
    float erosion_strength = 0.01f;
    
    // Biome
    bool generate_biomes = true;
    float sea_level = 0.0f;
};

/**
 * @brief Biome Types
 */
enum class BiomeType {
    DeepOcean,
    Ocean,
    Beach,
    Desert,
    Plains,
    Forest,
    Hills,
    Mountains,
    Snow,
    Swamp,
    Taiga
};

/**
 * @brief Terrain Generator
 */
class TerrainGenerator {
public:
    TerrainConfig config;
    
private:
    std::vector<TerrainVertex> vertices;
    std::vector<uint32_t> indices;
    std::vector<float> heightmap;
    
    std::mt19937 rng;
    
    // Cached noise permutation
    std::vector<int> perm_x;
    std::vector<int> perm_y;
    
public:
    TerrainGenerator(const TerrainConfig& cfg = TerrainConfig());
    ~TerrainGenerator() = default;

    // Generation
    void generate();
    void generate_heightmap();
    void generate_normals();
    void generate_biomes();
    
    // Noise functions
    float perlin_noise(float x, float y);
    float perlin_noise_2d(float x, float y);
    float simplex_noise(float x, float y);
    float value_noise(float x, float y);
    float fbm(float x, float y, int octaves, float persistence, float lacunarity);
    
    // Erosion
    void apply_thermal_erosion(float dt);
    void apply_hydraulic_erosion(float iterations);
    void apply_hydrophobic_erosion();
    
    // Terrain operations
    void smooth(float radius);
    void terrace(int levels);
    void erode_channel(const Vec3& start, const Vec3& end, float depth);
    
    // Height queries
    float get_height_at(float x, float z) const;
    float get_height_at(int x, int z) const;
    Vec3 get_normal_at(float x, float z) const;
    float get_biome_temperature(float x, float z) const;
    float get_biome_moisture(float x, float z) const;
    BiomeType get_biome_at(float x, float z) const;
    
    // Mesh generation
    void build_mesh();
    void update_mesh();
    
    // Data access
    const std::vector<TerrainVertex>& get_vertices() const { return vertices; }
    const std::vector<uint32_t>& get_indices() const { return indices; }
    const std::vector<float>& get_heightmap() const { return heightmap; }
    
    int get_width() const { return config.width; }
    int get_depth() const { return config.depth; }
    float get_world_size() const { return config.world_size; }
    
    // Coordinate conversion
    Vec2 world_to_grid(float x, float z) const;
    Vec2 grid_to_world(int x, int z) const;
    
    // Modification
    void modify_height(int x, int z, float delta);
    void set_height(int x, int z, float height);
    
private:
    void initialize_noise();
    int heightmap_index(int x, int z) const { return z * config.width + x; }
    bool in_bounds(int x, int z) const {
        return x >= 0 && x < config.width && z >= 0 && z < config.depth;
    }
    
    // Gradients for Perlin
    Vec2 gradient(int hash, float x, float y) const;
    float fade(float t) const { return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f); }
    float lerp(float a, float b, float t) const { return a + t * (b - a); }
    
    // Biome determination
    BiomeType determine_biome(float height, float moisture, float temperature) const;
    
    // Erosion helpers
    float calculate_slope(int x, int z) const;
    void apply_sediment_transport(int x, int z, float amount);
};

/**
 * @brief Terrain Chunk (for LOD)
 */
class TerrainChunk {
public:
    struct Config {
        int grid_size = 64;
        float chunk_world_size = 16.0f;
    } config;
    
private:
    int chunk_x, chunk_z;
    std::vector<TerrainVertex> vertices;
    std::vector<uint32_t> indices;
    std::shared_ptr<TerrainGenerator> generator;
    
    int lod_level;
    
public:
    TerrainChunk(int x, int z, const std::shared_ptr<TerrainGenerator>& gen, 
                 int lod = 0);
    ~TerrainChunk() = default;
    
    void generate();
    void update(const Vec3& camera_position);
    
    const std::vector<TerrainVertex>& get_vertices() const { return vertices; }
    const std::vector<uint32_t>& get_indices() const { return indices; }
    
    int get_chunk_x() const { return chunk_x; }
    int get_chunk_z() const { return chunk_z; }
    int get_lod() const { return lod_level; }
    
    Vec2 get_world_offset() const;
};

/**
 * @brief Terrain World
 */
class TerrainWorld {
public:
    struct Config {
        TerrainConfig terrain;
        int chunks_visible = 8;
        bool use_chunks = true;
    } config;
    
private:
    TerrainConfig terrain_config;
    std::shared_ptr<TerrainGenerator> generator;
    std::vector<std::shared_ptr<TerrainChunk>> chunks;
    
public:
    TerrainWorld(const Config& cfg = Config());
    ~TerrainWorld() = default;
    
    void initialize();
    void update(const Vec3& camera_position);
    
    std::shared_ptr<TerrainGenerator> get_generator() { return generator; }
    const std::shared_ptr<TerrainGenerator>& get_generator() const { return generator; }
    
    float get_height_at(float x, float z) const;
    Vec3 get_normal_at(float x, float z) const;
    BiomeType get_biome_at(float x, float z) const;
    
    const std::vector<std::shared_ptr<TerrainChunk>>& get_chunks() const { return chunks; }
};

} // namespace physics
} // namespace atlas
