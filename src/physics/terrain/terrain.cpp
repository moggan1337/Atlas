/**
 * @file terrain.cpp
 * @brief Terrain Generation Implementation
 */

#include "terrain.hpp"
#include <algorithm>
#include <cmath>

namespace atlas {
namespace physics {

TerrainGenerator::TerrainGenerator(const TerrainConfig& cfg)
    : config(cfg)
    , rng(cfg.seed)
{
    initialize_noise();
    generate();
}

void TerrainGenerator::initialize_noise() {
    // Initialize permutation tables for Perlin noise
    perm_x.resize(256);
    perm_y.resize(256);
    
    std::vector<int> p(256);
    for (int i = 0; i < 256; ++i) p[i] = i;
    
    // Shuffle
    for (int i = 255; i > 0; --i) {
        std::uniform_int_distribution<int> dist(0, i);
        int j = dist(rng);
        std::swap(p[i], p[j]);
    }
    
    // Double the permutation for overflow
    for (int i = 0; i < 256; ++i) {
        perm_x[i] = p[i];
        perm_y[i] = p[(i + 1) % 256];
    }
}

void TerrainGenerator::generate() {
    heightmap.resize(config.width * config.depth);
    vertices.resize(config.width * config.depth);
    
    generate_heightmap();
    generate_biomes();
    generate_normals();
    build_mesh();
}

void TerrainGenerator::generate_heightmap() {
    float half_width = config.world_size * 0.5f;
    
    for (int z = 0; z < config.depth; ++z) {
        for (int x = 0; x < config.width; ++x) {
            float world_x = (static_cast<float>(x) / config.width - 0.5f) * config.world_size;
            float world_z = (static_cast<float>(z) / config.depth - 0.5f) * config.world_size;
            
            // FBM noise
            float height = fbm(world_x * config.frequency, 
                              world_z * config.frequency,
                              config.octaves,
                              config.persistence,
                              config.lacunarity);
            
            // Normalize to 0-1
            height = (height + 1.0f) * 0.5f;
            
            // Apply power curve for more mountainous terrain
            height = std::pow(height, 1.5f);
            
            // Scale to max height
            height *= config.max_height;
            
            // Island mask (fade out at edges)
            float dist_x = static_cast<float>(x) / config.width - 0.5f;
            float dist_z = static_cast<float>(z) / config.depth - 0.5f;
            float dist = std::sqrt(dist_x * dist_x + dist_z * dist_z) * 2.0f;
            float mask = 1.0f - std::pow(dist, 3.0f);
            mask = std::max(0.0f, mask);
            
            height *= mask;
            
            heightmap[heightmap_index(x, z)] = height;
            
            // Update vertex
            TerrainVertex& v = vertices[heightmap_index(x, z)];
            v.height = height;
            v.position = Vec3(world_x, height, world_z);
            v.uv = Vec2(static_cast<float>(x) / config.width,
                       static_cast<float>(z) / config.depth);
        }
    }
    
    // Apply erosion
    if (config.enable_thermal_erosion) {
        for (int i = 0; i < config.erosion_iterations; ++i) {
            apply_thermal_erosion(config.erosion_strength);
        }
    }
}

void TerrainGenerator::generate_normals() {
    // Calculate normals from heightmap
    for (int z = 0; z < config.depth; ++z) {
        for (int x = 0; x < config.width; ++x) {
            float hL = get_height_at(x - 1, z);
            float hR = get_height_at(x + 1, z);
            float hD = get_height_at(x, z - 1);
            float hU = get_height_at(x, z + 1);
            
            Vec3 normal = Vec3(hL - hR, 2.0f, hD - hU).normalized();
            
            vertices[heightmap_index(x, z)].normal = normal;
        }
    }
}

void TerrainGenerator::generate_biomes() {
    for (int z = 0; z < config.depth; ++z) {
        for (int x = 0; x < config.width; ++x) {
            int idx = heightmap_index(x, z);
            TerrainVertex& v = vertices[idx];
            
            float world_x = (static_cast<float>(x) / config.width - 0.5f) * config.world_size;
            float world_z = (static_cast<float>(z) / config.depth - 0.5f) * config.world_size;
            
            // Temperature based on latitude and height
            float latitude = std::abs(world_z / (config.world_size * 0.5f));
            v.temperature = 1.0f - latitude - (v.height / config.max_height) * 0.3f;
            v.temperature = std::clamp(v.temperature, 0.0f, 1.0f);
            
            // Moisture based on noise
            float moisture_noise = fbm(world_x * 0.02f, world_z * 0.02f, 3, 0.5f, 2.0f);
            v.moisture = (moisture_noise + 1.0f) * 0.5f;
            
            // Determine biome
            v.biome_type = static_cast<int>(determine_biome(v.height, v.moisture, v.temperature));
        }
    }
}

void TerrainGenerator::build_mesh() {
    indices.clear();
    
    for (int z = 0; z < config.depth - 1; ++z) {
        for (int x = 0; x < config.width - 1; ++x) {
            int idx = heightmap_index(x, z);
            
            // Triangle 1
            indices.push_back(idx);
            indices.push_back(idx + 1);
            indices.push_back(idx + config.width);
            
            // Triangle 2
            indices.push_back(idx + 1);
            indices.push_back(idx + config.width + 1);
            indices.push_back(idx + config.width);
        }
    }
}

float TerrainGenerator::perlin_noise(float x, float y) {
    // Simplified 2D Perlin noise
    int X = static_cast<int>(std::floor(x)) & 255;
    int Y = static_cast<int>(std::floor(y)) & 255;
    
    x -= std::floor(x);
    y -= std::floor(y);
    
    float u = fade(x);
    float v = fade(y);
    
    int A = perm_x[X] + Y;
    int B = perm_x[X + 1] + Y;
    
    return lerp(
        lerp(gradient(perm_x[A], x, y), gradient(perm_x[B], x - 1.0f, y), u),
        lerp(gradient(perm_x[A + 1], x, y - 1.0f), gradient(perm_x[B + 1], x - 1.0f, y - 1.0f), u),
        v
    );
}

float TerrainGenerator::perlin_noise_2d(float x, float y) {
    return perlin_noise(x, y);
}

float TerrainGenerator::simplex_noise(float x, float y) {
    // Simplified simplex-like noise (approximation)
    return perlin_noise(x, y);
}

float TerrainGenerator::value_noise(float x, float y) {
    int X = static_cast<int>(std::floor(x)) & 255;
    int Y = static_cast<int>(std::floor(y)) & 255;
    
    x -= std::floor(x);
    y -= std::floor(y);
    
    float u = fade(x);
    float v = fade(y);
    
    int A = perm_x[X] + Y;
    int B = perm_x[X + 1] + Y;
    
    return lerp(
        static_cast<float>(perm_x[A]) / 255.0f,
        static_cast<float>(perm_x[B]) / 255.0f,
        u
    );
}

float TerrainGenerator::fbm(float x, float y, int octaves, 
                            float persistence, float lacunarity) {
    float total = 0.0f;
    float amplitude = 1.0f;
    float frequency = 1.0f;
    float max_value = 0.0f;
    
    for (int i = 0; i < octaves; ++i) {
        total += perlin_noise(x * frequency, y * frequency) * amplitude;
        max_value += amplitude;
        amplitude *= persistence;
        frequency *= lacunarity;
    }
    
    return total / max_value;
}

Vec2 TerrainGenerator::gradient(int hash, float x, float y) const {
    int h = hash & 7;
    float u = h < 4 ? x : y;
    float v = h < 4 ? y : x;
    return ((h & 1) ? -u : u) + ((h & 2) ? -2.0f * v : 2.0f * v) == Vec2(0, 0) ? 
           Vec2(0, 0) : Vec2(((h & 1) ? -u : u), ((h & 2) ? -2.0f * v : 2.0f * v));
}

void TerrainGenerator::apply_thermal_erosion(float dt) {
    std::vector<float> new_heightmap = heightmap;
    
    for (int z = 1; z < config.depth - 1; ++z) {
        for (int x = 1; x < config.width - 1; ++x) {
            float h = heightmap[heightmap_index(x, z)];
            
            // Check neighbors
            float hL = heightmap[heightmap_index(x - 1, z)];
            float hR = heightmap[heightmap_index(x + 1, z)];
            float hD = heightmap[heightmap_index(x, z - 1)];
            float hU = heightmap[heightmap_index(x, z + 1)];
            
            // Talus angle (slope threshold)
            float talus = 0.1f;
            
            // Calculate gradients
            float max_slope = 0.0f;
            float total_diff = 0.0f;
            
            if (hL - h > talus) {
                max_slope = std::max(max_slope, hL - h);
                total_diff += hL - h;
            }
            if (hR - h > talus) {
                max_slope = std::max(max_slope, hR - h);
                total_diff += hR - h;
            }
            if (hD - h > talus) {
                max_slope = std::max(max_slope, hD - h);
                total_diff += hD - h;
            }
            if (hU - h > talus) {
                max_slope = std::max(max_slope, hU - h);
                total_diff += hU - h;
            }
            
            if (total_diff > talus) {
                float amount = dt * (max_slope - talus) * 0.5f;
                new_heightmap[heightmap_index(x, z)] -= amount;
            }
        }
    }
    
    heightmap = new_heightmap;
    
    // Update vertex heights
    for (int z = 0; z < config.depth; ++z) {
        for (int x = 0; x < config.width; ++x) {
            int idx = heightmap_index(x, z);
            vertices[idx].height = heightmap[idx];
            vertices[idx].position.y = heightmap[idx];
        }
    }
    
    generate_normals();
}

void TerrainGenerator::apply_hydraulic_erosion(float iterations) {
    // Simplified hydraulic erosion
    for (float iter = 0; iter < iterations; ++iter) {
        // Random drop placement
        float drop_x = static_cast<float>(rand() % config.width);
        float drop_z = static_cast<float>(rand() % config.depth);
        
        float energy = 1.0f;
        float sediment = 0.0f;
        
        // Follow gradient
        for (int step = 0; step < 50; ++step) {
            int x = static_cast<int>(drop_x);
            int z = static_cast<int>(drop_z);
            
            if (!in_bounds(x, z)) break;
            
            float h = get_height_at(x, z);
            
            // Pick direction of steepest descent
            float hL = get_height_at(x - 1, z);
            float hR = get_height_at(x + 1, z);
            float hD = get_height_at(x, z - 1);
            float hU = get_height_at(x, z + 1);
            
            float min_h = std::min(hL, std::min(hR, std::min(hD, hU)));
            
            if (min_h >= h) break;  // Local minimum
            
            // Move to lower neighbor
            Vec2 dir(0, 0);
            if (hL == min_h) dir = Vec2(-1, 0);
            else if (hR == min_h) dir = Vec2(1, 0);
            else if (hD == min_h) dir = Vec2(0, -1);
            else if (hU == min_h) dir = Vec2(0, 1);
            
            drop_x += dir.x;
            drop_z += dir.y;
            
            // Carry sediment
            float capacity = std::max(0.0f, (min_h - h) * energy);
            if (sediment > capacity) {
                // Deposit sediment
                sediment -= capacity;
                int idx = heightmap_index(x, z);
                heightmap[idx] += sediment * config.erosion_strength;
            } else {
                // Erode
                sediment = capacity;
                int idx = heightmap_index(x, z);
                heightmap[idx] -= sediment * config.erosion_strength;
            }
            
            energy *= 0.99f;  // Dissipate
        }
    }
    
    // Update vertices
    for (int z = 0; z < config.depth; ++z) {
        for (int x = 0; x < config.width; ++x) {
            int idx = heightmap_index(x, z);
            vertices[idx].height = heightmap[idx];
            vertices[idx].position.y = heightmap[idx];
        }
    }
    
    generate_normals();
}

float TerrainGenerator::get_height_at(float x, float z) const {
    Vec2 grid = world_to_grid(x, z);
    int x0 = static_cast<int>(grid.x);
    int z0 = static_cast<int>(grid.y);
    int x1 = std::min(x0 + 1, config.width - 1);
    int z1 = std::min(z0 + 1, config.depth - 1);
    
    x0 = std::max(0, x0);
    z0 = std::max(0, z0);
    
    float fx = grid.x - x0;
    float fz = grid.y - z0;
    
    float h00 = get_height_at(x0, z0);
    float h10 = get_height_at(x1, z0);
    float h01 = get_height_at(x0, z1);
    float h11 = get_height_at(x1, z1);
    
    float h0 = h00 * (1.0f - fx) + h10 * fx;
    float h1 = h01 * (1.0f - fx) + h11 * fx;
    
    return h0 * (1.0f - fz) + h1 * fz;
}

float TerrainGenerator::get_height_at(int x, int z) const {
    if (!in_bounds(x, z)) return 0.0f;
    return heightmap[heightmap_index(x, z)];
}

Vec3 TerrainGenerator::get_normal_at(float x, float z) const {
    float hL = get_height_at(x - 1, z);
    float hR = get_height_at(x + 1, z);
    float hD = get_height_at(x, z - 1);
    float hU = get_height_at(x, z + 1);
    
    return Vec3(hL - hR, 2.0f, hD - hU).normalized();
}

float TerrainGenerator::get_biome_temperature(float x, float z) const {
    return vertices[heightmap_index(static_cast<int>(x), static_cast<int>(z))].temperature;
}

float TerrainGenerator::get_biome_moisture(float x, float z) const {
    return vertices[heightmap_index(static_cast<int>(x), static_cast<int>(z))].moisture;
}

BiomeType TerrainGenerator::determine_biome(float height, float moisture, float temperature) const {
    if (height < config.sea_level - 5.0f) return BiomeType::DeepOcean;
    if (height < config.sea_level) return BiomeType::Ocean;
    if (height < config.sea_level + 2.0f) return BiomeType::Beach;
    
    if (temperature < 0.2f) {
        if (height > config.max_height * 0.8f) return BiomeType::Snow;
        return BiomeType::Taiga;
    }
    
    if (temperature > 0.7f) {
        if (moisture < 0.3f) return BiomeType::Desert;
        if (height > config.max_height * 0.7f) return BiomeType::Mountains;
        return BiomeType::Plains;
    }
    
    if (moisture < 0.3f) return BiomeType::Desert;
    if (moisture > 0.7f && height < config.sea_level + 5.0f) return BiomeType::Swamp;
    if (height > config.max_height * 0.7f) return BiomeType::Mountains;
    if (height > config.max_height * 0.5f) return BiomeType::Hills;
    if (moisture > 0.5f) return BiomeType::Forest;
    
    return BiomeType::Plains;
}

BiomeType TerrainGenerator::get_biome_at(float x, float z) const {
    int idx = heightmap_index(static_cast<int>(x), static_cast<int>(z));
    return static_cast<BiomeType>(vertices[idx].biome_type);
}

Vec2 TerrainGenerator::world_to_grid(float x, float z) const {
    return Vec2(
        (x / config.world_size + 0.5f) * config.width,
        (z / config.world_size + 0.5f) * config.depth
    );
}

Vec2 TerrainGenerator::grid_to_world(int x, int z) const {
    return Vec2(
        (static_cast<float>(x) / config.width - 0.5f) * config.world_size,
        (static_cast<float>(z) / config.depth - 0.5f) * config.world_size
    );
}

void TerrainGenerator::modify_height(int x, int z, float delta) {
    if (!in_bounds(x, z)) return;
    int idx = heightmap_index(x, z);
    heightmap[idx] += delta;
    vertices[idx].height = heightmap[idx];
    vertices[idx].position.y = heightmap[idx];
}

void TerrainGenerator::set_height(int x, int z, float height) {
    if (!in_bounds(x, z)) return;
    int idx = heightmap_index(x, z);
    heightmap[idx] = height;
    vertices[idx].height = height;
    vertices[idx].position.y = height;
}

void TerrainGenerator::smooth(float radius) {
    std::vector<float> smoothed = heightmap;
    int r = static_cast<int>(radius);
    
    for (int z = r; z < config.depth - r; ++z) {
        for (int x = r; x < config.width - r; ++x) {
            float sum = 0.0f;
            int count = 0;
            
            for (int dz = -r; dz <= r; ++dz) {
                for (int dx = -r; dx <= r; ++dx) {
                    sum += heightmap[heightmap_index(x + dx, z + dz)];
                    ++count;
                }
            }
            
            smoothed[heightmap_index(x, z)] = sum / count;
        }
    }
    
    heightmap = smoothed;
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].height = heightmap[i];
        vertices[i].position.y = heightmap[i];
    }
    
    generate_normals();
}

void TerrainGenerator::terrace(int levels) {
    float level_height = config.max_height / levels;
    
    for (auto& h : heightmap) {
        h = std::round(h / level_height) * level_height;
    }
    
    for (auto& v : vertices) {
        v.height = std::round(v.height / level_height) * level_height;
        v.position.y = v.height;
    }
    
    generate_normals();
}

void TerrainGenerator::erode_channel(const Vec3& start, const Vec3& end, float depth) {
    Vec2 start_grid = world_to_grid(start.x, start.z);
    Vec2 end_grid = world_to_grid(end.x, end.z);
    
    // Simple line erosion
    int steps = 100;
    for (int i = 0; i <= steps; ++i) {
        float t = static_cast<float>(i) / steps;
        float x = start_grid.x + (end_grid.x - start_grid.x) * t;
        float z = start_grid.y + (end_grid.y - start_grid.y) * t;
        
        int ix = static_cast<int>(x);
        int iz = static_cast<int>(z);
        
        for (int dz = -2; dz <= 2; ++dz) {
            for (int dx = -2; dx <= 2; ++dx) {
                float dist = std::sqrt(dx * dx + dz * dz);
                float amount = depth * (1.0f - dist / 3.0f);
                modify_height(ix + dx, iz + dz, -amount);
            }
        }
    }
}

// ============== TerrainChunk Implementation ==============

TerrainChunk::TerrainChunk(int x, int z, const std::shared_ptr<TerrainGenerator>& gen, int lod)
    : chunk_x(x), chunk_z(z), generator(gen), lod_level(lod) 
{
    generate();
}

void TerrainChunk::generate() {
    // Generate chunk vertices based on parent generator
    float chunk_size = config.chunk_world_size;
    int res = config.grid_size;
    
    vertices.resize(res * res);
    indices.clear();
    
    for (int z = 0; z < res; ++z) {
        for (int x = 0; x < res; ++x) {
            float world_x = chunk_x * chunk_size + (static_cast<float>(x) / res) * chunk_size;
            float world_z = chunk_z * chunk_size + (static_cast<float>(z) / res) * chunk_size;
            
            TerrainVertex v;
            v.position = Vec3(world_x, generator->get_height_at(world_x, world_z), world_z);
            v.normal = generator->get_normal_at(world_x, world_z);
            v.uv = Vec2(static_cast<float>(x) / res, static_cast<float>(z) / res);
            v.height = v.position.y;
            
            vertices[z * res + x] = v;
        }
    }
    
    // Build indices
    for (int z = 0; z < res - 1; ++z) {
        for (int x = 0; x < res - 1; ++x) {
            int idx = z * res + x;
            indices.push_back(idx);
            indices.push_back(idx + 1);
            indices.push_back(idx + res);
            indices.push_back(idx + 1);
            indices.push_back(idx + res + 1);
            indices.push_back(idx + res);
        }
    }
}

void TerrainChunk::update(const Vec3& camera_position) {
    // LOD could be updated here based on distance
}

Vec2 TerrainChunk::get_world_offset() const {
    return Vec2(
        chunk_x * config.chunk_world_size,
        chunk_z * config.chunk_world_size
    );
}

// ============== TerrainWorld Implementation ==============

TerrainWorld::TerrainWorld(const Config& cfg)
    : config(cfg)
    , generator(std::make_shared<TerrainGenerator>(cfg.terrain))
{}

void TerrainWorld::initialize() {
    if (config.use_chunks) {
        int chunks = config.chunks_visible;
        for (int z = -chunks/2; z < chunks/2; ++z) {
            for (int x = -chunks/2; x < chunks/2; ++x) {
                chunks.push_back(std::make_shared<TerrainChunk>(x, z, generator));
            }
        }
    }
}

void TerrainWorld::update(const Vec3& camera_position) {
    if (config.use_chunks) {
        for (auto& chunk : chunks) {
            chunk->update(camera_position);
        }
    }
}

float TerrainWorld::get_height_at(float x, float z) const {
    return generator->get_height_at(x, z);
}

Vec3 TerrainWorld::get_normal_at(float x, float z) const {
    return generator->get_normal_at(x, z);
}

BiomeType TerrainWorld::get_biome_at(float x, float z) const {
    return generator->get_biome_at(x, z);
}

} // namespace physics
} // namespace atlas
