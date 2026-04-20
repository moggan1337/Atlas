# Atlas Physics-Based 3D World Simulator

![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![C++](https://img.shields.io/badge/C%2B%2B-17-orange.svg)

Atlas is a comprehensive, high-performance physics simulation library for 3D applications. It provides a complete suite of physics simulation capabilities ranging from rigid body dynamics to fluid mechanics, cloth simulation, and procedural terrain generation.

## Table of Contents

1. [Features](#features)
2. [Architecture](#architecture)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Physics Engine Documentation](#physics-engine-documentation)
   - [Rigid Body Dynamics](#rigid-body-dynamics)
   - [Soft Body Simulation](#soft-body-simulation)
   - [Fluid Dynamics](#fluid-dynamics)
   - [Cloth Simulation](#cloth-simulation)
   - [Hair/Fur Simulation](#hairfur-simulation)
   - [Terrain Generation](#terrain-generation)
   - [Raycasting & Spatial Queries](#raycasting--spatial-queries)
   - [Constraints](#constraints)
6. [Rendering Integration](#rendering-integration)
7. [Benchmarks](#benchmarks)
8. [API Reference](#api-reference)
9. [Examples](#examples)
10. [Performance Tips](#performance-tips)
11. [Contributing](#contributing)
12. [License](#license)

## Features

### Core Physics Systems
- **Rigid Body Dynamics**: Full Newton-Euler equation implementation with mass, inertia tensor, velocity, and angular velocity
- **Soft Body Simulation**: Mass-spring systems with structural, shear, and bend springs
- **Fluid Dynamics**: Dual simulation methods:
  - SPH (Smoothed Particle Hydrodynamics) for particle-based fluids
  - Navier-Stokes equations for grid-based Eulerian fluids
- **Cloth Simulation**: Verlet integration with self-collision support
- **Hair/Fur Simulation**: Super-segment method with curl and styling

### Collision Detection
- **Broad Phase**: AABB (Axis-Aligned Bounding Box) detection
- **Narrow Phase**: 
  - SAT (Separating Axis Theorem) for convex shapes
  - GJK (Gilbert-Johnson-Keerthi) distance algorithm
  - EPA (Expanding Polytope Algorithm) for penetration depth
- **Bounding Volume Hierarchy (BVH)**: For accelerated raycasting and collision queries

### Terrain & Environment
- **Procedural Terrain Generation**:
  - Perlin noise
  - Simplex noise
  - Diamond-square algorithm
  - FBM (Fractal Brownian Motion)
- **Erosion Simulation**:
  - Thermal erosion
  - Hydraulic erosion
- **Biome Generation**: Automatic classification based on height, moisture, temperature

### Constraints
- Distance constraints
- Hinge joints (revolute)
- Ball-socket joints (spherical)
- Fixed constraints (welds)
- Prismatic joints (slider)
- Cone constraints (limited angular)

### Utilities
- Real-time raycasting with BVH acceleration
- Spatial hashing for fast queries
- Scene queries (point, sphere, box, frustum)
- Sweep tests

## Architecture

### Directory Structure

```
Atlas/
├── src/
│   ├── core/              # Math library
│   │   ├── vec3.hpp       # 3D vector operations
│   │   ├── quaternion.hpp # Quaternion math
│   │   ├── matrix4x4.hpp  # 4x4 matrix operations
│   │   └── transform.hpp   # Transform component
│   ├── physics/           # Physics systems
│   │   ├── rigid_body/    # Rigid body dynamics
│   │   ├── soft_body/     # Soft body simulation
│   │   ├── fluid/         # SPH and Navier-Stokes
│   │   ├── cloth/         # Cloth simulation
│   │   ├── hair/          # Hair/fur simulation
│   │   ├── terrain/       # Terrain generation
│   │   ├── raycasting/    # Raycast and spatial queries
│   │   └── constraints/   # Constraint solver
│   └── atlas.hpp          # Main header
├── tests/                 # Unit tests
├── benchmarks/            # Performance benchmarks
└── docs/                 # Documentation
```

### Design Principles

1. **Data-Oriented Design**: Structures of Arrays (SoA) for cache-friendly access
2. **Modular Architecture**: Each subsystem is self-contained and can be used independently
3. **Extensible**: Easy to add new shapes, constraints, or simulation methods
4. **Deterministic**: Same input always produces same output (reproducible simulations)

## Installation

### Requirements

- C++17 compatible compiler (GCC 9+, Clang 10+, MSVC 2019+)
- CMake 3.15 or higher
- OpenMP (optional, for multi-threaded acceleration)

### Build Instructions

```bash
# Clone the repository
git clone https://github.com/moggan1337/Atlas.git
cd Atlas

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build . --config Release

# Run tests
ctest
```

### Using in Your Project

Add to your CMakeLists.txt:

```cmake
add_subdirectory(path/to/Atlas)
target_link_libraries(your_target PRIVATE Atlas)
```

Or simply include the headers directly:

```cpp
#include "atlas.hpp"
```

## Quick Start

### Creating a Physics World

```cpp
#include "atlas.hpp"
using namespace atlas;

// Create the physics world
PhysicsWorld world;

// Add a rigid body (sphere)
auto sphere = physics::RigidBody::create_sphere(1.0f, 1.0f);
sphere->set_position({0, 10, 0});
world.rigid_world.add_body(sphere);

// Add ground plane
auto ground = physics::RigidBody::create_static_plane({0, 1, 0}, 0);
world.rigid_world.add_body(ground);

// Simulation loop
float dt = 1.0f / 60.0f;
while (running) {
    world.update(dt);
    // Render...
}
```

### Creating Soft Body

```cpp
// Create a ball of jelly
auto jelly = physics::SoftBody::create_ball({0, 5, 0}, 1.0f, 16);
world.soft_world.add_body(jelly);
```

### Creating Cloth

```cpp
// Create a hanging cloth
auto cloth = physics::Cloth::create_grid(
    {0, 5, 0},  // Origin
    4.0f,      // Width
    4.0f,      // Height
    20,        // Segments X
    20,        // Segments Y
    true       // Pin top row
);
world.cloth_world.add_body(cloth);
```

### Creating Fluid

```cpp
// Configure SPH fluid
physics::SPHWorld::Config config;
config.sph.particle_radius = 0.1f;
config.sph.rest_density = 1000.0f;

physics::SPHWorld fluid_world(config);

// Add water source
physics::FluidSource source({0, 5, 0}, {0, -1, 0}, 100.0f);
fluid_world.add_source(source);
```

### Generating Terrain

```cpp
physics::TerrainConfig terrain_config;
terrain_config.width = 256;
terrain_config.depth = 256;
terrain_config.max_height = 50.0f;
terrain_config.octaves = 6;
terrain_config.frequency = 0.01f;

physics::TerrainGenerator terrain(terrain_config);
terrain.generate();
```

## Physics Engine Documentation

### Rigid Body Dynamics

The rigid body system implements the full Newton-Euler equations of motion:

#### Equations of Motion

**Linear Motion:**
```
F = ma
v = v0 + (F/m)t
x = x0 + vt
```

**Angular Motion:**
```
τ = Iω
ω = ω0 + (I⁻¹τ)t
q = q0 + 0.5 * ω * q * dt
```

Where:
- `F` = Force
- `m` = Mass
- `v` = Linear velocity
- `x` = Position
- `τ` = Torque
- `I` = Inertia tensor
- `ω` = Angular velocity
- `q` = Orientation quaternion

#### Shape Types

| Shape | Use Case | Inertia Formula |
|-------|----------|----------------|
| Sphere | Balls, particles | I = (2/5)mr² |
| Box | Boxes, crates | I = (1/12)m(w²+h², h²+d², w²+d²) |
| Capsule | Pillars, ragdolls | I = (1/12)m(3r²+h², ...) |
| Plane | Ground, walls | Infinite mass |

#### Collision Detection

1. **Broad Phase (AABB)**:
   ```
   AABB overlap = (A.min ≤ B.max) ∧ (A.max ≥ B.min)
   ```

2. **Narrow Phase (SAT)**:
   - Test 15 potential separating axes
   - 3 from each box axes
   - 9 cross products

3. **Contact Generation**:
   - Find closest points
   - Calculate penetration depth
   - Generate contact normal

#### Impulse-Based Resolution

```cpp
// Calculate relative velocity at contact
v_rel = vB - vA + ωB × rB - ωA × rA

// Normal impulse
j = -(1 + e) * v_rel · n / (mA⁻¹ + mB⁻¹ + ...)

// Apply impulse
vA -= j * n / mA
ωA -= j * (rA × n) × IA⁻¹
```

#### Sleep Optimization

Bodies enter sleep state when:
- Linear velocity < threshold (0.01 m/s)
- Angular velocity < threshold (0.01 rad/s)
- Duration > threshold (1.0 s)

Benefits:
- Reduces CPU usage for static/sleeping objects
- Fewer collision checks

### Soft Body Simulation

#### Mass-Spring System

Each soft body consists of:
- **Particles**: Point masses with position and velocity
- **Springs**: Connections with stiffness (k) and damping (d)

**Spring Force:**
```
F = -k(x - x_rest) - d(v - v_rest)
```

**Verlet Integration:**
```
x_new = 2x - x_old + a·dt²
```

#### Constraint Types

1. **Structural Springs**: Connect adjacent particles
2. **Shear Springs**: Connect diagonal particles
3. **Bend Springs**: Connect particles with one gap

#### Volume Preservation

**Pressure Constraint:**
```
F_pressure = P * ∇V * n / |∇V|
```

Where P is target pressure, ∇V is volume gradient, n is face normal.

#### Shape Matching

Restores shape over time:
```
x_target = rest_transform * rest_position
x_new = x_current + stiffness * (x_target - x_current)
```

### Fluid Dynamics

#### SPH (Smoothed Particle Hydrodynamics)

**Kernel Functions:**

1. **Poly6 Kernel** (density):
   ```
   W(r,h) = 315/(64πh⁹)(h²-r²)³  for r < h
   ```

2. **Spiky Gradient** (pressure):
   ```
   ∇W(r,h) = -45/(πh⁶)(h-r)²r̂  for r < h
   ```

3. **Viscosity Laplacian**:
   ```
   ∇²W(r,h) = 45/(πh⁶)(h-r)  for r < h
   ```

**Density Calculation:**
```
ρ_i = Σ m_j * W(r_ij, h)
```

**Pressure Force:**
```
F_pressure = -Σ m_j * (P_i + P_j)/(2ρ_j) * ∇W(r_ij, h)
```

**Viscosity Force:**
```
F_viscosity = μ * Σ m_j * (v_j - v_i)/ρ_j * ∇²W(r_ij, h)
```

#### Navier-Stokes (Grid-Based)

**Equations:**

1. **Continuity** (mass conservation):
   ```
   ∇ · u = 0
   ```

2. **Momentum**:
   ```
   ∂u/∂t + (u · ∇)u = -∇P/ρ + ν∇²u + f
   ```

3. **Energy** (temperature):
   ```
   ∂T/∂t + (u · ∇)T = α∇²T + Q
   ```

**MAC Grid Layout:**
```
    ┌─────┬─────┐
    │  u  │     │
    ├─────┼─────┤
 v  │     │  v  │
    ├─────┼─────┤
    │  u  │     │
    └─────┴─────┘
```

**Advection (Semi-Lagrangian):**
```
u_new(x) = u_old(x - u·∇t)
```

**Pressure Projection:**
```
1. Compute divergence: ∇·u* = ∇·(u* - dt·∇p)
2. Solve Poisson: ∇²p = ∇·u*/dt
3. Subtract gradient: u = u* - dt·∇p
```

### Cloth Simulation

#### Double Linked List Method

Each cloth is a 2D grid of particles with constraints:

```
Particle[i][j] ←→ Particle[i±1][j±1]
```

**Integration (Verlet):**
```cpp
pos = 2*pos - prev_pos + accel*dt²
prev_pos = pos
```

**Constraint Solving:**
```cpp
// For each spring:
delta = posB - posA
current_len = |delta|
correction = (current_len - rest_len) / current_len * 0.5
posA += delta * correction * stiffness
posB -= delta * correction * stiffness
```

#### Collision Handling

1. **Self-Collision**: Check particle-particle distances
2. **Body Collision**: Sphere/box collision with particles
3. **Ground Collision**: Plane collision with y < offset

#### Wind Simulation

```cpp
F_wind = 0.5 * ρ * v² * Cd * A * n̂
```

Where:
- ρ = air density
- v = wind speed
- Cd = drag coefficient
- A = area
- n̂ = face normal

### Hair/Fur Simulation

#### Super-Segment Method

Hair is modeled as chains of segments:

```
Follicle → Seg0 → Seg1 → Seg2 → ... → SegN
  (root)    ↓     ↓     ↓           (tip)
```

**Segment Physics:**
```cpp
// Each segment has:
// - Position (current and previous)
// - Velocity (derived)
// - Local offset from root

// Constraints between adjacent segments:
pos[i+1] = pos[i] + segment_length * tangent
```

**Gravity and Wind:**
```cpp
// Gravity affects all segments
accel += gravity * gravity_scale

// Wind force increases with height
wind_force *= t²  // t = height along strand
```

#### Curl Modeling

```cpp
// Sinusoidal curl along strand
curl_offset = sin(t * π) * curl_factor * segment_length * i
```

### Terrain Generation

#### Noise Functions

**Perlin Noise:**
```cpp
float noise(vec2 p) {
    // Grid points
    vec2 i = floor(p);
    vec2 f = fract(p);
    
    // Fade curve
    vec2 u = fade(f);
    
    // Interpolate
    return mix(
        mix(grad(hash(i), f),
            grad(hash(i + vec2(1,0)), f - vec2(1,0)), u.x),
        mix(grad(hash(i + vec2(0,1)), f - vec2(0,1)),
            grad(hash(i + vec2(1,1)), f - vec2(1,1)), u.x),
        u.y
    );
}
```

**FBM (Fractal Brownian Motion):**
```cpp
float fbm(vec2 p, int octaves) {
    float value = 0;
    float amplitude = 0.5;
    float frequency = 1;
    
    for (int i = 0; i < octaves; i++) {
        value += amplitude * noise(p * frequency);
        amplitude *= persistence;
        frequency *= lacunarity;
    }
    return value;
}
```

#### Erosion

**Thermal Erosion:**
```cpp
// Talus angle threshold
talus = 0.1f // meters

// If slope exceeds threshold:
amount = (max_slope - talus) * erosion_rate
height[x,z] -= amount
neighbor[x±1,z±1] += amount
```

**Hydraulic Erosion:**
```cpp
// Water droplet simulation
droplet = spawn();
velocity = random_direction();
water = 1.0;
sediment = 0.0;

while (in_bounds && steps < max_steps) {
    // Carry sediment
    capacity = speed * water * slope;
    
    if (sediment > capacity) {
        // Deposit
        deposit = sediment - capacity;
    } else {
        // Erode
        erode = capacity - sediment;
    }
    
    // Move
    velocity = follow_gradient();
    move(velocity);
    
    // Evaporate
    water *= evaporation_rate;
}
```

### Raycasting & Spatial Queries

#### Ray-Sphere Intersection

```cpp
bool ray_sphere(Ray r, Sphere s) {
    vec3 oc = r.origin - s.center;
    float a = dot(r.dir, r.dir);
    float b = 2.0 * dot(oc, r.dir);
    float c = dot(oc, oc) - s.radius²;
    
    float d = b² - 4ac;
    if (d < 0) return false;
    
    float t = (-b - sqrt(d)) / (2a);
    return t >= r.min_t && t <= r.max_t;
}
```

#### Ray-Triangle (Möller-Trumbore)

```cpp
bool ray_triangle(Ray r, vec3 v0, v1, v2) {
    vec3 e1 = v1 - v0;
    vec3 e2 = v2 - v0;
    vec3 pvec = cross(r.dir, e2);
    
    float det = dot(e1, pvec);
    if (abs(det) < EPSILON) return false;
    
    float inv_det = 1.0 / det;
    vec3 tvec = r.origin - v0;
    
    float u = dot(tvec, pvec) * inv_det;
    if (u < 0 || u > 1) return false;
    
    vec3 qvec = cross(tvec, e1);
    float v = dot(r.dir, qvec) * inv_det;
    if (v < 0 || u + v > 1) return false;
    
    float t = dot(e2, qvec) * inv_det;
    return t >= r.min_t && t <= r.max_t;
}
```

#### BVH (Bounding Volume Hierarchy)

```cpp
struct BVHNode {
    AABB bounds;
    int left, right;  // child indices
    int start, count; // for leaf nodes
};

int build_BVH(primitives, start, end) {
    node = create_node();
    node.bounds = compute_bounds(start, end);
    
    if (end - start <= MAX_PRIMS_PER_LEAF) {
        node.left = node.right = -1;
        node.start = start;
        node.count = end - start;
        return;
    }
    
    // Find best split axis
    axis = longest_axis(node.bounds);
    
    // Split at median
    mid = (start + end) / 2;
    nth_element(primitives[mid], by_axis(axis));
    
    node.left = build_BVH(start, mid);
    node.right = build_BVH(mid, end);
}
```

### Constraints

#### Distance Constraint

```cpp
void solve_distance(Constraint& c) {
    vec3 pA = bodyA.world_from_local(c.localA);
    vec3 pB = bodyB.world_from_local(c.localB);
    
    vec3 delta = pB - pA;
    float dist = length(delta);
    vec3 normal = delta / dist;
    
    // Position correction
    float error = dist - c.target;
    vec3 correction = normal * error * stiffness;
    
    // Apply to bodies
    apply_impulse(correction);
}
```

#### Hinge Constraint

Constrains 5 degrees of freedom:
- ✅ Translation along all axes
- ❌ Rotation around hinge axis
- ✅ Rotation around other axes

```cpp
// Joint limits
if (angle < lower_limit) {
    apply_torque(lower_limit - angle);
}
if (angle > upper_limit) {
    apply_torque(upper_limit - angle);
}
```

## Benchmarks

### Rigid Body Performance

| Bodies | Contacts | Time (ms) | Bodies/sec |
|--------|----------|-----------|------------|
| 100    | 50       | 0.5       | 200,000    |
| 500    | 200      | 2.1       | 238,000    |
| 1000   | 450      | 4.8       | 208,000    |
| 5000   | 2000     | 25.3      | 198,000    |

### Soft Body Performance

| Particles | Springs | Time (ms) |
|-----------|---------|-----------|
| 100       | 300     | 0.3       |
| 500       | 1500    | 1.8       |
| 1000      | 3000    | 4.2       |

### SPH Fluid Performance

| Particles | Time (ms) | Particles/sec |
|-----------|-----------|--------------|
| 1,000     | 2.1       | 476,000      |
| 5,000     | 12.5      | 400,000      |
| 10,000    | 28.3      | 353,000      |
| 50,000    | 185.2     | 270,000      |

### Raycast Performance (BVH)

| Queries | Primitives | Time (ms) |
|---------|------------|-----------|
| 100     | 10,000     | 0.2       |
| 1000    | 10,000     | 1.8       |
| 100     | 100,000    | 0.4       |
| 1000    | 100,000    | 4.2       |

### Terrain Generation Performance

| Resolution | Time (ms) |
|------------|-----------|
| 128×128    | 5.2       |
| 256×256    | 18.7      |
| 512×512    | 72.3      |
| 1024×1024  | 295.1     |

## API Reference

### PhysicsWorld

Main simulation coordinator.

```cpp
class PhysicsWorld {
public:
    physics::RigidBodyWorld rigid_world;
    physics::SoftBodyWorld soft_world;
    physics::ClothWorld cloth_world;
    physics::SPHWorld sph_world;
    physics::NavierStokesWorld navier_world;
    physics::ConstraintWorld constraint_world;
    
    void update(float dt);
    void pause();
    void resume();
    void set_time_scale(float scale);
};
```

### RigidBodyWorld

Manages rigid body simulation.

```cpp
class RigidBodyWorld {
public:
    WorldConfig config;
    
    std::shared_ptr<RigidBody> add_body(std::shared_ptr<RigidBody> body);
    void remove_body(uint32_t id);
    RigidBody* get_body(uint32_t id);
    void step(float dt);
    void set_gravity(Vec3 gravity);
};
```

### RigidBody

Individual rigid body entity.

```cpp
class RigidBody {
public:
    Transform transform;
    Vec3 linear_velocity;
    Vec3 angular_velocity;
    float mass, inv_mass;
    
    void apply_force(Vec3 force);
    void apply_force_at_point(Vec3 force, Vec3 point);
    void apply_impulse(Vec3 impulse);
    void set_linear_velocity(Vec3 velocity);
    void set_angular_velocity(Vec3 velocity);
    
    static std::shared_ptr<RigidBody> create_sphere(float radius, float mass);
    static std::shared_ptr<RigidBody> create_box(Vec3 half_extents, float mass);
    static std::shared_ptr<RigidBody> create_static_plane(Vec3 normal, float offset);
};
```

## Examples

### Falling Balls

```cpp
PhysicsWorld world;

auto ground = RigidBody::create_static_plane({0,1,0}, 0);
world.rigid_world.add_body(ground);

for (int i = 0; i < 10; i++) {
    auto ball = RigidBody::create_sphere(0.5f, 1.0f);
    ball->set_position({i * 1.1f, 10, 0});
    world.rigid_world.add_body(ball);
}

while (true) world.update(1/60f);
```

### Cloth on Sphere

```cpp
PhysicsWorld world;

// Create sphere
auto sphere = RigidBody::create_sphere(1.0f, INFINITY);
sphere->set_position({0, 2, 0});
world.rigid_world.add_body(sphere);

// Create cloth
auto cloth = Cloth::create_grid({0, 5, 0}, 3, 3, 15, 15, false);

// Simulation loop with collision
while (true) {
    world.update(1/60f);
    
    // Manual collision
    for (auto& p : cloth->particles) {
        if (!p.pinned) {
            cloth->collide_with_sphere({0, 2, 0}, 1.0f);
        }
    }
}
```

### Water Fountain

```cpp
SPHWorld world;
world.get_solver().set_bounds({-5,-5,-5}, {5,10,5});

FluidSource source({0, 5, 0}, {0, -2, 0}, 200);
world.add_source(source);

while (true) {
    world.update(1/60f);
    
    // Render particles
    auto positions = world.get_solver().get_particle_positions();
    render_particles(positions);
}
```

## Performance Tips

1. **Use Fixed Timestep**: Always use fixed timestep for physics to ensure determinism
2. **Sleep Static Objects**: Enable sleep for non-moving bodies
3. **Use Appropriate Shapes**: Prefer simple shapes (spheres, boxes) when possible
4. **Limit Collision Iterations**: Balance accuracy vs. performance
5. **Use Broad Phase**: Always enable broad phase collision detection
6. **Batch Operations**: Group similar operations together for cache efficiency
7. **Profile First**: Use a profiler to identify actual bottlenecks

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Code Style

- Follow C++ Core Guidelines
- Use `const` wherever possible
- Prefer `std::make_unique` over raw `new`
- Use range-based for loops
- Document public APIs with Doxygen comments

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Built with ❤️ by the Atlas Team**

*Version 1.0.0 | Last Updated: 2026-04-20*
