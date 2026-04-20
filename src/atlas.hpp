/**
 * @file atlas.hpp
 * @brief Atlas Physics Simulator - Main Header
 * 
 * A comprehensive physics simulation library featuring:
 * - Rigid body dynamics with Newton-Euler equations
 * - Soft body simulation using mass-spring systems
 * - Fluid dynamics (SPH and Navier-Stokes)
 * - Cloth simulation with Verlet integration
 * - Hair/fur simulation
 * - Terrain generation with procedural algorithms
 * - Real-time raycasting and spatial queries
 * - Collision detection (BVH, SAT)
 * - Constraint solver (distance, hinge, ball-socket)
 */

#pragma once

#include "core/vec3.hpp"
#include "core/quaternion.hpp"
#include "core/matrix4x4.hpp"
#include "core/vec2.hpp"
#include "core/transform.hpp"

#include "physics/rigid_body/rigid_body.hpp"
#include "physics/rigid_body/rigid_body_world.hpp"

#include "physics/soft_body/soft_body.hpp"

#include "physics/fluid/sph_solver.hpp"
#include "physics/fluid/navier_stokes.hpp"

#include "physics/cloth/cloth.hpp"

#include "physics/hair/hair.hpp"

#include "physics/terrain/terrain.hpp"

#include "physics/raycasting/raycasting.hpp"

#include "physics/constraints/constraints.hpp"

#include <memory>

namespace atlas {

/**
 * @brief Atlas Physics World
 * 
 * Combines all physics subsystems into a unified world.
 */
class PhysicsWorld {
public:
    physics::RigidBodyWorld rigid_world;
    physics::SoftBodyWorld soft_world;
    physics::ClothWorld cloth_world;
    physics::SPHWorld sph_world;
    physics::NavierStokesWorld navier_stokes_world;
    physics::ConstraintWorld constraint_world;
    
    bool paused;
    float time_scale;
    
    PhysicsWorld() : paused(false), time_scale(1.0f) {
        rigid_world.config.substeps = 4;
        soft_world.config.substeps = 2;
        cloth_world.config.substeps = 4;
    }
    
    void update(float dt) {
        if (paused) return;
        dt *= time_scale;
        
        // Fixed timestep with accumulator
        static float accumulator = 0.0f;
        const float fixed_dt = 1.0f / 120.0f;
        
        accumulator += dt;
        
        while (accumulator >= fixed_dt) {
            // Update rigid bodies
            rigid_world.step_fixed(fixed_dt);
            
            // Solve constraints
            constraint_world.solve(fixed_dt);
            
            // Update soft bodies
            soft_world.step_fixed(fixed_dt);
            
            // Update cloth
            cloth_world.step_fixed(fixed_dt);
            
            // Update SPH fluid
            sph_world.update(fixed_dt);
            
            // Update Navier-Stokes
            navier_stokes_world.step(fixed_dt);
            
            accumulator -= fixed_dt;
        }
    }
    
    void pause() { paused = true; }
    void resume() { paused = false; }
    void set_time_scale(float scale) { time_scale = scale; }
};

} // namespace atlas
