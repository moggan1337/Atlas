/**
 * @file constraints.hpp
 * @brief Physics Constraints
 * 
 * Implements various constraint types:
 * - Distance constraints
 * - Hinge constraints
 * - Ball-socket (spherical) joints
 * - Fixed constraints
 * - Prismatic (slider) constraints
 */

#pragma once

#include "../../core/vec3.hpp"
#include "../../core/quaternion.hpp"
#include "../rigid_body/rigid_body.hpp"
#include <memory>
#include <vector>

namespace atlas {
namespace physics {

using core::Vec3;
using core::Quaternion;

/**
 * @brief Constraint Types
 */
enum class ConstraintType {
    Distance,
    Hinge,
    BallSocket,
    Fixed,
    Prismatic,
    Cone,
    Twist
};

/**
 * @brief Base Constraint
 */
class Constraint {
public:
    RigidBody* body_a;
    RigidBody* body_b;
    ConstraintType type;
    bool enabled;
    
    Constraint(ConstraintType t) : type(t), enabled(true) {}
    virtual ~Constraint() = default;
    
    virtual void pre_solve(float dt) {}
    virtual void solve() = 0;
    virtual void post_solve() {}
};

/**
 * @brief Distance Constraint
 */
class DistanceConstraint : public Constraint {
public:
    Vec3 local_anchor_a;
    Vec3 local_anchor_b;
    float target_distance;
    float stiffness;
    float damping;
    
    DistanceConstraint()
        : Constraint(ConstraintType::Distance)
        , target_distance(0.0f)
        , stiffness(0.9f)
        , damping(0.1f) {}
    
    DistanceConstraint(RigidBody* a, RigidBody* b, const Vec3& anchor_a, 
                      const Vec3& anchor_b, float dist)
        : Constraint(ConstraintType::Distance)
        , local_anchor_a(anchor_a)
        , local_anchor_b(anchor_b)
        , target_distance(dist)
        , stiffness(0.9f)
        , damping(0.1f) 
    {
        body_a = a;
        body_b = b;
    }
    
    void solve() override;
};

/**
 * @brief Hinge Constraint (revolute joint)
 */
class HingeConstraint : public Constraint {
public:
    Vec3 local_anchor_a;
    Vec3 local_anchor_b;
    Vec3 local_axis_a;
    Vec3 local_axis_b;
    
    bool use_limits;
    float lower_limit;
    float upper_limit;
    
    // Motor
    bool use_motor;
    float motor_velocity;
    float motor_torque;
    
    HingeConstraint()
        : Constraint(ConstraintType::Hinge)
        , use_limits(false)
        , lower_limit(-M_PI)
        , upper_limit(M_PI)
        , use_motor(false)
        , motor_velocity(0.0f)
        , motor_torque(0.0f) {}
    
    HingeConstraint(RigidBody* a, RigidBody* b, const Vec3& anchor_a, const Vec3& anchor_b)
        : Constraint(ConstraintType::Hinge)
        , local_anchor_a(anchor_a)
        , local_anchor_b(anchor_b)
        , use_limits(false)
        , lower_limit(-M_PI)
        , upper_limit(M_PI)
        , use_motor(false)
        , motor_velocity(0.0f)
        , motor_torque(0.0f)
    {
        body_a = a;
        body_b = b;
        local_axis_a = Vec3::unit_y();
        local_axis_b = Vec3::unit_y();
    }
    
    void pre_solve(float dt) override;
    void solve() override;
    void post_solve() override;
};

/**
 * @brief Ball-Socket Constraint (spherical joint)
 */
class BallSocketConstraint : public Constraint {
public:
    Vec3 local_anchor_a;
    Vec3 local_anchor_b;
    
    BallSocketConstraint() : Constraint(ConstraintType::BallSocket) {}
    
    BallSocketConstraint(RigidBody* a, RigidBody* b, const Vec3& anchor_a, const Vec3& anchor_b)
        : Constraint(ConstraintType::BallSocket)
        , local_anchor_a(anchor_a)
        , local_anchor_b(anchor_b)
    {
        body_a = a;
        body_b = b;
    }
    
    void solve() override;
};

/**
 * @brief Fixed Constraint (welds bodies together)
 */
class FixedConstraint : public Constraint {
public:
    Vec3 local_anchor_a;
    Vec3 local_anchor_b;
    Quaternion local_rotation_a;
    Quaternion local_rotation_b;
    
    FixedConstraint() : Constraint(ConstraintType::Fixed) {}
    
    FixedConstraint(RigidBody* a, RigidBody* b, const Vec3& anchor_a, const Vec3& anchor_b)
        : Constraint(ConstraintType::Fixed)
        , local_anchor_a(anchor_a)
        , local_anchor_b(anchor_b)
    {
        body_a = a;
        body_b = b;
    }
    
    void solve() override;
};

/**
 * @brief Prismatic Constraint (slider joint)
 */
class PrismaticConstraint : public Constraint {
public:
    Vec3 local_anchor_a;
    Vec3 local_anchor_b;
    Vec3 local_axis_a;
    Vec3 local_axis_b;
    
    bool use_limits;
    float lower_limit;
    float upper_limit;
    
    // Motor
    bool use_motor;
    float motor_force;
    
    PrismaticConstraint()
        : Constraint(ConstraintType::Prismatic)
        , use_limits(false)
        , lower_limit(-10.0f)
        , upper_limit(10.0f)
        , use_motor(false)
        , motor_force(100.0f)
    {}
    
    PrismaticConstraint(RigidBody* a, RigidBody* b, const Vec3& anchor_a, 
                       const Vec3& anchor_b, const Vec3& axis)
        : Constraint(ConstraintType::Prismatic)
        , local_anchor_a(anchor_a)
        , local_anchor_b(anchor_b)
        , local_axis_a(axis)
        , local_axis_b(axis)
        , use_limits(false)
        , lower_limit(-10.0f)
        , upper_limit(10.0f)
        , use_motor(false)
        , motor_force(100.0f)
    {
        body_a = a;
        body_b = b;
    }
    
    void pre_solve(float dt) override;
    void solve() override;
};

/**
 * @brief Cone Constraint (limited angular movement)
 */
class ConeConstraint : public Constraint {
public:
    Vec3 local_anchor_a;
    Vec3 local_anchor_b;
    float cone_angle;
    float twist_angle;
    
    ConeConstraint()
        : Constraint(ConstraintType::Cone)
        , cone_angle(M_PI / 4.0f)
        , twist_angle(M_PI / 4.0f) {}
    
    ConeConstraint(RigidBody* a, RigidBody* b, const Vec3& anchor_a, 
                  const Vec3& anchor_b, float cone, float twist)
        : Constraint(ConstraintType::Cone)
        , local_anchor_a(anchor_a)
        , local_anchor_b(anchor_b)
        , cone_angle(cone)
        , twist_angle(twist)
    {
        body_a = a;
        body_b = b;
    }
    
    void solve() override;
};

/**
 * @brief Constraint Solver World
 */
class ConstraintWorld {
public:
    std::vector<std::shared_ptr<Constraint>> constraints;
    int iterations;
    float tolerance;
    
    ConstraintWorld() : iterations(10), tolerance(1e-4f) {}
    
    std::shared_ptr<DistanceConstraint> add_distance_constraint(
        RigidBody* a, RigidBody* b, 
        const Vec3& anchor_a, const Vec3& anchor_b,
        float distance = -1.0f);
    
    std::shared_ptr<HingeConstraint> add_hinge_constraint(
        RigidBody* a, RigidBody* b,
        const Vec3& anchor_a, const Vec3& anchor_b);
    
    std::shared_ptr<BallSocketConstraint> add_ball_socket_constraint(
        RigidBody* a, RigidBody* b,
        const Vec3& anchor_a, const Vec3& anchor_b);
    
    std::shared_ptr<FixedConstraint> add_fixed_constraint(
        RigidBody* a, RigidBody* b,
        const Vec3& anchor_a, const Vec3& anchor_b);
    
    std::shared_ptr<PrismaticConstraint> add_prismatic_constraint(
        RigidBody* a, RigidBody* b,
        const Vec3& anchor_a, const Vec3& anchor_b,
        const Vec3& axis);
    
    void remove_constraint(Constraint* constraint);
    void clear();
    
    void solve(float dt);
    
private:
    void pre_solve_all(float dt);
    void solve_iteration();
    void post_solve_all();
};

} // namespace physics
} // namespace atlas
