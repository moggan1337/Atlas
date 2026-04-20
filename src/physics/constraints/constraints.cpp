/**
 * @file constraints.cpp
 * @brief Constraint Implementation
 */

#include "constraints.hpp"
#include <algorithm>
#include <cmath>

namespace atlas {
namespace physics {

void DistanceConstraint::solve() {
    if (!enabled) return;
    if (!body_a->has_finite_mass() && !body_b->has_finite_mass()) return;
    
    // Get world anchors
    Vec3 world_anchor_a = body_a->transform.transform_point(local_anchor_a);
    Vec3 world_anchor_b = body_b->transform.transform_point(local_anchor_b);
    
    Vec3 delta = world_anchor_b - world_anchor_a;
    float current_distance = delta.magnitude();
    
    if (current_distance < 1e-6f) return;
    
    Vec3 normal = delta / current_distance;
    
    // Position correction
    float error = current_distance - target_distance;
    float correction = error * stiffness;
    
    Vec3 impulse = normal * correction;
    
    // Apply impulse
    float inv_mass_a = body_a->inv_mass;
    float inv_mass_b = body_b->inv_mass;
    float total_inv_mass = inv_mass_a + inv_mass_b;
    
    if (total_inv_mass < 1e-6f) return;
    
    if (body_a->has_finite_mass()) {
        body_a->transform.position += normal * (correction * inv_mass_a / total_inv_mass);
    }
    if (body_b->has_finite_mass()) {
        body_b->transform.position -= normal * (correction * inv_mass_b / total_inv_mass);
    }
    
    // Velocity damping
    Vec3 vel_a = body_a->get_world_point_velocity(world_anchor_a);
    Vec3 vel_b = body_b->get_world_point_velocity(world_anchor_b);
    Vec3 rel_vel = vel_b - vel_a;
    float vel_along_normal = rel_vel.dot(normal);
    
    Vec3 damp_impulse = normal * (-vel_along_normal * damping);
    
    if (body_a->has_finite_mass()) {
        body_a->apply_impulse_at_point(-damp_impulse, world_anchor_a);
    }
    if (body_b->has_finite_mass()) {
        body_b->apply_impulse_at_point(damp_impulse, world_anchor_b);
    }
}

void HingeConstraint::pre_solve(float dt) {
    // Apply motor torque if enabled
    if (use_motor && body_a->has_finite_mass()) {
        Vec3 axis = body_a->transform.rotate(local_axis_a);
        body_a->apply_torque(axis * motor_torque);
    }
}

void HingeConstraint::solve() {
    if (!enabled) return;
    
    // Get world anchors
    Vec3 world_anchor_a = body_a->transform.transform_point(local_anchor_a);
    Vec3 world_anchor_b = body_b->transform.transform_point(local_anchor_b);
    
    // Position constraint - anchors must coincide
    Vec3 pos_error = world_anchor_b - world_anchor_a;
    
    // Calculate relative velocity
    Vec3 vel_a = body_a->get_world_point_velocity(world_anchor_a);
    Vec3 vel_b = body_b->get_world_point_velocity(world_anchor_b);
    Vec3 rel_vel = vel_b - vel_a;
    
    // Apply position correction
    float pos_correction = 0.9f;
    float baumgarte = 0.2f;
    
    Vec3 correction = pos_error * pos_correction + rel_vel * baumgarte * dt;
    
    float inv_mass_a = body_a->inv_mass;
    float inv_mass_b = body_b->inv_mass;
    float total_mass = inv_mass_a + inv_mass_b;
    
    if (total_mass < 1e-6f) return;
    
    if (body_a->has_finite_mass()) {
        body_a->transform.position += correction * (inv_mass_a / total_mass);
    }
    if (body_b->has_finite_mass()) {
        body_b->transform.position -= correction * (inv_mass_b / total_mass);
    }
    
    // Angular constraint - axis alignment
    Vec3 world_axis_a = body_a->transform.rotate(local_axis_a);
    Vec3 world_axis_b = body_b->transform.rotate(local_axis_b);
    
    // Get perpendicular axes
    Vec3 perp_a = world_axis_a.cross(Vec3::unit_x()).normalized();
    if (perp_a.magnitude_squared() < 0.01f) {
        perp_a = world_axis_a.cross(Vec3::unit_y()).normalized();
    }
    Vec3 perp_b = world_axis_b.cross(Vec3::unit_x()).normalized();
    if (perp_b.magnitude_squared() < 0.01f) {
        perp_b = world_axis_b.cross(Vec3::unit_y()).normalized();
    }
    
    // Constrain rotation around hinge axis
    float axis_error_a = world_axis_a.dot(perp_b);
    float axis_error_b = world_axis_b.dot(perp_a);
    
    // Angular correction
    if (body_a->has_finite_mass()) {
        Vec3 angular_correction = perp_a * (axis_error_a * baumgarte);
        body_a->angular_velocity -= angular_correction.component_mul(body_a->inv_inertia) * 0.1f;
    }
    if (body_b->has_finite_mass()) {
        Vec3 angular_correction = perp_b * (axis_error_b * baumgarte);
        body_b->angular_velocity -= angular_correction.component_mul(body_b->inv_inertia) * 0.1f;
    }
}

void HingeConstraint::post_solve() {
    // Limit angular movement if enabled
    if (!use_limits) return;
    
    // This would check and enforce joint limits
}

void BallSocketConstraint::solve() {
    if (!enabled) return;
    
    Vec3 world_anchor_a = body_a->transform.transform_point(local_anchor_a);
    Vec3 world_anchor_b = body_b->transform.transform_point(local_anchor_b);
    
    Vec3 pos_error = world_anchor_b - world_anchor_a;
    float error_mag = pos_error.magnitude();
    
    if (error_mag < 1e-6f) return;
    
    Vec3 correction = pos_error * 0.9f;
    
    float inv_mass_a = body_a->inv_mass;
    float inv_mass_b = body_b->inv_mass;
    float total_mass = inv_mass_a + inv_mass_b;
    
    if (total_mass < 1e-6f) return;
    
    if (body_a->has_finite_mass()) {
        body_a->transform.position += correction * (inv_mass_a / total_mass);
    }
    if (body_b->has_finite_mass()) {
        body_b->transform.position -= correction * (inv_mass_b / total_mass);
    }
    
    // Velocity damping
    Vec3 vel_a = body_a->get_world_point_velocity(world_anchor_a);
    Vec3 vel_b = body_b->get_world_point_velocity(world_anchor_b);
    Vec3 rel_vel = vel_b - vel_a;
    
    Vec3 damping = -rel_vel * 0.1f;
    
    if (body_a->has_finite_mass()) {
        body_a->apply_impulse_at_point(-damping, world_anchor_a);
    }
    if (body_b->has_finite_mass()) {
        body_b->apply_impulse_at_point(damping, world_anchor_b);
    }
}

void FixedConstraint::solve() {
    if (!enabled) return;
    
    Vec3 world_anchor_a = body_a->transform.transform_point(local_anchor_a);
    Vec3 world_anchor_b = body_b->transform.transform_point(local_anchor_b);
    
    // Position correction
    Vec3 pos_error = world_anchor_b - world_anchor_a;
    
    float inv_mass_a = body_a->inv_mass;
    float inv_mass_b = body_b->inv_mass;
    float total_mass = inv_mass_a + inv_mass_b;
    
    if (total_mass < 1e-6f) return;
    
    Vec3 correction = pos_error * 0.9f;
    
    if (body_a->has_finite_mass()) {
        body_a->transform.position += correction * (inv_mass_a / total_mass);
    }
    if (body_b->has_finite_mass()) {
        body_b->transform.position -= correction * (inv_mass_b / total_mass);
    }
    
    // Rotation correction - align orientations
    Quaternion q_a = body_a->transform.rotation * local_rotation_a;
    Quaternion q_b = body_b->transform.rotation * local_rotation_b;
    
    Quaternion q_diff = q_b * q_a.inverse();
    
    // Extract angle-axis
    float angle = std::acos(std::min(1.0f, q_diff.w)) * 2.0f;
    if (angle > 1e-6f) {
        Vec3 axis = Vec3(q_diff.x, q_diff.y, q_diff.z);
        if (axis.magnitude_squared() > 1e-6f) {
            axis = axis.normalized();
            
            // Apply angular correction
            if (body_a->has_finite_mass()) {
                Vec3 correction = axis * (angle * 0.1f);
                body_a->angular_velocity -= correction.component_mul(body_a->inv_inertia);
            }
            if (body_b->has_finite_mass()) {
                Vec3 correction = axis * (angle * 0.1f);
                body_b->angular_velocity += correction.component_mul(body_b->inv_inertia);
            }
        }
    }
}

void PrismaticConstraint::pre_solve(float dt) {
    if (use_motor && body_a->has_finite_mass()) {
        Vec3 axis = body_a->transform.rotate(local_axis_a);
        body_a->apply_force(axis * motor_force, 
                           body_a->transform.transform_point(local_anchor_a));
    }
}

void PrismaticConstraint::solve() {
    if (!enabled) return;
    
    Vec3 world_anchor_a = body_a->transform.transform_point(local_anchor_a);
    Vec3 world_anchor_b = body_b->transform.transform_point(local_anchor_b);
    Vec3 world_axis_a = body_a->transform.rotate(local_axis_a);
    
    // Position constraint along axis
    Vec3 pos_error = world_anchor_b - world_anchor_a;
    float along_axis = pos_error.dot(world_axis_a);
    float perp_error = (pos_error - world_axis_a * along_axis).magnitude();
    
    // Constrain perpendicular movement
    Vec3 perp_correction = (pos_error - world_axis_a * along_axis) * 0.9f;
    
    float inv_mass_a = body_a->inv_mass;
    float inv_mass_b = body_b->inv_mass;
    float total_mass = inv_mass_a + inv_mass_b;
    
    if (total_mass < 1e-6f) return;
    
    if (body_a->has_finite_mass()) {
        body_a->transform.position += perp_correction * (inv_mass_a / total_mass);
    }
    if (body_b->has_finite_mass()) {
        body_b->transform.position -= perp_correction * (inv_mass_b / total_mass);
    }
    
    // Angular constraints - constrain rotation around axis
    Vec3 world_axis_b = body_b->transform.rotate(local_axis_b);
    
    // Constrain axis misalignment
    Vec3 axis_error = world_axis_a.cross(world_axis_b);
    
    if (body_a->has_finite_mass()) {
        body_a->angular_velocity -= axis_error.component_mul(body_a->inv_inertia) * 0.1f;
    }
    if (body_b->has_finite_mass()) {
        body_b->angular_velocity += axis_error.component_mul(body_b->inv_inertia) * 0.1f;
    }
}

void ConeConstraint::solve() {
    if (!enabled) return;
    
    Vec3 world_anchor_a = body_a->transform.transform_point(local_anchor_a);
    Vec3 world_anchor_b = body_b->transform.transform_point(local_anchor_b);
    
    // Position constraint
    Vec3 pos_error = world_anchor_b - world_anchor_a;
    float error_mag = pos_error.magnitude();
    
    if (error_mag > 1e-4f) {
        Vec3 correction = pos_error * 0.9f;
        
        float inv_mass_a = body_a->inv_mass;
        float inv_mass_b = body_b->inv_mass;
        float total_mass = inv_mass_a + inv_mass_b;
        
        if (total_mass > 1e-6f) {
            if (body_a->has_finite_mass()) {
                body_a->transform.position += correction * (inv_mass_a / total_mass);
            }
            if (body_b->has_finite_mass()) {
                body_b->transform.position -= correction * (inv_mass_b / total_mass);
            }
        }
    }
    
    // Cone constraint - limit Y-axis rotation
    Vec3 local_y_a = body_a->transform.rotate(Vec3::unit_y());
    Vec3 local_y_b = body_b->transform.rotate(Vec3::unit_y());
    
    float dot = local_y_a.dot(local_y_b);
    float angle = std::acos(std::clamp(dot, -1.0f, 1.0f));
    
    if (angle > cone_angle) {
        Vec3 axis = local_y_a.cross(local_y_b);
        if (axis.magnitude_squared() > 1e-6f) {
            axis.normalize();
            float excess = angle - cone_angle;
            
            if (body_a->has_finite_mass()) {
                Vec3 correction = axis * (excess * 0.1f);
                body_a->angular_velocity -= correction.component_mul(body_a->inv_inertia);
            }
            if (body_b->has_finite_mass()) {
                Vec3 correction = axis * (excess * 0.1f);
                body_b->angular_velocity += correction.component_mul(body_b->inv_inertia);
            }
        }
    }
}

// ============== ConstraintWorld Implementation ==============

std::shared_ptr<DistanceConstraint> ConstraintWorld::add_distance_constraint(
    RigidBody* a, RigidBody* b, 
    const Vec3& anchor_a, const Vec3& anchor_b,
    float distance) 
{
    auto constraint = std::make_shared<DistanceConstraint>(a, b, anchor_a, anchor_b, distance);
    
    if (distance < 0.0f) {
        Vec3 world_a = a->transform.transform_point(anchor_a);
        Vec3 world_b = b->transform.transform_point(anchor_b);
        constraint->target_distance = (world_b - world_a).magnitude();
    }
    
    constraints.push_back(constraint);
    return constraint;
}

std::shared_ptr<HingeConstraint> ConstraintWorld::add_hinge_constraint(
    RigidBody* a, RigidBody* b,
    const Vec3& anchor_a, const Vec3& anchor_b) 
{
    auto constraint = std::make_shared<HingeConstraint>(a, b, anchor_a, anchor_b);
    constraints.push_back(constraint);
    return constraint;
}

std::shared_ptr<BallSocketConstraint> ConstraintWorld::add_ball_socket_constraint(
    RigidBody* a, RigidBody* b,
    const Vec3& anchor_a, const Vec3& anchor_b) 
{
    auto constraint = std::make_shared<BallSocketConstraint>(a, b, anchor_a, anchor_b);
    constraints.push_back(constraint);
    return constraint;
}

std::shared_ptr<FixedConstraint> ConstraintWorld::add_fixed_constraint(
    RigidBody* a, RigidBody* b,
    const Vec3& anchor_a, const Vec3& anchor_b) 
{
    auto constraint = std::make_shared<FixedConstraint>(a, b, anchor_a, anchor_b);
    constraints.push_back(constraint);
    return constraint;
}

std::shared_ptr<PrismaticConstraint> ConstraintWorld::add_prismatic_constraint(
    RigidBody* a, RigidBody* b,
    const Vec3& anchor_a, const Vec3& anchor_b,
    const Vec3& axis) 
{
    auto constraint = std::make_shared<PrismaticConstraint>(a, b, anchor_a, anchor_b, axis);
    constraints.push_back(constraint);
    return constraint;
}

void ConstraintWorld::remove_constraint(Constraint* constraint) {
    constraints.erase(
        std::remove_if(constraints.begin(), constraints.end(),
            [constraint](const std::shared_ptr<Constraint>& c) { return c.get() == constraint; }),
        constraints.end()
    );
}

void ConstraintWorld::clear() {
    constraints.clear();
}

void ConstraintWorld::solve(float dt) {
    pre_solve_all(dt);
    
    for (int i = 0; i < iterations; ++i) {
        solve_iteration();
    }
    
    post_solve_all();
}

void ConstraintWorld::pre_solve_all(float dt) {
    for (auto& c : constraints) {
        if (c->enabled) {
            c->pre_solve(dt);
        }
    }
}

void ConstraintWorld::solve_iteration() {
    for (auto& c : constraints) {
        if (c->enabled) {
            c->solve();
        }
    }
}

void ConstraintWorld::post_solve_all() {
    for (auto& c : constraints) {
        if (c->enabled) {
            c->post_solve();
        }
    }
}

} // namespace physics
} // namespace atlas
