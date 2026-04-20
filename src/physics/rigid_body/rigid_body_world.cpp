/**
 * @file rigid_body_world.cpp
 * @brief Rigid Body World Implementation
 */

#include "rigid_body_world.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace atlas {
namespace physics {

RigidBodyWorld::RigidBodyWorld(const WorldConfig& cfg)
    : config(cfg) {}

std::shared_ptr<RigidBody> RigidBodyWorld::add_body(const std::shared_ptr<RigidBody>& body) {
    body->id = next_body_id++;
    bodies.push_back(body);
    body_id_map[body->id] = bodies.size() - 1;
    return body;
}

void RigidBodyWorld::remove_body(uint32_t id) {
    auto it = body_id_map.find(id);
    if (it != body_id_map.end()) {
        size_t index = it->second;
        bodies.erase(bodies.begin() + index);
        // Update map indices
        for (auto& pair : body_id_map) {
            if (pair.second > index) {
                pair.second--;
            }
        }
        body_id_map.erase(it);
    }
}

void RigidBodyWorld::remove_body(const std::shared_ptr<RigidBody>& body) {
    remove_body(body->id);
}

RigidBody* RigidBodyWorld::get_body(uint32_t id) {
    auto it = body_id_map.find(id);
    if (it != body_id_map.end()) {
        return bodies[it->second].get();
    }
    return nullptr;
}

void RigidBodyWorld::clear() {
    bodies.clear();
    body_id_map.clear();
    manifolds.clear();
    next_body_id = 0;
}

void RigidBodyWorld::set_gravity(const Vec3& gravity) {
    config.gravity = gravity;
}

void RigidBodyWorld::set_collision_callback(
    std::function<void(RigidBody&, RigidBody&, const ContactManifold&)> callback
) {
    collision_callback = callback;
}

void RigidBodyWorld::step(float dt) {
    dt = std::min(dt, config.max_dt);
    
    float sub_dt = dt / static_cast<float>(config.substeps);
    for (int i = 0; i < config.substeps; ++i) {
        step_fixed(sub_dt);
    }
}

void RigidBodyWorld::step_fixed(float dt) {
    // Apply gravity and external forces
    for (auto& body : bodies) {
        if (body->has_finite_mass() && !body->is_sleeping) {
            body->apply_gravity(config.gravity);
        }
    }

    // Integrate velocities
    integrate_bodies(dt);

    // Broad phase collision detection
    broadphase();

    // Narrow phase collision detection
    narrowphase();

    // Warm start contacts
    warm_start_contacts();

    // Solve velocity constraints
    solve_velocity_constraints(dt);

    // Integrate positions
    for (auto& body : bodies) {
        if (body->has_finite_mass() && !body->is_sleeping) {
            body->integrate_positions(dt);
        }
    }

    // Solve position constraints ( Baumgarte stabilization )
    solve_position_constraints(dt);

    // Clear forces
    for (auto& body : bodies) {
        body->clear_accumulators();
    }

    // Update sleep states
    if (config.allow_sleep) {
        update_sleep_states(dt);
    }

    // Clamp velocities
    for (auto& body : bodies) {
        clamp_velocity(body.get());
    }
}

void RigidBodyWorld::integrate_bodies(float dt) {
    for (auto& body : bodies) {
        if (body->has_finite_mass() && !body->is_sleeping) {
            body->integrate_velocities(dt);
        }
    }
}

void RigidBodyWorld::broadphase() {
    broadphase_pairs.clear();
    
    // Simple O(n²) broadphase with AABB test
    // For large worlds, use SAP or BVH
    for (size_t i = 0; i < bodies.size(); ++i) {
        for (size_t j = i + 1; j < bodies.size(); ++j) {
            RigidBody* a = bodies[i].get();
            RigidBody* b = bodies[j].get();
            
            // Skip if both are static/kinematic
            if (a->state == BodyState::Static && b->state == BodyState::Static) {
                continue;
            }
            
            // Check collision filtering
            if ((a->collision_mask & b->collision_group) == 0 ||
                (b->collision_mask & a->collision_group) == 0) {
                continue;
            }
            
            // AABB intersection test
            if (aabb_intersect(a, b)) {
                broadphase_pairs.push_back({a, b});
            }
        }
    }
}

bool RigidBodyWorld::aabb_intersect(const RigidBody* a, const RigidBody* b) const {
    Vec3 a_min = a->get_aabb_min();
    Vec3 a_max = a->get_aabb_max();
    Vec3 b_min = b->get_aabb_min();
    Vec3 b_max = b->get_aabb_max();

    return (a_min.x <= b_max.x && a_max.x >= b_min.x) &&
           (a_min.y <= b_max.y && a_max.y >= b_min.y) &&
           (a_min.z <= b_max.z && a_max.z >= b_min.z);
}

void RigidBodyWorld::narrowphase() {
    manifolds.clear();
    
    for (auto& pair : broadphase_pairs) {
        ContactManifold manifold;
        if (generate_contacts(pair.body_a, pair.body_b, manifold)) {
            manifolds.push_back(manifold);
            
            // Call collision callback
            if (collision_callback) {
                collision_callback(*pair.body_a, *pair.body_b, manifold);
            }
        }
    }
}

bool RigidBodyWorld::generate_contacts(RigidBody* a, RigidBody* b, ContactManifold& manifold) {
    manifold.body_a = a;
    manifold.body_b = b;
    manifold.contacts.clear();

    ShapeType type_a = a->shape ? a->shape->get_type() : ShapeType::Sphere;
    ShapeType type_b = b->shape ? b->shape->get_type() : ShapeType::Sphere;

    // Dispatch to appropriate collision function
    if (type_a == ShapeType::Sphere && type_b == ShapeType::Sphere) {
        return sphere_vs_sphere(a, b, manifold);
    } else if (type_a == ShapeType::Box && type_b == ShapeType::Box) {
        return box_vs_box(a, b, manifold);
    } else if (type_a == ShapeType::Sphere && type_b == ShapeType::Box) {
        return sphere_vs_box(a, b, manifold);
    } else if (type_a == ShapeType::Box && type_b == ShapeType::Sphere) {
        return sphere_vs_box(b, a, manifold);
    } else if (type_a == ShapeType::Sphere && type_b == ShapeType::Plane) {
        return sphere_vs_plane(a, b, manifold);
    } else if (type_a == ShapeType::Box && type_b == ShapeType::Plane) {
        return box_vs_plane(a, b, manifold);
    } else if (type_b == ShapeType::Plane && type_a == ShapeType::Sphere) {
        return sphere_vs_plane(b, a, manifold);
    } else if (type_b == ShapeType::Plane && type_a == ShapeType::Box) {
        return box_vs_plane(b, a, manifold);
    }

    return false;
}

bool RigidBodyWorld::sphere_vs_sphere(RigidBody* a, RigidBody* b, ContactManifold& manifold) {
    SphereShape* sphere_a = static_cast<SphereShape*>(a->shape.get());
    SphereShape* sphere_b = static_cast<SphereShape*>(b->shape.get());

    Vec3 pos_a = a->transform.position;
    Vec3 pos_b = b->transform.position;
    
    Vec3 delta = pos_b - pos_a;
    float dist_sq = delta.magnitude_squared();
    float radius_sum = sphere_a->radius + sphere_b->radius;

    if (dist_sq > radius_sum * radius_sum) {
        return false;
    }

    float dist = std::sqrt(dist_sq);
    
    ContactPoint contact;
    if (dist < 1e-6f) {
        // Spheres are at same position
        contact.position = pos_a;
        contact.normal = Vec3::unit_y();
        contact.penetration = radius_sum;
    } else {
        contact.normal = delta / dist;
        contact.position = pos_a + contact.normal * sphere_a->radius;
        contact.penetration = radius_sum - dist;
    }

    // Combine material properties
    contact.restitution = std::min(a->material.restitution, b->material.restitution);
    contact.friction = std::sqrt(a->material.friction * b->material.friction);
    contact.normal_impulse = Vec3::zero();
    contact.tangent_impulse = Vec3::zero();

    manifold.contacts.push_back(contact);
    manifold.contact_count = 1;
    manifold.normal = contact.normal;

    return true;
}

bool RigidBodyWorld::box_vs_box(RigidBody* a, RigidBody* b, ContactManifold& manifold) {
    BoxShape* box_a = static_cast<BoxShape*>(a->shape.get());
    BoxShape* box_b = static_cast<BoxShape*>(b->shape.get());

    // Use SAT for box-box collision
    std::vector<Vec3> axes = get_box_axes(a);
    std::vector<Vec3> axes_b = get_box_axes(b);
    axes.insert(axes.end(), axes_b.begin(), axes_b.end());

    Vec3 normal;
    float penetration;

    if (!sat_test(axes, a, b, normal, penetration)) {
        return false;
    }

    // Find contact points (simplified - just one)
    ContactPoint contact;
    contact.normal = normal;
    contact.penetration = penetration;
    contact.restitution = std::min(a->material.restitution, b->material.restitution);
    contact.friction = std::sqrt(a->material.friction * b->material.friction);
    contact.position = (a->transform.position + b->transform.position) * 0.5f;
    contact.normal_impulse = Vec3::zero();
    contact.tangent_impulse = Vec3::zero();

    manifold.contacts.push_back(contact);
    manifold.contact_count = 1;

    return true;
}

bool RigidBodyWorld::sphere_vs_box(RigidBody* sphere, RigidBody* box, ContactManifold& manifold) {
    SphereShape* sphere_shape = static_cast<SphereShape*>(sphere->shape.get());
    BoxShape* box_shape = static_cast<BoxShape*>(box->shape.get());

    // Transform sphere center to box local space
    Vec3 local_sphere = sphere->transform.world_to_local(sphere->transform.position);
    Vec3 local_extents = box_shape->half_extents;

    // Find closest point on box to sphere center
    Vec3 closest = local_sphere;
    closest = closest.component_mul(local_extents);
    closest = closest.component_div(local_extents);
    closest = closest.component_mul(local_extents);

    Vec3 delta = local_sphere - closest;
    float dist_sq = delta.magnitude_squared();

    if (dist_sq > sphere_shape->radius * sphere_shape->radius) {
        return false;
    }

    float dist = std::sqrt(dist_sq);
    
    ContactPoint contact;
    if (dist < 1e-6f) {
        // Sphere center is inside box
        Vec3 dist_to_surface = local_extents.component_mul(Vec3::one()) - local_sphere.abs();
        float min_dist = std::min(dist_to_surface.x, 
                                 std::min(dist_to_surface.y, dist_to_surface.z));
        
        if (min_dist == dist_to_surface.x) {
            contact.normal = Vec3(local_sphere.x > 0 ? 1.0f : -1.0f, 0, 0);
        } else if (min_dist == dist_to_surface.y) {
            contact.normal = Vec3(0, local_sphere.y > 0 ? 1.0f : -1.0f, 0);
        } else {
            contact.normal = Vec3(0, 0, local_sphere.z > 0 ? 1.0f : -1.0f);
        }
        
        contact.penetration = sphere_shape->radius + min_dist;
        closest = local_sphere + contact.normal * (local_extents.component_mul(Vec3::one()) - contact.normal * min_dist);
    } else {
        contact.normal = delta / dist;
        contact.penetration = sphere_shape->radius - dist;
        closest = local_sphere - contact.normal * dist;
    }

    // Transform to world space
    contact.position = box->transform.local_to_world(closest);
    contact.normal = box->transform.local_to_world_direction(contact.normal);
    contact.restitution = std::min(sphere->material.restitution, box->material.restitution);
    contact.friction = std::sqrt(sphere->material.friction * box->material.friction);
    contact.normal_impulse = Vec3::zero();
    contact.tangent_impulse = Vec3::zero();

    manifold.contacts.push_back(contact);
    manifold.contact_count = 1;

    return true;
}

bool RigidBodyWorld::sphere_vs_plane(RigidBody* sphere, RigidBody* plane, ContactManifold& manifold) {
    SphereShape* sphere_shape = static_cast<SphereShape*>(sphere->shape.get());
    PlaneShape* plane_shape = static_cast<PlaneShape*>(plane->shape.get());

    Vec3 sphere_pos = sphere->transform.position;
    Vec3 plane_normal = plane_shape->normal;

    // Distance from sphere center to plane
    float dist = (sphere_pos - plane_normal * plane_shape->offset).dot(plane_normal);

    if (dist > sphere_shape->radius) {
        return false;
    }

    ContactPoint contact;
    contact.position = sphere_pos - plane_normal * dist;
    contact.normal = dist > 0.0f ? plane_normal : -plane_normal;
    contact.penetration = sphere_shape->radius - dist;
    contact.restitution = std::min(sphere->material.restitution, plane->material.restitution);
    contact.friction = std::sqrt(sphere->material.friction * plane->material.friction);
    contact.normal_impulse = Vec3::zero();
    contact.tangent_impulse = Vec3::zero();

    manifold.contacts.push_back(contact);
    manifold.contact_count = 1;
    manifold.normal = contact.normal;

    return true;
}

bool RigidBodyWorld::box_vs_plane(RigidBody* box, RigidBody* plane, ContactManifold& manifold) {
    BoxShape* box_shape = static_cast<BoxShape*>(box->shape.get());
    PlaneShape* plane_shape = static_cast<PlaneShape*>(plane->shape.get());

    // Get box corners in world space
    Vec3 corners[8] = {
        Vec3(-box_shape->half_extents.x, -box_shape->half_extents.y, -box_shape->half_extents.z),
        Vec3(box_shape->half_extents.x, -box_shape->half_extents.y, -box_shape->half_extents.z),
        Vec3(-box_shape->half_extents.x, box_shape->half_extents.y, -box_shape->half_extents.z),
        Vec3(box_shape->half_extents.x, box_shape->half_extents.y, -box_shape->half_extents.z),
        Vec3(-box_shape->half_extents.x, -box_shape->half_extents.y, box_shape->half_extents.z),
        Vec3(box_shape->half_extents.x, -box_shape->half_extents.y, box_shape->half_extents.z),
        Vec3(-box_shape->half_extents.x, box_shape->half_extents.y, box_shape->half_extents.z),
        Vec3(box_shape->half_extents.x, box_shape->half_extents.y, box_shape->half_extents.z)
    };

    Vec3 plane_normal = plane_shape->normal;
    Vec3 max_corner = box->transform.position;
    float max_penetration = -INFINITY;
    Vec3 contact_point = box->transform.position;

    bool has_contact = false;

    for (int i = 0; i < 8; ++i) {
        Vec3 world_corner = box->transform.transform_point(corners[i]);
        float dist = (world_corner - plane_normal * plane_shape->offset).dot(plane_normal);

        if (dist < 0.0f) {
            has_contact = true;
            if (dist > max_penetration) {
                max_penetration = dist;
                contact_point = world_corner;
            }
        }
    }

    if (!has_contact) {
        return false;
    }

    ContactPoint contact;
    contact.position = contact_point;
    contact.normal = max_penetration > 0.0f ? plane_normal : -plane_normal;
    contact.penetration = -max_penetration;
    contact.restitution = std::min(box->material.restitution, plane->material.restitution);
    contact.friction = std::sqrt(box->material.friction * plane->material.friction);
    contact.normal_impulse = Vec3::zero();
    contact.tangent_impulse = Vec3::zero();

    manifold.contacts.push_back(contact);
    manifold.contact_count = 1;
    manifold.normal = contact.normal;

    return true;
}

std::vector<Vec3> RigidBodyWorld::get_box_axes(RigidBody* box) {
    std::vector<Vec3> axes;
    axes.push_back(box->transform.get_right());
    axes.push_back(box->transform.get_up());
    axes.push_back(box->transform.get_forward());
    return axes;
}

bool RigidBodyWorld::sat_test(const std::vector<Vec3>& axes, RigidBody* a, RigidBody* b,
                             Vec3& normal, float& penetration) {
    normal = Vec3::zero();
    penetration = 0.0f;
    bool first_axis = true;

    for (const Vec3& axis : axes) {
        if (axis.magnitude_squared() < 1e-6f) continue;

        Vec3 axis_norm = axis.normalized();

        float min_a = INFINITY, max_a = -INFINITY;
        float min_b = INFINITY, max_b = -INFINITY;

        // Get AABB on axis for both boxes
        BoxShape* box_a = static_cast<BoxShape*>(a->shape.get());
        BoxShape* box_b = static_cast<BoxShape*>(b->shape.get());

        Vec3 corners_a[8] = {
            Vec3(-box_a->half_extents.x, -box_a->half_extents.y, -box_a->half_extents.z),
            Vec3(box_a->half_extents.x, -box_a->half_extents.y, -box_a->half_extents.z),
            Vec3(-box_a->half_extents.x, box_a->half_extents.y, -box_a->half_extents.z),
            Vec3(box_a->half_extents.x, box_a->half_extents.y, -box_a->half_extents.z),
            Vec3(-box_a->half_extents.x, -box_a->half_extents.y, box_a->half_extents.z),
            Vec3(box_a->half_extents.x, -box_a->half_extents.y, box_a->half_extents.z),
            Vec3(-box_a->half_extents.x, box_a->half_extents.y, box_a->half_extents.z),
            Vec3(box_a->half_extents.x, box_a->half_extents.y, box_a->half_extents.z)
        };

        Vec3 corners_b[8] = {
            Vec3(-box_b->half_extents.x, -box_b->half_extents.y, -box_b->half_extents.z),
            Vec3(box_b->half_extents.x, -box_b->half_extents.y, -box_b->half_extents.z),
            Vec3(-box_b->half_extents.x, box_b->half_extents.y, -box_b->half_extents.z),
            Vec3(box_b->half_extents.x, box_b->half_extents.y, -box_b->half_extents.z),
            Vec3(-box_b->half_extents.x, -box_b->half_extents.y, box_b->half_extents.z),
            Vec3(box_b->half_extents.x, -box_b->half_extents.y, box_b->half_extents.z),
            Vec3(-box_b->half_extents.x, box_b->half_extents.y, box_b->half_extents.z),
            Vec3(box_b->half_extents.x, box_b->half_extents.y, box_b->half_extents.z)
        };

        for (int i = 0; i < 8; ++i) {
            Vec3 world_a = a->transform.transform_point(corners_a[i]);
            Vec3 world_b = b->transform.transform_point(corners_b[i]);

            float proj_a = world_a.dot(axis_norm);
            float proj_b = world_b.dot(axis_norm);

            min_a = std::min(min_a, proj_a);
            max_a = std::max(max_a, proj_a);
            min_b = std::min(min_b, proj_b);
            max_b = std::max(max_b, proj_b);
        }

        // Check for gap
        if (max_a < min_b || max_b < min_a) {
            return false;
        }

        float axis_penetration = std::min(max_a - min_b, max_b - min_a);
        if (first_axis || axis_penetration < penetration) {
            penetration = axis_penetration;
            normal = axis_norm;
            // Ensure normal points from a to b
            Vec3 d = b->transform.position - a->transform.position;
            if (d.dot(normal) < 0.0f) {
                normal = -normal;
            }
        }

        first_axis = false;
    }

    return true;
}

bool RigidBodyWorld::gjk_collision(RigidBody* a, RigidBody* b, Vec3& simplex_out) {
    // Simplified GJK - full implementation would be more extensive
    Vec3 direction = b->transform.position - a->transform.position;
    if (direction.magnitude_squared() < 1e-6f) {
        direction = Vec3::unit_x();
    }

    simplex_out = direction;
    return true; // Simplified
}

void RigidBodyWorld::solve_contacts(float dt) {
    for (int iter = 0; iter < config.solver.max_iterations; ++iter) {
        for (auto& manifold : manifolds) {
            for (auto& contact : manifold.contacts) {
                resolve_contact(manifold, dt);
            }
        }
    }
}

void RigidBodyWorld::solve_velocity_constraints(float dt) {
    for (auto& manifold : manifolds) {
        for (auto& contact : manifold.contacts) {
            apply_contact_impulse(manifold, contact, dt);
        }
    }
}

void RigidBodyWorld::solve_position_constraints(float dt) {
    for (auto& manifold : manifolds) {
        apply_position_correction(manifold);
    }
}

void RigidBodyWorld::resolve_contact(ContactManifold& manifold, float dt) {
    RigidBody* a = manifold.body_a;
    RigidBody* b = manifold.body_b;

    if (!a->has_finite_mass() && !b->has_finite_mass()) return;

    for (auto& contact : manifold.contacts) {
        // Calculate relative velocity at contact point
        Vec3 r_a = contact.position - a->transform.position;
        Vec3 r_b = contact.position - b->transform.position;

        Vec3 vel_a = a->linear_velocity + a->angular_velocity.cross(r_a);
        Vec3 vel_b = b->linear_velocity + b->angular_velocity.cross(r_b);
        Vec3 rel_vel = vel_b - vel_a;

        float vel_along_normal = rel_vel.dot(contact.normal);

        // Don't resolve if velocities are separating
        if (vel_along_normal > 0.0f) continue;

        // Calculate impulse scalar
        float e = contact.restitution;
        float j = -(1.0f + e) * vel_along_normal;

        // Calculate mass contributions
        float inv_mass_sum = a->inv_mass + b->inv_mass;

        // Add rotational contribution
        Vec3 cross_a = r_a.cross(contact.normal);
        Vec3 cross_b = r_b.cross(contact.normal);
        
        float rot_factor_a = cross_a.dot(cross_a.component_mul(a->inv_inertia));
        float rot_factor_b = cross_b.dot(cross_b.component_mul(b->inv_inertia));
        
        j /= inv_mass_sum + rot_factor_a + rot_factor_b;

        // Apply normal impulse
        Vec3 impulse = contact.normal * j;
        a->apply_impulse_at_point(-impulse, contact.position);
        b->apply_impulse_at_point(impulse, contact.position);

        contact.normal_impulse += impulse;

        // Friction impulse
        Vec3 tangent = rel_vel - contact.normal * vel_along_normal;
        float tangent_mag = tangent.magnitude();
        
        if (tangent_mag > 1e-6f) {
            tangent = tangent / tangent_mag;
            
            float jt = -tangent_mag / inv_mass_sum;
            
            // Coulomb friction
            float mu = contact.friction;
            jt = std::clamp(jt, -mu * j, mu * j);

            Vec3 friction_impulse = tangent * jt;
            a->apply_impulse_at_point(-friction_impulse, contact.position);
            b->apply_impulse_at_point(friction_impulse, contact.position);

            contact.tangent_impulse += friction_impulse;
        }
    }
}

void RigidBodyWorld::apply_contact_impulse(ContactManifold& manifold, ContactPoint& contact, float dt) {
    // Velocity-based impulse application (simplified)
}

void RigidBodyWorld::warm_start_contacts() {
    for (auto& manifold : manifolds) {
        RigidBody* a = manifold.body_a;
        RigidBody* b = manifold.body_b;

        for (auto& contact : manifold.contacts) {
            Vec3 r_a = contact.position - a->transform.position;
            Vec3 r_b = contact.position - b->transform.position;

            // Apply cached impulses
            Vec3 total_impulse = contact.normal_impulse + contact.tangent_impulse;
            
            a->linear_velocity -= total_impulse * a->inv_mass;
            a->angular_velocity -= r_a.cross(total_impulse).component_mul(a->inv_inertia);
            
            b->linear_velocity += total_impulse * b->inv_mass;
            b->angular_velocity += r_b.cross(total_impulse).component_mul(b->inv_inertia);
        }
    }
}

void RigidBodyWorld::apply_position_correction(ContactManifold& manifold) {
    RigidBody* a = manifold.body_a;
    RigidBody* b = manifold.body_b;

    float total_penetration = 0.0f;
    Vec3 total_correction = Vec3::zero();

    for (auto& contact : manifold.contacts) {
        if (contact.penetration > config.solver.slop) {
            float correction = config.solver.baumgarte_scale * 
                              (contact.penetration - config.solver.slop);
            total_correction += contact.normal * correction;
            total_penetration += contact.penetration;
        }
    }

    if (total_penetration > 1e-6f) {
        Vec3 correction = total_correction / (a->inv_mass + b->inv_mass);
        
        if (a->has_finite_mass()) {
            a->transform.position -= correction * a->inv_mass;
        }
        if (b->has_finite_mass()) {
            b->transform.position += correction * b->inv_mass;
        }
    }
}

void RigidBodyWorld::update_sleep_states(float dt) {
    for (auto& body : bodies) {
        body->update_sleep_state(dt);
    }
}

void RigidBodyWorld::clamp_velocity(RigidBody* body) {
    float speed = body->linear_velocity.magnitude();
    if (speed > config.solver.max_velocity) {
        body->linear_velocity *= config.solver.max_velocity / speed;
    }

    float ang_speed = body->angular_velocity.magnitude();
    if (ang_speed > config.solver.max_angular_velocity) {
        body->angular_velocity *= config.solver.max_angular_velocity / ang_speed;
    }
}

void RigidBodyWorld::query_point(const Vec3& point, std::vector<RigidBody*>& results) {
    for (auto& body : bodies) {
        Vec3 min_point = body->get_aabb_min();
        Vec3 max_point = body->get_aabb_max();

        if (point.x >= min_point.x && point.x <= max_point.x &&
            point.y >= min_point.y && point.y <= max_point.y &&
            point.z >= min_point.z && point.z <= max_point.z) {
            results.push_back(body.get());
        }
    }
}

void RigidBodyWorld::query_aabb(const Vec3& min, const Vec3& max, std::vector<RigidBody*>& results) {
    for (auto& body : bodies) {
        Vec3 body_min = body->get_aabb_min();
        Vec3 body_max = body->get_aabb_max();

        if (min.x <= body_max.x && max.x >= body_min.x &&
            min.y <= body_max.y && max.y >= body_min.y &&
            min.z <= body_max.z && max.z >= body_min.z) {
            results.push_back(body.get());
        }
    }
}

void RigidBodyWorld::query_ray(const Vec3& origin, const Vec3& direction, float max_distance,
                                std::vector<RaycastHit>& hits) {
    for (auto& body : bodies) {
        Vec3 hit_point, hit_normal;
        float hit_fraction;

        if (body->shape && body->shape->raycast(origin, direction, max_distance,
                                                 hit_point, hit_normal, hit_fraction)) {
            RaycastHit hit;
            hit.body = body.get();
            hit.point = hit_point;
            hit.normal = hit_normal;
            hit.fraction = hit_fraction;
            hits.push_back(hit);
        }
    }

    // Sort by fraction
    std::sort(hits.begin(), hits.end(), [](const RaycastHit& a, const RaycastHit& b) {
        return a.fraction < b.fraction;
    });
}

std::vector<std::pair<Vec3, Vec3>> RigidBodyWorld::get_debug_lines() {
    std::vector<std::pair<Vec3, Vec3>> lines;

    // Draw contact manifolds
    for (auto& manifold : manifolds) {
        for (auto& contact : manifold.contacts) {
            lines.push_back({contact.position, contact.position + manifold.normal * 0.2f});
        }
    }

    return lines;
}

std::vector<Vec3> RigidBodyWorld::get_debug_points() {
    std::vector<Vec3> points;

    for (auto& manifold : manifolds) {
        for (auto& contact : manifold.contacts) {
            points.push_back(contact.position);
        }
    }

    return points;
}

} // namespace physics
} // namespace atlas
