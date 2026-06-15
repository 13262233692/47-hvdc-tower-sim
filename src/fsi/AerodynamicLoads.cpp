#include "fsi/AerodynamicLoads.hpp"
#include "common/MathUtils.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>

namespace hvdc {
namespace fsi {

AerodynamicLoads::AerodynamicLoads(fem::FEModel* struct_model,
                                   fvm::FluidSolver* fluid_solver,
                                   const AerodynamicConfig& cfg)
    : struct_model_(struct_model), fluid_solver_(fluid_solver), config_(cfg)
{
    if (struct_model_) {
        conductor_ids_ = struct_model_->conductor_element_ids();
        tower_ids_ = struct_model_->tower_element_ids();
    }
}

Real AerodynamicLoads::compute_power_law_wind(Real height) const {
    Real U_ref = math::vec3_norm(Vec3{config_.wind_direction[0], 
                                       config_.wind_direction[1], 
                                       config_.wind_direction[2]});
    if (U_ref < EPS) return 0.0;
    if (height < 1.0) height = 1.0;
    if (config_.reference_height < EPS) return U_ref;
    return U_ref * std::pow(height / config_.reference_height, 
                             config_.power_law_exponent);
}

Vec3 AerodynamicLoads::compute_wind_velocity_at_point(const Vec3& point, Real time) const {
    Vec3 wind_dir_norm = math::vec3_normalize(
        {config_.wind_direction[0], config_.wind_direction[1], config_.wind_direction[2]});
    
    Real speed;
    if (config_.use_power_law_wind_profile) {
        speed = compute_power_law_wind(point[2]);
    } else {
        speed = math::vec3_norm(Vec3{config_.wind_direction[0], 
                                      config_.wind_direction[1], 
                                      config_.wind_direction[2]});
    }
    
    Vec3 velocity = math::vec3_scale(wind_dir_norm, speed);
    
    if (config_.gust_amplitude > EPS) {
        Real gust = config_.gust_amplitude * 
                    std::sin(2.0 * PI * config_.gust_frequency * time);
        Vec3 gust_vec = math::vec3_scale(wind_dir_norm, gust);
        math::vec3_add_inplace(velocity, gust_vec);
    }
    
    if (config_.compute_buffeting) {
        Real phase = time * 2.0 * PI;
        Vec3 buffeting = {
            config_.turbulence_intensity * speed * 0.1 * std::sin(phase * 3.7 + point[0]),
            config_.turbulence_intensity * speed * 0.05 * std::sin(phase * 5.2 + point[1]),
            config_.turbulence_intensity * speed * 0.03 * std::sin(phase * 7.1 + point[2])
        };
        math::vec3_add_inplace(velocity, buffeting);
    }
    
    return velocity;
}

Real AerodynamicLoads::compute_drag_coefficient(
    Real Reynolds, Real relative_roughness, Real ice_thickness_ratio) const
{
    Real Cd = config_.drag_coefficient_cylinder;
    
    if (ice_thickness_ratio > EPS) {
        Cd *= config_.ice_drag_coefficient_increase;
    }
    
    if (Reynolds < EPS) return Cd;
    
    if (Reynolds < 1.0e3) {
        Cd = 10.0 / std::max(Reynolds, 1.0);
    } else if (Reynolds < 2.0e5) {
        Cd = 1.2;
    } else if (Reynolds < 3.0e5) {
        Real t = (Reynolds - 2.0e5) / 1.0e5;
        Cd = 1.2 - t * 0.8;
    } else if (Reynolds < 3.5e6) {
        Real t = std::min(1.0, (Reynolds - 3.0e5) / 3.2e6);
        Cd = 0.4 + t * 0.3;
    } else {
        Cd = 0.7;
    }
    
    if (relative_roughness > EPS) {
        Cd += 0.3 * std::min(1.0, relative_roughness * 100.0);
    }
    
    return Cd;
}

Real AerodynamicLoads::compute_lift_coefficient(
    Real angle_of_attack, Real Reynolds, bool galloping_possible) const
{
    Real alpha = angle_of_attack;
    Real Cl = config_.lift_coefficient_cylinder;
    
    if (std::fabs(alpha) < 0.001) return Cl;
    
    Cl = 2.0 * PI * alpha;
    
    if (galloping_possible) {
        Real abs_alpha = std::fabs(alpha);
        if (abs_alpha > 0.1 && abs_alpha < 1.2) {
            Real icicle_effect = -3.0 * alpha * 
                std::exp(-8.0 * std::pow(alpha - 0.6, 2));
            Cl += icicle_effect;
        }
    }
    
    return Cl;
}

void AerodynamicLoads::resolve_force_components(
    Real Cd, Real Cl, Real dynamic_pressure,
    Real diameter, Real length,
    const Vec3& flow_dir, const Vec3& lift_dir,
    Vec3& force) const
{
    Real A_frontal = diameter * length;
    Real drag_mag = Cd * dynamic_pressure * A_frontal;
    Real lift_mag = Cl * dynamic_pressure * A_frontal;
    
    force = math::vec3_add(
        math::vec3_scale(flow_dir, drag_mag),
        math::vec3_scale(lift_dir, lift_mag)
    );
}

void AerodynamicLoads::compute_face_loads(
    const fvm::FluidState& state,
    const std::vector<std::array<Vec3, 3>>& grad_U,
    const std::vector<Vec3>& grad_p,
    std::vector<FaceAeroLoad>& face_loads,
    const std::vector<Index>& interface_faces) const
{
    face_loads.assign(fluid_solver_->grid()->num_faces(), FaceAeroLoad{});
    
    Vec3 U_inf = {config_.wind_direction[0], 
                   config_.wind_direction[1], 
                   config_.wind_direction[2]};
    Real U_inf_mag = math::vec3_norm(U_inf);
    Real q_inf = 0.5 * config_.rho_air * U_inf_mag * U_inf_mag;
    
    for (Index fi : interface_faces) {
        if (fi >= fluid_solver_->grid()->num_faces()) continue;
        const auto& face = fluid_solver_->grid()->face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0 || ci >= state.size()) continue;
        
        FaceAeroLoad load;
        Vec3 n = face.normal();
        Real A = face.area();
        
        load.pressure = state.p.gauge_pressure(ci);
        load.pressure_force = math::vec3_scale(n, -load.pressure * A);
        
        load.c_p = (q_inf > EPS) ? load.pressure / q_inf : 0.0;
        
        if (ci < static_cast<Index>(grad_U.size())) {
            Vec3 shear = {0, 0, 0};
            Real mu_eff = state.mu[ci];
            if (config_.use_rain_effects) {
                mu_eff *= (1.0 + 0.1 * config_.rain_intensity);
            }
            const auto& grad = grad_U[ci];
            for (Index i = 0; i < 3; ++i) {
                for (Index j = 0; j < 3; ++j) {
                    Real tau_ij = mu_eff * (grad[i][j] + grad[j][i]);
                    shear[i] += tau_ij * n[j];
                }
            }
            load.viscous_force = math::vec3_scale(shear, -A);
            load.shear_stress_mag = math::vec3_norm(shear);
            load.c_f = (q_inf > EPS) ? load.shear_stress_mag / q_inf : 0.0;
        }
        
        load.total_force = math::vec3_add(load.pressure_force, load.viscous_force);
        
        Vec3 U_local = state.U[ci];
        Real U_mag = math::vec3_norm(U_local);
        Real L = face.area() > EPS ? std::sqrt(face.area()) : 1.0;
        load.reynolds = (state.nu[ci] > EPS) ? U_mag * L / state.nu[ci] : 0.0;
        load.mach = (340.0 > EPS) ? U_mag / 340.0 : 0.0;
        
        face_loads[fi] = load;
    }
}

void AerodynamicLoads::compute_element_loads(
    const Vector& structural_displacement,
    const Vector& structural_velocity,
    std::vector<ElementAeroLoad>& elem_loads,
    Real time) const
{
    Index n_elems = struct_model_->num_elements();
    elem_loads.assign(n_elems, ElementAeroLoad{});
    
    for (Index ei = 0; ei < n_elems; ++ei) {
        const auto* elem = struct_model_->element(ei);
        if (!elem) continue;
        
        const auto* n1 = elem->node(0);
        const auto* n2 = elem->node(1);
        if (!n1 || !n2) continue;
        
        Vec3 x1 = n1->displaced_coords(structural_displacement);
        Vec3 x2 = n2->displaced_coords(structural_displacement);
        Vec3 x_mid = math::vec3_scale(math::vec3_add(x1, x2), 0.5);
        
        Vec3 v1 = {0, 0, 0};
        Vec3 v2 = {0, 0, 0};
        for (Index d = 0; d < 3; ++d) {
            Index gd1 = n1->dof_start() + d;
            Index gd2 = n2->dof_start() + d;
            if (gd1 < structural_velocity.size()) v1[d] = structural_velocity[gd1];
            if (gd2 < structural_velocity.size()) v2[d] = structural_velocity[gd2];
        }
        Vec3 v_mid = math::vec3_scale(math::vec3_add(v1, v2), 0.5);
        
        Vec3 wind_vel = compute_wind_velocity_at_point(x_mid, time);
        Vec3 rel_vel = math::vec3_sub(wind_vel, v_mid);
        Real rel_speed = math::vec3_norm(rel_vel);
        Vec3 flow_dir = (rel_speed > EPS) ? math::vec3_scale(rel_vel, 1.0/rel_speed) : Vec3{1,0,0};
        
        Vec3 elem_dir = math::vec3_normalize(math::vec3_sub(x2, x1));
        Vec3 vertical = {0, 0, 1.0};
        Vec3 lift_dir_candidate = math::vec3_cross(flow_dir, elem_dir);
        Vec3 lift_dir;
        if (math::vec3_norm2(lift_dir_candidate) > EPS) {
            lift_dir = math::vec3_normalize(lift_dir_candidate);
        } else {
            lift_dir = vertical;
        }
        
        Real length = elem->length_initial();
        if (auto* cond = dynamic_cast<const fem::ConductorElement*>(elem)) {
            length = cond->catenary_length();
        }
        
        Real diameter = 0.0;
        Real ice_thickness = 0.0;
        Real ice_ratio = 0.0;
        
        if (auto* beam = dynamic_cast<const fem::BeamElementNL*>(elem)) {
            if (beam->section()) {
                diameter = std::max(beam->section()->h, beam->section()->b);
                ice_thickness = beam->section()->thickness_ice;
            }
        } else if (auto* truss = dynamic_cast<const fem::TrussElementNL*>(elem)) {
            if (truss->section()) {
                diameter = truss->section()->D;
                ice_thickness = truss->section()->thickness_ice;
            }
        } else if (auto* cond = dynamic_cast<const fem::ConductorElement*>(elem)) {
            if (cond->section()) {
                diameter = std::max(cond->section()->h, cond->section()->wind_drag_diameter());
                ice_thickness = cond->section()->thickness_ice;
            }
        }
        
        if (diameter < EPS) diameter = 0.02;
        Real D_eff = diameter;
        if (ice_thickness > EPS) D_eff = diameter + 2.0 * ice_thickness;
        ice_ratio = (diameter > EPS) ? ice_thickness / diameter : 0.0;
        
        Real nu_air = config_.mu_air / config_.rho_air;
        Real Re = (nu_air > EPS) ? rel_speed * D_eff / nu_air : 0.0;
        Real roughness = (D_eff > EPS) ? config_.terrain_roughness / D_eff : 0.0;
        
        Real Cd = compute_drag_coefficient(Re, roughness, ice_ratio);
        Real alpha = find_angle_of_attack(rel_vel, elem_dir, vertical);
        Real Cl = compute_lift_coefficient(alpha, Re, ice_thickness > EPS);
        
        Vec3 vortex_force = {0, 0, 0};
        Real vortex_freq = 0.0;
        if (config_.compute_vortex_induced_vibration && rel_speed > EPS) {
            compute_vortex_shedding_force(rel_speed, D_eff, length, time, vortex_force, vortex_freq);
        }
        
        Vec3 buffeting_force = {0, 0, 0};
        if (config_.compute_buffeting && rel_speed > EPS) {
            buffeting_force = compute_buffeting_force(rel_vel, time,
                config_.turbulence_length_scale, config_.turbulence_intensity, length);
        }
        
        Real q = 0.5 * config_.rho_air * rel_speed * rel_speed;
        Vec3 force_per_len;
        resolve_force_components(Cd, Cl, q, D_eff, 1.0, flow_dir, lift_dir, force_per_len);
        
        math::vec3_add_inplace(force_per_len, math::vec3_scale(vortex_force, 1.0 / std::max(length, EPS)));
        math::vec3_add_inplace(force_per_len, math::vec3_scale(buffeting_force, 1.0 / std::max(length, EPS)));
        
        ElementAeroLoad elem_load;
        elem_load.element_id = ei;
        elem_load.force_per_unit_length = force_per_len;
        elem_load.drag_coeff = Cd;
        elem_load.lift_coeff = Cl;
        elem_load.reference_diameter = D_eff;
        elem_load.wind_velocity = wind_vel;
        elem_load.relative_speed = rel_speed;
        
        Vec3 moment_arm = math::vec3_cross(
            math::vec3_sub(x_mid, x1),
            force_per_len
        );
        elem_load.moment_per_unit_length = moment_arm;
        
        elem_loads[ei] = elem_load;
    }
}

Real AerodynamicLoads::find_angle_of_attack(
    const Vec3& flow_velocity, const Vec3& element_direction, const Vec3& element_normal) const
{
    Vec3 flow_dir = math::vec3_normalize(flow_velocity);
    Vec3 elem_t = math::vec3_normalize(element_direction);
    
    Real dot = math::vec3_dot(flow_dir, elem_t);
    Real cos_alpha = std::fabs(dot);
    if (cos_alpha > 1.0 - EPS) return 0.0;
    
    Real sin_alpha = std::sqrt(1.0 - cos_alpha * cos_alpha);
    
    Vec3 cross_vec = math::vec3_cross(flow_dir, elem_t);
    Real sign = math::vec3_dot(cross_vec, element_normal) > 0.0 ? 1.0 : -1.0;
    
    return sign * std::asin(sin_alpha);
}

void AerodynamicLoads::compute_conductor_aerodynamic_loads(
    const std::vector<Index>& conductor_element_ids,
    const Vector& structural_displacement,
    const Vector& structural_velocity,
    std::vector<Vec3>& nodal_loads,
    Real time) const
{
    Index nnodes = struct_model_->num_nodes();
    nodal_loads.assign(nnodes, {0, 0, 0});
    
    std::vector<ElementAeroLoad> elem_loads;
    compute_element_loads(structural_displacement, structural_velocity,
                          elem_loads, time);
    
    for (Index ei : conductor_element_ids) {
        if (ei >= static_cast<Index>(elem_loads.size())) continue;
        const auto& load = elem_loads[ei];
        const auto* elem = struct_model_->element(ei);
        if (!elem) continue;
        
        Real L = load.element_id >= 0 ? elem->length_initial() : 1.0;
        Vec3 total_force = math::vec3_scale(load.force_per_unit_length, L);
        Vec3 half_force = math::vec3_scale(total_force, 0.5);
        
        const auto* n1 = elem->node(0);
        const auto* n2 = elem->node(1);
        if (n1 && n1->id() < nnodes) {
            math::vec3_add_inplace(nodal_loads[n1->id()], half_force);
        }
        if (n2 && n2->id() < nnodes) {
            math::vec3_add_inplace(nodal_loads[n2->id()], half_force);
        }
    }
}

void AerodynamicLoads::compute_tower_aerodynamic_loads(
    const std::vector<Index>& tower_element_ids,
    const Vector& structural_displacement,
    const Vector& structural_velocity,
    std::vector<Vec3>& nodal_loads,
    Real time) const
{
    compute_conductor_aerodynamic_loads(tower_element_ids, structural_displacement,
                                         structural_velocity, nodal_loads, time);
}

void AerodynamicLoads::compute_equivalent_nodal_loads(
    const fvm::FluidState& state,
    const std::vector<std::array<Vec3, 3>>& grad_U,
    const std::vector<Vec3>& grad_p,
    Vector& struct_nodal_loads,
    Real time) const
{
    Index ndofs = struct_model_->total_dofs();
    if (struct_nodal_loads.size() != ndofs) {
        struct_nodal_loads.resize(ndofs);
    }
    struct_nodal_loads.zero();
    
    Vector disp_zeros(ndofs);
    disp_zeros.zero();
    Vector vel_zeros(ndofs);
    vel_zeros.zero();
    
    std::vector<ElementAeroLoad> elem_loads;
    compute_element_loads(disp_zeros, vel_zeros, elem_loads, time);
    
    for (Index ei = 0; ei < struct_model_->num_elements(); ++ei) {
        if (ei >= static_cast<Index>(elem_loads.size())) continue;
        const auto& load = elem_loads[ei];
        const auto* elem = struct_model_->element(ei);
        if (!elem) continue;
        
        Real L = elem->length_initial();
        Vec3 total_force = math::vec3_scale(load.force_per_unit_length, L);
        
        for (Index ni = 0; ni < elem->num_nodes(); ++ni) {
            const auto* node = elem->node(ni);
            if (!node) continue;
            Index start = node->dof_start();
            for (Index d = 0; d < 3; ++d) {
                Index gdof = start + d;
                if (gdof < ndofs && node->is_dof_active(d)) {
                    struct_nodal_loads.add(gdof, total_force[d] * 0.5);
                }
            }
        }
    }
    
    (void)state; (void)grad_U; (void)grad_p;
}

void AerodynamicLoads::compute_vortex_shedding_force(
    Real flow_speed, Real diameter, Real length,
    Real time, Vec3& force, Real& vortex_freq) const
{
    force = {0, 0, 0};
    vortex_freq = 0.0;
    if (flow_speed < EPS || diameter < EPS) return;
    
    vortex_freq = config_.strouhal_number * flow_speed / diameter;
    Real omega = 2.0 * PI * vortex_freq;
    
    Real phase = omega * time;
    Real q = 0.5 * config_.rho_air * flow_speed * flow_speed;
    Real A = diameter * length;
    
    Real lift_amplitude = 0.3 * q * A;
    Real drag_amplitude = 0.05 * q * A;
    
    Vec3 crosswind = {0.0, 1.0, 0.0};
    Vec3 inline_dir = {1.0, 0.0, 0.0};
    
    force = math::vec3_add(
        math::vec3_scale(crosswind, lift_amplitude * std::sin(phase)),
        math::vec3_scale(inline_dir, drag_amplitude * std::sin(2.0 * phase))
    );
}

Vec3 AerodynamicLoads::compute_buffeting_force(
    const Vec3& mean_velocity, Real time,
    Real length_scale, Real turbulence_intensity,
    Real element_length) const
{
    Real U_mag = math::vec3_norm(mean_velocity);
    if (U_mag < EPS || length_scale < EPS) return {0,0,0};
    
    Vec3 dir = math::vec3_scale(mean_velocity, 1.0 / U_mag);
    Vec3 side = math::vec3_normalize(math::vec3_cross(dir, {0,0,1}));
    Vec3 up = math::vec3_normalize(math::vec3_cross(dir, side));
    
    Real q = 0.5 * config_.rho_air * U_mag * U_mag;
    Real A = 1.0 * element_length;
    
    Vec3 gusts = {
        turbulence_intensity * std::sin(2.0 * PI * U_mag * time / length_scale),
        0.8 * turbulence_intensity * std::sin(2.0 * PI * U_mag * time / (length_scale * 0.7) + 1.3),
        0.5 * turbulence_intensity * std::sin(2.0 * PI * U_mag * time / (length_scale * 1.5) + 2.7)
    };
    
    Vec3 force;
    for (Index i = 0; i < 3; ++i) force[i] = 0.0;
    math::vec3_add_inplace(force, math::vec3_scale(dir, 2.0 * q * A * gusts[0]));
    math::vec3_add_inplace(force, math::vec3_scale(side, 1.5 * q * A * gusts[1]));
    math::vec3_add_inplace(force, math::vec3_scale(up, 1.0 * q * A * gusts[2]));
    
    return force;
}

void AerodynamicLoads::compute_total_forces_and_moments(
    const std::vector<ElementAeroLoad>& loads,
    Vec3& total_force, Vec3& total_moment,
    const Vec3& reference_point) const
{
    total_force = {0, 0, 0};
    total_moment = {0, 0, 0};
    
    for (const auto& load : loads) {
        const auto* elem = struct_model_->element(load.element_id);
        if (!elem) continue;
        const auto* n1 = elem->node(0);
        const auto* n2 = elem->node(1);
        if (!n1 || !n2) continue;
        
        Vec3 x1 = n1->coords();
        Vec3 x2 = n2->coords();
        Vec3 x_mid = math::vec3_scale(math::vec3_add(x1, x2), 0.5);
        Real L = elem->length_initial();
        
        Vec3 F = math::vec3_scale(load.force_per_unit_length, L);
        Vec3 M = math::vec3_add(
            math::vec3_cross(math::vec3_sub(x_mid, reference_point), F),
            math::vec3_scale(load.moment_per_unit_length, L)
        );
        
        math::vec3_add_inplace(total_force, F);
        math::vec3_add_inplace(total_moment, M);
    }
}

Real AerodynamicLoads::compute_RANS_pressure_coefficient(
    Real pressure, const Vec3& velocity, const Vec3& free_stream_velocity) const
{
    Real q_inf = 0.5 * config_.rho_air * math::vec3_norm2(free_stream_velocity);
    if (q_inf < EPS) return 0.0;
    return pressure / q_inf;
}

Real AerodynamicLoads::compute_friction_coefficient(
    Real wall_shear_stress, const Vec3& free_stream_velocity) const
{
    Real q_inf = 0.5 * config_.rho_air * math::vec3_norm2(free_stream_velocity);
    if (q_inf < EPS) return 0.0;
    return wall_shear_stress / q_inf;
}

void AerodynamicLoads::apply_ice_shape_effects(
    Real ice_thickness, Real diameter,
    Real& Cd_multiplier, Real& Cl_offset, Real& effective_diameter) const
{
    Cd_multiplier = 1.0;
    Cl_offset = 0.0;
    effective_diameter = diameter;
    
    if (ice_thickness < EPS) return;
    
    Real ratio = ice_thickness / std::max(diameter, EPS);
    
    Cd_multiplier = 1.0 + 2.5 * (1.0 - std::exp(-ratio * 8.0));
    Cl_offset = -1.5 * ratio * std::exp(-2.0 * std::pow(ratio - 0.15, 2));
    effective_diameter = diameter + 2.0 * ice_thickness;
}

Vec3 AerodynamicLoads::rotate_2D(const Vec3& vec, Real angle) const {
    Real c = std::cos(angle), s = std::sin(angle);
    return {c * vec[0] - s * vec[1], s * vec[0] + c * vec[1], vec[2]};
}

void AerodynamicLoads::compute_crosswind_forces(
    Real speed, Real diameter, Real length,
    Real angle_of_attack, Real Re,
    Vec3& lift, Vec3& moment) const
{
    lift = {0, 0, 0};
    moment = {0, 0, 0};
    if (speed < EPS || diameter < EPS) return;
    
    Real q = 0.5 * config_.rho_air * speed * speed;
    Real A = diameter * length;
    Real Cl = compute_lift_coefficient(angle_of_attack, Re, true);
    
    Vec3 lift_dir = {0, 1, 0};
    lift = math::vec3_scale(lift_dir, Cl * q * A);
    
    moment = {0, 0, 0.1 * Cl * q * A * diameter};
}

} // namespace fsi
} // namespace hvdc
