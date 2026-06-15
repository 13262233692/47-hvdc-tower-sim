#include "fsi/LoadTransfer.hpp"
#include "common/MathUtils.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>

namespace hvdc {
namespace fsi {

LoadTransfer::LoadTransfer(fem::FEModel* struct_model,
                           fvm::Grid* fluid_grid,
                           InterfaceMapper* mapper,
                           AerodynamicLoads* aero)
    : struct_model_(struct_model), fluid_grid_(fluid_grid),
      mapper_(mapper), aero_(aero)
{
}

void LoadTransfer::transfer_direct_nearest(
    const std::vector<Vec3>& source,
    Vector& target,
    bool force_to_disp) const
{
    if (!mapper_ || !struct_model_) return;
    
    Index ndofs = struct_model_->total_dofs();
    if (target.size() != ndofs) target.resize(ndofs);
    target.zero();
    
    const auto& fluid_to_struct = mapper_->fluid_to_struct_map();
    
    for (Index fi = 0; fi < fluid_to_struct.size(); ++fi) {
        if (fi >= static_cast<Index>(source.size())) continue;
        const Vec3& src_val = source[fi];
        
        for (const auto& pair : fluid_to_struct[fi]) {
            const auto* node = struct_model_->get_node(pair.struct_node_id);
            if (!node) continue;
            Index start = node->dof_start();
            for (Index d = 0; d < 3; ++d) {
                Index gdof = start + d;
                if (gdof < ndofs && node->is_dof_active(d)) {
                    target.add(gdof, pair.weight * src_val[d]);
                }
            }
        }
    }
    
    (void)force_to_disp;
}

void LoadTransfer::transfer_conservative(
    const std::vector<Vec3>& source,
    Vector& target,
    const std::vector<Real>& source_areas,
    const std::vector<Real>& target_areas) const
{
    if (!mapper_ || !struct_model_) return;
    
    Index ndofs = struct_model_->total_dofs();
    if (target.size() != ndofs) target.resize(ndofs);
    target.zero();
    
    Vec3 total_source = {0, 0, 0};
    for (Index i = 0; i < static_cast<Index>(source.size()); ++i) {
        Real A = (i < static_cast<Index>(source_areas.size())) ? source_areas[i] : 1.0;
        for (Index d = 0; d < 3; ++d) total_source[d] += source[i][d] * A;
    }
    
    const auto& fluid_to_struct = mapper_->fluid_to_struct_map();
    
    for (Index fi = 0; fi < fluid_to_struct.size(); ++fi) {
        if (fi >= static_cast<Index>(source.size())) continue;
        Real A_src = (fi < static_cast<Index>(source_areas.size())) ? source_areas[fi] : 1.0;
        
        Real weight_sum = 0.0;
        for (const auto& pair : fluid_to_struct[fi]) weight_sum += pair.weight;
        if (weight_sum < EPS) continue;
        
        for (const auto& pair : fluid_to_struct[fi]) {
            const auto* node = struct_model_->get_node(pair.struct_node_id);
            if (!node) continue;
            Real A_tgt = (pair.struct_node_id < static_cast<Index>(target_areas.size())) 
                         ? target_areas[pair.struct_node_id] : 1.0;
            Real w = pair.weight / weight_sum * A_src / std::max(A_tgt, EPS);
            Index start = node->dof_start();
            for (Index d = 0; d < 3; ++d) {
                Index gdof = start + d;
                if (gdof < ndofs && node->is_dof_active(d)) {
                    target.add(gdof, w * source[fi][d]);
                }
            }
        }
    }
}

void LoadTransfer::transfer_lumped(
    const std::vector<Vec3>& distributed_per_len,
    Vector& target,
    bool consistent) const
{
    if (!struct_model_) return;
    
    Index ndofs = struct_model_->total_dofs();
    if (target.size() != ndofs) target.resize(ndofs);
    target.zero();
    
    for (Index ei = 0; ei < struct_model_->num_elements(); ++ei) {
        const auto* elem = struct_model_->element(ei);
        if (!elem || ei >= static_cast<Index>(distributed_per_len.size())) continue;
        
        Vec3 q = distributed_per_len[ei];
        Real L = elem->length_initial();
        Vec3 total = math::vec3_scale(q, L);
        
        if (consistent) {
            for (Index ni = 0; ni < elem->num_nodes(); ++ni) {
                const auto* node = elem->node(ni);
                if (!node) continue;
                Index start = node->dof_start();
                Vec3 node_force;
                if (elem->num_nodes() == 2) {
                    for (Index d = 0; d < 3; ++d) {
                        node_force[d] = total[d] * 0.5;
                    }
                } else {
                    for (Index d = 0; d < 3; ++d) {
                        node_force[d] = total[d] / elem->num_nodes();
                    }
                }
                for (Index d = 0; d < 3; ++d) {
                    Index gdof = start + d;
                    if (gdof < ndofs && node->is_dof_active(d)) {
                        target.add(gdof, node_force[d]);
                    }
                }
            }
        } else {
            for (Index ni = 0; ni < elem->num_nodes(); ++ni) {
                const auto* node = elem->node(ni);
                if (!node) continue;
                Index start = node->dof_start();
                Real w = 1.0 / elem->num_nodes();
                for (Index d = 0; d < 3; ++d) {
                    Index gdof = start + d;
                    if (gdof < ndofs && node->is_dof_active(d)) {
                        target.add(gdof, total[d] * w);
                    }
                }
            }
        }
    }
}

void LoadTransfer::transfer_fluid_to_structural(
    const std::vector<Vec3>& fluid_face_forces,
    Vector& structural_node_forces) const
{
    if (load_scheme_ == LoadTransferScheme::DirectNearest ||
        load_scheme_ == LoadTransferScheme::ConservativeInterpolation)
    {
        if (mapper_) {
            mapper_->map_fluid_force_to_struct(fluid_face_forces, structural_node_forces);
        } else {
            transfer_direct_nearest(fluid_face_forces, structural_node_forces, true);
        }
    }
}

void LoadTransfer::transfer_fluid_pressure_to_structural(
    const Vector& fluid_pressure_cell,
    Vector& structural_node_forces) const
{
    if (!fluid_grid_ || !mapper_) return;
    
    std::vector<Vec3> face_forces(fluid_grid_->num_faces(), {0, 0, 0});
    const auto& fluid_iface_faces = mapper_->fluid_interface_faces();
    
    for (Index fi : fluid_iface_faces) {
        if (fi >= fluid_grid_->num_faces()) continue;
        const auto& face = fluid_grid_->face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        Real p = 0.0;
        if (ci >= 0 && ci < fluid_pressure_cell.size()) {
            p = fluid_pressure_cell[ci];
        }
        Vec3 force = math::vec3_scale(face.normal(), -p * face.area());
        face_forces[fi] = force;
    }
    
    transfer_fluid_to_structural(face_forces, structural_node_forces);
}

void LoadTransfer::transfer_structural_to_fluid(
    const Vector& structural_displacement,
    std::vector<Vec3>& fluid_node_displacement) const
{
    if (!mapper_) return;
    
    std::vector<Vec3> face_disp;
    mapper_->map_struct_displacement_to_fluid(structural_displacement, face_disp);
    
    fluid_node_displacement.assign(fluid_grid_->num_nodes(), {0, 0, 0});
    
    const auto& struct_iface_nodes = mapper_->struct_interface_nodes();
    const auto& fluid_iface_faces = mapper_->fluid_interface_faces();
    
    std::vector<Real> weights(fluid_grid_->num_nodes(), 0.0);
    
    for (Index fi : fluid_iface_faces) {
        if (fi >= static_cast<Index>(face_disp.size())) continue;
        const auto& face = fluid_grid_->face(fi);
        for (Index ni = 0; ni < face.num_nodes(); ++ni) {
            Index fn = face.node_id(ni);
            if (fn < fluid_grid_->num_nodes()) {
                for (Index d = 0; d < 3; ++d) {
                    fluid_node_displacement[fn][d] += face_disp[fi][d];
                }
                weights[fn] += 1.0;
            }
        }
    }
    
    for (Index ni = 0; ni < fluid_grid_->num_nodes(); ++ni) {
        if (weights[ni] > EPS) {
            for (Index d = 0; d < 3; ++d) {
                fluid_node_displacement[ni][d] /= weights[ni];
            }
        }
    }
    
    (void)struct_iface_nodes;
}

void LoadTransfer::transfer_structural_velocity_to_fluid(
    const Vector& structural_velocity,
    std::vector<Vec3>& fluid_face_velocity) const
{
    if (!mapper_) return;
    mapper_->map_struct_velocity_to_fluid(structural_velocity, fluid_face_velocity);
}

void LoadTransfer::transfer_face_aero_to_nodal(
    const std::vector<FaceAeroLoad>& face_loads,
    Vector& structural_node_forces) const
{
    if (!fluid_grid_) return;
    
    std::vector<Vec3> face_forces(face_loads.size());
    for (Index i = 0; i < static_cast<Index>(face_loads.size()); ++i) {
        face_forces[i] = face_loads[i].total_force;
    }
    transfer_fluid_to_structural(face_forces, structural_node_forces);
}

void LoadTransfer::transfer_element_aero_to_nodal(
    const std::vector<ElementAeroLoad>& elem_loads,
    Vector& structural_node_forces) const
{
    if (!struct_model_) return;
    
    std::vector<Vec3> distributed(elem_loads.size());
    for (Index i = 0; i < static_cast<Index>(elem_loads.size()); ++i) {
        distributed[i] = elem_loads[i].force_per_unit_length;
    }
    
    bool consistent = (load_scheme_ == LoadTransferScheme::ConsistentNodal ||
                       load_scheme_ == LoadTransferScheme::Mortar);
    transfer_lumped(distributed, structural_node_forces, consistent);
}

void LoadTransfer::compute_consistent_nodal_forces(
    const std::vector<Vec3>& distributed_load_per_length,
    Vector& structural_node_forces) const
{
    transfer_lumped(distributed_load_per_length, structural_node_forces, true);
}

void LoadTransfer::compute_lumped_nodal_forces(
    const std::vector<Vec3>& distributed_load_per_length,
    Vector& structural_node_forces) const
{
    transfer_lumped(distributed_load_per_length, structural_node_forces, false);
}

void LoadTransfer::compute_mortar_integration_masses(
    const std::vector<Index>& interface_faces,
    std::vector<Real>& mass_source,
    std::vector<Real>& mass_target) const
{
    if (!fluid_grid_ || !mapper_) return;
    
    mass_source.assign(interface_faces.size(), 0.0);
    const auto& struct_iface_nodes = mapper_->struct_interface_nodes();
    mass_target.assign(struct_iface_nodes.size(), 0.0);
    
    for (Index i = 0; i < static_cast<Index>(interface_faces.size()); ++i) {
        Index fi = interface_faces[i];
        if (fi < fluid_grid_->num_faces()) {
            mass_source[i] = fluid_grid_->face(fi).area();
        }
    }
    
    const auto& fluid_to_struct = mapper_->fluid_to_struct_map();
    for (Index si = 0; si < fluid_to_struct.size(); ++si) {
        for (const auto& pair : fluid_to_struct[si]) {
            if (pair.struct_node_id < static_cast<Index>(mass_target.size())) {
                mass_target[pair.struct_node_id] += pair.weight * mass_source[si];
            }
        }
    }
}

void LoadTransfer::redistribute_forces_for_conservation(
    const Vector& input_forces,
    Vector& output_forces,
    Real desired_total_force_per_direction[3]) const
{
    if (!struct_model_) return;
    
    output_forces = input_forces;
    
    Real actual_total[3] = {0, 0, 0};
    const auto& struct_iface = mapper_->struct_interface_nodes();
    
    for (Index sid : struct_iface) {
        const auto* node = struct_model_->get_node(sid);
        if (!node) continue;
        Index start = node->dof_start();
        for (Index d = 0; d < 3; ++d) {
            Index gdof = start + d;
            if (gdof < input_forces.size()) actual_total[d] += input_forces[gdof];
        }
    }
    
    Real delta[3];
    for (Index d = 0; d < 3; ++d) delta[d] = desired_total_force_per_direction[d] - actual_total[d];
    
    Index n = static_cast<Index>(struct_iface.size());
    if (n < 1) return;
    
    for (Index sid : struct_iface) {
        const auto* node = struct_model_->get_node(sid);
        if (!node) continue;
        Index start = node->dof_start();
        for (Index d = 0; d < 3; ++d) {
            Index gdof = start + d;
            if (gdof < output_forces.size() && node->is_dof_active(d)) {
                output_forces.add(gdof, delta[d] / n);
            }
        }
    }
}

void LoadTransfer::transfer_temperature_from_structural_to_fluid(
    const Vector& structural_temp,
    Vector& fluid_wall_temperature) const
{
    if (!mapper_ || !fluid_grid_) return;
    
    fluid_wall_temperature.assign(fluid_grid_->num_cells(), 293.15);
    const auto& fluid_iface_faces = mapper_->fluid_interface_faces();
    const auto& fluid_to_struct = mapper_->fluid_to_struct_map();
    
    for (Index fi = 0; fi < fluid_iface_faces.size(); ++fi) {
        Index ff = fluid_iface_faces[fi];
        if (ff >= fluid_grid_->num_faces()) continue;
        Index ci = fluid_grid_->face(ff).owner();
        if (ci < 0 || ci >= fluid_grid_->num_cells()) continue;
        
        Real T = 0.0, wsum = 0.0;
        if (fi < fluid_to_struct.size()) {
            for (const auto& pair : fluid_to_struct[fi]) {
                const auto* node = struct_model_->get_node(pair.struct_node_id);
                if (!node) continue;
                Index start = node->dof_start();
                Index temp_dof = start + 6;
                if (temp_dof < structural_temp.size()) {
                    T += pair.weight * structural_temp[temp_dof];
                    wsum += pair.weight;
                }
            }
        }
        if (wsum > EPS) fluid_wall_temperature[ci] = T / wsum;
    }
}

void LoadTransfer::transfer_heat_flux_from_fluid_to_structural(
    const Vector& fluid_heat_flux,
    Vector& structural_heat_load) const
{
    if (!mapper_ || !struct_model_) return;
    
    Index ndofs = struct_model_->total_dofs();
    if (structural_heat_load.size() != ndofs) structural_heat_load.resize(ndofs);
    structural_heat_load.zero();
    
    const auto& fluid_iface_faces = mapper_->fluid_interface_faces();
    const auto& fluid_to_struct = mapper_->fluid_to_struct_map();
    
    for (Index fi = 0; fi < fluid_iface_faces.size(); ++fi) {
        Index ff = fluid_iface_faces[fi];
        if (ff >= fluid_grid_->num_faces()) continue;
        Index ci = fluid_grid_->face(ff).owner();
        if (ci < 0 || ci >= fluid_heat_flux.size()) continue;
        Real q = fluid_heat_flux[ci] * fluid_grid_->face(ff).area();
        
        if (fi >= fluid_to_struct.size()) continue;
        for (const auto& pair : fluid_to_struct[fi]) {
            const auto* node = struct_model_->get_node(pair.struct_node_id);
            if (!node) continue;
            Index start = node->dof_start();
            Index heat_dof = start + 6;
            if (heat_dof < ndofs) {
                structural_heat_load.add(heat_dof, pair.weight * q);
            }
        }
    }
}

Real LoadTransfer::check_force_conservation(
    const std::vector<Vec3>& fluid_forces,
    const Vector& struct_forces) const
{
    Vec3 total_fluid = {0, 0, 0};
    for (const auto& f : fluid_forces) math::vec3_add_inplace(total_fluid, f);
    
    Vec3 total_struct = {0, 0, 0};
    if (mapper_) {
        for (Index sid : mapper_->struct_interface_nodes()) {
            const auto* node = struct_model_->get_node(sid);
            if (!node) continue;
            Index start = node->dof_start();
            for (Index d = 0; d < 3; ++d) {
                Index gdof = start + d;
                if (gdof < struct_forces.size()) total_struct[d] += struct_forces[gdof];
            }
        }
    }
    
    Real diff = math::vec3_norm(math::vec3_sub(total_fluid, total_struct));
    Real mag = math::vec3_norm(total_fluid);
    return (mag > EPS) ? diff / mag : 0.0;
}

Real LoadTransfer::check_energy_conservation(
    const Vector& struct_vel,
    const Vector& struct_force,
    const std::vector<Vec3>& fluid_wall_vel,
    const std::vector<Vec3>& fluid_force) const
{
    Real P_struct = 0.0;
    Index n_struct = std::min(struct_vel.size(), struct_force.size());
    for (Index i = 0; i < n_struct; ++i) P_struct += struct_vel[i] * struct_force[i];
    
    Real P_fluid = 0.0;
    Index n_fluid = std::min(fluid_wall_vel.size(), fluid_force.size());
    for (Index i = 0; i < n_fluid; ++i) {
        for (Index d = 0; d < 3; ++d) {
            P_fluid += fluid_wall_vel[i][d] * fluid_force[i][d];
        }
    }
    
    Real denom = std::max(std::fabs(P_struct), std::fabs(P_fluid));
    if (denom < EPS) return 0.0;
    return std::fabs(P_struct + P_fluid) / denom;
}

} // namespace fsi
} // namespace hvdc
