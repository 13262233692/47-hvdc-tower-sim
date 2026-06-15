#include "fsi/InterfaceMapper.hpp"
#include "common/MathUtils.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <queue>

namespace hvdc {
namespace fsi {

InterfaceMapper::InterfaceMapper(fem::FEModel* struct_model,
                                 fvm::Grid* fluid_grid,
                                 MappingMethod method)
    : struct_model_(struct_model), fluid_grid_(fluid_grid), method_(method)
{
}

void InterfaceMapper::set_structural_interface_nodes(const std::vector<Index>& node_ids) {
    struct_interface_nodes_ = node_ids;
    struct_positions_.resize(node_ids.size());
    for (Index i = 0; i < static_cast<Index>(node_ids.size()); ++i) {
        const auto* node = struct_model_->get_node(node_ids[i]);
        if (node) struct_positions_[i] = node->coords();
    }
}

void InterfaceMapper::set_fluid_interface_faces(const std::vector<Index>& face_ids) {
    fluid_interface_faces_ = face_ids;
    fluid_face_positions_.resize(face_ids.size());
    fluid_face_normals_.resize(face_ids.size());
    fluid_face_areas_.resize(face_ids.size());
    for (Index i = 0; i < static_cast<Index>(face_ids.size()); ++i) {
        const auto& face = fluid_grid_->face(face_ids[i]);
        fluid_face_positions_[i] = face.centroid();
        fluid_face_normals_[i] = face.normal();
        fluid_face_areas_[i] = face.area();
    }
}

void InterfaceMapper::set_fluid_interface_patch(Index boundary_patch_id) {
    auto range = fluid_grid_->boundary_patch_range(boundary_patch_id);
    if (range.first < 0) return;
    std::vector<Index> face_ids;
    for (Index fi = range.first; fi < range.second; ++fi) {
        face_ids.push_back(fi);
    }
    set_fluid_interface_faces(face_ids);
}

Real InterfaceMapper::rbf_gaussian(Real r, Real eps) const {
    return std::exp(-(eps * r) * (eps * r));
}

Real InterfaceMapper::rbf_multiquadric(Real r, Real eps) const {
    return std::sqrt(1.0 + (eps * r) * (eps * r));
}

Real InterfaceMapper::rbf_tps(Real r) const {
    if (r < EPS) return 0.0;
    return r * r * std::log(r);
}

Real InterfaceMapper::compute_influence_radius() const {
    if (struct_positions_.empty() || fluid_face_positions_.empty()) return 0.0;
    
    Real min_dist = INF;
    for (const auto& sp : struct_positions_) {
        for (const auto& fp : fluid_face_positions_) {
            Real d = math::vec3_norm(math::vec3_sub(sp, fp));
            if (d > EPS) min_dist = std::min(min_dist, d);
        }
    }
    return min_dist;
}

void InterfaceMapper::normalize_weights(std::vector<InterfacePair>& pairs) const {
    if (pairs.empty()) return;
    Real total = 0.0;
    for (const auto& p : pairs) total += p.weight;
    if (total < EPS) {
        if (!pairs.empty()) pairs[0].weight = 1.0;
        return;
    }
    for (auto& p : pairs) p.weight /= total;
}

void InterfaceMapper::compute_nearest_pairs(
    const std::vector<Vec3>& source,
    const std::vector<Vec3>& target,
    std::vector<std::vector<InterfacePair>>& mapping,
    Real max_radius,
    Index max_neighbors) const
{
    Real r2 = max_radius * max_radius;
    mapping.assign(target.size(), {});
    
    for (Index ti = 0; ti < static_cast<Index>(target.size()); ++ti) {
        std::priority_queue<std::pair<Real, Index>> pq;
        
        for (Index si = 0; si < static_cast<Index>(source.size()); ++si) {
            Real d2 = math::vec3_norm2(math::vec3_sub(target[ti], source[si]));
            if (d2 <= r2) {
                pq.emplace(-d2, si);
                if (static_cast<Index>(pq.size()) > max_neighbors) pq.pop();
            }
        }
        
        while (!pq.empty()) {
            auto [neg_d2, si] = pq.top();
            pq.pop();
            Real d = std::sqrt(-neg_d2);
            InterfacePair p;
            p.struct_node_id = si;
            p.distance = d;
            p.weight = (method_ == MappingMethod::NearestNeighbor) ? 1.0 
                       : rbf_gaussian(d, rbf_eps_);
            mapping[ti].push_back(p);
        }
        
        std::reverse(mapping[ti].begin(), mapping[ti].end());
        normalize_weights(mapping[ti]);
    }
}

void InterfaceMapper::build_mapping(Real search_radius_factor) {
    if (!struct_model_ || !fluid_grid_) return;
    
    if (struct_interface_nodes_.empty()) {
        auto& nodes = struct_model_->interface_nodes();
        std::vector<Index> ids;
        for (const auto* n : nodes) ids.push_back(n->id());
        set_structural_interface_nodes(ids);
    }
    
    if (fluid_interface_faces_.empty()) {
        const auto& cell_ids = fluid_grid_->interface_cell_ids();
        std::vector<Index> face_ids;
        for (Index ci : cell_ids) {
            const auto& cell = fluid_grid_->cell(ci);
            for (Index fi = 0; fi < cell.num_faces() && fi < cell.face_ids().size(); ++fi) {
                Index fid = cell.face_ids()[fi];
                if (fid >= 0) face_ids.push_back(fid);
            }
        }
        set_fluid_interface_faces(face_ids);
    }
    
    Real infl_rad = compute_influence_radius();
    if (infl_rad < EPS) infl_rad = 1.0;
    Real search_radius = infl_rad * search_radius_factor;
    
    compute_nearest_pairs(fluid_face_positions_, struct_positions_,
                          struct_to_fluid_, search_radius, 10);
    
    compute_nearest_pairs(struct_positions_, fluid_face_positions_,
                          fluid_to_struct_, search_radius, 10);
    
    for (Index si = 0; si < struct_to_fluid_.size(); ++si) {
        for (auto& p : struct_to_fluid_[si]) {
            if (p.struct_node_id >= 0 && p.struct_node_id < static_cast<Index>(struct_interface_nodes_.size())) {
                p.struct_node_id = struct_interface_nodes_[si];
            }
            if (!fluid_interface_faces_.empty()) {
                // p.fluid_face_id = fluid_interface_faces_[...];
            }
        }
    }
    
    HVDC_LOG_INFO("FSI interface mapping built: struct nodes="
                  << struct_interface_nodes_.size()
                  << ", fluid faces=" << fluid_interface_faces_.size()
                  << ", method=" << static_cast<Index>(method_));
}

void InterfaceMapper::build_conservative_mapping() {
    if (!struct_model_ || !fluid_grid_) return;
    build_mapping(3.0);
}

void InterfaceMapper::build_rbf_mapping(Real support_radius, Real shape_param) {
    rbf_eps_ = shape_param;
    build_mapping(support_radius);
}

void InterfaceMapper::map_struct_displacement_to_fluid(
    const Vector& struct_disp,
    std::vector<Vec3>& fluid_face_displacement) const
{
    fluid_face_displacement.assign(fluid_grid_->num_faces(), {0, 0, 0});
    
    for (Index fi = 0; fi < fluid_interface_faces_.size(); ++fi) {
        Index ff_id = fluid_interface_faces_[fi];
        Vec3 disp = {0, 0, 0};
        
        for (const auto& pair : fluid_to_struct_[fi]) {
            if (pair.struct_node_id < 0) continue;
            const auto* node = struct_model_->get_node(pair.struct_node_id);
            if (!node) continue;
            Index start = node->dof_start();
            Vec3 nd = {0, 0, 0};
            for (Index d = 0; d < 3; ++d) {
                Index gdof = start + d;
                if (gdof < struct_disp.size()) nd[d] = struct_disp[gdof];
            }
            for (Index d = 0; d < 3; ++d) disp[d] += pair.weight * nd[d];
        }
        
        if (ff_id < static_cast<Index>(fluid_face_displacement.size())) {
            fluid_face_displacement[ff_id] = disp;
        }
    }
}

void InterfaceMapper::map_struct_velocity_to_fluid(
    const Vector& struct_vel,
    std::vector<Vec3>& fluid_face_velocity) const
{
    map_struct_displacement_to_fluid(struct_vel, fluid_face_velocity);
}

void InterfaceMapper::map_fluid_force_to_struct(
    const std::vector<Vec3>& fluid_face_force,
    Vector& struct_node_force) const
{
    Index ndofs = struct_model_->total_dofs();
    if (struct_node_force.size() != ndofs) {
        struct_node_force.resize(ndofs);
    }
    struct_node_force.zero();
    
    for (Index fi = 0; fi < fluid_interface_faces_.size(); ++fi) {
        Index ff_id = fluid_interface_faces_[fi];
        if (ff_id >= static_cast<Index>(fluid_face_force.size())) continue;
        const Vec3& ff = fluid_face_force[ff_id];
        
        for (const auto& pair : fluid_to_struct_[fi]) {
            if (pair.struct_node_id < 0) continue;
            const auto* node = struct_model_->get_node(pair.struct_node_id);
            if (!node) continue;
            Index start = node->dof_start();
            for (Index d = 0; d < 3; ++d) {
                Index gdof = start + d;
                if (gdof < ndofs && node->is_dof_active(d)) {
                    struct_node_force.add(gdof, pair.weight * ff[d]);
                }
            }
        }
    }
}

void InterfaceMapper::map_fluid_pressure_to_struct(
    const std::vector<Real>& fluid_pressure,
    Vector& struct_node_force) const
{
    std::vector<Vec3> face_forces(fluid_grid_->num_faces(), {0, 0, 0});
    
    for (Index fi = 0; fi < fluid_interface_faces_.size(); ++fi) {
        Index ff_id = fluid_interface_faces_[fi];
        if (ff_id >= fluid_grid_->num_faces()) continue;
        const auto& face = fluid_grid_->face(ff_id);
        Real p = 0.0;
        if (ff_id < static_cast<Index>(fluid_pressure.size())) {
            p = fluid_pressure[ff_id];
        } else {
            Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
            p = 0.0;
        }
        Vec3 force = math::vec3_scale(face.normal(), -p * face.area());
        face_forces[ff_id] = force;
    }
    
    map_fluid_force_to_struct(face_forces, struct_node_force);
}

Real InterfaceMapper::compute_conservation_error(
    const std::vector<Vec3>& face_forces,
    const Vector& struct_forces) const
{
    Vec3 sum_fluid = {0, 0, 0};
    for (Index fi = 0; fi < fluid_interface_faces_.size(); ++fi) {
        Index ff_id = fluid_interface_faces_[fi];
        if (ff_id < static_cast<Index>(face_forces.size())) {
            math::vec3_add_inplace(sum_fluid, face_forces[ff_id]);
        }
    }
    
    Vec3 sum_struct = {0, 0, 0};
    for (Index ni : struct_interface_nodes_) {
        const auto* node = struct_model_->get_node(ni);
        if (!node) continue;
        Index start = node->dof_start();
        for (Index d = 0; d < 3; ++d) {
            Index gdof = start + d;
            if (gdof < struct_forces.size()) {
                sum_struct[d] += struct_forces[gdof];
            }
        }
    }
    
    Real fluid_mag = math::vec3_norm(sum_fluid);
    if (fluid_mag < EPS) return 0.0;
    return math::vec3_norm(math::vec3_sub(sum_fluid, sum_struct)) / fluid_mag;
}

void InterfaceMapper::deform_fluid_mesh(
    const Vector& struct_disp,
    std::vector<Vec3>& new_node_coords,
    Real spring_stiffness) const
{
    Index nf_nodes = fluid_grid_->num_nodes();
    new_node_coords = fluid_grid_->nodes();
    
    std::vector<Vec3> disp(nf_nodes, {0, 0, 0});
    std::vector<Real> weights(nf_nodes, 0.0);
    
    for (Index ni = 0; ni < struct_interface_nodes_.size(); ++ni) {
        Index sn_id = struct_interface_nodes_[ni];
        const auto* node = struct_model_->get_node(sn_id);
        if (!node) continue;
        Index start = node->dof_start();
        Vec3 nd = {0, 0, 0};
        for (Index d = 0; d < 3; ++d) {
            Index gdof = start + d;
            if (gdof < struct_disp.size()) nd[d] = struct_disp[gdof];
        }
        
        Vec3 sp = struct_positions_[ni];
        Real r_max = compute_influence_radius() * 5.0;
        Real r2_max = r_max * r_max;
        
        for (Index fni = 0; fni < nf_nodes; ++fni) {
            Vec3 diff = math::vec3_sub(new_node_coords[fni], sp);
            Real d2 = math::vec3_norm2(diff);
            if (d2 < r2_max) {
                Real d = std::sqrt(d2);
                Real w = std::exp(-spring_stiffness * d2 / (r2_max * 0.1));
                for (Index di = 0; di < 3; ++di) {
                    disp[fni][di] += w * nd[di];
                    weights[fni] += w;
                }
            }
        }
    }
    
    for (Index fni = 0; fni < nf_nodes; ++fni) {
        if (weights[fni] > EPS) {
            for (Index di = 0; di < 3; ++di) {
                new_node_coords[fni][di] += disp[fni][di] / weights[fni];
            }
        }
    }
}

} // namespace fsi
} // namespace hvdc
