#include "fem/Assembler.hpp"
#include "common/MathUtils.hpp"
#include "common/MPIManager.hpp"
#include "common/Timer.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>

extern hvdc::Timer g_timer;

namespace hvdc {
namespace fem {

Assembler::Assembler(FEModel* model)
    : model_(model)
{
}

void Assembler::initialize() {
    if (!model_) {
        HVDC_LOG_ERROR("Assembler::initialize: model is null");
        return;
    }
    
    total_dofs_ = model_->total_dofs();
    
    Index n_elems = model_->num_elements();
    element_dof_map_.resize(n_elems);
    
    for (Index ei = 0; ei < n_elems; ++ei) {
        const Element* elem = model_->element(ei);
        gather_element_dofs(elem, element_dof_map_[ei]);
    }
    
    auto& mpi = MPIManager::instance();
    Index s, e;
    mpi.partition_range(n_elems, s, e);
    local_elem_range_ = {s, e};
    
    initialized_ = true;
    
    HVDC_LOG_INFO("Assembler initialized: total_dofs=" << total_dofs_
                  << ", elements=" << n_elems
                  << ", local elements: " << s << "-" << e);
}

void Assembler::gather_element_dofs(const Element* elem,
                                      IndexVec& global_dofs) const
{
    global_dofs.clear();
    if (!elem) return;
    
    for (Index ni = 0; ni < elem->num_nodes(); ++ni) {
        const Node* n = elem->node(ni);
        if (!n) continue;
        Index start = n->dof_start();
        for (Index d = 0; d < n->num_dofs(); ++d) {
            global_dofs.push_back(start + d);
        }
    }
}

void Assembler::scatter_12x12_to_global(
    SparseMatrix& K_global,
    const Mat12x12& K_elem,
    const IndexVec& dof_ids,
    Real scale)
{
    Index ndofs = std::min<Index>(12, static_cast<Index>(dof_ids.size()));
    for (Index i = 0; i < ndofs; ++i) {
        Index gi = dof_ids[i];
        if (gi < 0) continue;
        for (Index j = 0; j < ndofs; ++j) {
            Index gj = dof_ids[j];
            if (gj < 0) continue;
            Real val = K_elem[i][j] * scale;
            if (std::fabs(val) > EPS) {
                K_global.add(gi, gj, val);
            }
        }
    }
}

void Assembler::scatter_12vec_to_global(
    Vector& F_global,
    const Vec12& F_elem,
    const IndexVec& dof_ids,
    Real scale)
{
    Index ndofs = std::min<Index>(12, static_cast<Index>(dof_ids.size()));
    for (Index i = 0; i < ndofs; ++i) {
        Index gi = dof_ids[i];
        if (gi >= 0 && gi < F_global.size()) {
            F_global.add(gi, F_elem[i] * scale);
        }
    }
}

void Assembler::process_single_element(
    const Element* elem,
    const Vector& displacement,
    bool include_geometric,
    Mat12x12& K_T,
    Vec12& F_int,
    Vec12& F_ext)
{
    Vec12 d_elem{};
    for (Index ni = 0; ni < elem->num_nodes(); ++ni) {
        Index start = elem->node(ni)->dof_start();
        for (Index d = 0; d < 6; ++d) {
            Index gdof = start + d;
            if (gdof < displacement.size()) {
                d_elem[ni * 6 + d] = displacement[gdof];
            }
        }
    }
    
    elem->tangent_stiffness_matrix(K_T, d_elem, include_geometric);
    F_int = elem->internal_forces(d_elem);
    F_ext = elem->equivalent_nodal_loads();
}

void Assembler::process_element_range(
    Index elem_start,
    Index elem_end,
    const Vector& displacement,
    bool include_geometric,
    SparseMatrix& K_T_local,
    Vector& F_int_local,
    Vector& F_ext_local,
    bool compute_K,
    bool compute_F)
{
    for (Index ei = elem_start; ei < elem_end; ++ei) {
        const Element* elem = model_->element(ei);
        if (!elem) continue;
        
        Mat12x12 K_T{};
        Vec12 F_int{}, F_ext{};
        const auto& dof_ids = element_dof_map_[ei];
        
        if (compute_K && compute_F) {
            process_single_element(elem, displacement, include_geometric,
                                    K_T, F_int, F_ext);
            scatter_12x12_to_global(K_T_local, K_T, dof_ids);
            scatter_12vec_to_global(F_int_local, F_int, dof_ids);
            scatter_12vec_to_global(F_ext_local, F_ext, dof_ids);
        } else if (compute_F) {
            Vec12 d_elem{};
            for (Index ni = 0; ni < elem->num_nodes(); ++ni) {
                Index start = elem->node(ni)->dof_start();
                for (Index d = 0; d < 6; ++d) {
                    Index gdof = start + d;
                    if (gdof < displacement.size()) {
                        d_elem[ni * 6 + d] = displacement[gdof];
                    }
                }
            }
            F_int = elem->internal_forces(d_elem);
            F_ext = elem->equivalent_nodal_loads();
            scatter_12vec_to_global(F_int_local, F_int, dof_ids);
            scatter_12vec_to_global(F_ext_local, F_ext, dof_ids);
        } else if (compute_K) {
            Vec12 d_elem{};
            for (Index ni = 0; ni < elem->num_nodes(); ++ni) {
                Index start = elem->node(ni)->dof_start();
                for (Index d = 0; d < 6; ++d) {
                    Index gdof = start + d;
                    if (gdof < displacement.size()) {
                        d_elem[ni * 6 + d] = displacement[gdof];
                    }
                }
            }
            elem->tangent_stiffness_matrix(K_T, d_elem, include_geometric);
            scatter_12x12_to_global(K_T_local, K_T, dof_ids);
        }
    }
}

void Assembler::assemble_tangent_stiffness(
    const Vector& displacement,
    AssemblyResult& result,
    bool include_geometric,
    bool include_mass,
    bool use_parallel)
{
    if (!initialized_) initialize();
    if (!initialized_) return;
    
    Timer local_timer;
    local_timer.start("assemble");
    
    Index n_dofs = model_->total_dofs();
    Index nnz_est = estimate_total_nnz();
    
    result.K_T.reset(n_dofs, n_dofs);
    result.F_int.resize(n_dofs);
    result.F_int.zero();
    result.F_ext.resize(n_dofs);
    result.F_ext.zero();
    
    Index elem_start = 0;
    Index elem_end = model_->num_elements();
    
    if (use_parallel) {
        elem_start = local_elem_range_[0];
        elem_end = local_elem_range_[1];
    }
    
    process_element_range(elem_start, elem_end, displacement, include_geometric,
                           result.K_T, result.F_int, result.F_ext,
                           true, true);
    
    result.K_T.finalize();
    
    if (include_mass) {
        result.M.reset(n_dofs, n_dofs);
        assemble_mass_matrix(result.M, true, use_parallel);
    }
    
    result.nnz_K = result.K_T.nnz();
    result.assembly_time = local_timer.elapsed("assemble");
    
    HVDC_LOG_INFO("Tangent stiffness assembled: nnz=" << result.nnz_K
                  << ", time=" << result.assembly_time << "s");
}

void Assembler::assemble_internal_forces(
    const Vector& displacement,
    Vector& F_int,
    bool use_parallel)
{
    if (!initialized_) initialize();
    
    Index n_dofs = model_->total_dofs();
    if (F_int.size() != n_dofs) {
        F_int.resize(n_dofs);
    }
    F_int.zero();
    
    SparseMatrix dummy_K;
    Vector dummy_F_ext(n_dofs);
    dummy_F_ext.zero();
    
    Index elem_start = 0;
    Index elem_end = model_->num_elements();
    if (use_parallel) {
        elem_start = local_elem_range_[0];
        elem_end = local_elem_range_[1];
    }
    
    process_element_range(elem_start, elem_end, displacement, false,
                           dummy_K, F_int, dummy_F_ext, false, true);
}

void Assembler::assemble_external_loads(
    Vector& F_ext,
    bool include_gravity,
    bool include_equivalent)
{
    Index n_dofs = model_->total_dofs();
    if (F_ext.size() != n_dofs) {
        F_ext.resize(n_dofs);
    }
    F_ext.zero();
    
    Index n_elems = model_->num_elements();
    for (Index ei = 0; ei < n_elems; ++ei) {
        const Element* elem = model_->element(ei);
        if (!elem) continue;
        
        Vec12 F_elem{};
        bool has_load = false;
        
        if (include_equivalent) {
            Vec12 F_eq = elem->equivalent_nodal_loads();
            for (Index i = 0; i < 12; ++i) F_elem[i] += F_eq[i];
            has_load = true;
        }
        
        if (include_gravity) {
            Vec12 F_g = elem->gravity_loads();
            for (Index i = 0; i < 12; ++i) F_elem[i] += F_g[i];
            has_load = true;
        }
        
        if (has_load) {
            scatter_12vec_to_global(F_ext, F_elem, element_dof_map_[ei]);
        }
    }
}

void Assembler::assemble_mass_matrix(
    SparseMatrix& M,
    bool lumped,
    bool use_parallel)
{
    if (!initialized_) initialize();
    
    Index n_dofs = model_->total_dofs();
    if (M.rows() != n_dofs || M.cols() != n_dofs) {
        M.reset(n_dofs, n_dofs);
    }
    
    Index n_elems = model_->num_elements();
    Index s = 0, e = n_elems;
    if (use_parallel) {
        s = local_elem_range_[0];
        e = local_elem_range_[1];
    }
    
    for (Index ei = s; ei < e; ++ei) {
        const Element* elem = model_->element(ei);
        if (!elem) continue;
        
        Mat12x12 M_elem{};
        elem->mass_matrix(M_elem);
        
        if (lumped) {
            Mat12x12 M_lumped{};
            for (Index i = 0; i < 12; ++i) {
                Real row_sum = 0.0;
                for (Index j = 0; j < 12; ++j) row_sum += std::fabs(M_elem[i][j]);
                M_lumped[i][i] = row_sum > EPS ? row_sum : M_elem[i][i];
            }
            scatter_12x12_to_global(M, M_lumped, element_dof_map_[ei]);
        } else {
            scatter_12x12_to_global(M, M_elem, element_dof_map_[ei]);
        }
    }
    
    M.finalize();
}

void Assembler::compute_residual(
    const Vector& displacement,
    const Vector& F_extra,
    Vector& residual,
    AssemblyResult* cached_result)
{
    AssemblyResult local_result;
    AssemblyResult* result = cached_result ? cached_result : &local_result;
    
    if (!cached_result) {
        assemble_tangent_stiffness(displacement, *result, true, false, true);
    }
    
    Vector F_ext_total(result->F_ext);
    if (F_extra.size() > 0) {
        F_ext_total += F_extra;
    }
    
    residual = displacement;
    residual.resize(model_->total_dofs());
    for (Index i = 0; i < model_->total_dofs(); ++i) {
        residual[i] = F_ext_total[i] - result->F_int[i];
    }
}

void Assembler::apply_dirichlet_bc(
    SparseMatrix& K,
    Vector& F,
    const Vector& displacement,
    Real penalty_factor)
{
    Index n_dofs = model_->total_dofs();
    IndexVec cdofs = model_->constrained_dofs();
    
    for (Index gdof : cdofs) {
        Index node_id = gdof / 6;
        Index local_dof = gdof % 6;
        const Node* node = model_->get_node(node_id);
        if (!node) continue;
        
        Real bc_val = node->bc_value(local_dof);
        Real diff = displacement[gdof] - bc_val;
        
        IndexVec cols = {gdof};
        Vec vals = {penalty_factor};
        K.set_row(gdof, cols, vals);
        
        F.set(gdof, penalty_factor * bc_val);
        (void)diff;
    }
}

void Assembler::apply_dof_constraints(
    SparseMatrix& K,
    Vector& F,
    const IndexVec& constrained_dofs,
    const Vec& bc_values,
    Real penalty_factor)
{
    for (Index idx = 0; idx < static_cast<Index>(constrained_dofs.size()); ++idx) {
        Index gdof = constrained_dofs[idx];
        Real bc_val = (idx < static_cast<Index>(bc_values.size())) ? bc_values[idx] : 0.0;
        
        Real old_diag = K(gdof, gdof);
        K.set(gdof, gdof, old_diag + penalty_factor);
        F.set(gdof, F[gdof] + penalty_factor * bc_val);
    }
}

void Assembler::estimate_nnz_pattern(
    IndexVec& d_nnz,
    IndexVec& o_nnz,
    Index nnz_per_row_estimate)
{
    Index n_dofs = model_->total_dofs();
    auto& mpi = MPIManager::instance();
    
    Index local_start, local_end;
    mpi.partition_range(n_dofs, local_start, local_end);
    
    Index local_rows = local_end - local_start;
    d_nnz.assign(local_rows, nnz_per_row_estimate);
    o_nnz.assign(local_rows, nnz_per_row_estimate / 3);
}

Index Assembler::estimate_total_nnz() const {
    Index nnz = 0;
    Index n_elems = model_->num_elements();
    for (Index ei = 0; ei < n_elems; ++ei) {
        Index dofs_per_elem = static_cast<Index>(element_dof_map_[ei].size());
        nnz += dofs_per_elem * dofs_per_elem;
    }
    return nnz;
}

void Assembler::assemble_tangent_stiffness(
    const Vector& displacement,
    SparseMatrix& K_out,
    FEModel* model)
{
    if (model && !model_) model_ = model;
    AssemblyResult result;
    assemble_tangent_stiffness(displacement, result, true, false, true);
    K_out = std::move(result.K_T);
}

void Assembler::assemble_mass_matrix(
    SparseMatrix& M,
    FEModel* model,
    bool lumped)
{
    if (model && !model_) model_ = model;
    assemble_mass_matrix(M, lumped, true);
}

void Assembler::assemble_damping_matrix(
    SparseMatrix& C,
    FEModel* model,
    Real alpha,
    Real beta)
{
    if (model && !model_) model_ = model;
    if (!model_ || !initialized_) initialize();
    Index ndofs = model_->total_dofs();
    
    SparseMatrix M, K;
    assemble_mass_matrix(M, model_, true);
    AssemblyResult res;
    Vector zero_u(ndofs); zero_u.zero();
    assemble_tangent_stiffness(zero_u, res, true, false, false);
    K = std::move(res.K_T);
    
    C = SparseMatrix(ndofs, ndofs, std::min(M.nnz() + K.nnz(), ndofs * 60));
    for (Index r = 0; r < ndofs; ++r) {
        for (Index j = M.row_start(r); j < M.row_start(r + 1); ++j) {
            Index c = M.col(j);
            C.add(r, c, alpha * M.value(j));
        }
        for (Index j = K.row_start(r); j < K.row_start(r + 1); ++j) {
            Index c = K.col(j);
            C.add(r, c, beta * K.value(j));
        }
    }
    C.finalize();
}

void Assembler::multiply_stiffness_vector(
    const SparseMatrix& K,
    const Vector& x,
    Vector& y,
    Real scale_K,
    Real scale_y)
{
    Index n = K.rows();
    if (y.size() != n) y.resize(n);
    if (std::fabs(scale_y) < EPS) {
        y.zero();
    } else {
        for (Index i = 0; i < n; ++i) y[i] *= scale_y;
    }
    for (Index r = 0; r < n; ++r) {
        Real sum = 0.0;
        for (Index j = K.row_start(r); j < K.row_start(r + 1); ++j) {
            Index c = K.col(j);
            if (c < x.size()) sum += K.value(j) * x[c];
        }
        y.add(r, scale_K * sum);
    }
}

void Assembler::apply_dirichlet_bc(
    SparseMatrix& K,
    Vector& F,
    const Vector& displacement,
    const FEModel& model,
    Real penalty_factor)
{
    (void)displacement;
    Index ndofs = model.total_dofs();
    const auto& cdofs = model.constrained_dofs();
    Vec bc_vals(cdofs.size(), 0.0);
    for (Index i = 0; i < cdofs.size(); ++i) {
        const Node* n = nullptr;
        for (Index ni = 0; ni < model.num_nodes() && !n; ++ni) {
            const Node* cn = model.node(ni);
            if (cn && cn->dof_start() <= cdofs[i] && cdofs[i] < cn->dof_start() + cn->num_dofs()) {
                Index local = cdofs[i] - cn->dof_start();
                if (local < 7) bc_vals[i] = cn->bc_value(local);
                n = cn;
            }
        }
    }
    apply_dof_constraints(K, F, cdofs, bc_vals, penalty_factor);
}

} // namespace fem
} // namespace hvdc
