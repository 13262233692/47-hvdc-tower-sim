#ifndef HVDC_FEM_ASSEMBLER_HPP
#define HVDC_FEM_ASSEMBLER_HPP

#include "common/Types.hpp"
#include "common/Vector.hpp"
#include "common/SparseMatrix.hpp"
#ifdef HVDC_USE_PETSC
#include "common/PETScManager.hpp"
#endif
#include "fem/FEModel.hpp"
#include "fem/Element.hpp"
#include "fem/BeamElementNL.hpp"
#include "fem/TrussElementNL.hpp"
#include "fem/ConductorElement.hpp"

#include <vector>
#include <memory>

namespace hvdc {
namespace fem {

struct AssemblyResult {
    SparseMatrix K_T;
    SparseMatrix K_Mat;
    SparseMatrix K_Geo;
    SparseMatrix M;
    Vector F_int;
    Vector F_ext;
    Vector F_gravity;
    Vector F_coupling;
    Vector Residual;
    Index nnz_K = 0;
    Index nnz_M = 0;
    Real assembly_time = 0.0;
    Real element_process_time = 0.0;
};

class Assembler {
public:
    Assembler() = default;
    explicit Assembler(FEModel* model);
    
    void set_model(FEModel* model) { model_ = model; }
    FEModel* model() { return model_; }
    const FEModel* model() const { return model_; }
    
    void initialize();
    
    void assemble_tangent_stiffness(
        const Vector& displacement,
        AssemblyResult& result,
        bool include_geometric = true,
        bool include_mass = false,
        bool use_parallel = true);
    
    void assemble_internal_forces(
        const Vector& displacement,
        Vector& F_int,
        bool use_parallel = true);
    
    void assemble_external_loads(
        Vector& F_ext,
        bool include_gravity = true,
        bool include_equivalent = true);
    
    void assemble_mass_matrix(
        SparseMatrix& M,
        bool lumped = true,
        bool use_parallel = true);
    
    void compute_residual(
        const Vector& displacement,
        const Vector& F_extra,
        Vector& residual,
        AssemblyResult* cached_result = nullptr);
    
    void apply_dirichlet_bc(
        SparseMatrix& K,
        Vector& F,
        const Vector& displacement,
        Real penalty_factor = 1.0e18);
    
    void apply_dirichlet_bc(
        SparseMatrix& K,
        Vector& F,
        const Vector& displacement,
        const FEModel& model,
        Real penalty_factor = 1.0e18);
    
    void apply_dof_constraints(
        SparseMatrix& K,
        Vector& F,
        const IndexVec& constrained_dofs,
        const Vec& bc_values,
        Real penalty_factor = 1.0e18);
    
    void estimate_nnz_pattern(
        IndexVec& d_nnz,
        IndexVec& o_nnz,
        Index nnz_per_row_estimate = 60);
    
    Index estimate_total_nnz() const;
    
    void assemble_tangent_stiffness(
        const Vector& displacement,
        SparseMatrix& K_out,
        FEModel* model = nullptr);
    
    void assemble_mass_matrix(
        SparseMatrix& M,
        FEModel* model,
        bool lumped = true);
    
    void assemble_damping_matrix(
        SparseMatrix& C,
        FEModel* model,
        Real alpha = 0.05,
        Real beta = 0.01);
    
    void multiply_stiffness_vector(
        const SparseMatrix& K,
        const Vector& x,
        Vector& y,
        Real scale_K = 1.0,
        Real scale_y = 0.0);

protected:
    FEModel* model_ = nullptr;
    
    Index total_dofs_ = 0;
    bool initialized_ = false;
    
    std::vector<IndexVec> element_dof_map_;
    IndexVec local_elem_range_;
    
    void gather_element_dofs(const Element* elem,
                              IndexVec& global_dofs) const;
    
    void scatter_12x12_to_global(
        SparseMatrix& K_global,
        const Mat12x12& K_elem,
        const IndexVec& dof_ids,
        Real scale = 1.0);
    
    void scatter_12vec_to_global(
        Vector& F_global,
        const Vec12& F_elem,
        const IndexVec& dof_ids,
        Real scale = 1.0);
    
    void process_single_element(
        const Element* elem,
        const Vector& displacement,
        bool include_geometric,
        Mat12x12& K_T,
        Vec12& F_int,
        Vec12& F_ext);
    
    void process_element_range(
        Index elem_start,
        Index elem_end,
        const Vector& displacement,
        bool include_geometric,
        SparseMatrix& K_T_local,
        Vector& F_int_local,
        Vector& F_ext_local,
        bool compute_K,
        bool compute_F);
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_ASSEMBLER_HPP
