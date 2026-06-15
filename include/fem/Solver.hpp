#ifndef HVDC_FEM_SOLVER_HPP
#define HVDC_FEM_SOLVER_HPP

#include "common/Types.hpp"
#include "common/Vector.hpp"
#include "common/SparseMatrix.hpp"
#ifdef HVDC_USE_PETSC
#include "common/PETScManager.hpp"
#endif
#include "fem/FEModel.hpp"
#include "fem/Assembler.hpp"

#include <string>
#include <memory>

namespace hvdc {
namespace fem {

struct SolverConfig {
    Real tol_abs = 1.0e-8;
    Real tol_rel = 1.0e-6;
    Index max_iters = 100;
    Index max_line_search = 10;
    Real line_search_alpha = 0.5;
    Real line_search_beta = 0.8;
    bool use_displacement_control = false;
    Real load_increment = 1.0;
    bool predict_solution = true;
    bool update_stiffness_each_iter = true;
    std::string ksp_type = "gmres";
    std::string pc_type = "gamg";
    Real krylov_rtol = 1.0e-6;
    Real krylov_atol = 1.0e-50;
    Index krylov_max_iters = 20000;
};

struct SolverResult {
    Vector displacement;
    Vector velocity;
    Vector acceleration;
    Vector reaction_forces;
    Real residual_norm = 0.0;
    Real initial_residual_norm = 0.0;
    Index newton_iters = 0;
    Index total_linear_iters = 0;
    bool converged = false;
    std::string converged_reason;
    Real solve_time = 0.0;
    Real assemble_time = 0.0;
    std::vector<Real> residual_history;
    std::vector<Index> linear_iter_history;
};

class Solver {
public:
    Solver() = default;
    explicit Solver(FEModel* model, Assembler* assembler = nullptr);
    
    void set_model(FEModel* model) { model_ = model; }
    void set_assembler(Assembler* assem) { assembler_ = assem; }
    void set_config(const SolverConfig& config) { config_ = config; }
    
    SolverConfig& config() { return config_; }
    const SolverConfig& config() const { return config_; }
    
    void initialize();
    
    bool solve_static_nonlinear(
        const Vector& F_ext,
        const Vector& u_guess,
        SolverResult& result);
    
    bool solve_static_linear(
        const Vector& F_ext,
        SolverResult& result);
    
    bool solve_transient_newmark(
        Real dt,
        const Vector& F_ext,
        const Vector& u_prev,
        const Vector& v_prev,
        const Vector& a_prev,
        SolverResult& result,
        Real beta = 0.25,
        Real gamma = 0.5);
    
    bool solve_modal(
        Index num_modes,
        std::vector<Real>& frequencies,
        std::vector<Vector>& mode_shapes);
    
    Vector compute_reactions(const Vector& displacement);
    
    void compute_strains_stresses(
        const Vector& displacement,
        std::vector<Real>& elem_strain,
        std::vector<Real>& elem_stress,
        std::vector<Vec3>& elem_force);
    
    bool solve_linear_system(
        SparseMatrix& K,
        Vector& solution,
        const Vector& rhs);
    
    bool apply_newmark_time_integration(
        SparseMatrix& K,
        SparseMatrix& M,
        SparseMatrix& C,
        const Vector& F_ext,
        Vector& u,
        Vector& v,
        Vector& a,
        Real dt,
        Real beta = 0.25,
        Real gamma = 0.5);

protected:
    FEModel* model_ = nullptr;
    Assembler* assembler_ = nullptr;
    std::unique_ptr<Assembler> owned_assembler_;
    SolverConfig config_;
    bool initialized_ = false;
    
#ifdef HVDC_USE_PETSC
    PETScSolver petsc_solver_;
#endif
    
    bool solve_linear_system_impl(
        const SparseMatrix& K,
        const Vector& rhs,
        Vector& solution);
    
    bool newton_raphson_iteration(
        const Vector& F_ext,
        const Vector& F_extra,
        SolverResult& result,
        bool is_transient = false,
        Real dt = 0.0,
        const Vector* u_prev = nullptr,
        const Vector* v_prev = nullptr,
        const Vector* a_prev = nullptr,
        Real beta = 0.25,
        Real gamma = 0.5);
    
    Real line_search(
        Vector& du,
        const Vector& residual,
        const Vector& F_ext_total,
        AssemblyResult& assembly);
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_SOLVER_HPP
