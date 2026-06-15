#ifndef HVDC_FVM_FLUIDSOLVER_HPP
#define HVDC_FVM_FLUIDSOLVER_HPP

#include "common/Types.hpp"
#include "common/Vector.hpp"
#include "common/SparseMatrix.hpp"
#ifdef HVDC_USE_PETSC
#include "common/PETScManager.hpp"
#endif
#include "fvm/Grid.hpp"
#include "fvm/Field.hpp"
#include "fvm/Discretization.hpp"
#include "fvm/NSComposable.hpp"
#include "fvm/BoundaryCondition.hpp"

#include <memory>
#include <string>
#include <vector>

namespace hvdc {
namespace fvm {

enum class SolverAlgorithm : UInt8 {
    SIMPLE = 0,
    SIMPLEC = 1,
    PISO = 2,
    PIMPLE = 3,
    SteadySIMPLE = 4
};

struct FluidSolverConfig {
    SolverAlgorithm algorithm = SolverAlgorithm::SIMPLE;
    Index max_correctors = 6;
    Index max_non_orthogonal_correctors = 0;
    Index pimple_outer_iters = 3;
    Real residual_tol = 1.0e-4;
    Real residual_rel_tol = 1.0e-3;
    Real pressure_tol = 1.0e-3;
    
    Index max_linear_iters_U = 2000;
    Index max_linear_iters_p = 4000;
    Index max_linear_iters_k = 2000;
    Index max_linear_iters_eps = 2000;
    
    Real linear_rtol_U = 1.0e-5;
    Real linear_rtol_p = 1.0e-7;
    Real linear_rtol_turb = 1.0e-5;
    
    bool enable_turbulence = false;
    Index turbulence_start_iter = 5;
    Index turbulence_update_interval = 1;
    
    std::string linear_solver_U = "PBiCGStab";
    std::string linear_solver_p = "GAMG";
    std::string linear_solver_turb = "PBiCGStab";
    
    std::string preconditioner_U = "DILU";
    std::string preconditioner_p = "GAMG";
    std::string preconditioner_turb = "DILU";
    
    Real adaptive_dt_cfl = 0.5;
    Real adaptive_dt_max = 1.0;
    Real adaptive_dt_min = 1.0e-6;
    bool use_adaptive_dt = false;
    
    bool compute_forces_every_step = true;
    Index output_interval = 10;
};

struct FluidIterationStats {
    Real initial_res_U[3] = {0, 0, 0};
    Real final_res_U[3] = {0, 0, 0};
    Real final_res_p = 0;
    Real initial_res_p = 0;
    Real mass_imbalance = 0;
    Real cfl = 0;
    Index total_linear_iters = 0;
    Index corrector_iters = 0;
    Real solve_time = 0;
    std::vector<Real> p_residual_history;
    std::vector<Real> u_residual_history;
};

struct FluidSolutionResult {
    bool converged = false;
    std::string reason;
    Index total_iters = 0;
    Index total_linear_iters = 0;
    Real total_time = 0;
    FluidIterationStats final_stats;
    std::vector<FluidIterationStats> iteration_history;
    std::vector<Vec3> aerodynamic_forces;
};

class FluidSolver {
public:
    FluidSolver() = default;
    FluidSolver(Grid* grid, const FluidSolverConfig& cfg,
                const FlowConfig& flow_cfg = FlowConfig(),
                const DiscretizationConfig& disc_cfg = DiscretizationConfig());
    
    void set_grid(Grid* grid) { grid_ = grid; ns_.set_grid(grid); }
    void set_config(const FluidSolverConfig& cfg) { config_ = cfg; }
    void set_boundary_conditions(BoundaryConditionManager* bcm) { bc_mgr_ = bcm; }
    
    const FluidSolverConfig& config() const { return config_; }
    BoundaryConditionManager* boundary_manager() { return bc_mgr_; }
    Grid* grid() { return grid_; }
    
    FluidState& state() { return state_; }
    const FluidState& state() const { return state_; }
    FluidState& state_prev() { return state_prev_; }
    
    NSComposable& ns() { return ns_; }
    const NSComposable& ns() const { return ns_; }
    
    void initialize();
    void reset_solution();
    
    void set_free_stream_velocity(const Vec3& U_inf);
    void initialize_fields();
    FluidState& current_state() { return state_; }
    const FluidState& current_state() const { return state_; }
    
    bool solve_steady(Index max_iters = 1000, Real tol = 1.0e-6);
    bool solve_steady(FluidSolutionResult& result, Index max_iters = 1000);
    bool solve_transient_step(Real dt, Index max_iters = 50, Real tol = 1.0e-5);
    bool solve_transient_step(
        Real dt,
        FluidSolutionResult& result,
        const FluidState* state_old = nullptr);
    
    FluidIterationStats perform_simple_iteration(Index iter_count = 0);
    FluidIterationStats perform_piso_iteration(Real dt, Index n_correctors = 2);
    
    void save_state();
    void restore_state();
    
    void advance_time(Real dt);
    
    Vec3 compute_aerodynamic_force_on_patch(const std::string& patch_name);
    std::vector<std::pair<std::string, Vec3>> compute_all_boundary_forces();
    
    std::vector<Real> get_face_pressure_on_patch(const std::string& patch_name) const;
    std::vector<Vec3> get_face_velocity_on_patch(const std::string& patch_name) const;
    
    void update_interface_velocity(const std::vector<Index>& face_ids,
                                    const std::vector<Vec3>& velocities);
    
    Real compute_suggested_dt() const;
    
    void print_iteration_report(const FluidIterationStats& stats, Index iter) const;

protected:
    Grid* grid_ = nullptr;
    FluidSolverConfig config_;
    BoundaryConditionManager* bc_mgr_ = nullptr;
    
    NSComposable ns_;
    FluidState state_;
    FluidState state_prev_;
    FluidState state_old_;
    Vec3 U_inf_{0.0, 0.0, 0.0};
    
    std::vector<Vec3> face_velocity_;
    std::vector<Real> face_mass_flux_;
    std::vector<Real> face_pressure_;
    std::vector<std::array<Vec3, 3>> grad_U_;
    std::vector<Vec3> grad_p_;
    
    SparseMatrix A_U[3];
    SparseMatrix A_p;
    SparseMatrix A_k, A_eps;
    Vector b_U[3];
    Vector b_p;
    Vector b_k, b_eps;
    Vector diag_inv_Ap[3];
    VectorField HbyA_;
    
#ifdef HVDC_USE_PETSC
    PETScSolver petsc_solver_U_;
    PETScSolver petsc_solver_p_;
    PETScSolver petsc_solver_turb_;
#endif
    
    bool initialized_ = false;
    Index current_iter_ = 0;
    
    void solve_momentum_predictor(Index component);
    void solve_pressure_correction();
    void solve_turbulence_equations(Real dt);
    
    void apply_all_boundary_conditions();
    void update_gradients();
    
    void under_relax_momentum();
    void under_relax_pressure();
    void under_relax_turbulence();
    
    bool check_convergence(const FluidIterationStats& stats);
    
    Real solve_linear_system(SparseMatrix& A, Vector& x, const Vector& b,
                              const std::string& solver_name,
                              const std::string& prec_name,
                              Index max_iters, Real rtol,
                              Real atol = 1.0e-50, Index* iters = nullptr);
};

} // namespace fvm
} // namespace hvdc

#endif // HVDC_FVM_FLUIDSOLVER_HPP
