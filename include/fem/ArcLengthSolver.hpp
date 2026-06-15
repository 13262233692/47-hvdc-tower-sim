#ifndef HVDC_FEM_ARC_LENGTH_SOLVER_HPP
#define HVDC_FEM_ARC_LENGTH_SOLVER_HPP

#include "common/Types.hpp"
#include "common/Vector.hpp"
#include "common/SparseMatrix.hpp"
#include "fem/FEModel.hpp"
#include "fem/Assembler.hpp"

#include <string>
#include <memory>
#include <vector>

namespace hvdc {
namespace fem {

enum class ArcLengthMethod {
    CRISFIELD_CYLINDRICAL,
    CRISFIELD_SPHERICAL,
    RIKS_WEMPENER,
    DISPLACEMENT_CONTROL
};

struct ArcLengthConfig {
    Real tol_abs = 1.0e-8;
    Real tol_rel = 1.0e-6;
    Index max_newton_iters = 50;
    Index max_path_steps = 2000;
    
    Real initial_arc_length = 0.1;
    Real min_arc_length = 1.0e-6;
    Real max_arc_length = 10.0;
    
    Index ideal_newton_iters = 4;
    Real arc_length_growth_factor = 1.5;
    Real arc_length_shrink_factor = 0.5;
    
    Real load_scaling_alpha = 1.0;
    
    ArcLengthMethod method = ArcLengthMethod::CRISFIELD_CYLINDRICAL;
    
    bool auto_switch_direction = true;
    Real singularity_det_threshold = 1.0e-14;
    
    Index max_line_search = 5;
    Real line_search_alpha = 0.5;
    Real line_search_beta = 0.8;
    
    bool compute_eigenvalue_at_critical = true;
    bool save_history = true;
    
    std::string ksp_type = "gmres";
    std::string pc_type = "gamg";
    Real krylov_rtol = 1.0e-8;
    Real krylov_atol = 1.0e-50;
    Index krylov_max_iters = 20000;
    
    Real target_load_factor = -1.0;
    Real target_displacement_norm = -1.0;
    Real max_load_factor = -1.0;
    Real max_displacement_norm = -1.0;
};

struct ArcLengthStepResult {
    Index step_id = 0;
    Real load_factor = 0.0;
    Real delta_load_factor = 0.0;
    Real displacement_norm = 0.0;
    Real delta_displacement_norm = 0.0;
    Real arc_length = 0.0;
    
    Vector displacement;
    Vector residual;
    Real residual_norm = 0.0;
    Real initial_residual_norm = 0.0;
    Index newton_iters = 0;
    Index total_linear_iters = 0;
    
    bool converged = false;
    std::string converged_reason;
    bool is_critical_point = false;
    Real stiffness_det_sign = 1.0;
    Real min_eigenvalue_estimate = 0.0;
    
    bool path_follow_failed = false;
    std::string fail_reason;
    
    Real step_time = 0.0;
    std::vector<Real> newton_residual_history;
};

struct ArcLengthResult {
    bool success = false;
    std::string fail_reason;
    
    std::vector<ArcLengthStepResult> steps;
    
    Real total_time = 0.0;
    Index total_steps = 0;
    Index total_newton_iters = 0;
    Index total_linear_iters = 0;
    
    Index critical_points_detected = 0;
    Real collapse_load_factor = -1.0;
    Real collapse_displacement_norm = -1.0;
    
    std::vector<Real> load_factor_history;
    std::vector<Real> displacement_norm_history;
    std::vector<Real> arc_length_history;
    std::vector<Index> newton_iter_history;
    
    Real final_load_factor = 0.0;
    Real final_displacement_norm = 0.0;
    Vector final_displacement;
};

class ArcLengthSolver {
public:
    ArcLengthSolver() = default;
    explicit ArcLengthSolver(FEModel* model, Assembler* assembler = nullptr);
    
    void set_model(FEModel* model) { model_ = model; }
    void set_assembler(Assembler* assem) { assembler_ = assem; }
    void set_config(const ArcLengthConfig& config) { config_ = config; }
    
    ArcLengthConfig& config() { return config_; }
    const ArcLengthConfig& config() const { return config_; }
    
    void initialize();
    
    bool solve_path(
        const Vector& F_reference,
        const Vector& u_initial,
        ArcLengthResult& result);
    
    bool solve_single_step(
        const Vector& F_reference,
        Vector& u_current,
        Real& lambda_current,
        Real& delta_lambda_prev,
        Vector& delta_u_prev,
        Real current_arc_length,
        ArcLengthStepResult& step_result);
    
protected:
    FEModel* model_ = nullptr;
    Assembler* assembler_ = nullptr;
    std::unique_ptr<Assembler> owned_assembler_;
    ArcLengthConfig config_;
    bool initialized_ = false;
    
    Real F_ref_norm_sq_ = 0.0;
    
    bool solve_linear_system_impl(
        const SparseMatrix& K,
        const Vector& rhs,
        Vector& solution);
    
    bool crisfield_iteration(
        const Vector& F_ref,
        Vector& u,
        Real& lambda,
        const Vector& delta_u_prev,
        Real delta_lambda_prev,
        Real arc_length,
        ArcLengthStepResult& result);
    
    Real compute_constraint_equation(
        const Vector& delta_u,
        Real delta_lambda,
        Real arc_length) const;
    
    Real compute_constraint_derivative_du(
        const Vector& delta_u,
        const Vector& du,
        Real delta_lambda,
        Real arc_length) const;
    
    bool detect_critical_point(
        const SparseMatrix& K_T,
        Real& det_sign,
        Real& min_eig_estimate);
    
    Real adaptive_arc_length(
        Real current_arc_length,
        Index newton_iters_used) const;
    
    Real arc_length_line_search(
        const Vector& delta_u_trial,
        Real delta_lambda_trial,
        const Vector& residual,
        Real constraint_residual,
        const Vector& F_ref,
        const Vector& u_current,
        Real lambda_current,
        AssemblyResult& assembly);
    
    void save_step_history(
        ArcLengthResult& result,
        const ArcLengthStepResult& step) const;
    
    bool check_termination_criteria(
        const ArcLengthStepResult& step) const;
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_ARC_LENGTH_SOLVER_HPP
