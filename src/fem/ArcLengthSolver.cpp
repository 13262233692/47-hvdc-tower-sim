#include "fem/ArcLengthSolver.hpp"
#include "fem/Solver.hpp"
#include "common/MathUtils.hpp"
#include "common/Timer.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>

extern hvdc::Timer g_timer;

namespace hvdc {
namespace fem {

ArcLengthSolver::ArcLengthSolver(FEModel* model, Assembler* assembler)
    : model_(model), assembler_(assembler)
{
    if (!assembler_ && model_) {
        owned_assembler_ = std::make_unique<Assembler>(model_);
        assembler_ = owned_assembler_.get();
    }
}

void ArcLengthSolver::initialize() {
    if (!model_) return;
    if (!assembler_) {
        owned_assembler_ = std::make_unique<Assembler>(model_);
        assembler_ = owned_assembler_.get();
    }
    assembler_->initialize();
    initialized_ = true;
    
    HVDC_LOG_INFO("ArcLengthSolver initialized: method=" 
                  << static_cast<int>(config_.method));
}

bool ArcLengthSolver::solve_linear_system_impl(
    const SparseMatrix& K,
    const Vector& rhs,
    Vector& solution)
{
    Index n = rhs.size();
    solution.resize(n);
    solution.zero();
    
    const auto& rp = K.row_ptr();
    const auto& ci = K.col_idx();
    const auto& v = K.values();
    
    std::vector<Real> diag(n, 0.0);
    for (Index i = 0; i < n; ++i) {
        for (Index j = rp[i]; j < rp[i + 1]; ++j) {
            if (ci[j] == i) {
                diag[i] = v[j];
                break;
            }
        }
    }
    
    Vector r = rhs;
    for (Index iter = 0; iter < config_.krylov_max_iters; ++iter) {
        Real r_norm = r.norm();
        if (r_norm < config_.krylov_atol) return true;
        if (iter > 0 && r_norm < config_.krylov_rtol * rhs.norm()) return true;
        
        Vector p(n);
        for (Index i = 0; i < n; ++i) {
            p[i] = (std::fabs(diag[i]) > EPS) ? r[i] / diag[i] : r[i];
        }
        
        Vector Kp(n);
        K.multiply(p, Kp);
        
        Real pKp = p.dot(Kp);
        if (std::fabs(pKp) < EPS) {
            for (Index i = 0; i < n; ++i) {
                solution[i] += p[i];
            }
            break;
        }
        
        Real alpha = p.dot(r) / pKp;
        for (Index i = 0; i < n; ++i) {
            solution[i] += alpha * p[i];
            r[i] -= alpha * Kp[i];
        }
    }
    
    return true;
}

Real ArcLengthSolver::compute_constraint_equation(
    const Vector& delta_u,
    Real delta_lambda,
    Real arc_length) const
{
    Real du_sq = delta_u.dot(delta_u);
    Real dl_sq = delta_lambda * delta_lambda * F_ref_norm_sq_;
    Real alpha_sq = config_.load_scaling_alpha * config_.load_scaling_alpha;
    
    switch (config_.method) {
        case ArcLengthMethod::CRISFIELD_CYLINDRICAL:
            return du_sq + alpha_sq * dl_sq - arc_length * arc_length;
            
        case ArcLengthMethod::CRISFIELD_SPHERICAL:
            return du_sq + alpha_sq * dl_sq - arc_length * arc_length;
            
        case ArcLengthMethod::RIKS_WEMPENER:
            return du_sq + alpha_sq * dl_sq - arc_length * arc_length;
            
        case ArcLengthMethod::DISPLACEMENT_CONTROL:
            return du_sq - arc_length * arc_length;
            
        default:
            return du_sq + alpha_sq * dl_sq - arc_length * arc_length;
    }
}

Real ArcLengthSolver::compute_constraint_derivative_du(
    const Vector& delta_u,
    const Vector& du,
    Real delta_lambda,
    Real arc_length) const
{
    (void)arc_length;
    Real d_du = 2.0 * delta_u.dot(du);
    
    if (config_.method == ArcLengthMethod::DISPLACEMENT_CONTROL) {
        return d_du;
    }
    
    Real alpha_sq = config_.load_scaling_alpha * config_.load_scaling_alpha;
    return d_du + 2.0 * alpha_sq * delta_lambda * F_ref_norm_sq_;
}

bool ArcLengthSolver::detect_critical_point(
    const SparseMatrix& K_T,
    Real& det_sign,
    Real& min_eig_estimate)
{
    const auto& rp = K_T.row_ptr();
    const auto& ci = K_T.col_idx();
    const auto& v = K_T.values();
    
    Index n = K_T.rows();
    
    Real prod_sign = 1.0;
    Real min_diag = std::numeric_limits<Real>::max();
    Real max_diag = -std::numeric_limits<Real>::max();
    
    for (Index i = 0; i < n; ++i) {
        Real diag_val = 0.0;
        for (Index j = rp[i]; j < rp[i + 1]; ++j) {
            if (ci[j] == i) {
                diag_val = v[j];
                break;
            }
        }
        if (diag_val < 0.0) {
            prod_sign *= -1.0;
        }
        min_diag = std::min(min_diag, diag_val);
        max_diag = std::max(max_diag, diag_val);
    }
    
    det_sign = prod_sign;
    min_eig_estimate = min_diag;
    
    if (min_diag < config_.singularity_det_threshold && 
        max_diag > config_.singularity_det_threshold) {
        return true;
    }
    
    return false;
}

Real ArcLengthSolver::adaptive_arc_length(
    Real current_arc_length,
    Index newton_iters_used) const
{
    Real ratio = static_cast<Real>(config_.ideal_newton_iters) 
                 / static_cast<Real>(std::max<Index>(newton_iters_used, 1));
    
    Real factor = std::pow(
        (ratio > 1.0) ? config_.arc_length_growth_factor 
                      : config_.arc_length_shrink_factor,
        std::fabs(ratio - 1.0));
    
    Real new_arc_length = current_arc_length * factor;
    
    new_arc_length = std::max(new_arc_length, config_.min_arc_length);
    new_arc_length = std::min(new_arc_length, config_.max_arc_length);
    
    return new_arc_length;
}

Real ArcLengthSolver::arc_length_line_search(
    const Vector& delta_u_trial,
    Real delta_lambda_trial,
    const Vector& residual,
    Real constraint_residual,
    const Vector& F_ref,
    const Vector& u_current,
    Real lambda_current,
    AssemblyResult&)
{
    (void)constraint_residual;
    (void)F_ref;
    
    Real alpha = 1.0;
    Index n = u_current.size();
    
    Real initial_merit = 0.5 * residual.dot(residual);
    if (initial_merit < EPS) return alpha;
    
    Vector u_test(n);
    Vector F_int_test(n);
    Vector resid_test(n);
    
    for (Index ls_iter = 0; ls_iter < config_.max_line_search; ++ls_iter) {
        for (Index j = 0; j < n; ++j) {
            u_test[j] = u_current[j] + alpha * delta_u_trial[j];
        }
        
        assembler_->assemble_internal_forces(u_test, F_int_test, false);
        
        Real lambda_test = lambda_current + alpha * delta_lambda_trial;
        for (Index j = 0; j < n; ++j) {
            resid_test[j] = lambda_test * F_ref[j] - F_int_test[j];
        }
        
        for (Index gdof : model_->constrained_dofs()) {
            resid_test.set(gdof, 0.0);
        }
        
        Real test_merit = 0.5 * resid_test.dot(resid_test);
        
        Real target = initial_merit - config_.line_search_beta * alpha * initial_merit;
        
        if (test_merit <= target + EPS) {
            return alpha;
        }
        
        alpha *= config_.line_search_alpha;
    }
    
    return alpha;
}

bool ArcLengthSolver::crisfield_iteration(
    const Vector& F_ref,
    Vector& u,
    Real& lambda,
    const Vector& delta_u_prev,
    Real delta_lambda_prev,
    Real arc_length,
    ArcLengthStepResult& result)
{
    Index n = model_->total_dofs();
    
    AssemblyResult assembly;
    
    Vector delta_u = delta_u_prev;
    Real delta_lambda = delta_lambda_prev;
    
    Vector u_prev_step = u;
    Real lambda_prev_step = lambda;
    
    Vector residual(n);
    Vector scaled_F_ref(n);
    for (Index i = 0; i < n; ++i) {
        scaled_F_ref[i] = lambda * F_ref[i];
    }
    
    assembler_->assemble_tangent_stiffness(u, assembly, true, false, true);
    
    for (Index i = 0; i < n; ++i) {
        residual[i] = lambda * F_ref[i] - assembly.F_int[i];
    }
    for (Index gdof : model_->constrained_dofs()) {
        residual.set(gdof, 0.0);
    }
    
    result.initial_residual_norm = residual.norm();
    result.residual_norm = result.initial_residual_norm;
    
    if (result.initial_residual_norm < config_.tol_abs) {
        result.converged = true;
        result.converged_reason = "initial_residual";
        return true;
    }
    
    assembler_->apply_dirichlet_bc(assembly.K_T, residual, u);
    
    Vector du_R(n);
    solve_linear_system_impl(assembly.K_T, residual, du_R);
    
    for (Index gdof : model_->constrained_dofs()) {
        du_R.set(gdof, 0.0);
    }
    
    Vector du_F(n);
    solve_linear_system_impl(assembly.K_T, F_ref, du_F);
    
    for (Index gdof : model_->constrained_dofs()) {
        du_F.set(gdof, 0.0);
    }
    
    Real alpha_sq = config_.load_scaling_alpha * config_.load_scaling_alpha;
    
    bool direction_switched = false;
    if (config_.auto_switch_direction && delta_u_prev.size() > 0) {
        Real dot_product = du_R.dot(delta_u_prev);
        if (dot_product < 0.0) {
            for (Index i = 0; i < n; ++i) {
                du_R[i] = -du_R[i];
            }
            direction_switched = true;
        }
    }
    
    result.converged = false;
    Index iter;
    
    for (iter = 0; iter < config_.max_newton_iters; ++iter) {
        result.newton_residual_history.push_back(result.residual_norm);
        
        assembler_->apply_dirichlet_bc(assembly.K_T, residual, u);
        
        du_R.zero();
        solve_linear_system_impl(assembly.K_T, residual, du_R);
        
        for (Index gdof : model_->constrained_dofs()) {
            du_R.set(gdof, 0.0);
        }
        
        du_F.zero();
        solve_linear_system_impl(assembly.K_T, F_ref, du_F);
        
        for (Index gdof : model_->constrained_dofs()) {
            du_F.set(gdof, 0.0);
        }
        
        Real a = 2.0 * delta_u.dot(du_R);
        Real b = 2.0 * delta_u.dot(du_F) 
                 + 2.0 * alpha_sq * delta_lambda * F_ref_norm_sq_;
        
        Real g = compute_constraint_equation(delta_u, delta_lambda, arc_length);
        
        Real delta_lambda_new = 0.0;
        if (std::fabs(b) > EPS) {
            delta_lambda_new = -(g + a) / b;
        } else {
            delta_lambda_new = 0.0;
        }
        
        if (direction_switched && iter == 0) {
            delta_lambda_new = -delta_lambda_new;
        }
        
        Vector delta_u_trial(n);
        for (Index i = 0; i < n; ++i) {
            delta_u_trial[i] = du_R[i] + delta_lambda_new * du_F[i];
        }
        
        Real line_alpha = arc_length_line_search(
            delta_u_trial, delta_lambda_new, 
            residual, g, F_ref, u, lambda, assembly);
        
        for (Index i = 0; i < n; ++i) {
            delta_u[i] += line_alpha * delta_u_trial[i];
            u[i] += line_alpha * delta_u_trial[i];
        }
        delta_lambda += line_alpha * delta_lambda_new;
        lambda = lambda_prev_step + delta_lambda;
        
        assembler_->assemble_tangent_stiffness(u, assembly, true, false, true);
        
        for (Index i = 0; i < n; ++i) {
            residual[i] = lambda * F_ref[i] - assembly.F_int[i];
        }
        for (Index gdof : model_->constrained_dofs()) {
            residual.set(gdof, 0.0);
        }
        
        result.residual_norm = residual.norm();
        
        Real rel = result.initial_residual_norm > EPS 
                   ? result.residual_norm / result.initial_residual_norm
                   : result.residual_norm;
        
        HVDC_LOG_DEBUG("Crisfield iter " << iter 
                       << ": ||R||=" << result.residual_norm
                       << " (rel=" << rel << ")"
                       << ", λ=" << lambda
                       << ", Δλ=" << delta_lambda);
        
        if (result.residual_norm < config_.tol_abs || rel < config_.tol_rel) {
            result.converged = true;
            result.converged_reason = "tolerance";
            break;
        }
    }
    
    result.newton_iters = iter;
    result.residual = residual;
    result.displacement = u;
    result.delta_load_factor = delta_lambda;
    result.load_factor = lambda;
    result.delta_displacement_norm = delta_u.norm();
    result.displacement_norm = u.norm();
    
    Real det_sign;
    Real min_eig;
    result.is_critical_point = detect_critical_point(assembly.K_T, det_sign, min_eig);
    result.stiffness_det_sign = det_sign;
    result.min_eigenvalue_estimate = min_eig;
    
    if (result.converged) {
        HVDC_LOG_INFO("Crisfield converged after " << iter 
                      << " iters: ||R||=" << result.residual_norm
                      << ", λ=" << lambda
                      << ", critical=" << (result.is_critical_point ? "YES" : "no"));
    } else {
        HVDC_LOG_WARNING("Crisfield did NOT converge after " << iter
                         << " iters: ||R||=" << result.residual_norm);
    }
    
    return result.converged;
}

bool ArcLengthSolver::solve_single_step(
    const Vector& F_reference,
    Vector& u_current,
    Real& lambda_current,
    Real& delta_lambda_prev,
    Vector& delta_u_prev,
    Real current_arc_length,
    ArcLengthStepResult& step_result)
{
    step_result.arc_length = current_arc_length;
    
    return crisfield_iteration(
        F_reference,
        u_current,
        lambda_current,
        delta_u_prev,
        delta_lambda_prev,
        current_arc_length,
        step_result);
}

void ArcLengthSolver::save_step_history(
    ArcLengthResult& result,
    const ArcLengthStepResult& step) const
{
    if (!config_.save_history) return;
    
    result.steps.push_back(step);
    result.load_factor_history.push_back(step.load_factor);
    result.displacement_norm_history.push_back(step.displacement_norm);
    result.arc_length_history.push_back(step.arc_length);
    result.newton_iter_history.push_back(step.newton_iters);
    
    result.total_steps++;
    result.total_newton_iters += step.newton_iters;
    
    if (step.is_critical_point) {
        result.critical_points_detected++;
        if (result.collapse_load_factor < 0.0 || 
            std::fabs(step.load_factor) < std::fabs(result.collapse_load_factor)) {
            result.collapse_load_factor = step.load_factor;
            result.collapse_displacement_norm = step.displacement_norm;
        }
    }
}

bool ArcLengthSolver::check_termination_criteria(
    const ArcLengthStepResult& step) const
{
    if (config_.target_load_factor > 0.0) {
        if (std::fabs(step.load_factor) >= config_.target_load_factor) {
            return true;
        }
    }
    
    if (config_.target_displacement_norm > 0.0) {
        if (step.displacement_norm >= config_.target_displacement_norm) {
            return true;
        }
    }
    
    if (config_.max_load_factor > 0.0) {
        if (std::fabs(step.load_factor) >= config_.max_load_factor) {
            return true;
        }
    }
    
    if (config_.max_displacement_norm > 0.0) {
        if (step.displacement_norm >= config_.max_displacement_norm) {
            return true;
        }
    }
    
    return false;
}

bool ArcLengthSolver::solve_path(
    const Vector& F_reference,
    const Vector& u_initial,
    ArcLengthResult& result)
{
    if (!initialized_) {
        HVDC_LOG_ERROR("ArcLengthSolver not initialized");
        result.success = false;
        result.fail_reason = "not_initialized";
        return false;
    }
    
    Index n = model_->total_dofs();
    
    F_ref_norm_sq_ = F_reference.dot(F_reference);
    if (F_ref_norm_sq_ < EPS) {
        F_ref_norm_sq_ = 1.0;
    }
    
    Vector u_current = u_initial;
    Real lambda_current = 0.0;
    Real current_arc_length = config_.initial_arc_length;
    
    Vector delta_u_prev(n);
    delta_u_prev.zero();
    Real delta_lambda_prev = 0.0;
    
    Vector initial_F_int(n);
    assembler_->assemble_internal_forces(u_current, initial_F_int, false);
    
    {
        AssemblyResult assembly;
        assembler_->assemble_tangent_stiffness(u_current, assembly, true, false, true);
        
        Vector initial_rhs(n);
        for (Index i = 0; i < n; ++i) {
            initial_rhs[i] = F_reference[i];
        }
        
        assembler_->apply_dirichlet_bc(assembly.K_T, initial_rhs, u_current);
        
        Vector delta_u_initial(n);
        delta_u_initial.zero();
        solve_linear_system_impl(assembly.K_T, initial_rhs, delta_u_initial);
        
        for (Index gdof : model_->constrained_dofs()) {
            delta_u_initial.set(gdof, 0.0);
        }
        
        Real du_norm = delta_u_initial.norm();
        if (du_norm > EPS) {
            Real scale = current_arc_length / du_norm;
            for (Index i = 0; i < n; ++i) {
                delta_u_prev[i] = delta_u_initial[i] * scale;
            }
            delta_lambda_prev = scale;
            
            for (Index i = 0; i < n; ++i) {
                u_current[i] += delta_u_prev[i];
            }
            lambda_current = delta_lambda_prev;
        }
    }
    
    Timer path_timer;
    path_timer.tic();
    
    Index failed_attempts = 0;
    const Index max_failed_attempts = 10;
    
    for (Index step_id = 0; step_id < config_.max_path_steps; ++step_id) {
        ArcLengthStepResult step_result;
        step_result.step_id = step_id;
        
        Timer step_timer;
        step_timer.tic();
        
        Vector u_saved = u_current;
        Real lambda_saved = lambda_current;
        
        bool step_ok = crisfield_iteration(
            F_reference,
            u_current,
            lambda_current,
            delta_u_prev,
            delta_lambda_prev,
            current_arc_length,
            step_result);
        
        step_result.step_time = step_timer.toc();
        
        if (step_ok) {
            delta_u_prev.zero();
            for (Index i = 0; i < n; ++i) {
                delta_u_prev[i] = u_current[i] - u_saved[i];
            }
            delta_lambda_prev = lambda_current - lambda_saved;
            
            save_step_history(result, step_result);
            
            HVDC_LOG_INFO("Step " << step_id 
                          << ": λ=" << step_result.load_factor
                          << ", ||u||=" << step_result.displacement_norm
                          << ", Δs=" << current_arc_length
                          << ", iters=" << step_result.newton_iters
                          << ", critical=" << (step_result.is_critical_point ? "YES" : "no"));
            
            failed_attempts = 0;
            
            current_arc_length = adaptive_arc_length(
                current_arc_length, step_result.newton_iters);
            
            if (check_termination_criteria(step_result)) {
                result.success = true;
                HVDC_LOG_INFO("Path following completed: reached target criteria");
                break;
            }
            
            if (step_result.is_critical_point && 
                config_.auto_switch_direction) {
                delta_lambda_prev = -delta_lambda_prev;
                for (Index i = 0; i < n; ++i) {
                    delta_u_prev[i] = -delta_u_prev[i];
                }
                HVDC_LOG_INFO("Critical point detected: switching path direction");
            }
            
        } else {
            failed_attempts++;
            HVDC_LOG_WARNING("Step " << step_id 
                             << " failed, retrying with Δs=" 
                             << current_arc_length * config_.arc_length_shrink_factor);
            
            u_current = u_saved;
            lambda_current = lambda_saved;
            current_arc_length *= config_.arc_length_shrink_factor;
            
            if (current_arc_length < config_.min_arc_length || 
                failed_attempts >= max_failed_attempts) {
                result.success = false;
                result.fail_reason = "path_divergence";
                HVDC_LOG_ERROR("Path following diverged after " 
                               << step_id << " steps");
                break;
            }
        }
    }
    
    result.total_time = path_timer.toc();
    
    if (result.steps.size() > 0) {
        const auto& last = result.steps.back();
        result.final_load_factor = last.load_factor;
        result.final_displacement_norm = last.displacement_norm;
        result.final_displacement = last.displacement;
    }
    
    if (result.success) {
        HVDC_LOG_INFO("Arc-length path following completed successfully: "
                      << result.total_steps << " steps, "
                      << result.total_newton_iters << " Newton iters, "
                      << "time=" << result.total_time << "s");
    }
    
    return result.success;
}

} // namespace fem
} // namespace hvdc
