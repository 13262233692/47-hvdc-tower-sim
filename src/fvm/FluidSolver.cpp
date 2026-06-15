#include "fvm/FluidSolver.hpp"
#include "common/MathUtils.hpp"
#include "common/MPIManager.hpp"
#include "common/Timer.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>

extern hvdc::Timer g_timer;

namespace hvdc {
namespace fvm {

FluidSolver::FluidSolver(Grid* grid, const FluidSolverConfig& cfg,
                         const FlowConfig& flow_cfg,
                         const DiscretizationConfig& disc_cfg)
    : grid_(grid), config_(cfg), ns_(grid, flow_cfg, disc_cfg)
{
}

void FluidSolver::initialize() {
    if (!grid_) {
        HVDC_LOG_ERROR("FluidSolver::initialize: grid is null");
        return;
    }
    
    ns_.initialize_state(state_);
    state_prev_ = state_;
    state_old_ = state_;
    
    Index nc = grid_->num_cells();
    Index nf = grid_->num_faces();
    
    face_velocity_.resize(nf);
    face_mass_flux_.resize(nf, 0.0);
    face_pressure_.resize(nf, 0.0);
    grad_U_.resize(nc);
    grad_p_.resize(nc);
    HbyA_.resize(nc);
    
    for (Index comp = 0; comp < 3; ++comp) {
        diag_inv_Ap[comp].resize(nc);
    }
    
#ifdef HVDC_USE_PETSC
    petsc_solver_U_.create();
    petsc_solver_p_.create();
    petsc_solver_turb_.create();
#endif
    
    ns_.set_flow_config(ns_.flow_config());
    
    apply_all_boundary_conditions();
    update_gradients();
    
    for (Index fi = 0; fi < nf; ++fi) {
        const Face& face = grid_->face(fi);
        Index P = face.owner();
        Index N = face.neighbor();
        Vec3 Sf = math::vec3_scale(face.normal(), face.area());
        if (face.is_on_boundary()) {
            Index ci = (P >= 0) ? P : N;
            face_velocity_[fi] = state_.U[ci];
            face_mass_flux_[fi] = state_.rho[ci] * math::vec3_dot(state_.U[ci], Sf);
            face_pressure_[fi] = state_.p[ci];
        } else {
            Real gc = face.weight();
            for (Index k = 0; k < 3; ++k) {
                face_velocity_[fi][k] = gc * state_.U[P][k] + (1.0 - gc) * state_.U[N][k];
            }
            Real rho_f = gc * state_.rho[P] + (1.0 - gc) * state_.rho[N];
            face_mass_flux_[fi] = rho_f * math::vec3_dot(face_velocity_[fi], Sf);
            face_pressure_[fi] = gc * state_.p[P] + (1.0 - gc) * state_.p[N];
        }
    }
    
    if (bc_mgr_) {
        bc_mgr_->apply_all_pressure_bc(*grid_, state_, state_.p, face_pressure_);
        bc_mgr_->apply_all_velocity_bc(*grid_, state_, state_.U, face_velocity_);
    }
    
    state_.update_thermophysical_properties();
    
    initialized_ = true;
    current_iter_ = 0;
    
    GridStats stats;
    grid_->compute_stats(stats);
    HVDC_LOG_INFO("Fluid solver initialized: " 
                  << stats.num_cells << " cells, "
                  << stats.num_faces << " faces, "
                  << "ref rho=" << state_.rho_ref
                  << ", ref U=" << math::vec3_norm(state_.U[0]) << "m/s");
}

void FluidSolver::reset_solution() {
    state_ = state_prev_;
    current_iter_ = 0;
    apply_all_boundary_conditions();
}

void FluidSolver::save_state() {
    state_old_ = state_;
}

void FluidSolver::restore_state() {
    state_ = state_old_;
}

void FluidSolver::advance_time(Real dt) {
    (void)dt;
    state_prev_ = state_;
}

void FluidSolver::apply_all_boundary_conditions() {
    if (!bc_mgr_) return;
    bc_mgr_->apply_all_velocity_bc(*grid_, state_, state_.U, face_velocity_);
    bc_mgr_->apply_all_pressure_bc(*grid_, state_, state_.p, face_pressure_);
    bc_mgr_->apply_all_turbulence_bc(*grid_, state_);
    bc_mgr_->apply_all_temperature_bc(*grid_, state_.T);
}

void FluidSolver::update_gradients() {
    ns_.disc().compute_vector_gradient(state_.U, face_velocity_, grad_U_);
    
    Index nc = grid_->num_cells();
    for (Index ci = 0; ci < nc; ++ci) {
        grad_p_[ci] = state_.p.pressure_gradient(ci, *grid_);
    }
}

Real FluidSolver::solve_linear_system(SparseMatrix& A, Vector& x, const Vector& b,
                                       const std::string& solver_name,
                                       const std::string& prec_name,
                                       Index max_iters, Real rtol,
                                       Real atol, Index* iters)
{
#ifdef HVDC_USE_PETSC
    PETScSolver* solver = nullptr;
    if (solver_name == config_.linear_solver_p) {
        solver = &petsc_solver_p_;
    } else if (solver_name == config_.linear_solver_turb) {
        solver = &petsc_solver_turb_;
    } else {
        solver = &petsc_solver_U_;
    }
    
    solver->set_ksp_type(solver_name.empty() ? "gmres" : solver_name);
    solver->set_preconditioner(prec_name.empty() ? "bjacobi" : prec_name);
    solver->set_tolerances(rtol, atol, 1.0e5, max_iters);
    
    PETScMatrix A_petsc;
    PETScVector b_petsc(b.size()), x_petsc(x.size());
    A_petsc.from_sparse_matrix(A);
    b_petsc.from_vector(b);
    x_petsc.from_vector(x);
    
    solver->set_operator(A_petsc);
    solver->solve(b_petsc, x_petsc);
    x_petsc.to_vector(x);
    
    if (iters) *iters = solver->iterations();
    return solver->residual_norm();
#else
    Index n = b.size();
    if (x.size() != n) x.resize(n);
    
    Vector r = b;
    A.multiply(x, r);
    for (Index i = 0; i < n; ++i) r[i] = b[i] - r[i];
    
    Vector diag_inv(n);
    for (Index i = 0; i < n; ++i) {
        Real d = A(i, i);
        diag_inv[i] = (std::fabs(d) > EPS) ? 1.0 / d : 1.0;
    }
    
    Vector p(n), Ap(n);
    for (Index i = 0; i < n; ++i) p[i] = diag_inv[i] * r[i];
    
    Real r_norm0 = r.dot(r);
    if (r_norm0 < atol * atol) {
        if (iters) *iters = 0;
        return 0.0;
    }
    
    Index iter;
    for (iter = 0; iter < max_iters; ++iter) {
        A.multiply(p, Ap);
        Real pAp = p.dot(Ap);
        if (std::fabs(pAp) < EPS) break;
        
        Real r_dot_z = 0.0;
        for (Index i = 0; i < n; ++i) r_dot_z += r[i] * diag_inv[i] * r[i];
        Real alpha = r_dot_z / pAp;
        
        for (Index i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        
        Real r_norm = r.dot(r);
        Real rel = r_norm / r_norm0;
        if (rel < rtol * rtol || r_norm < atol * atol) break;
        
        Real new_r_dot_z = 0.0;
        for (Index i = 0; i < n; ++i) new_r_dot_z += r[i] * diag_inv[i] * r[i];
        Real beta = new_r_dot_z / std::max(r_dot_z, EPS);
        
        for (Index i = 0; i < n; ++i) {
            p[i] = diag_inv[i] * r[i] + beta * p[i];
        }
    }
    
    if (iters) *iters = iter;
    return std::sqrt(r.dot(r) / r_norm0);
#endif
}

void FluidSolver::solve_momentum_predictor(Index component) {
    Index nc = grid_->num_cells();
    
    const FluidState* prev_state = 
        (config_.algorithm == SolverAlgorithm::SIMPLE ||
         config_.algorithm == SolverAlgorithm::SteadySIMPLE) ? nullptr : &state_prev_;
    
    Real dt = 1.0;
    if (config_.algorithm != SolverAlgorithm::SteadySIMPLE &&
        config_.algorithm != SolverAlgorithm::SIMPLE) {
        dt = compute_suggested_dt();
    }
    
    ns_.assemble_momentum_predictor(
        state_, prev_state, dt, component,
        A_U[component], b_U[component], HbyA_,
        diag_inv_Ap[component], grad_p_);
    
    if (bc_mgr_) {
        bc_mgr_->update_linear_system_bc(*grid_, "U", component,
                                          A_U[component], b_U[component], state_.U);
    }
    
    ScalarField U_comp = state_.U.component(component);
    Vector x_vec(U_comp.raw());
    
    Index iters;
    Real res = solve_linear_system(
        A_U[component], x_vec, b_U[component],
        config_.linear_solver_U, config_.preconditioner_U,
        config_.max_linear_iters_U, config_.linear_rtol_U,
        1.0e-50, &iters);
    
    state_.U.set_component(component, ScalarField(x_vec.raw()));
    
    ns_.disc().under_relax_vector_field(state_.U, state_prev_.U, config_.disc_config().relaxation_U);
    
    (void)res;
}

void FluidSolver::solve_pressure_correction() {
    ns_.assemble_pressure_correction(
        state_, A_p, diag_inv_Ap[0], HbyA_,
        A_p, b_p, face_velocity_, face_mass_flux_, grad_p_);
    
    ScalarField dp(grid_->num_cells(), 0.0, FieldLocation::CellCenter, "dp");
    Vector x_vec(dp.raw());
    
    Index iters;
    Real res = solve_linear_system(
        A_p, x_vec, b_p,
        config_.linear_solver_p, config_.preconditioner_p,
        config_.max_linear_iters_p, config_.linear_rtol_p,
        1.0e-50, &iters);
    
    dp = ScalarField(x_vec.raw());
    
    ns_.correct_velocity_and_pressure(
        state_, A_p, diag_inv_Ap[0], HbyA_, dp,
        face_velocity_, face_mass_flux_, grad_p_);
    
    ns_.disc().under_relax_field(state_.p, state_prev_.p, config_.disc_config().relaxation_p);
    
    (void)res;
}

void FluidSolver::solve_turbulence_equations(Real dt) {
    if (!config_.enable_turbulence) return;
    if (current_iter_ < config_.turbulence_start_iter) return;
    
    ns_.assemble_turbulence_equations(
        state_, &state_prev_, dt, grad_U_,
        A_k, b_k, A_eps, b_eps);
    
    ScalarField k_new = state_.turb.k;
    Vector x_vec(k_new.raw());
    Index iters;
    solve_linear_system(A_k, x_vec, b_k,
        config_.linear_solver_turb, config_.preconditioner_turb,
        config_.max_linear_iters_k, config_.linear_rtol_turb,
        1.0e-50, &iters);
    k_new = ScalarField(x_vec.raw());
    
    ScalarField eps_new = state_.turb.epsilon;
    Vector x_vec_eps(eps_new.raw());
    solve_linear_system(A_eps, x_vec_eps, b_eps,
        config_.linear_solver_turb, config_.preconditioner_turb,
        config_.max_linear_iters_eps, config_.linear_rtol_turb,
        1.0e-50, &iters);
    eps_new = ScalarField(x_vec_eps.raw());
    
    ns_.disc().under_relax_field(k_new, state_.turb.k, config_.disc_config().relaxation_k);
    ns_.disc().under_relax_field(eps_new, state_.turb.epsilon, config_.disc_config().relaxation_eps);
    
    ns_.update_turbulence_from_solution(state_, k_new, eps_new);
    
    if (bc_mgr_) bc_mgr_->apply_all_turbulence_bc(*grid_, state_);
}

void FluidSolver::under_relax_momentum() {
}

void FluidSolver::under_relax_pressure() {
}

void FluidSolver::under_relax_turbulence() {
}

FluidIterationStats FluidSolver::perform_simple_iteration(Index iter_count) {
    FluidIterationStats stats;
    current_iter_++;
    
    Timer iter_timer;
    iter_timer.start("iter");
    
    for (Index comp = 0; comp < 3; ++comp) {
        solve_momentum_predictor(comp);
        
        Vector resid(grid_->num_cells());
        A_U[comp].multiply(Vector(state_.U.component(comp).raw()), resid);
        for (Index i = 0; i < resid.size(); ++i) resid[i] -= b_U[comp][i];
        stats.final_res_U[comp] = resid.norm();
        if (iter_count == 0) stats.initial_res_U[comp] = stats.final_res_U[comp];
    }
    
    apply_all_boundary_conditions();
    update_gradients();
    
    solve_pressure_correction();
    
    Vector p_resid(grid_->num_cells());
    A_p.multiply(Vector(state_.p.raw()), p_resid);
    for (Index i = 0; i < p_resid.size(); ++i) p_resid[i] -= b_p[i];
    stats.final_res_p = p_resid.norm();
    stats.initial_res_p = (iter_count == 0) ? stats.final_res_p : stats.initial_res_p;
    
    Real dt = (config_.algorithm == SolverAlgorithm::SteadySIMPLE) ? 1.0 : compute_suggested_dt();
    solve_turbulence_equations(dt);
    
    apply_all_boundary_conditions();
    update_gradients();
    
    state_.update_thermophysical_properties();
    
    stats.mass_imbalance = ns_.compute_mass_imbalance(face_mass_flux_);
    stats.cfl = ns_.compute_cfl_number(state_, face_mass_flux_, dt);
    stats.corrector_iters = config_.max_correctors;
    stats.p_residual_history.push_back(stats.final_res_p);
    stats.u_residual_history.push_back(
        std::max({stats.final_res_U[0], stats.final_res_U[1], stats.final_res_U[2]}));
    
    iter_timer.stop("iter");
    stats.solve_time = iter_timer.elapsed("iter");
    
    return stats;
}

FluidIterationStats FluidSolver::perform_piso_iteration(Real dt, Index n_correctors) {
    FluidIterationStats stats;
    current_iter_++;
    
    Timer iter_timer;
    iter_timer.start("iter");
    
    for (Index comp = 0; comp < 3; ++comp) {
        solve_momentum_predictor(comp);
    }
    apply_all_boundary_conditions();
    
    for (Index corr = 0; corr < n_correctors; ++corr) {
        update_gradients();
        solve_pressure_correction();
        apply_all_boundary_conditions();
        
        for (Index comp = 0; comp < 3; ++comp) {
            solve_momentum_predictor(comp);
        }
        apply_all_boundary_conditions();
    }
    
    solve_turbulence_equations(dt);
    apply_all_boundary_conditions();
    update_gradients();
    state_.update_thermophysical_properties();
    
    stats.mass_imbalance = ns_.compute_mass_imbalance(face_mass_flux_);
    stats.cfl = ns_.compute_cfl_number(state_, face_mass_flux_, dt);
    
    iter_timer.stop("iter");
    stats.solve_time = iter_timer.elapsed("iter");
    stats.corrector_iters = n_correctors;
    
    return stats;
}

bool FluidSolver::solve_steady(FluidSolutionResult& result, Index max_iters) {
    if (!initialized_) initialize();
    
    Timer solve_timer;
    solve_timer.start("solve_steady");
    
    result.converged = false;
    result.iteration_history.reserve(max_iters);
    
    for (Index iter = 0; iter < max_iters; ++iter) {
        FluidIterationStats stats;
        
        if (config_.algorithm == SolverAlgorithm::SIMPLE ||
            config_.algorithm == SolverAlgorithm::SteadySIMPLE) {
            stats = perform_simple_iteration(iter);
        } else {
            Real dt = config_.use_adaptive_dt ? compute_suggested_dt() : 1.0;
            stats = perform_piso_iteration(dt, config_.max_correctors);
        }
        
        result.iteration_history.push_back(stats);
        result.total_linear_iters += stats.total_linear_iters;
        
        if (iter % config_.output_interval == 0 || iter == max_iters - 1) {
            print_iteration_report(stats, iter);
        }
        
        if (check_convergence(stats)) {
            result.converged = true;
            result.reason = "tolerance";
            result.total_iters = iter + 1;
            break;
        }
        
        if (config_.compute_forces_every_step) {
            auto forces = compute_all_boundary_forces();
            result.aerodynamic_forces.clear();
            for (const auto& f : forces) result.aerodynamic_forces.push_back(f.second);
        }
    }
    
    solve_timer.stop("solve_steady");
    result.total_time = solve_timer.elapsed("solve_steady");
    
    if (result.iteration_history.size() > 0) {
        result.final_stats = result.iteration_history.back();
    }
    
    if (!result.converged) {
        result.reason = "max_iters";
        result.total_iters = max_iters;
    }
    
    HVDC_LOG_INFO("Steady fluid solve " << (result.converged ? "converged" : "did NOT converge")
                  << " after " << result.total_iters << " iterations in "
                  << result.total_time << "s");
    
    return result.converged;
}

bool FluidSolver::solve_transient_step(
    Real dt,
    FluidSolutionResult& result,
    const FluidState* state_old)
{
    if (!initialized_) initialize();
    if (state_old) state_old_ = *state_old;
    
    save_state();
    state_prev_ = state_old ? *state_old : state_old_;
    
    if (config_.algorithm == SolverAlgorithm::PIMPLE) {
        for (Index outer = 0; outer < config_.pimple_outer_iters; ++outer) {
            perform_piso_iteration(dt, config_.max_correctors);
        }
    } else {
        perform_piso_iteration(dt, config_.max_correctors);
    }
    
    result.converged = true;
    result.reason = "transient_step";
    result.total_iters = 1;
    
    return true;
}

bool FluidSolver::check_convergence(const FluidIterationStats& stats) {
    Real max_U_res = *std::max_element(stats.final_res_U, stats.final_res_U + 3);
    
    bool abs_conv = (max_U_res < config_.residual_tol) &&
                     (stats.final_res_p < config_.pressure_tol);
    
    Real rel_factor = 1.0;
    for (Index i = 0; i < 3; ++i) {
        if (stats.initial_res_U[i] > EPS) {
            rel_factor = std::min(rel_factor, stats.final_res_U[i] / stats.initial_res_U[i]);
        }
    }
    bool rel_conv = rel_factor < config_.residual_rel_tol;
    
    return abs_conv || rel_conv;
}

void FluidSolver::print_iteration_report(const FluidIterationStats& stats, Index iter) const {
    auto& mpi = MPIManager::instance();
    if (!mpi.is_root()) return;
    
    auto flow_stats = ns_.compute_flow_stats(state_);
    
    std::cout << "  Iter " << std::setw(5) << iter
              << ": Ux=" << std::scientific << std::setw(8) << std::setprecision(2) << stats.final_res_U[0]
              << " Uy=" << std::scientific << std::setw(8) << std::setprecision(2) << stats.final_res_U[1]
              << " Uz=" << std::scientific << std::setw(8) << std::setprecision(2) << stats.final_res_U[2]
              << " p=" << std::scientific << std::setw(8) << std::setprecision(2) << stats.final_res_p
              << " CFL=" << std::fixed << std::setw(5) << std::setprecision(2) << stats.cfl
              << " Um=" << std::fixed << std::setw(6) << std::setprecision(2) << flow_stats[1]
              << " |dm|=" << std::scientific << std::setw(8) << std::setprecision(2) << stats.mass_imbalance
              << std::endl;
}

Vec3 FluidSolver::compute_aerodynamic_force_on_patch(const std::string& patch_name) {
    if (!bc_mgr_) return {0, 0, 0};
    BoundaryCondition* bc = bc_mgr_->bc_by_name(patch_name);
    if (!bc) return {0, 0, 0};
    return ns_.compute_total_force_on_boundary(bc->patch_id(), state_, grad_p_, grad_U_);
}

std::vector<std::pair<std::string, Vec3>> FluidSolver::compute_all_boundary_forces() {
    std::vector<std::pair<std::string, Vec3>> forces;
    if (!bc_mgr_) return forces;
    
    for (Index i = 0; i < bc_mgr_->num_bcs(); ++i) {
        const auto* bc = bc_mgr_->bc(i);
        Vec3 f = ns_.compute_total_force_on_boundary(bc->patch_id(), state_, grad_p_, grad_U_);
        forces.emplace_back(bc->name(), f);
    }
    return forces;
}

std::vector<Real> FluidSolver::get_face_pressure_on_patch(const std::string& patch_name) const {
    std::vector<Real> pressures;
    if (!bc_mgr_) return pressures;
    const BoundaryCondition* bc = bc_mgr_->bc_by_name(patch_name);
    if (!bc) return pressures;
    
    auto range = grid_->boundary_patch_range(bc->patch_id());
    if (range.first < 0) return pressures;
    
    pressures.reserve(range.second - range.first);
    for (Index fi = range.first; fi < range.second; ++fi) {
        if (fi < static_cast<Index>(face_pressure_.size())) {
            pressures.push_back(face_pressure_[fi]);
        }
    }
    return pressures;
}

std::vector<Vec3> FluidSolver::get_face_velocity_on_patch(const std::string& patch_name) const {
    std::vector<Vec3> velocities;
    if (!bc_mgr_) return velocities;
    const BoundaryCondition* bc = bc_mgr_->bc_by_name(patch_name);
    if (!bc) return velocities;
    
    auto range = grid_->boundary_patch_range(bc->patch_id());
    if (range.first < 0) return velocities;
    
    velocities.reserve(range.second - range.first);
    for (Index fi = range.first; fi < range.second; ++fi) {
        if (fi < static_cast<Index>(face_velocity_.size())) {
            velocities.push_back(face_velocity_[fi]);
        }
    }
    return velocities;
}

void FluidSolver::update_interface_velocity(const std::vector<Index>& face_ids,
                                              const std::vector<Vec3>& velocities) {
    for (Index i = 0; i < static_cast<Index>(face_ids.size()); ++i) {
        Index fi = face_ids[i];
        if (fi < static_cast<Index>(face_velocity_.size()) && 
            i < static_cast<Index>(velocities.size())) {
            face_velocity_[fi] = velocities[i];
            const Face& face = grid_->face(fi);
            Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
            if (ci >= 0 && ci < state_.U.size()) {
                state_.U[ci] = velocities[i];
            }
        }
    }
}

Real FluidSolver::compute_suggested_dt() const {
    Real cfl = ns_.compute_cfl_number(state_, face_mass_flux_, 1.0);
    if (cfl < EPS) return config_.adaptive_dt_max;
    
    Real dt = config_.adaptive_dt_cfl / cfl;
    return std::clamp(dt, config_.adaptive_dt_min, config_.adaptive_dt_max);
}

void FluidSolver::set_free_stream_velocity(const Vec3& U_inf) {
    if (!grid_) return;
    U_inf_ = U_inf;
    if (!initialized_) {
        state_.U.resize(grid_->num_cells());
        state_.p.resize(grid_->num_cells());
        for (Index ci = 0; ci < grid_->num_cells(); ++ci) {
            state_.U[ci] = U_inf;
            state_.p[ci] = 0.0;
        }
    }
}

void FluidSolver::initialize_fields() {
    if (!grid_) return;
    Index nc = grid_->num_cells();
    Index nf = grid_->num_faces();
    
    state_.U.assign(nc, U_inf_);
    state_.p.assign(nc, 0.0);
    state_.p_ref = 101325.0;
    state_.rho.assign(nc, 1.225);
    state_.mu.assign(nc, 1.789e-5);
    state_.nu.assign(nc, state_.mu[0] / state_.rho[0]);
    state_.turb.k.assign(nc, 0.01);
    state_.turb.epsilon.assign(nc, 0.1);
    state_.T.assign(nc, 293.15);
    
    face_velocity_.assign(nf, Vec3{U_inf_[0], U_inf_[1], U_inf_[2]});
    face_mass_flux_.assign(nf, 0.0);
    face_pressure_.assign(nf, 0.0);
    grad_U_.assign(nc, std::array<Vec3,3>{Vec3{0,0,0}, Vec3{0,0,0}, Vec3{0,0,0}});
    grad_p_.assign(nc, Vec3{0,0,0});
    
    initialized_ = true;
}

bool FluidSolver::solve_steady(Index max_iters, Real tol) {
    FluidSolutionResult result;
    return solve_steady(result, max_iters);
    (void)tol;
}

bool FluidSolver::solve_transient_step(Real dt, Index max_iters, Real tol) {
    FluidSolutionResult result;
    return solve_transient_step(dt, result, nullptr);
    (void)max_iters;
    (void)tol;
}

} // namespace fvm
} // namespace hvdc
