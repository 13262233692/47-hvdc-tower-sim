#include "fsi/CouplingManager.hpp"
#include "common/MathUtils.hpp"
#include "common/Logger.hpp"
#include "common/PETScManager.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace hvdc {
namespace fsi {

CouplingManager::CouplingManager(
    fem::FEModel* struct_model,
    fem::Assembler* assembler,
    fem::Solver* struct_solver,
    fvm::FluidSolver* fluid_solver,
    InterfaceMapper* mapper,
    AerodynamicLoads* aero,
    LoadTransfer* transfer,
    const SimulationConfig& cfg)
    : struct_model_(struct_model), assembler_(assembler),
      struct_solver_(struct_solver), fluid_solver_(fluid_solver),
      mapper_(mapper), aero_(aero), transfer_(transfer),
      config_(cfg), dt_(cfg.dt),
      max_iter_(cfg.max_coupling_iters), tol_(cfg.coupling_tol)
{
    timer_.register_section("Total");
    timer_.register_section("Structural");
    timer_.register_section("Fluid");
    timer_.register_section("Coupling");
    timer_.register_section("Transfer");
}

void CouplingManager::do_one_sided_coupling(
    Vector& struct_disp, Vector& struct_vel, Vector& struct_acc,
    fvm::FluidState& fluid_state,
    Real time,
    Vector& interface_forces,
    Vector& interface_displacements,
    bool compute_geometry_updates)
{
    timer_.start("Transfer");
    if (compute_geometry_updates && fluid_solver_ && mapper_) {
        update_fluid_mesh_boundary_motion(struct_disp);
    }
    timer_.stop("Transfer");
    
    timer_.start("Fluid");
    if (fluid_solver_) {
        if (config_.analysis == AnalysisType::StaticLinear ||
            config_.analysis == AnalysisType::StaticNonlinear) {
            fluid_solver_->solve_steady(1000, 1.0e-6);
        } else {
            fluid_solver_->solve_transient_step(dt_, 50, 1.0e-5);
        }
        fluid_state = fluid_solver_->current_state();
    }
    timer_.stop("Fluid");
    
    timer_.start("Transfer");
    Vector aero_forces(struct_model_ ? struct_model_->total_dofs() : 6);
    if (struct_model_) aero_forces.zero();
    
    if (aero_ && struct_model_) {
        std::vector<ElementAeroLoad> elem_loads;
        aero_->compute_element_loads(struct_disp, struct_vel, elem_loads, time);
        if (transfer_) {
            transfer_->transfer_element_aero_to_nodal(elem_loads, aero_forces);
        }
    }
    timer_.stop("Transfer");
    
    timer_.start("Structural");
    reassemble_structure_with_aero_loads(struct_disp, struct_vel, struct_acc,
                                          aero_forces, time);
    timer_.stop("Structural");
    
    interface_forces = aero_forces;
    compute_updated_interface_displacements(struct_disp, interface_displacements);
}

void CouplingManager::compute_updated_interface_displacements(
    const Vector& struct_disp,
    Vector& interface_disp_vec) const
{
    if (!mapper_ || !struct_model_) return;
    
    Index n_iface = mapper_->num_struct_interface_nodes();
    interface_disp_vec.assign(n_iface * 3, 0.0);
    
    Index idx = 0;
    for (Index nid : mapper_->struct_interface_nodes()) {
        const auto* node = struct_model_->get_node(nid);
        if (!node) { idx += 3; continue; }
        Index start = node->dof_start();
        for (Index d = 0; d < 3; ++d) {
            Index gdof = start + d;
            if (gdof < struct_disp.size()) {
                interface_disp_vec[idx + d] = struct_disp[gdof];
            }
        }
        idx += 3;
    }
}

void CouplingManager::reassemble_structure_with_aero_loads(
    Vector& struct_disp, Vector& struct_vel, Vector& struct_acc,
    const Vector& aero_forces,
    Real time)
{
    (void)struct_vel; (void)struct_acc; (void)time;
    if (!struct_solver_ || !assembler_ || !struct_model_) return;
    
    SparseMatrix K;
    Vector F_ext(struct_model_->total_dofs());
    assembler_->assemble_tangent_stiffness(struct_disp, K, struct_model_);
    
    F_ext.zero();
    Index min_sz = std::min(aero_forces.size(), F_ext.size());
    for (Index i = 0; i < min_sz; ++i) {
        F_ext.add(i, aero_forces[i]);
    }
    
    struct_model_->apply_body_force(F_ext);
    struct_model_->apply_ice_load(F_ext);
    
    assembler_->apply_dirichlet_bc(K, F_ext, struct_disp, *struct_model_);
    
    bool is_dynamic = (config_.analysis_type == AnalysisType::TransientLinear ||
                       config_.analysis_type == AnalysisType::TransientNonlinear ||
                       config_.analysis == AnalysisType::TransientLinear ||
                       config_.analysis == AnalysisType::TransientNonlinear);
    if (is_dynamic) {
        SparseMatrix M, C;
        assembler_->assemble_mass_matrix(M, struct_model_);
        assembler_->assemble_damping_matrix(C, struct_model_, 0.05, 0.01);
        
        Vector F_effective = F_ext;
        struct_solver_->apply_newmark_time_integration(
            K, M, C, F_effective, struct_disp, struct_vel, struct_acc,
            dt_, config_.beta, config_.gamma);
    } else {
        Vector du(struct_model_->total_dofs());
        du.zero();
        Vector residual = F_ext;
        assembler_->multiply_stiffness_vector(K, struct_disp, residual, -1.0, 1.0);
        
        const Vector& rhs = residual;
        struct_solver_->solve_linear_system(K, du, rhs);
        for (Index i = 0; i < struct_disp.size(); ++i) {
            struct_disp.add(i, du[i]);
        }
    }
}

void CouplingManager::update_fluid_mesh_boundary_motion(const Vector& struct_disp) {
    if (!fluid_solver_ || !mapper_ || !struct_model_) return;
    
    std::vector<Vec3> face_disp_vec;
    mapper_->map_struct_displacement_to_fluid(struct_disp, face_disp_vec);
    
    if (!face_disp_vec.empty()) {
        auto* fvm_grid = fluid_solver_->grid();
        std::vector<Vec3> velocities;
        if (dt_ > EPS) {
            velocities.reserve(face_disp_vec.size());
            for (const auto& v : face_disp_vec) {
                velocities.push_back(Vec3{v[0]/dt_, v[1]/dt_, v[2]/dt_});
            }
        }
        const auto& face_ids = mapper_->fluid_interface_faces();
        for (Index k = 0; k < static_cast<Index>(face_ids.size()); ++k) {
            Index fi = face_ids[k];
            if (fi >= fvm_grid->num_faces()) continue;
            auto& bc = fvm_grid->mutable_face_bc(fi);
            if (bc && bc->type() == fvm::BCType::InterfaceFSI) {
                auto* fsi_bc = dynamic_cast<fvm::BCInterfaceFSI*>(bc.get());
                if (fsi_bc && k < static_cast<Index>(velocities.size())) {
                    std::vector<Vec3> single_vel(1, velocities[k]);
                    fsi_bc->set_interface_velocity(single_vel);
                    std::vector<Vec3> single_disp(1, face_disp_vec[k]);
                    fsi_bc->set_interface_displacement(single_disp);
                }
            }
        }
    }
}

CouplingConvergenceData CouplingManager::perform_fsi_step(
    Vector& struct_disp,
    Vector& struct_vel,
    Vector& struct_acc,
    fvm::FluidState& fluid_state,
    Real time)
{
    switch (method_) {
        case CouplingIterationMethod::Explicit:
            return perform_explicit_step(struct_disp, struct_vel, struct_acc, fluid_state, time);
        case CouplingIterationMethod::Implicit_Aitken:
            return perform_implicit_aitken_step(struct_disp, struct_vel, struct_acc, fluid_state, time);
        case CouplingIterationMethod::Implicit_Anderson:
            return perform_implicit_anderson_step(struct_disp, struct_vel, struct_acc, fluid_state, time);
        case CouplingIterationMethod::Implicit_IQN_ILS:
            return perform_implicit_iqn_ils_step(struct_disp, struct_vel, struct_acc, fluid_state, time);
        case CouplingIterationMethod::Implicit_Broyden:
            return perform_implicit_broyden_step(struct_disp, struct_vel, struct_acc, fluid_state, time);
        default:
            return perform_explicit_step(struct_disp, struct_vel, struct_acc, fluid_state, time);
    }
}

CouplingConvergenceData CouplingManager::perform_explicit_step(
    Vector& struct_disp,
    Vector& struct_vel,
    Vector& struct_acc,
    fvm::FluidState& fluid_state,
    Real time)
{
    CouplingConvergenceData data;
    Vector interface_forces, interface_disps;
    
    timer_.start("Coupling");
    do_one_sided_coupling(struct_disp, struct_vel, struct_acc,
                          fluid_state, time,
                          interface_forces, interface_disps, true);
    timer_.stop("Coupling");
    
    data.iterations = 1;
    data.converged = true;
    data.force_norm = interface_forces.norm();
    data.disp_norm = interface_disps.norm();
    return data;
}

CouplingConvergenceData CouplingManager::perform_implicit_aitken_step(
    Vector& struct_disp,
    Vector& struct_vel,
    Vector& struct_acc,
    fvm::FluidState& fluid_state,
    Real time)
{
    CouplingConvergenceData data;
    Vector x_old = struct_disp;
    Vector interface_forces, interface_disps;
    Real omega = aitken_omega_;
    std::vector<Vector> x_hist;
    x_hist.reserve(max_iter_ + 1);
    
    for (Index iter = 0; iter < max_iter_; ++iter) {
        data.iterations = iter + 1;
        Vector x_iter = struct_disp;
        
        timer_.start("Coupling");
        do_one_sided_coupling(struct_disp, struct_vel, struct_acc,
                              fluid_state, time,
                              interface_forces, interface_disps, iter > 0);
        timer_.stop("Coupling");
        
        compute_interface_residual(x_old, struct_disp, x_iter, interface_forces, data);
        
        if (check_convergence(data)) {
            data.converged = true;
            aitken_omega_ = omega;
            return data;
        }
        
        x_hist.push_back(struct_disp);
        apply_aitken_relaxation(x_old, struct_disp, struct_disp, omega, x_hist);
        
        x_old = x_old;
        data.disp_norm_prev = data.disp_norm;
        data.force_norm_prev = data.force_norm;
    }
    
    HVDC_LOG_WARNING("Aitken coupling did not converge in " << max_iter_
                   << " iterations, disp_res=" << data.disp_norm);
    data.converged = false;
    return data;
}

CouplingConvergenceData CouplingManager::perform_implicit_anderson_step(
    Vector& struct_disp,
    Vector& struct_vel,
    Vector& struct_acc,
    fvm::FluidState& fluid_state,
    Real time)
{
    CouplingConvergenceData data;
    const Index m = 5;
    
    Vector disp_old = struct_disp;
    Vector force_old(struct_model_ ? struct_model_->total_dofs() : 6);
    if (struct_model_) force_old.zero();
    
    std::vector<Vector> x_hist;
    std::vector<Vector> r_hist;
    
    for (Index iter = 0; iter < max_iter_; ++iter) {
        data.iterations = iter + 1;
        Vector interface_forces, interface_disps;
        Vector x_iter = struct_disp;
        
        timer_.start("Coupling");
        do_one_sided_coupling(struct_disp, struct_vel, struct_acc,
                              fluid_state, time,
                              interface_forces, interface_disps, iter > 0);
        timer_.stop("Coupling");
        
        compute_interface_residual(disp_old, struct_disp, force_old, interface_forces, data);
        
        if (check_convergence(data)) {
            data.converged = true;
            return data;
        }
        
        Vector r = struct_disp;
        r -= x_iter;
        x_hist.push_back(struct_disp);
        r_hist.push_back(r);
        
        if (x_hist.size() > m + 1) {
            x_hist.erase(x_hist.begin());
            r_hist.erase(r_hist.begin());
        }
        
        if (iter >= 1) {
            apply_anderson_mixing(x_hist, r_hist, struct_disp, m);
        }
        
        disp_old = struct_disp;
        force_old = interface_forces;
        data.disp_norm_prev = data.disp_norm;
        data.force_norm_prev = data.force_norm;
    }
    
    HVDC_LOG_WARNING("Anderson coupling did not converge in " << max_iter_
                   << " iterations");
    data.converged = false;
    return data;
}

CouplingConvergenceData CouplingManager::perform_implicit_iqn_ils_step(
    Vector& struct_disp,
    Vector& struct_vel,
    Vector& struct_acc,
    fvm::FluidState& fluid_state,
    Real time)
{
    CouplingConvergenceData data;
    
    std::vector<Vector> disp_hist;
    std::vector<Vector> force_hist;
    
    Vector disp_old = struct_disp;
    Vector interface_forces(struct_model_ ? struct_model_->total_dofs() : 6);
    if (struct_model_) interface_forces.zero();
    
    for (Index iter = 0; iter < max_iter_; ++iter) {
        data.iterations = iter + 1;
        Vector interface_disps;
        Vector x_trial = struct_disp;
        
        timer_.start("Coupling");
        do_one_sided_coupling(struct_disp, struct_vel, struct_acc,
                              fluid_state, time,
                              interface_forces, interface_disps, iter > 0);
        timer_.stop("Coupling");
        
        compute_interface_residual(disp_old, struct_disp, disp_hist.size() > 0 ? force_hist.back() : interface_forces, interface_forces, data);
        
        if (check_convergence(data)) {
            data.converged = true;
            return data;
        }
        
        disp_hist.push_back(struct_disp);
        force_hist.push_back(interface_forces);
        
        if (iter >= 1) {
            apply_iqn_ils(disp_hist, force_hist, struct_disp);
        }
        
        disp_old = disp_hist[disp_hist.size() - 1];
        data.disp_norm_prev = data.disp_norm;
        data.force_norm_prev = data.force_norm;
    }
    
    HVDC_LOG_WARNING("IQN-ILS coupling did not converge in " << max_iter_
                   << " iterations");
    data.converged = false;
    return data;
}

CouplingConvergenceData CouplingManager::perform_implicit_broyden_step(
    Vector& struct_disp,
    Vector& struct_vel,
    Vector& struct_acc,
    fvm::FluidState& fluid_state,
    Real time)
{
    CouplingConvergenceData data;
    
    Index n = struct_disp.size();
    SparseMatrix J_inv;
    J_inv.reset(n, n);
    for (Index i = 0; i < n; ++i) {
        J_inv.set(i, i, 1.0);
    }
    
    Vector prev_x = struct_disp;
    Vector prev_F(n);
    prev_F.zero();
    
    for (Index iter = 0; iter < max_iter_; ++iter) {
        data.iterations = iter + 1;
        Vector interface_forces, interface_disps;
        
        timer_.start("Coupling");
        do_one_sided_coupling(struct_disp, struct_vel, struct_acc,
                              fluid_state, time,
                              interface_forces, interface_disps, iter > 0);
        timer_.stop("Coupling");
        
        compute_interface_residual(prev_x, struct_disp, prev_F, interface_forces, data);
        
        if (check_convergence(data)) {
            data.converged = true;
            return data;
        }
        
        Vector delta_x = struct_disp;
        delta_x -= prev_x;
        Vector delta_F = interface_forces;
        delta_F -= prev_F;
        
        update_broyden_jacobian(J_inv, delta_x, delta_F, 1.0);
        
        Vector correction(n);
        correction.zero();
        for (Index i = 0; i < n; ++i) {
            correction.add(i, -0.1 * interface_forces[i]);
        }
        J_inv.multiply(correction, struct_disp);
        for (Index i = 0; i < n; ++i) {
            struct_disp.add(i, prev_x[i]);
        }
        
        prev_x = struct_disp;
        prev_F = interface_forces;
        data.disp_norm_prev = data.disp_norm;
        data.force_norm_prev = data.force_norm;
    }
    
    HVDC_LOG_WARNING("Broyden coupling did not converge in " << max_iter_
                   << " iterations");
    data.converged = false;
    return data;
}

void CouplingManager::compute_interface_residual(
    const Vector& disp_old, const Vector& disp_new,
    const Vector& force_old, const Vector& force_new,
    CouplingConvergenceData& data) const
{
    data.disp_norm = compute_normalized_residual(disp_new, disp_old);
    data.force_norm = compute_normalized_residual(force_new, force_old);
    data.residual_norm = std::max(data.disp_norm, data.force_norm);
    
    Real tx = 0, ty = 0, tz = 0;
    Index n = force_new.size() / 3;
    for (Index i = 0; i < n; ++i) {
        if (i*3 < force_new.size()) tx += force_new[i*3];
        if (i*3+1 < force_new.size()) ty += force_new[i*3+1];
        if (i*3+2 < force_new.size()) tz += force_new[i*3+2];
    }
    data.total_force_x = tx;
    data.total_force_y = ty;
    data.total_force_z = tz;
}

bool CouplingManager::check_convergence(const CouplingConvergenceData& data) const {
    switch (conv_criterion_) {
        case FSIConvergenceCriterion::Displacement:
            return data.disp_norm < tol_;
        case FSIConvergenceCriterion::Force:
            return data.force_norm < tol_;
        case FSIConvergenceCriterion::DisplacementAndForce:
            return data.disp_norm < tol_ && data.force_norm < tol_;
        case FSIConvergenceCriterion::ResidualNorm:
            return data.residual_norm < tol_;
        default:
            return data.disp_norm < tol_;
    }
}

void CouplingManager::apply_aitken_relaxation(
    const Vector& x_old, const Vector& x_new,
    Vector& x_relaxed, Real& omega,
    const std::vector<Vector>& x_history) const
{
    (void)x_history;
    Index n = std::min({x_old.size(), x_new.size(), x_relaxed.size()});
    
    Vector dx_new = x_new;
    dx_new -= x_old;
    
    Real denom = dx_new.dot(dx_new);
    if (denom < EPS) {
        omega = std::min(1.0, 2.0 * omega);
    } else {
        if (x_history.size() >= 2) {
            Vector dx_prev = x_history[x_history.size() - 1];
            dx_prev -= x_history[x_history.size() - 2];
            Vector diff = dx_new;
            diff -= dx_prev;
            omega = omega * (1.0 - dx_new.dot(diff) / denom);
        }
        omega = std::max(0.001, std::min(2.0, omega));
    }
    
    for (Index i = 0; i < n; ++i) {
        x_relaxed[i] = x_old[i] + omega * (x_new[i] - x_old[i]);
    }
}

void CouplingManager::apply_anderson_mixing(
    const std::vector<Vector>& x_history,
    const std::vector<Vector>& r_history,
    Vector& x_new,
    Index m) const
{
    if (x_history.size() < 2) return;
    
    Index n = x_new.size();
    Index k = static_cast<Index>(x_history.size()) - 1;
    Index mk = std::min(m, k);
    
    if (mk < 1) return;
    
    Index start = k - mk;
    std::vector<std::vector<Real>> mat(mk, std::vector<Real>(mk, 0.0));
    std::vector<Real> rhs(mk, 0.0);
    
    for (Index i = 0; i < mk; ++i) {
        Index gi = start + i;
        if (gi >= r_history.size()) break;
        for (Index j = 0; j < mk; ++j) {
            Index gj = start + j;
            if (gj >= r_history.size()) break;
            Vector ri = r_history[gi + 1];
            ri -= r_history[start];
            Vector rj = r_history[gj + 1];
            rj -= r_history[start];
            mat[i][j] = ri.dot(rj);
        }
        Vector ri = r_history[gi + 1];
        ri -= r_history[start];
        Vector rk = r_history[k];
        rk -= r_history[start];
        rhs[i] = ri.dot(r_history[k]);
    }
    
    std::vector<Real> alpha(mk, 0.0);
    alpha[mk - 1] = 1.0;
    
    x_new = x_history.back();
    for (Index i = 0; i < mk; ++i) {
        Index gi = start + i + 1;
        if (gi >= x_history.size()) continue;
        Real a = alpha[i] / std::max<Real>(mk, 1.0);
        for (Index j = 0; j < n; ++j) {
            x_new[j] += a * (x_history[gi][j] - x_history[start][j]);
        }
    }
}

void CouplingManager::apply_iqn_ils(
    const std::vector<Vector>& displacement_history,
    const std::vector<Vector>& force_history,
    Vector& next_displacement) const
{
    if (displacement_history.size() < 3 || force_history.size() < 3) return;
    
    Index n = next_displacement.size();
    Index m = static_cast<Index>(displacement_history.size()) - 1;
    
    std::vector<Vector> delta_V;
    std::vector<Vector> delta_W;
    delta_V.reserve(m);
    delta_W.reserve(m);
    
    for (Index k = 1; k < static_cast<Index>(displacement_history.size()); ++k) {
        Vector dV = displacement_history[k];
        dV -= displacement_history[k-1];
        Vector dW = force_history[k];
        dW -= force_history[k-1];
        delta_V.push_back(std::move(dV));
        delta_W.push_back(std::move(dW));
    }
    
    next_displacement = displacement_history.back();
}

void CouplingManager::update_broyden_jacobian(
    SparseMatrix& J_approx,
    const Vector& delta_x, const Vector& delta_F,
    Real gamma) const
{
    (void)J_approx; (void)delta_x; (void)delta_F; (void)gamma;
}

Real CouplingManager::compute_normalized_residual(
    const Vector& current, const Vector& previous) const
{
    if (current.size() == 0) return 0.0;
    
    Index n = std::min(current.size(), previous.size());
    if (n == 0) return 0.0;
    
    Real num = 0.0;
    Real den = 0.0;
    for (Index i = 0; i < n; ++i) {
        Real diff = current[i] - previous[i];
        num += diff * diff;
        den += current[i] * current[i] + 1.0;
    }
    return std::sqrt(num / den);
}

void CouplingManager::run_coupled_simulation(
    Vector& struct_displacement,
    Vector& struct_velocity,
    Vector& struct_acceleration,
    fvm::FluidState& fluid_state,
    std::function<void(Index, Real, const Vector&, const fvm::FluidState&)> callback)
{
    Real t = config_.t_start;
    Real t_end = config_.t_end;
    Real dt = dt_;
    Index step = 0;
    
    HVDC_LOG_INFO("Starting coupled FSI simulation: t=[" << t << ", " << t_end 
                 << "], dt=" << dt << ", method=" << static_cast<Index>(method_));
    
    history_.displacement_history.clear();
    history_.force_history.clear();
    history_.time_history.clear();
    
    while (t < t_end - EPS) {
        if (t + dt > t_end) dt = t_end - t;
        dt_ = dt;
        
        HVDC_LOG_INFO("FSI step " << step << ": t=" << t << ", dt=" << dt);
        
        CouplingConvergenceData data = perform_fsi_step(
            struct_displacement, struct_velocity, struct_acceleration,
            fluid_state, t);
        
        history_.displacement_history.push_back(struct_displacement);
        history_.time_history.push_back(t + dt);
        
        HVDC_LOG_INFO("  Step " << step << " completed: iters=" << data.iterations
                     << ", converged=" << data.converged
                     << ", disp_res=" << data.disp_norm
                     << ", force_res=" << data.force_norm);
        
        if (callback) {
            callback(step, t + dt, struct_displacement, fluid_state);
        }
        
        if (data.converged) {
            Real factor = 1.0;
            if (data.iterations <= max_iter_ / 5) factor = 1.3;
            else if (data.iterations <= max_iter_ / 3) factor = 1.1;
            else if (data.iterations >= max_iter_ * 4 / 5) factor = 0.7;
            else if (data.iterations >= max_iter_ * 2 / 3) factor = 0.85;
            dt = std::min(config_.dt_max, std::max(config_.dt_min, dt * factor));
        } else {
            dt = std::max(config_.dt_min, dt * 0.5);
        }
        
        t += dt_;
        step++;
    }
    
    HVDC_LOG_INFO("Coupled FSI simulation completed: total_steps=" << step);
}

void CouplingManager::compute_coupling_diagnostics(
    const Vector& struct_disp,
    const fvm::FluidState& fluid_state,
    Real time) const
{
    (void)struct_disp; (void)fluid_state; (void)time;
}

Real CouplingManager::estimate_interface_added_mass(
    const Vector& trial_disp,
    const fvm::FluidState& fluid_state,
    Real dt) const
{
    (void)trial_disp; (void)fluid_state; (void)dt;
    return 0.0;
}

} // namespace fsi
} // namespace hvdc
