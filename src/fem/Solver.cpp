#include "fem/Solver.hpp"
#include "common/MathUtils.hpp"
#include "common/MPIManager.hpp"
#include "common/Timer.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>

extern hvdc::Timer g_timer;

namespace hvdc {
namespace fem {

Solver::Solver(FEModel* model, Assembler* assembler)
    : model_(model), assembler_(assembler)
{
    if (!assembler_ && model_) {
        owned_assembler_ = std::make_unique<Assembler>(model_);
        assembler_ = owned_assembler_.get();
    }
}

void Solver::initialize() {
    if (!model_) return;
    if (!assembler_) {
        owned_assembler_ = std::make_unique<Assembler>(model_);
        assembler_ = owned_assembler_.get();
    }
    assembler_->initialize();
    initialized_ = true;
    
#ifdef HVDC_USE_PETSC
    petsc_solver_.create();
    petsc_solver_.set_ksp_type(config_.ksp_type);
    petsc_solver_.set_preconditioner(config_.pc_type);
    petsc_solver_.set_tolerances(config_.krylov_rtol, config_.krylov_atol,
                                  1.0e5, config_.krylov_max_iters);
#endif
    
    HVDC_LOG_INFO("FEM Solver initialized: type=" << config_.ksp_type
                  << ", pc=" << config_.pc_type);
}

bool Solver::solve_linear_system_impl(
    const SparseMatrix& K,
    const Vector& rhs,
    Vector& solution)
{
#ifdef HVDC_USE_PETSC
    PETScMatrix A_petsc;
    PETScVector b_petsc(static_cast<Index>(rhs.size()));
    PETScVector x_petsc(static_cast<Index>(solution.size()));
    
    A_petsc.from_sparse_matrix(K);
    b_petsc.from_vector(rhs);
    x_petsc.zero();
    
    petsc_solver_.set_operator(A_petsc);
    bool ok = petsc_solver_.solve(b_petsc, x_petsc);
    x_petsc.to_vector(solution);
    
    return ok;
#else
    Index n = rhs.size();
    solution.resize(n);
    solution.zero();
    
    const auto& rp = K.row_ptr();
    const auto& ci = K.col_idx();
    const auto& v = K.values();
    
    std::vector<Real> diag(n);
    for (Index i = 0; i < n; ++i) {
        for (Index j = rp[i]; j < rp[i+1]; ++j) {
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
        if (iter > 0 && r_norm / config_.krylov_atol < config_.krylov_rtol) return true;
        
        Vector p(n);
        for (Index i = 0; i < n; ++i) {
            p[i] = (std::fabs(diag[i]) > EPS) ? r[i] / diag[i] : r[i];
        }
        
        Vector Kp(n);
        K.multiply(p, Kp);
        
        Real pKp = p.dot(Kp);
        if (std::fabs(pKp) < EPS) break;
        
        Real alpha = p.dot(r) / pKp;
        for (Index i = 0; i < n; ++i) {
            solution[i] += alpha * p[i];
            r[i] -= alpha * Kp[i];
        }
    }
    
    return true;
#endif
}

bool Solver::solve_static_nonlinear(
    const Vector& F_ext,
    const Vector& u_guess,
    SolverResult& result)
{
    if (!initialized_) initialize();
    
    Timer solver_timer;
    solver_timer.start("solve");
    
    Index n = model_->total_dofs();
    result.displacement = u_guess;
    if (result.displacement.size() != n) {
        result.displacement.resize(n);
        model_->gather_displacement_vec(result.displacement.raw());
    }
    result.velocity.resize(n);
    result.velocity.zero();
    result.acceleration.resize(n);
    result.acceleration.zero();
    result.reaction_forces.resize(n);
    result.reaction_forces.zero();
    result.residual_history.clear();
    result.linear_iter_history.clear();
    
    Vector F_zero_extra(n);
    F_zero_extra.zero();
    
    bool ok = newton_raphson_iteration(F_ext, F_zero_extra, result);
    
    solver_timer.stop("solve");
    result.solve_time = solver_timer.elapsed("solve");
    
    if (ok) {
        result.reaction_forces = compute_reactions(result.displacement);
    }
    
    return ok;
}

bool Solver::solve_static_linear(
    const Vector& F_ext,
    SolverResult& result)
{
    if (!initialized_) initialize();
    
    Timer timer;
    timer.start("solve_linear");
    
    Index n = model_->total_dofs();
    result.displacement.resize(n);
    result.displacement.zero();
    result.velocity.resize(n);
    result.velocity.zero();
    result.acceleration.resize(n);
    result.acceleration.zero();
    result.reaction_forces.resize(n);
    result.reaction_forces.zero();
    
    model_->gather_displacement_vec(result.displacement.raw());
    
    AssemblyResult assembly;
    assembler_->assemble_tangent_stiffness(result.displacement, assembly, false, false, true);
    assembler_->apply_dirichlet_bc(assembly.K_T, result.reaction_forces, result.displacement);
    
    Vector rhs = F_ext;
    for (Index gdof : model_->constrained_dofs()) {
        rhs.set(gdof, result.reaction_forces[gdof]);
    }
    
    bool ok = solve_linear_system_impl(assembly.K_T, rhs, result.displacement);
    
    timer.stop("solve_linear");
    result.solve_time = timer.elapsed("solve_linear");
    result.converged = ok;
    result.newton_iters = 1;
    
    if (ok) {
        result.reaction_forces = compute_reactions(result.displacement);
    }
    
    return ok;
}

bool Solver::newton_raphson_iteration(
    const Vector& F_ext,
    const Vector& F_extra,
    SolverResult& result,
    bool is_transient,
    Real dt,
    const Vector* u_prev,
    const Vector* v_prev,
    const Vector* a_prev,
    Real beta,
    Real gamma)
{
    Index n = model_->total_dofs();
    
    AssemblyResult assembly;
    assembler_->assemble_tangent_stiffness(result.displacement, assembly, true, is_transient, true);
    
    Vector F_total = F_ext;
    if (F_extra.size() > 0) F_total += F_extra;
    
    Vector residual(n);
    for (Index i = 0; i < n; ++i) {
        residual[i] = F_total[i] - assembly.F_int[i];
    }
    
    Vector M_times_a(n);
    if (is_transient && u_prev && v_prev && a_prev && !assembly.M.values().empty()) {
        Vector a_pred = *a_prev;
        assembly.M.multiply(a_pred, M_times_a);
        for (Index i = 0; i < n; ++i) {
            residual[i] -= M_times_a[i];
        }
        
        SparseMatrix Keff = assembly.K_T;
        for (Index i = 0; i < n; ++i) {
            for (Index j = assembly.K_T.row_ptr()[i]; j < assembly.K_T.row_ptr()[i+1]; ++j) {
                Index c = assembly.K_T.col_idx()[j];
                Real m_val = assembly.M(i, c);
                Keff.add(i, c, m_val / (beta * dt * dt));
            }
        }
        (void)gamma;
        assembly.K_T = Keff;
    }
    
    for (Index gdof : model_->constrained_dofs()) {
        residual.set(gdof, 0.0);
    }
    
    result.initial_residual_norm = residual.norm();
    result.residual_norm = result.initial_residual_norm;
    
    if (result.initial_residual_norm < config_.tol_abs) {
        result.converged = true;
        result.converged_reason = "initial_residual";
        HVDC_LOG_INFO("Newton converged initially: ||R||=" << result.initial_residual_norm);
        return true;
    }
    
    result.converged = false;
    Index iter;
    
    for (iter = 0; iter < config_.max_iters; ++iter) {
        result.residual_history.push_back(result.residual_norm);
        
        assembler_->apply_dirichlet_bc(assembly.K_T, residual, result.displacement);
        
        Vector du(n);
        du.zero();
        bool solve_ok = solve_linear_system_impl(assembly.K_T, residual, du);
        
        if (!solve_ok) {
            HVDC_LOG_WARNING("Linear solver failed at NR iter " << iter);
        }
        
        Real alpha = line_search(du, residual, F_total, assembly);
        
        for (Index i = 0; i < n; ++i) {
            result.displacement.add(i, alpha * du[i]);
        }
        
        assembler_->assemble_tangent_stiffness(result.displacement, assembly, true, is_transient, true);
        
        for (Index i = 0; i < n; ++i) {
            residual[i] = F_total[i] - assembly.F_int[i];
        }
        
        if (is_transient && !assembly.M.values().empty()) {
            assembly.M.multiply(result.acceleration, M_times_a);
            for (Index i = 0; i < n; ++i) {
                residual[i] -= M_times_a[i];
            }
        }
        
        for (Index gdof : model_->constrained_dofs()) {
            residual.set(gdof, 0.0);
        }
        
        result.residual_norm = residual.norm();
        
        Real rel = result.initial_residual_norm > EPS 
                   ? result.residual_norm / result.initial_residual_norm
                   : result.residual_norm;
        
        HVDC_LOG_DEBUG("NR iter " << iter 
                       << ": ||R||=" << result.residual_norm
                       << " (rel=" << rel << ")"
                       << ", alpha=" << alpha);
        
        if (result.residual_norm < config_.tol_abs || rel < config_.tol_rel) {
            result.converged = true;
            result.converged_reason = "tolerance";
            break;
        }
    }
    
    result.newton_iters = iter;
    
    if (result.converged) {
        HVDC_LOG_INFO("Newton-Raphson converged after " << iter 
                      << " iters: ||R||=" << result.residual_norm);
    } else {
        HVDC_LOG_WARNING("Newton-Raphson did NOT converge after " << iter
                         << " iters: ||R||=" << result.residual_norm);
    }
    
    return result.converged;
}

Real Solver::line_search(
    Vector& du,
    const Vector& residual,
    const Vector& F_ext_total,
    AssemblyResult& assembly)
{
    Real alpha = 1.0;
    Real initial_merit = 0.5 * residual.dot(residual);
    
    Index n = static_cast<Index>(du.size());
    Vector u_test(n);
    Vector resid_test(n);
    
    for (Index i = 0; i < config_.max_line_search; ++i) {
        u_test = assembly.F_int;
        for (Index j = 0; j < n; ++j) {
            u_test[j] += alpha * du[j];
        }
        
        Vector dummy_F(n);
        dummy_F.zero();
        
        for (Index j = 0; j < n; ++j) {
            resid_test[j] = F_ext_total[j] - u_test[j];
        }
        Real merit = 0.5 * resid_test.dot(resid_test);
        
        if (merit < initial_merit * (1.0 - 1.0e-4 * alpha)) {
            return alpha;
        }
        
        alpha *= config_.line_search_beta;
    }
    
    return alpha * config_.line_search_beta;
}

bool Solver::solve_transient_newmark(
    Real dt,
    const Vector& F_ext,
    const Vector& u_prev,
    const Vector& v_prev,
    const Vector& a_prev,
    SolverResult& result,
    Real beta,
    Real gamma)
{
    if (!initialized_) initialize();
    
    Index n = model_->total_dofs();
    result.displacement = u_prev;
    result.velocity = v_prev;
    result.acceleration = a_prev;
    
    if (result.displacement.size() != n) result.displacement.resize(n);
    if (result.velocity.size() != n) result.velocity.resize(n);
    if (result.acceleration.size() != n) result.acceleration.resize(n);
    
    Vector u_pred(n), v_pred(n);
    for (Index i = 0; i < n; ++i) {
        u_pred[i] = u_prev[i] + dt * v_prev[i] + 0.5 * dt * dt * (1.0 - 2.0 * beta) * a_prev[i];
        v_pred[i] = v_prev[i] + dt * (1.0 - gamma) * a_prev[i];
    }
    
    Vector effective_F = F_ext;
    AssemblyResult mass_assembly;
    assembler_->assemble_mass_matrix(mass_assembly.M, true, true);
    
    Vector a_prev_contrib(n);
    mass_assembly.M.multiply(a_prev, a_prev_contrib);
    for (Index i = 0; i < n; ++i) {
        effective_F[i] += a_prev_contrib[i] * (0.5 - beta) / (beta * dt * dt)
                        + mass_assembly.M(i, i) * (v_prev[i] / (beta * dt)
                          + (0.5 / beta - 1.0) * a_prev[i]);
    }
    
    Vector F_extra(n);
    F_extra.zero();
    for (Index i = 0; i < n; ++i) {
        F_extra[i] = mass_assembly.M(i, i) * u_pred[i] / (beta * dt * dt);
    }
    
    result.displacement = u_pred;
    bool ok = newton_raphson_iteration(F_ext, F_extra, result, true, dt,
                                        &u_prev, &v_prev, &a_prev, beta, gamma);
    
    if (ok) {
        for (Index i = 0; i < n; ++i) {
            result.acceleration[i] = (result.displacement[i] - u_pred[i]) / (beta * dt * dt);
            result.velocity[i] = v_pred[i] + gamma * dt * result.acceleration[i];
        }
        result.reaction_forces = compute_reactions(result.displacement);
    }
    
    return ok;
}

bool Solver::solve_modal(
    Index num_modes,
    std::vector<Real>& frequencies,
    std::vector<Vector>& mode_shapes)
{
#ifdef HVDC_USE_PETSC
    if (!initialized_) initialize();
    
    Index n = model_->total_dofs();
    num_modes = std::min(num_modes, static_cast<Index>(n / 10));
    num_modes = std::max<Index>(1, num_modes);
    
    AssemblyResult assembly;
    Vector zero_disp(n);
    zero_disp.zero();
    
    assembler_->assemble_tangent_stiffness(zero_disp, assembly, false, true, true);
    
    PETScMatrix K_petsc, M_petsc;
    K_petsc.from_sparse_matrix(assembly.K_T);
    M_petsc.from_sparse_matrix(assembly.M);
    
    EPS eps;
    EPSCreate(PETSC_COMM_WORLD, &eps);
    EPSSetOperators(eps, K_petsc.petsc_mat(), M_petsc.petsc_mat());
    EPSSetProblemType(eps, EPS_GHEP);
    EPSSetType(eps, EPSKRYLOVSCHUR);
    EPSSetDimensions(eps, num_modes, PETSC_DEFAULT, PETSC_DEFAULT);
    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
    EPSSetFromOptions(eps);
    EPSSolve(eps);
    
    Index nconv;
    EPSGetConverged(eps, &nconv);
    
    frequencies.resize(nconv);
    mode_shapes.resize(nconv, Vector(n));
    
    for (Index i = 0; i < nconv; ++i) {
        PetscScalar kr, ki;
        PETScVector x_petsc(n);
        EPSGetEigenpair(eps, i, &kr, &ki, x_petsc.petsc_vec(), nullptr);
        Real lambda = PetscRealPart(kr);
        frequencies[i] = (lambda > EPS) ? std::sqrt(lambda) / (2.0 * PI) : 0.0;
        x_petsc.to_vector(mode_shapes[i]);
    }
    
    EPSDestroy(&eps);
    
    HVDC_LOG_INFO("Modal analysis: " << nconv << " modes extracted");
    return nconv > 0;
#else
    frequencies.clear();
    mode_shapes.clear();
    return false;
#endif
}

Vector Solver::compute_reactions(const Vector& displacement) {
    Index n = model_->total_dofs();
    Vector reactions(n);
    reactions.zero();
    
    AssemblyResult assembly;
    assembler_->assemble_internal_forces(displacement, assembly.F_int, true);
    
    for (Index gdof : model_->constrained_dofs()) {
        reactions.set(gdof, assembly.F_int[gdof]);
    }
    
    return reactions;
}

void Solver::compute_strains_stresses(
    const Vector& displacement,
    std::vector<Real>& elem_strain,
    std::vector<Real>& elem_stress,
    std::vector<Vec3>& elem_force)
{
    Index n_elems = model_->num_elements();
    elem_strain.assign(n_elems, 0.0);
    elem_stress.assign(n_elems, 0.0);
    elem_force.assign(n_elems, {0.0, 0.0, 0.0});
    
    for (Index ei = 0; ei < n_elems; ++ei) {
        Element* elem = model_->element(ei);
        if (!elem) continue;
        
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
        
        if (auto* beam = dynamic_cast<BeamElementNL*>(elem)) {
            Real N, Vy, Vz, T, My, Mz;
            Mat3x3 T_mat;
            Vec3 dx;
            beam->compute_transformation(T_mat, dx, d_elem);
            Mat3x3 Tt = math::mat3_transpose(T_mat);
            
            Vec12 d_local{};
            for (Index bi = 0; bi < 4; ++bi) {
                Vec3 gb = {d_elem[bi*3], d_elem[bi*3+1], d_elem[bi*3+2]};
                Vec3 lb = math::mat3_mul_vec3(Tt, gb);
                d_local[bi*3] = lb[0]; d_local[bi*3+1] = lb[1]; d_local[bi*3+2] = lb[2];
            }
            
            beam->compute_local_forces(d_local, d_local, N, Vy, Vz, T, My, Mz);
            
            Real A = beam->section() ? beam->section()->A : 1.0;
            elem_strain[ei] = N / (beam->material()->E * A);
            elem_stress[ei] = N / A;
            elem_force[ei] = {N, Vy, Vz};
            
        } else if (auto* truss = dynamic_cast<TrussElementNL*>(elem)) {
            Real strain, stress, force;
            truss->compute_axial_response(d_elem, strain, stress, force);
            elem_strain[ei] = strain;
            elem_stress[ei] = stress;
            elem_force[ei] = {force, 0.0, 0.0};
            
        } else if (auto* cond = dynamic_cast<ConductorElement*>(elem)) {
            elem_strain[ei] = cond->axial_strain(d_elem);
            elem_stress[ei] = cond->axial_stress(d_elem);
            elem_force[ei] = {cond->current_tension(d_elem), 0.0, 0.0};
        }
    }
}

bool Solver::solve_linear_system(
    SparseMatrix& K,
    Vector& solution,
    const Vector& rhs)
{
    return solve_linear_system_impl(K, rhs, solution);
}

bool Solver::apply_newmark_time_integration(
    SparseMatrix& K,
    SparseMatrix& M,
    SparseMatrix& C,
    const Vector& F_ext,
    Vector& u,
    Vector& v,
    Vector& a,
    Real dt,
    Real beta,
    Real gamma)
{
    Index ndofs = u.size();
    if (v.size() != ndofs) v.resize(ndofs);
    if (a.size() != ndofs) a.resize(ndofs);
    
    Real a0 = 1.0 / (beta * dt * dt);
    Real a1 = gamma / (beta * dt);
    Real a2 = 1.0 / (beta * dt);
    Real a3 = 1.0 / (2.0 * beta) - 1.0;
    Real a4 = gamma / beta - 1.0;
    Real a5 = dt / 2.0 * (gamma / beta - 2.0);
    Real a6 = dt * (1.0 - gamma);
    Real a7 = gamma * dt;
    
    Vector u_pred = u;
    for (Index i = 0; i < ndofs; ++i) {
        u_pred[i] += dt * v[i] + 0.5 * dt * dt * (1.0 - 2.0 * beta) * a[i];
    }
    Vector a_new(ndofs), v_new(ndofs);
    for (Index i = 0; i < ndofs; ++i) {
        a_new[i] = a0 * (u_pred[i] - u[i]) - a2 * v[i] - a3 * a[i];
        v_new[i] = v[i] + a6 * a[i] + a7 * a_new[i];
    }
    
    SparseMatrix K_eff = K;
    if (M.rows() == ndofs) {
        for (Index r = 0; r < ndofs; ++r) {
            for (Index j = M.row_start(r); j < M.row_start(r + 1); ++j) {
                Index c = M.col(j);
                K_eff.add(r, c, a0 * M.value(j));
            }
        }
    }
    if (C.rows() == ndofs) {
        for (Index r = 0; r < ndofs; ++r) {
            for (Index j = C.row_start(r); j < C.row_start(r + 1); ++j) {
                Index c = C.col(j);
                K_eff.add(r, c, a1 * C.value(j));
            }
        }
    }
    K_eff.finalize();
    
    Vector F_eff = F_ext;
    if (M.rows() == ndofs) {
        Vector Ma(ndofs);
        for (Index r = 0; r < ndofs; ++r) {
            Real sum = 0.0;
            for (Index j = M.row_start(r); j < M.row_start(r + 1); ++j) {
                Index c = M.col(j);
                if (c < ndofs) sum += M.value(j) * a_new[c];
            }
            Ma[r] = sum;
        }
        for (Index i = 0; i < ndofs; ++i) F_eff[i] += Ma[i];
    }
    if (C.rows() == ndofs) {
        Vector Cv(ndofs);
        for (Index r = 0; r < ndofs; ++r) {
            Real sum = 0.0;
            for (Index j = C.row_start(r); j < C.row_start(r + 1); ++j) {
                Index c = C.col(j);
                if (c < ndofs) sum += C.value(j) * v_new[c];
            }
            Cv[r] = sum;
        }
        for (Index i = 0; i < ndofs; ++i) F_eff[i] += Cv[i];
    }
    
    Vector du(ndofs);
    bool ok = solve_linear_system(K_eff, du, F_eff);
    if (!ok) return false;
    
    for (Index i = 0; i < ndofs; ++i) {
        u[i] = u_pred[i] + du[i];
        Real du_val = du[i];
        a[i] = a0 * du_val - a2 * v[i] - a3 * a[i];
        v[i] = v[i] + a6 * a[i] + a7 * (a_new[i] + a0 * du_val);
    }
    
    return true;
}

} // namespace fem
} // namespace hvdc
