#include "fvm/NSComposable.hpp"
#include "common/MathUtils.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>

namespace hvdc {
namespace fvm {

NSComposable::NSComposable(Grid* grid, const FlowConfig& flow_cfg,
                             const DiscretizationConfig& disc_cfg)
    : grid_(grid), flow_cfg_(flow_cfg), disc_(grid, disc_cfg)
{
}

void NSComposable::initialize_state(FluidState& state) const {
    Index nc = grid_->num_cells();
    state.resize(nc);
    state.rho_ref = flow_cfg_.rho_ref;
    state.mu_ref = flow_cfg_.mu_ref;
    state.T_ref = flow_cfg_.T_ref;
    state.p_ref = flow_cfg_.p_ref;
    state.initialize_defaults(nc);
}

void NSComposable::compute_effective_viscosity(
    const FluidState& state,
    ScalarField& nu_eff,
    ScalarField& Gamma_U,
    ScalarField& Gamma_k,
    ScalarField& Gamma_eps) const
{
    Index nc = state.size();
    nu_eff.resize(nc);
    Gamma_U.resize(nc);
    Gamma_k.resize(nc);
    Gamma_eps.resize(nc);
    
    for (Index ci = 0; ci < nc; ++ci) {
        Real nu_lam = state.nu[ci];
        Real nu_turb = (flow_cfg_.enable_turbulence) ? state.turb.nut(ci) : 0.0;
        nu_eff[ci] = nu_lam + nu_turb;
        Gamma_U[ci] = state.rho[ci] * nu_eff[ci];
        if (flow_cfg_.enable_turbulence) {
            Gamma_k[ci] = state.rho[ci] * (nu_lam + nu_turb / flow_cfg_.sigma_k);
            Gamma_eps[ci] = state.rho[ci] * (nu_lam + nu_turb / flow_cfg_.sigma_eps);
        } else {
            Gamma_k[ci] = Gamma_U[ci];
            Gamma_eps[ci] = Gamma_U[ci];
        }
    }
}

void NSComposable::assemble_momentum_predictor(
    FluidState& state,
    const FluidState* state_prev,
    Real dt,
    Index component,
    SparseMatrix& A_U,
    Vector& b_U,
    VectorField& HbyA,
    Vector& diag_inv_Ap,
    VectorField& grad_p) const
{
    Index nc = state.size();
    grad_p.resize(nc);
    
    ScalarField p_field = state.p;
    std::vector<Vec3> grad_p_vec(nc);
    for (Index ci = 0; ci < nc; ++ci) {
        grad_p_vec[ci] = state.p.pressure_gradient(ci, *grid_);
        grad_p[ci] = grad_p_vec[ci];
    }
    
    ScalarField nu_eff, Gamma_U, Gamma_k, Gamma_eps;
    compute_effective_viscosity(state, nu_eff, Gamma_U, Gamma_k, Gamma_eps);
    
    ScalarField U_comp = state.U.component(component);
    std::vector<std::array<Vec3, 3>> grad_U;
    disc_.compute_vector_gradient(state.U, {}, grad_U);
    
    std::vector<Real> face_phi(grid_->num_faces(), 0.0);
    std::vector<Vec3> grad_vecp(nc);
    for (Index ci = 0; ci < nc; ++ci) grad_vecp[ci] = grad_p_vec[ci];
    
    disc_.discretize_momentum_equation(
        state.U, p_field, state.rho, nu_eff, grad_U, grad_vecp,
        component, A_U, b_U,
        state_prev ? state.rho[0] : 0.0,
        state_prev ? &state_prev->U : nullptr, dt);
    
    diag_inv_Ap.resize(nc);
    HbyA.resize(nc);
    for (Index ci = 0; ci < nc; ++ci) {
        Real diag = A_U(ci, ci);
        diag_inv_Ap[ci] = (std::fabs(diag) > EPS) ? 1.0 / diag : 0.0;
    }
    
    for (Index ci = 0; ci < nc; ++ci) {
        Real H_U = b_U[ci];
        for (Index j = A_U.row_ptr()[ci]; j < A_U.row_ptr()[ci+1]; ++j) {
            Index c = A_U.col_idx()[j];
            if (c != ci) {
                H_U -= A_U.values()[j] * state.U[c][component];
            }
        }
        HbyA[ci][component] = H_U * diag_inv_Ap[ci];
    }
}

void NSComposable::assemble_pressure_correction(
    FluidState& state,
    const SparseMatrix& Ap,
    const Vector& diag_inv_Ap,
    const VectorField& HbyA,
    SparseMatrix& A_p,
    Vector& b_p,
    VectorField& face_velocity,
    std::vector<Real>& face_mass_flux,
    const VectorField& grad_p) const
{
    (void)Ap;
    disc_.discretize_pressure_equation(
        state.U, state.rho, Ap, diag_inv_Ap, A_p, b_p, &HbyA);
    
    disc_.compute_face_velocity_rhie_chow(
        state.U, state.p, Ap, diag_inv_Ap,
        face_velocity, face_mass_flux, grad_p);
}

void NSComposable::correct_velocity_and_pressure(
    FluidState& state,
    const SparseMatrix& Ap,
    const Vector& diag_inv_Ap,
    const VectorField& HbyA,
    const ScalarField& dp,
    VectorField& face_velocity,
    std::vector<Real>& face_mass_flux,
    const VectorField& grad_p) const
{
    Index nc = state.size();
    
    VectorField dp_grad(nc);
    std::vector<Real> face_dp(grid_->num_faces(), 0.0);
    std::vector<Vec3> dp_grad_vec(nc);
    disc_.compute_scalar_gradient(dp, face_dp, dp_grad_vec);
    
    VectorField dp_correction = HbyA;
    for (Index ci = 0; ci < nc; ++ci) {
        dp_correction[ci] = dp_grad_vec[ci];
    }
    
    disc_.correct_velocity(
        state.U, Ap, diag_inv_Ap, state.rho, grad_p, dp_correction);
    
    for (Index ci = 0; ci < nc; ++ci) {
        state.p[ci] += dp[ci];
    }
    
    disc_.compute_face_velocity_rhie_chow(
        state.U, state.p, Ap, diag_inv_Ap,
        face_velocity, face_mass_flux, grad_p);
}

ScalarField NSComposable::compute_kinetic_energy_dissipation(
    const FluidState& state,
    const std::vector<std::array<Vec3, 3>>& grad_U) const
{
    Index nc = state.size();
    ScalarField epsilon(nc, 0.0);
    
    for (Index ci = 0; ci < nc; ++ci) {
        Mat3x3 S;
        compute_strain_rate_tensor(grad_U[ci], S);
        Real mag_S = compute_strain_rate_magnitude(S);
        Real nut = state.turb.nut(ci);
        Real k_val = state.turb.k[ci];
        if (k_val > EPS && flow_cfg_.enable_turbulence) {
            epsilon[ci] = flow_cfg_.C_mu * k_val * k_val * k_val 
                          / (k_val + EPS);
        }
        epsilon[ci] += 2.0 * state.nu[ci] * mag_S * mag_S;
    }
    return epsilon;
}

void NSComposable::assemble_turbulence_equations(
    FluidState& state,
    const FluidState* state_prev,
    Real dt,
    const std::vector<std::array<Vec3, 3>>& grad_U,
    SparseMatrix& A_k, Vector& b_k,
    SparseMatrix& A_eps, Vector& b_eps) const
{
    if (!flow_cfg_.enable_turbulence) return;
    
    Index nc = state.size();
    ScalarField P_k = compute_production_k(state, grad_U);
    
    ScalarField nu_eff, Gamma_U, Gamma_k, Gamma_eps;
    compute_effective_viscosity(state, nu_eff, Gamma_U, Gamma_k, Gamma_eps);
    
    std::vector<Vec3> grad_k(nc);
    {
        std::vector<Real> face_k(grid_->num_faces(), 0.0);
        disc_.compute_scalar_gradient(state.turb.k, face_k, grad_k);
    }
    
    disc_.discretize_convection_diffusion_scalar(
        state.turb.k, state.rho, Gamma_k, grad_k, A_k, b_k,
        state_prev ? state.rho[0] : 0.0,
        state_prev ? &state_prev->turb.k : nullptr, dt, &state.U);
    
    for (Index ci = 0; ci < nc; ++ci) {
        Real k_val = std::max(state.turb.k[ci], 1.0e-12);
        Real rho = state.rho[ci];
        Real Vp = grid_->cell(ci).volume();
        
        b_k.add(ci, rho * P_k[ci] * Vp);
        Real epsilon = std::max(state.turb.epsilon[ci], 1.0e-12);
        A_k.add(ci, ci, rho * epsilon / k_val * Vp);
    }
    
    std::vector<Vec3> grad_eps(nc);
    {
        std::vector<Real> face_eps(grid_->num_faces(), 0.0);
        disc_.compute_scalar_gradient(state.turb.epsilon, face_eps, grad_eps);
    }
    
    disc_.discretize_convection_diffusion_scalar(
        state.turb.epsilon, state.rho, Gamma_eps, grad_eps, A_eps, b_eps,
        state_prev ? state.rho[0] : 0.0,
        state_prev ? &state_prev->turb.epsilon : nullptr, dt, &state.U);
    
    for (Index ci = 0; ci < nc; ++ci) {
        Real k_val = std::max(state.turb.k[ci], 1.0e-12);
        Real eps_val = std::max(state.turb.epsilon[ci], 1.0e-12);
        Real rho = state.rho[ci];
        Real Vp = grid_->cell(ci).volume();
        
        b_eps.add(ci, flow_cfg_.C1_eps * rho * eps_val / k_val * P_k[ci] * Vp);
        A_eps.add(ci, ci, flow_cfg_.C2_eps * rho * eps_val / k_val * Vp);
    }
}

void NSComposable::update_turbulence_from_solution(
    FluidState& state,
    const ScalarField& k_new,
    const ScalarField& eps_new) const
{
    Index nc = state.size();
    for (Index ci = 0; ci < nc; ++ci) {
        state.turb.k[ci] = std::max(k_new[ci], 1.0e-12);
        state.turb.epsilon[ci] = std::max(eps_new[ci], 1.0e-12);
    }
    state.turb.compute_nut_from_k_epsilon(flow_cfg_.C_mu);
}

Real NSComposable::compute_cfl_number(
    const FluidState& state,
    const std::vector<Real>& face_flux,
    Real dt) const
{
    Index nc = grid_->num_cells();
    Real cfl_max = 0.0;
    
    for (Index ci = 0; ci < nc; ++ci) {
        const Cell& cell = grid_->cell(ci);
        Real sum_in = 0.0;
        Real sum_out = 0.0;
        
        for (Index fi = 0; fi < cell.num_faces() && fi < cell.face_ids().size(); ++fi) {
            Index face_id = cell.face_ids()[fi];
            if (face_id < 0) continue;
            const Face& face = grid_->face(face_id);
            Real flux = face_flux[face_id];
            if (face.owner() == ci) {
                if (flux > 0) sum_out += flux;
                else sum_in -= flux;
            } else {
                if (flux > 0) sum_in += flux;
                else sum_out -= flux;
            }
        }
        
        Real V = std::max(cell.volume(), EPS);
        Real cfl = std::max(sum_in, sum_out) * dt / (state.rho[ci] * V);
        cfl_max = std::max(cfl_max, cfl);
    }
    
    return cfl_max;
}

Vec3 NSComposable::compute_total_force_on_boundary(
    Index patch_id,
    const FluidState& state,
    const VectorField& grad_p,
    const std::vector<std::array<Vec3, 3>>& grad_U) const
{
    Vec3 force = {0.0, 0.0, 0.0};
    auto range = grid_->boundary_patch_range(patch_id);
    if (range.first < 0) return force;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid_->face(fi);
        Vec3 nA = math::vec3_scale(face.normal(), face.area());
        Index ci = face.owner() >= 0 ? face.owner() : face.neighbor();
        if (ci < 0) continue;
        
        Vec3 n = face.normal();
        Real p = state.p[ci] - state.p.reference_pressure();
        Vec3 pressure_force = math::vec3_scale(nA, p);
        math::vec3_add_inplace(force, pressure_force);
        
        if (ci < static_cast<Index>(grad_U.size()) && flow_cfg_.enable_turbulence) {
            const auto& grad = grad_U[ci];
            Vec3 tau_n = {0, 0, 0};
            for (Index i = 0; i < 3; ++i) {
                for (Index j = 0; j < 3; ++j) {
                    Real tau_ij = state.mu[ci] * (grad[i][j] + grad[j][i]);
                    if (flow_cfg_.enable_turbulence) {
                        tau_ij += state.rho[ci] * state.turb.nut(ci) * (grad[i][j] + grad[j][i]);
                        if (i == j) tau_ij -= 2.0/3.0 * state.rho[ci] * state.turb.k[ci];
                    }
                    tau_n[i] += tau_ij * n[j];
                }
            }
            math::vec3_add_inplace(force, math::vec3_scale(tau_n, -face.area()));
        }
    }
    (void)grad_p;
    return force;
}

Real NSComposable::compute_mass_imbalance(const std::vector<Real>& face_flux) const {
    Index nc = grid_->num_cells();
    Real max_imbalance = 0.0;
    
    for (Index ci = 0; ci < nc; ++ci) {
        const Cell& cell = grid_->cell(ci);
        Real sum = 0.0;
        for (Index fi = 0; fi < cell.num_faces() && fi < cell.face_ids().size(); ++fi) {
            Index face_id = cell.face_ids()[fi];
            if (face_id < 0) continue;
            const Face& face = grid_->face(face_id);
            if (face.owner() == ci) sum += face_flux[face_id];
            else if (face.neighbor() == ci) sum -= face_flux[face_id];
        }
        max_imbalance = std::max(max_imbalance, std::fabs(sum));
    }
    return max_imbalance;
}

std::array<Real, 6> NSComposable::compute_flow_stats(const FluidState& state) const {
    Index nc = state.size();
    Real U_max = 0.0, p_min = INF, p_max = -INF;
    Real U_avg = 0.0, p_avg = 0.0, k_max = 0.0;
    Real total_vol = 0.0;
    
    for (Index ci = 0; ci < nc; ++ci) {
        Real V = grid_->cell(ci).volume();
        if (V <= 0) continue;
        total_vol += V;
        Real speed = state.U.speed(ci);
        U_max = std::max(U_max, speed);
        U_avg += speed * V;
        p_min = std::min(p_min, state.p[ci]);
        p_max = std::max(p_max, state.p[ci]);
        p_avg += state.p[ci] * V;
        if (flow_cfg_.enable_turbulence) {
            k_max = std::max(k_max, state.turb.k[ci]);
        }
    }
    
    if (total_vol > EPS) {
        U_avg /= total_vol;
        p_avg /= total_vol;
    }
    
    return {U_max, U_avg, p_min, p_max, p_avg, k_max};
}

ScalarField NSComposable::compute_production_k(
    const FluidState& state,
    const std::vector<std::array<Vec3, 3>>& grad_U) const
{
    Index nc = state.size();
    ScalarField P_k(nc, 0.0);
    
    for (Index ci = 0; ci < nc; ++ci) {
        Mat3x3 S;
        compute_strain_rate_tensor(grad_U[ci], S);
        Real mag_S = compute_strain_rate_magnitude(S);
        P_k[ci] = state.turb.nut(ci) * 2.0 * mag_S * mag_S;
    }
    return P_k;
}

void NSComposable::compute_strain_rate_tensor(
    const std::array<Vec3, 3>& grad,
    Mat3x3& S) const
{
    for (Index i = 0; i < 3; ++i) {
        for (Index j = 0; j < 3; ++j) {
            S[i][j] = 0.5 * (grad[i][j] + grad[j][i]);
        }
    }
}

Real NSComposable::compute_strain_rate_magnitude(const Mat3x3& S) const {
    Real sum = 0.0;
    for (Index i = 0; i < 3; ++i) {
        for (Index j = 0; j < 3; ++j) {
            sum += 2.0 * S[i][j] * S[i][j];
        }
    }
    return std::sqrt(std::max(sum, 0.0));
}

} // namespace fvm
} // namespace hvdc
