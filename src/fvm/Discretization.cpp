#include "fvm/Discretization.hpp"
#include "common/MathUtils.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>

namespace hvdc {
namespace fvm {

Discretization::Discretization(Grid* grid, const DiscretizationConfig& cfg)
    : grid_(grid), config_(cfg)
{
}

Real Discretization::interpolate_central(
    Index face_id,
    const ScalarField& phi) const
{
    const Face& face = grid_->face(face_id);
    Index P = face.owner();
    Index N = face.neighbor();
    
    if (face.is_on_boundary()) {
        if (P >= 0) return phi[P];
        if (N >= 0) return phi[N];
        return 0.0;
    }
    
    Real gc = face.weight();
    return gc * phi[P] + (1.0 - gc) * phi[N];
}

Vec3 Discretization::interpolate_central_vector(
    Index face_id,
    const VectorField& U) const
{
    const Face& face = grid_->face(face_id);
    Index P = face.owner();
    Index N = face.neighbor();
    
    if (face.is_on_boundary()) {
        Index ci = (P >= 0) ? P : N;
        return (ci >= 0) ? U[ci] : Vec3{0,0,0};
    }
    
    Real gc = face.weight();
    Vec3 result;
    for (Index i = 0; i < 3; ++i) {
        result[i] = gc * U[P][i] + (1.0 - gc) * U[N][i];
    }
    return result;
}

Real Discretization::interpolate_second_order_upwind(
    Index face_id,
    const ScalarField& phi,
    const std::vector<Vec3>& grad,
    Real m_dot_sign) const
{
    const Face& face = grid_->face(face_id);
    Index P = face.owner();
    Index N = face.neighbor();
    
    if (face.is_on_boundary()) {
        return (m_dot_sign >= 0) ? phi[P] : phi[N >= 0 ? N : P];
    }
    
    Index upwind = (m_dot_sign >= 0) ? P : N;
    Index downwind = (m_dot_sign >= 0) ? N : P;
    (void)downwind;
    
    const Cell& cell = grid_->cell(upwind);
    const Vec3& d = face.delta();
    Real f_norm = face.delta_mag();
    if (f_norm < EPS) return phi[upwind];
    
    Vec3 d_unit = math::vec3_scale(d, 1.0 / f_norm);
    Vec3 face_vec = math::vec3_sub(face.centroid(), cell.centroid());
    Real proj = math::vec3_dot(face_vec, d_unit);
    
    if (upwind < static_cast<Index>(grad.size())) {
        Real grad_component = math::vec3_dot(grad[upwind], face_vec);
        return phi[upwind] + grad_component;
    }
    return phi[upwind];
}

void Discretization::interpolate_face_scalar(
    const ScalarField& cell_field,
    const std::vector<Vec3>& cell_grad,
    Index face_id,
    Real& face_value,
    Real& face_upwind_value) const
{
    const Face& face = grid_->face(face_id);
    face_value = interpolate_central(face_id, cell_field);
    face_upwind_value = face_value;
    
    if (face.is_on_boundary()) {
        return;
    }
    
    Index P = face.owner();
    Index N = face.neighbor();
    
    face_upwind_value = cell_field[P];
    (void)N; (void)cell_grad;
}

void Discretization::interpolate_face_vector(
    const VectorField& cell_field,
    const std::vector<std::array<Vec3, 3>>& cell_grad,
    Index face_id,
    Vec3& face_value,
    Vec3& face_upwind_value) const
{
    (void)cell_grad;
    face_value = interpolate_central_vector(face_id, cell_field);
    const Face& face = grid_->face(face_id);
    Index P = face.owner();
    if (P >= 0) face_upwind_value = cell_field[P];
    else face_upwind_value = face_value;
}

void Discretization::compute_scalar_gradient(
    const ScalarField& phi,
    const std::vector<Real>& face_phi,
    std::vector<Vec3>& grad) const
{
    Index nc = grid_->num_cells();
    grad.assign(nc, math::vec3_zero());
    
    for (Index ci = 0; ci < nc; ++ci) {
        const Cell& cell = grid_->cell(ci);
        Real V = cell.volume();
        if (V < EPS) continue;
        
        Vec3 sum = math::vec3_zero();
        for (Index fi = 0; fi < cell.num_faces() && fi < cell.face_ids().size(); ++fi) {
            Index face_id = cell.face_ids()[fi];
            if (face_id < 0) continue;
            const Face& face = grid_->face(face_id);
            
            Real phi_f;
            if (face_id < static_cast<Index>(face_phi.size())) {
                phi_f = face_phi[face_id];
            } else {
                phi_f = interpolate_central(face_id, phi);
            }
            
            Vec3 nA = math::vec3_scale(face.normal(), face.area());
            math::vec3_add_inplace(sum, math::vec3_scale(nA, phi_f));
        }
        
        grad[ci] = math::vec3_scale(sum, 1.0 / V);
    }
}

void Discretization::compute_vector_gradient(
    const VectorField& U,
    const std::vector<Vec3>& face_U,
    std::vector<std::array<Vec3, 3>>& grad) const
{
    Index nc = grid_->num_cells();
    grad.assign(nc, {});
    for (auto& g : grad) g.fill(math::vec3_zero());
    
    for (Index ci = 0; ci < nc; ++ci) {
        const Cell& cell = grid_->cell(ci);
        Real V = cell.volume();
        if (V < EPS) continue;
        
        std::array<Vec3, 3> sums;
        sums.fill(math::vec3_zero());
        
        for (Index fi = 0; fi < cell.num_faces() && fi < cell.face_ids().size(); ++fi) {
            Index face_id = cell.face_ids()[fi];
            if (face_id < 0) continue;
            const Face& face = grid_->face(face_id);
            
            Vec3 Uf;
            if (face_id < static_cast<Index>(face_U.size())) {
                Uf = face_U[face_id];
            } else {
                Uf = interpolate_central_vector(face_id, U);
            }
            
            Vec3 nA = math::vec3_scale(face.normal(), face.area());
            for (Index comp = 0; comp < 3; ++comp) {
                math::vec3_add_inplace(sums[comp], math::vec3_scale(nA, Uf[comp]));
            }
        }
        
        for (Index comp = 0; comp < 3; ++comp) {
            grad[ci][comp] = math::vec3_scale(sums[comp], 1.0 / V);
        }
    }
}

Real Discretization::face_mass_flux(
    Index face_id,
    const Vec3& face_U,
    Real face_rho) const
{
    const Face& face = grid_->face(face_id);
    Vec3 nA = math::vec3_scale(face.normal(), face.area());
    return face_rho * math::vec3_dot(face_U, nA);
}

void Discretization::discretize_convection_diffusion_scalar(
    const ScalarField& phi,
    const ScalarField& rho,
    const ScalarField& Gamma,
    const std::vector<Vec3>& grad_phi,
    SparseMatrix& A,
    Vector& b,
    Real transient_rho,
    const ScalarField* phi_prev,
    Real dt,
    const VectorField* U) const
{
    Index nc = grid_->num_cells();
    A.reset(nc, nc);
    b.resize(nc);
    b.zero();
    
    std::vector<Real> face_phi(grid_->num_faces());
    std::vector<Real> face_rho(grid_->num_faces());
    std::vector<Real> face_Gamma(grid_->num_faces());
    std::vector<Real> face_flux(grid_->num_faces(), 0.0);
    
    for (Index fi = 0; fi < grid_->num_faces(); ++fi) {
        const Face& face = grid_->face(fi);
        Index P = face.owner();
        Index N = face.neighbor();
        Real gc = face.weight();
        
        if (face.is_on_boundary()) {
            Index ci = (P >= 0) ? P : N;
            face_phi[fi] = phi[ci];
            face_rho[fi] = rho[ci];
            face_Gamma[fi] = Gamma[ci];
        } else {
            face_rho[fi] = gc * rho[P] + (1.0 - gc) * rho[N];
            face_Gamma[fi] = gc * Gamma[P] + (1.0 - gc) * Gamma[N];
            face_phi[fi] = interpolate_central(fi, phi);
        }
        
        if (U) {
            Vec3 Uf = interpolate_central_vector(fi, *U);
            face_flux[fi] = face_mass_flux(fi, Uf, face_rho[fi]);
        }
    }
    
    for (Index ci = 0; ci < nc; ++ci) {
        const Cell& cell = grid_->cell(ci);
        Real Vp = cell.volume();
        
        Real a_p = 0.0;
        Real b_p = 0.0;
        
        if (transient_rho > 0.0 && phi_prev && dt > EPS) {
            a_p += transient_rho * Vp / dt;
            b_p += transient_rho * Vp / dt * (*phi_prev)[ci];
        }
        
        for (Index fi = 0; fi < cell.num_faces() && fi < cell.face_ids().size(); ++fi) {
            Index face_id = cell.face_ids()[fi];
            if (face_id < 0) continue;
            
            const Face& face = grid_->face(face_id);
            Vec3 Sf = math::vec3_scale(face.normal(), face.area());
            Real Af = face.area();
            Real dmag = face.delta_mag();
            if (dmag < EPS) continue;
            
            Index P = face.owner();
            Index N = face.neighbor();
            
            Real m_dot = face_flux[face_id];
            Real Gam_f = face_Gamma[face_id];
            Real D_f = Gam_f * Af / dmag;
            
            Real Ff = m_dot;
            Real a_nb = 0.0;
            
            if (face.is_on_boundary()) {
                Real phi_bc = face_phi[face_id];
                if (m_dot >= 0) {
                    b_p += std::max(Ff, 0.0) * phi_bc + D_f * phi_bc;
                    a_p += std::max(Ff, 0.0) + D_f;
                } else {
                    a_p += D_f;
                }
            } else {
                switch (config_.conv_scheme) {
                    case ConvectionScheme::FirstOrderUpwind:
                        a_nb = D_f + std::max(-Ff, 0.0);
                        a_p += D_f + std::max(Ff, 0.0);
                        break;
                        
                    case ConvectionScheme::SecondOrderUpwind:
                    default:
                        Real blend = config_.upwind_blending;
                        Real central_coeff = 0.5;
                        Real phi_up_P, phi_up_N;
                        
                        if (!grad_phi.empty()) {
                            phi_up_P = interpolate_second_order_upwind(
                                face_id, phi, grad_phi, 1.0);
                            phi_up_N = interpolate_second_order_upwind(
                                face_id, phi, grad_phi, -1.0);
                        } else {
                            phi_up_P = phi[P];
                            phi_up_N = phi[N];
                        }
                        
                        a_nb = D_f * central_coeff + std::max(-Ff, 0.0);
                        a_p = a_p + D_f * central_coeff + std::max(Ff, 0.0);
                        
                        if (std::fabs(m_dot) > EPS && blend < 1.0) {
                            Real extra = (m_dot > 0) 
                                ? m_dot * (phi_up_P - (1.0 - blend) * phi[P])
                                : -m_dot * (phi_up_N - (1.0 - blend) * phi[N]);
                            b_p -= extra * (1.0 - blend);
                        }
                        break;
                }
                
                if (ci == P) {
                    if (N >= 0) A.add(P, N, -a_nb);
                    a_p += 0;
                } else if (ci == N) {
                    if (P >= 0) A.add(N, P, -a_nb);
                    a_p += 0;
                }
            }
        }
        
        if (a_p > EPS) {
            A.add(ci, ci, a_p);
        }
        b.add(ci, b_p);
    }
    
    A.finalize();
}

void Discretization::discretize_momentum_equation(
    const VelocityField& U,
    const ScalarField& p,
    const ScalarField& rho,
    const ScalarField& nu_eff,
    const std::vector<std::array<Vec3, 3>>& grad_U,
    const std::vector<Vec3>& grad_p,
    Index component,
    SparseMatrix& A,
    Vector& b,
    Real transient_rho,
    const VelocityField* U_prev,
    Real dt) const
{
    Index nc = grid_->num_cells();
    A.reset(nc, nc);
    b.resize(nc);
    b.zero();
    
    std::vector<Real> face_rho(grid_->num_faces());
    std::vector<Real> face_nu(grid_->num_faces());
    std::vector<Real> face_flux(grid_->num_faces(), 0.0);
    std::vector<Vec3> face_vel(grid_->num_faces());
    
    for (Index fi = 0; fi < grid_->num_faces(); ++fi) {
        const Face& face = grid_->face(fi);
        Index P = face.owner();
        Index N = face.neighbor();
        Real gc = face.weight();
        
        if (face.is_on_boundary()) {
            Index ci = (P >= 0) ? P : N;
            face_rho[fi] = rho[ci];
            face_nu[fi] = nu_eff[ci];
            face_vel[fi] = U[ci];
        } else {
            face_rho[fi] = gc * rho[P] + (1.0 - gc) * rho[N];
            face_nu[fi] = gc * nu_eff[P] + (1.0 - gc) * nu_eff[N];
            face_vel[fi] = interpolate_central_vector(fi, U);
        }
        
        face_flux[fi] = face_mass_flux(fi, face_vel[fi], face_rho[fi]);
    }
    
    for (Index ci = 0; ci < nc; ++ci) {
        const Cell& cell = grid_->cell(ci);
        Real Vp = cell.volume();
        
        Real a_p = 0.0;
        Real b_p = 0.0;
        
        if (transient_rho > 0.0 && U_prev && dt > EPS) {
            a_p += rho[ci] * Vp / dt;
            b_p += rho[ci] * Vp / dt * (*U_prev)[ci][component];
        }
        
        Vec3 grad_p_ci = (ci < static_cast<Index>(grad_p.size())) 
                          ? grad_p[ci] : Vec3{0,0,0};
        b_p -= Vp * grad_p_ci[component];
        
        for (Index fi = 0; fi < cell.num_faces() && fi < cell.face_ids().size(); ++fi) {
            Index face_id = cell.face_ids()[fi];
            if (face_id < 0) continue;
            
            const Face& face = grid_->face(face_id);
            Vec3 Sf = math::vec3_scale(face.normal(), face.area());
            Real Af = face.area();
            Real dmag = face.delta_mag();
            if (dmag < EPS) continue;
            
            Index P = face.owner();
            Index N = face.neighbor();
            
            Real m_dot = face_flux[face_id];
            Real nu_f = face_nu[face_id];
            Real rho_f = face_rho[fi];
            Real Gamma_f = rho_f * nu_f;
            Real D_f = Gamma_f * Af / dmag;
            
            if (face.is_on_boundary()) {
                Real U_bc = face_vel[face_id][component];
                
                a_p += std::max(m_dot, 0.0) + D_f;
                b_p += (std::max(-m_dot, 0.0) + D_f) * U_bc;
                
            } else {
                Real Ff = m_dot;
                Real a_nb = D_f + std::max(-Ff, 0.0);
                a_p += D_f + std::max(Ff, 0.0);
                
                if (ci == P && N >= 0) {
                    A.add(P, N, -a_nb);
                } else if (ci == N && P >= 0) {
                    A.add(N, P, -a_nb);
                }
                
                if (!grad_U.empty() && ci < static_cast<Index>(grad_U.size())) {
                    const auto& grad = grad_U[ci];
                    Vec3 grad_comp = grad[component];
                    Vec3 nonortho = math::vec3_sub(Sf, 
                        math::vec3_scale(face.delta(), 
                            Af * Af / std::max(math::vec3_dot(face.delta(), Sf), EPS)));
                    b_p -= Gamma_f * math::vec3_dot(grad_comp, nonortho);
                }
            }
        }
        
        if (a_p > EPS) {
            A.add(ci, ci, a_p);
        }
        b.add(ci, b_p);
    }
    
    A.finalize();
}

void Discretization::discretize_pressure_equation(
    const VelocityField& U,
    const ScalarField& rho,
    const SparseMatrix& Ap,
    const Vector& diag_inv_Ap,
    SparseMatrix& A_p,
    Vector& b_p,
    const VectorField* HbyA) const
{
    Index nc = grid_->num_cells();
    A_p.reset(nc, nc);
    b_p.resize(nc);
    b_p.zero();
    
    for (Index ci = 0; ci < nc; ++ci) {
        const Cell& cell = grid_->cell(ci);
        Real Vp = cell.volume();
        Real sum_flux = 0.0;
        Real diag_ap = (ci < diag_inv_Ap.size()) ? diag_inv_Ap[ci] : 1.0;
        
        for (Index fi = 0; fi < cell.num_faces() && fi < cell.face_ids().size(); ++fi) {
            Index face_id = cell.face_ids()[fi];
            if (face_id < 0) continue;
            const Face& face = grid_->face(face_id);
            Vec3 Sf = math::vec3_scale(face.normal(), face.area());
            Real Af = face.area();
            Real dmag = face.delta_mag();
            if (dmag < EPS) continue;
            
            Index P = face.owner();
            Index N = face.neighbor();
            
            Real gc = face.weight();
            Real r_P = (P >= 0 && P < diag_inv_Ap.size()) ? diag_inv_Ap[P] : diag_ap;
            Real r_N = (N >= 0 && N < diag_inv_Ap.size()) ? diag_inv_Ap[N] : diag_ap;
            Real r_f = gc * r_P + (1.0 - gc) * r_N;
            Real rho_f = gc * rho[P] + (1.0 - gc) * rho[N];
            
            Real a_nb = rho_f * Af * Af * r_f / std::max(dmag, EPS);
            
            Vec3 Uf = interpolate_central_vector(face_id, U);
            Real flux = rho_f * math::vec3_dot(Uf, Sf);
            
            if (HbyA) {
                Vec3 HA_f = interpolate_central_vector(face_id, *HbyA);
                flux = rho_f * math::vec3_dot(HA_f, Sf);
            }
            
            if (face.is_on_boundary()) {
                A_p.add(ci, ci, a_nb);
                b_p.add(ci, flux);
            } else {
                A_p.add(ci, ci, a_nb);
                if (ci == P && N >= 0) A_p.add(P, N, -a_nb);
                if (ci == N && P >= 0) A_p.add(N, P, -a_nb);
                if (ci == P) b_p.add(ci, flux);
                if (ci == N) b_p.add(ci, -flux);
            }
        }
        
        if (cell.is_on_boundary()) {
            A_p.add(ci, ci, 1.0e-10);
        }
        (void)Vp;
    }
    
    A_p.finalize();
}

void Discretization::compute_face_velocity_rhie_chow(
    const VelocityField& U,
    const ScalarField& p,
    const SparseMatrix& Ap,
    const Vector& diag_inv_Ap,
    std::vector<Vec3>& face_velocity,
    std::vector<Real>& face_mass_flux_out,
    const std::vector<Vec3>& grad_p) const
{
    Index nf = grid_->num_faces();
    face_velocity.resize(nf);
    face_mass_flux_out.assign(nf, 0.0);
    
    for (Index fi = 0; fi < nf; ++fi) {
        const Face& face = grid_->face(fi);
        Index P = face.owner();
        Index N = face.neighbor();
        Real gc = face.weight();
        Vec3 Sf = math::vec3_scale(face.normal(), face.area());
        
        if (face.is_on_boundary()) {
            Index ci = (P >= 0) ? P : N;
            face_velocity[fi] = U[ci];
            Real rho_f = rho[ci];
            face_mass_flux_out[fi] = rho_f * math::vec3_dot(U[ci], Sf);
            continue;
        }
        
        Real r_P = (P >= 0 && P < diag_inv_Ap.size()) ? diag_inv_Ap[P] : 1.0;
        Real r_N = (N >= 0 && N < diag_inv_Ap.size()) ? diag_inv_Ap[N] : 1.0;
        Real gc_inv = 1.0 - gc;
        
        Vec3 Uf = interpolate_central_vector(fi, U);
        
        Vec3 gp_P = (P >= 0 && P < static_cast<Index>(grad_p.size())) ? grad_p[P] : Vec3{0,0,0};
        Vec3 gp_N = (N >= 0 && N < static_cast<Index>(grad_p.size())) ? grad_p[N] : Vec3{0,0,0};
        
        Vec3 gp_f;
        for (Index i = 0; i < 3; ++i) {
            gp_f[i] = gc * gp_P[i] + gc_inv * gp_N[i];
        }
        
        Vec3 delta_p = math::vec3_sub(face.delta(), 
            math::vec3_scale(face.delta(), 0.0));
        Real dmag = face.delta_mag();
        Real p_diff = (p[N] - p[P]) / std::max(dmag, EPS);
        
        Vec3 correction = math::vec3_zero();
        Real r_f = gc * r_P + gc_inv * r_N;
        for (Index i = 0; i < 3; ++i) {
            correction[i] = -r_f * (gp_f[i] - delta_p * face.normal()[i]);
        }
        
        for (Index i = 0; i < 3; ++i) {
            face_velocity[fi][i] = Uf[i] + correction[i];
        }
        
        Real rho_f = gc * rho[P] + gc_inv * rho[N];
        face_mass_flux_out[fi] = rho_f * math::vec3_dot(face_velocity[fi], Sf);
    }
}

void Discretization::apply_relaxation(
    SparseMatrix& A,
    Vector& b,
    Vector& field,
    Real relax_factor) const
{
    if (relax_factor >= 1.0) return;
    
    Index n = field.size();
    for (Index i = 0; i < n; ++i) {
        Real diag = A(i, i);
        if (diag > EPS) {
            A.set(i, i, diag / relax_factor);
            Real off_diag_contrib = diag / relax_factor * (1.0 - relax_factor) * field[i];
            b.add(i, off_diag_contrib);
        }
    }
}

void Discretization::correct_velocity(
    VelocityField& U,
    const SparseMatrix& Ap,
    const Vector& diag_inv_Ap,
    const ScalarField& rho,
    const VectorField& grad_p,
    const VectorField& dp_correction) const
{
    Index nc = U.size();
    for (Index ci = 0; ci < nc; ++ci) {
        Real rA = (ci < diag_inv_Ap.size()) ? diag_inv_Ap[ci] : 1.0;
        if (ci < static_cast<Index>(dp_correction.size())) {
            for (Index comp = 0; comp < 3; ++comp) {
                U[ci][comp] -= rA * dp_correction[ci][comp];
            }
        } else if (ci < static_cast<Index>(grad_p.size())) {
            for (Index comp = 0; comp < 3; ++comp) {
                U[ci][comp] -= rA * grad_p[ci][comp];
            }
        }
        (void)rho; (void)Ap;
    }
}

void Discretization::under_relax_field(
    ScalarField& field,
    const ScalarField& prev_field,
    Real relax_factor) const
{
    Index n = field.size();
    for (Index i = 0; i < n; ++i) {
        field[i] = prev_field[i] + relax_factor * (field[i] - prev_field[i]);
    }
}

void Discretization::under_relax_vector_field(
    VectorField& field,
    const VectorField& prev_field,
    Real relax_factor) const
{
    Index n = field.size();
    for (Index i = 0; i < n; ++i) {
        for (Index j = 0; j < 3; ++j) {
            field[i][j] = prev_field[i][j] + relax_factor * (field[i][j] - prev_field[i][j]);
        }
    }
}

Real Discretization::compute_convective_weight(
    Index face_id,
    Real m_dot,
    Real Pe) const
{
    (void)face_id;
    if (std::fabs(Pe) < 0.1) return 0.5;
    if (Pe > 10.0) return 1.0;
    if (Pe < -10.0) return 0.0;
    return 0.5 + 0.5 * std::tanh(0.1 * Pe);
}

void Discretization::compute_TVD_limiter(
    Real phi_r,
    ConvectionScheme scheme,
    Real& psi) const
{
    psi = 0.0;
    switch (scheme) {
        case ConvectionScheme::MinMod:
            psi = std::max(0.0, std::min(1.0, phi_r));
            break;
        case ConvectionScheme::MUSCL:
            psi = std::max(0.0, std::min(2.0, std::min(2.0 * phi_r, 0.5 * (1.0 + phi_r))));
            break;
        case ConvectionScheme::QUICK:
            if (phi_r > 0.0) {
                psi = std::min(8.0 / 3.0 * phi_r, 
                    std::min(0.75 + 0.25 * phi_r, 2.0));
            }
            break;
        default:
            psi = 1.0;
            break;
    }
}

} // namespace fvm
} // namespace hvdc
