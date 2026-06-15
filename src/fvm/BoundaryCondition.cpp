#include "fvm/BoundaryCondition.hpp"
#include "fvm/Grid.hpp"
#include "common/MathUtils.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>

namespace hvdc {
namespace fvm {

BoundaryCondition::BoundaryCondition(const std::string& name, Index patch_id, BCType type)
    : name_(name), patch_id_(patch_id), type_(type)
{
}

void BoundaryCondition::apply_velocity(
    const Grid& grid,
    const FluidState& state,
    VelocityField& U,
    std::vector<Vec3>& face_velocity) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0 || fi >= static_cast<Index>(face_velocity.size())) continue;
        
        Vec3 n = face.normal();
        
        switch (type_) {
            case BCType::InletVelocity: {
                face_velocity[fi] = U_inlet_;
                if (ci < U.size()) U[ci] = U_inlet_;
                break;
            }
            case BCType::WallNoSlip:
            case BCType::WallMoving: {
                face_velocity[fi] = {0.0, 0.0, 0.0};
                if (ci < U.size()) U[ci] = {0.0, 0.0, 0.0};
                break;
            }
            case BCType::Symmetry: {
                if (ci < U.size()) {
                    Real Un = math::vec3_dot(U[ci], n);
                    Vec3 Ut = math::vec3_sub(U[ci], math::vec3_scale(n, Un));
                    face_velocity[fi] = Ut;
                    U[ci] = Ut;
                }
                break;
            }
            case BCType::OutletPressure:
            case BCType::OutletFlow:
            case BCType::ZeroGradient: {
                if (ci < U.size()) face_velocity[fi] = U[ci];
                break;
            }
            default:
                if (ci < U.size()) face_velocity[fi] = U[ci];
                break;
        }
        (void)state;
    }
}

void BoundaryCondition::apply_pressure(
    const Grid& grid,
    const FluidState& state,
    PressureField& p,
    std::vector<Real>& face_pressure) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0) continue;
        if (fi >= static_cast<Index>(face_pressure.size())) continue;
        
        switch (type_) {
            case BCType::OutletPressure:
            case BCType::InletPressure: {
                face_pressure[fi] = p_value_;
                break;
            }
            case BCType::InletTotalPressure: {
                face_pressure[fi] = p0_;
                if (ci < state.rho.size() && ci < state.U.size()) {
                    Real q = 0.5 * state.rho[ci] * state.U.speed_squared(ci);
                    p[ci] = p0_ - q;
                }
                break;
            }
            case BCType::InletVelocity:
            case BCType::WallNoSlip:
            case BCType::Symmetry:
            case BCType::ZeroGradient: {
                face_pressure[fi] = p[ci];
                break;
            }
            default:
                face_pressure[fi] = p[ci];
                break;
        }
    }
}

void BoundaryCondition::apply_turbulence(
    const Grid& grid,
    FluidState& state) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0) continue;
        
        switch (type_) {
            case BCType::InletVelocity:
            case BCType::InletPressure:
            case BCType::InletTotalPressure: {
                state.turb.k[ci] = std::max(k_inlet_, 1.0e-12);
                state.turb.epsilon[ci] = std::max(eps_inlet_, 1.0e-12);
                if (omega_inlet_ > EPS) {
                    state.turb.omega[ci] = std::max(omega_inlet_, 1.0e-12);
                }
                break;
            }
            case BCType::WallNoSlip: {
                state.turb.k[ci] = EPS_MEDIUM;
                state.turb.epsilon[ci] = EPS_MEDIUM;
                break;
            }
            default:
                break;
        }
    }
}

void BoundaryCondition::apply_temperature(
    const Grid& grid,
    ScalarField& T) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0) continue;
        
        if (type_ == BCType::InletVelocity || 
            type_ == BCType::InletPressure ||
            type_ == BCType::WallNoSlip) {
            T[ci] = T_value_;
        }
    }
}

void BoundaryCondition::update_matrix_rhs_for_bc(
    const Grid& grid,
    Index component,
    SparseMatrix& A,
    Vector& b,
    const VectorField& field) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0 || ci >= A.rows()) continue;
        
        Real bc_val = 0.0;
        bool apply_fixed = false;
        
        switch (type_) {
            case BCType::InletVelocity:
                bc_val = U_inlet_[component];
                apply_fixed = true;
                break;
            case BCType::WallNoSlip:
            case BCType::WallMoving:
                bc_val = 0.0;
                apply_fixed = true;
                break;
            case BCType::Symmetry: {
                Vec3 n = face.normal();
                bc_val = field[ci][component] * (1.0 - n[component] * n[component]);
                break;
            }
            default:
                break;
        }
        
        if (apply_fixed) {
            Real old_diag = A(ci, ci);
            Real factor = old_diag * 1.0e5;
            A.set(ci, ci, old_diag + factor);
            b.set(ci, b[ci] + factor * bc_val);
        }
    }
}

void BCPressureOutlet::apply_pressure(
    const Grid& grid,
    const FluidState& state,
    PressureField& p,
    std::vector<Real>& face_pressure) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        if (fi < static_cast<Index>(face_pressure.size())) {
            face_pressure[fi] = p_value_;
        }
    }
    (void)state; (void)p;
}

void BCInletVelocity::apply_velocity(
    const Grid& grid,
    const FluidState& state,
    VelocityField& U,
    std::vector<Vec3>& face_velocity) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0) continue;
        
        if (fi < static_cast<Index>(face_velocity.size())) {
            face_velocity[fi] = U_inlet_;
        }
        if (ci < U.size()) {
            U[ci] = U_inlet_;
        }
    }
    (void)state;
}

void BCInletVelocity::apply_turbulence(
    const Grid& grid,
    FluidState& state) const
{
    BoundaryCondition::apply_turbulence(grid, state);
}

void BCInletVelocity::compute_turbulence_from_velocity() {
    Real I = 0.01;
    Real U_mag = math::vec3_norm(U_inlet_);
    k_inlet_ = 1.5 * (U_mag * I) * (U_mag * I);
    Real L = 1.0;
    Real C_mu = 0.09;
    eps_inlet_ = std::pow(C_mu, 0.75) * std::pow(k_inlet_, 1.5) / L;
    omega_inlet_ = eps_inlet_ / (C_mu * k_inlet_);
}

void BCWallNoSlip::apply_velocity(
    const Grid& grid,
    const FluidState& state,
    VelocityField& U,
    std::vector<Vec3>& face_velocity) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0) continue;
        
        if (fi < static_cast<Index>(face_velocity.size())) {
            face_velocity[fi] = wall_velocity_;
        }
        if (ci < U.size()) {
            U[ci] = wall_velocity_;
        }
    }
    (void)state;
}

Real BCWallNoSlip::compute_y_plus(Real u_tau, Real y, Real nu) const {
    return u_tau * y / std::max(nu, EPS);
}

Real BCWallNoSlip::compute_wall_shear_stress(
    Real U_tangential, Real y, Real nu,
    Real& u_tau, Real& y_plus) const
{
    if (U_tangential < EPS || y < EPS) {
        u_tau = 0.0;
        y_plus = 0.0;
        return 0.0;
    }
    
    Real C_mu = 0.09;
    Real kappa = 0.41;
    Real E = 9.8;
    
    u_tau = std::sqrt(C_mu * 0.01);
    y_plus = compute_y_plus(u_tau, y, nu);
    
    if (y_plus < 11.0) {
        u_tau = std::sqrt(nu * U_tangential / y);
    } else {
        u_tau = kappa * U_tangential / std::log(E * y_plus);
    }
    
    y_plus = compute_y_plus(u_tau, y, nu);
    return u_tau * u_tau;
}

void BCSymmetry::apply_velocity(
    const Grid& grid,
    const FluidState& state,
    VelocityField& U,
    std::vector<Vec3>& face_velocity) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0 || ci >= U.size()) continue;
        
        Vec3 n = face.normal();
        Real Un = math::vec3_dot(U[ci], n);
        Vec3 Ut = math::vec3_sub(U[ci], math::vec3_scale(n, Un));
        if (fi < static_cast<Index>(face_velocity.size())) {
            face_velocity[fi] = Ut;
        }
        U[ci] = Ut;
    }
    (void)state;
}

void BCFarField::apply_velocity(
    const Grid& grid,
    const FluidState& state,
    VelocityField& U,
    std::vector<Vec3>& face_velocity) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0) continue;
        
        Vec3 n = face.normal();
        Real flux = 0.0;
        if (ci < state.U.size()) {
            flux = math::vec3_dot(state.U[ci], n);
        }
        
        Vec3 U_set = (flux >= 0) ? (ci < U.size() ? U[ci] : U_inlet_) : U_inlet_;
        if (fi < static_cast<Index>(face_velocity.size())) {
            face_velocity[fi] = U_set;
        }
        if (ci < U.size()) {
            U[ci] = U_set;
        }
    }
}

void BCFarField::apply_pressure(
    const Grid& grid,
    const FluidState& state,
    PressureField& p,
    std::vector<Real>& face_pressure) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0 || fi >= static_cast<Index>(face_pressure.size())) continue;
        
        Vec3 n = face.normal();
        Real flux = (ci < state.U.size()) ? math::vec3_dot(state.U[ci], n) : 0.0;
        if (flux >= 0 && ci < p.size()) {
            face_pressure[fi] = p[ci];
        } else {
            face_pressure[fi] = p_value_;
        }
    }
}

void BCInterfaceFSI::apply_velocity(
    const Grid& grid,
    const FluidState& state,
    VelocityField& U,
    std::vector<Vec3>& face_velocity) const
{
    auto range = grid.boundary_patch_range(patch_id_);
    if (range.first < 0) return;
    
    Index face_offset = range.first;
    
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& face = grid.face(fi);
        Index ci = (face.owner() >= 0) ? face.owner() : face.neighbor();
        if (ci < 0) continue;
        
        Vec3 wall_v = {0, 0, 0};
        Index local_face_idx = fi - face_offset;
        if (local_face_idx < static_cast<Index>(interface_velocity_.size())) {
            wall_v = interface_velocity_[local_face_idx];
        }
        
        if (fi < static_cast<Index>(face_velocity.size())) {
            face_velocity[fi] = wall_v;
        }
        if (ci < U.size() && is_moving_wall_) {
            U[ci] = wall_v;
        }
    }
    (void)state;
}

Index BoundaryConditionManager::add_bc(std::unique_ptr<BoundaryCondition> bc) {
    Index id = static_cast<Index>(bcs_.size());
    name_map_[bc->name()] = id;
    patch_map_[bc->patch_id()] = id;
    bcs_.push_back(std::move(bc));
    return id;
}

BoundaryCondition* BoundaryConditionManager::bc_by_name(const std::string& name) {
    auto it = name_map_.find(name);
    return (it != name_map_.end()) ? bcs_[it->second].get() : nullptr;
}

BoundaryCondition* BoundaryConditionManager::bc_by_patch(Index patch_id) {
    auto it = patch_map_.find(patch_id);
    return (it != patch_map_.end()) ? bcs_[it->second].get() : nullptr;
}

void BoundaryConditionManager::apply_all_velocity_bc(
    const Grid& grid,
    const FluidState& state,
    VelocityField& U,
    std::vector<Vec3>& face_velocity) const
{
    for (const auto& bc : bcs_) {
        bc->apply_velocity(grid, state, U, face_velocity);
    }
}

void BoundaryConditionManager::apply_all_pressure_bc(
    const Grid& grid,
    const FluidState& state,
    PressureField& p,
    std::vector<Real>& face_pressure) const
{
    for (const auto& bc : bcs_) {
        bc->apply_pressure(grid, state, p, face_pressure);
    }
}

void BoundaryConditionManager::apply_all_turbulence_bc(
    const Grid& grid,
    FluidState& state) const
{
    for (const auto& bc : bcs_) {
        bc->apply_turbulence(grid, state);
    }
}

void BoundaryConditionManager::apply_all_temperature_bc(
    const Grid& grid,
    ScalarField& T) const
{
    for (const auto& bc : bcs_) {
        bc->apply_temperature(grid, T);
    }
}

void BoundaryConditionManager::update_linear_system_bc(
    const Grid& grid,
    const std::string& field_name,
    Index component,
    SparseMatrix& A,
    Vector& b,
    const VectorField& field) const
{
    (void)field_name;
    for (const auto& bc : bcs_) {
        bc->update_matrix_rhs_for_bc(grid, component, A, b, field);
    }
}

} // namespace fvm
} // namespace hvdc
