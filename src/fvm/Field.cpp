#include "fvm/Field.hpp"
#include "fvm/Grid.hpp"
#include "common/MathUtils.hpp"
#include <algorithm>
#include <cmath>

namespace hvdc {
namespace fvm {

template <typename T>
Field<T>::Field(Index size, FieldLocation loc, const std::string& name)
    : data_(size), name_(name), location_(loc) {}

template <typename T>
Field<T>::Field(Index size, T value, FieldLocation loc, const std::string& name)
    : data_(size, value), name_(name), location_(loc) {}

template class Field<Real>;
template class Field<Vec3>;
template class Field<Mat3x3>;

ScalarField VelocityField::component(Index i) const {
    ScalarField comp(size(), 0.0, FieldLocation::CellCenter, name_ + "_" + std::to_string(i));
    for (Index ci = 0; ci < size(); ++ci) {
        comp[ci] = data_[ci][i];
    }
    return comp;
}

void VelocityField::set_component(Index i, const ScalarField& field) {
    Index n = std::min(size(), field.size());
    for (Index ci = 0; ci < n; ++ci) {
        data_[ci][i] = field[ci];
    }
}

Real VelocityField::speed(Index cell_id) const {
    return std::sqrt(speed_squared(cell_id));
}

Real VelocityField::speed_squared(Index cell_id) const {
    return math::vec3_norm2(data_[cell_id]);
}

ScalarField VelocityField::speed_field() const {
    ScalarField sp(size(), 0.0, FieldLocation::CellCenter, "speed");
    for (Index ci = 0; ci < size(); ++ci) {
        sp[ci] = speed(ci);
    }
    return sp;
}

ScalarField VelocityField::kinetic_energy() const {
    ScalarField ke(size(), 0.0, FieldLocation::CellCenter, "k_kinetic");
    for (Index ci = 0; ci < size(); ++ci) {
        ke[ci] = 0.5 * speed_squared(ci);
    }
    return ke;
}

Vec3 VelocityField::face_interpolated(Index face_id, Real weight,
                                       const Vec3& owner_val, const Vec3& neighbor_val) const
{
    (void)face_id;
    Vec3 result;
    for (Index i = 0; i < 3; ++i) {
        result[i] = weight * owner_val[i] + (1.0 - weight) * neighbor_val[i];
    }
    return result;
}

Vec3 PressureField::pressure_gradient(Index cell_id, const Grid& grid) const {
    Vec3 grad = {0, 0, 0};
    const Cell& cell = grid.cell(cell_id);
    Real V = cell.volume();
    if (V < EPS) return grad;
    
    for (Index fi = 0; fi < cell.num_faces(); ++fi) {
        Index face_id = cell.face_id(fi);
        if (face_id < 0) continue;
        const Face& face = grid.face(face_id);
        
        Vec3 n = face.normal();
        Real A = face.area();
        Index owner = face.owner();
        Index neighbor = face.neighbor();
        
        Real p_face;
        if (face.is_on_boundary()) {
            if (owner == cell_id) {
                p_face = data_[owner];
            } else {
                p_face = (neighbor >= 0) ? data_[neighbor] : data_[owner];
            }
        } else {
            Real w = face.weight();
            p_face = w * data_[owner] + (1.0 - w) * data_[neighbor];
        }
        
        math::vec3_add_inplace(grad, math::vec3_scale(n, p_face * A));
    }
    
    math::vec3_scale_inplace(grad, 1.0 / V);
    return grad;
}

TurbulenceField::TurbulenceField(Index n_cells)
    : k(n_cells, 0.0, FieldLocation::CellCenter, "k"),
      epsilon(n_cells, 0.0, FieldLocation::CellCenter, "epsilon"),
      omega(n_cells, 0.0, FieldLocation::CellCenter, "omega"),
      nu_t(n_cells, 0.0, FieldLocation::CellCenter, "nu_t"),
      mut(n_cells, 0.0, FieldLocation::CellCenter, "mut")
{}

Real TurbulenceField::Reynolds_stress_xx(Index cell_id, Real Ux_grad) const {
    return -2.0 / 3.0 * k[cell_id] + 2.0 * nu_t[cell_id] * Ux_grad;
}

void TurbulenceField::compute_nut_from_k_epsilon(Real C_mu) {
    Index n = k.size();
    for (Index ci = 0; ci < n; ++ci) {
        Real eps = std::max(epsilon[ci], 1.0e-12);
        nu_t[ci] = C_mu * k[ci] * k[ci] / eps;
    }
}

void TurbulenceField::compute_nut_from_k_omega() {
    Index n = k.size();
    for (Index ci = 0; ci < n; ++ci) {
        Real om = std::max(omega[ci], 1.0e-12);
        nu_t[ci] = k[ci] / om;
    }
}

void FluidState::resize(Index n) {
    U.resize(n);
    p.resize(n);
    turb = TurbulenceField(n);
    T.resize(n);
    rho.resize(n);
    mu.resize(n);
    nu.resize(n);
}

void FluidState::initialize_defaults(Index n_cells) {
    resize(n_cells);
    U.fill({0.0, 0.0, 0.0});
    p.fill(p_ref);
    T.fill(T_ref);
    rho.fill(rho_ref);
    mu.fill(mu_ref);
    Real nu_ref = mu_ref / rho_ref;
    nu.fill(nu_ref);
    
    Real I = 0.01;
    Real k_inlet = 1.5 * (20.0 * I) * (20.0 * I);
    Real L = 1.0;
    Real C_mu = 0.09;
    Real eps_inlet = std::pow(C_mu, 0.75) * std::pow(k_inlet, 1.5) / L;
    
    turb.k.fill(k_inlet);
    turb.epsilon.fill(eps_inlet);
    turb.omega.fill(eps_inlet / (C_mu * k_inlet));
    turb.compute_nut_from_k_epsilon();
}

void FluidState::update_thermophysical_properties() {
    Index n = T.size();
    for (Index ci = 0; ci < n; ++ci) {
        rho[ci] = p_ref / (287.058 * std::max(T[ci], 200.0));
        Real T_ratio = T[ci] / T_ref;
        mu[ci] = mu_ref * std::pow(T_ratio, 1.5) * (T_ref + 110.4) / (T[ci] + 110.4);
        if (rho[ci] > EPS) {
            nu[ci] = mu[ci] / rho[ci];
        }
    }
}

} // namespace fvm
} // namespace hvdc
