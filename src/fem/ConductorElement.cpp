#include "fem/ConductorElement.hpp"
#include "common/MathUtils.hpp"
#include <cmath>
#include <algorithm>

namespace hvdc {
namespace fem {

ConductorElement::ConductorElement()
    : T0_(0.0), sag_(0.0), DAF_(1.0)
{
    type_ = ElementType::ConductorCable;
}

ConductorElement::ConductorElement(Index id, Node* n1, Node* n2,
                                    const Material* mat, const BeamSection* sec,
                                    Real initial_tension, Real span_sag)
    : material_(mat), section_(sec), T0_(initial_tension), sag_(span_sag), DAF_(1.0)
{
    type_ = ElementType::ConductorCable;
    id_ = id;
    nodes_[0] = n1;
    nodes_[1] = n2;
    initialize();
}

void ConductorElement::initialize() {
    if (!nodes_[0] || !nodes_[1]) return;
    
    Vec3 x1 = nodes_[0]->coords();
    Vec3 x2 = nodes_[1]->coords();
    Vec3 dx = math::vec3_sub(x2, x1);
    L0_ = math::vec3_norm(dx);
    
    if (material_) material_id_ = material_->id;
    if (section_) section_id_ = section_->id;
    
    if (T0_ > 0.0 && section_ && material_ && sag_ <= 0.0) {
        if (L0_ > 0.0 && section_->A > 0.0) {
            Real w = material_->rho * section_->A * GRAVITY;
            Real H = T0_;
            Real L_span = L0_;
            sag_ = w * L_span * L_span / (8.0 * H);
        }
    }
}

void ConductorElement::stiffness_matrix(Mat12x12& K_out) const {
    Vec12 zero{};
    tangent_stiffness_matrix(K_out, zero, false);
}

void ConductorElement::tangent_stiffness_matrix(
    Mat12x12& K_T_out, const Vec12& displacement,
    bool include_geometric) const
{
    math::mat_set_zero(K_T_out);
    if (!material_ || !section_) return;
    
    Vec3 x1 = nodes_[0]->displaced_coords({displacement[0], displacement[1], displacement[2]});
    Vec3 x2 = nodes_[1]->displaced_coords({displacement[6], displacement[7], displacement[8]});
    Vec3 dx = math::vec3_sub(x2, x1);
    Real L_chord = math::vec3_norm(dx);
    
    if (L_chord < EPS) return;
    
    Vec3 e = math::vec3_scale(dx, 1.0 / L_chord);
    
    Real L_cat = catenary_length();
    Real delta_L = (L_chord - L0_) / L0_;
    Real eps_total = (L_cat - L0_) / L0_;
    if (L0_ > EPS) {
        eps_total += (delta_L > 0) ? delta_L * 0.5 : delta_L;
    }
    
    if (delta_T_ != 0.0) {
        eps_total -= material_->alpha * delta_T_;
    }
    
    Real T_current;
    if (section_->A > EPS) {
        Real sigma = material_->E * eps_total;
        T_current = sigma * section_->A + T0_;
    } else {
        T_current = T0_;
    }
    T_current = std::max(T_current, 0.0);
    
    Real E_eff = compute_effective_modulus(T_current);
    Real EA_eff = E_eff * section_->A;
    Real L_eff = std::max(L_chord, L_cat);
    Real EA_L_eff = EA_eff / L_eff;
    Real T_over_L = include_geometric ? (T_current / L_eff) : 0.0;
    
    compute_truss_equivalent_stiffness(K_T_out, e, EA_L_eff, T_over_L, include_geometric);
}

void ConductorElement::mass_matrix(Mat12x12& M_out) const {
    math::mat_set_zero(M_out);
    if (!material_ || !section_) return;
    
    Real L_cat = catenary_length();
    Real m_total = material_->rho * section_->A * L_cat;
    if (section_->A_ice > 0.0) {
        m_total += 917.0 * section_->A_ice * L_cat;
    }
    
    for (Index bi = 0; bi < 2; ++bi) {
        for (Index i = 0; i < 3; ++i) {
            M_out[bi * 6 + i][bi * 6 + i] = 0.5 * m_total;
        }
    }
}

Vec12 ConductorElement::internal_forces(const Vec12& displacement) const {
    Vec12 f_int{};
    if (!material_ || !section_) return f_int;
    
    Vec3 x1 = nodes_[0]->displaced_coords({displacement[0], displacement[1], displacement[2]});
    Vec3 x2 = nodes_[1]->displaced_coords({displacement[6], displacement[7], displacement[8]});
    Vec3 dx = math::vec3_sub(x2, x1);
    Real L_chord = math::vec3_norm(dx);
    
    if (L_chord < EPS) return f_int;
    
    Vec3 e = math::vec3_scale(dx, 1.0 / L_chord);
    Real T = current_tension(displacement);
    
    Vec3 f1 = math::vec3_scale(e, -T);
    Vec3 f2 = math::vec3_scale(e, T);
    
    for (Index i = 0; i < 3; ++i) {
        f_int[i]     = f1[i];
        f_int[6 + i] = f2[i];
    }
    
    return f_int;
}

Vec12 ConductorElement::equivalent_nodal_loads() const {
    Vec12 loads{};
    if (!section_ || !material_) return loads;
    
    Real L_cat = catenary_length();
    Real w = material_->rho * section_->A * GRAVITY;
    if (section_->A_ice > 0.0) {
        w += 917.0 * section_->A_ice * GRAVITY;
    }
    Real W_total = w * L_cat;
    
    if (sag_ > 0.0 && L0_ > 0.0) {
        Real L = L0_;
        Real h_sag = sag_;
        Real V = w * L / 2.0;
        Real H = w * L * L / (8.0 * h_sag);
        
        Vec3 x1 = nodes_[0]->coords();
        Vec3 x2 = nodes_[1]->coords();
        Vec3 dx = math::vec3_sub(x2, x1);
        Vec3 e_horiz = math::vec3_normalize(Vec3{dx[0], dx[1], 0.0});
        if (math::vec3_norm2(e_horiz) < EPS) {
            e_horiz = {1.0, 0.0, 0.0};
        }
        
        Vec3 e_vert = {0.0, 0.0, -1.0};
        
        Vec3 F1 = math::vec3_add(math::vec3_scale(e_horiz, -H),
                                  math::vec3_scale(e_vert, -V));
        Vec3 F2 = math::vec3_add(math::vec3_scale(e_horiz, H),
                                  math::vec3_scale(e_vert, -V));
        
        for (Index i = 0; i < 3; ++i) {
            loads[i]     = F1[i];
            loads[6 + i] = F2[i];
        }
    } else {
        for (Index bi = 0; bi < 2; ++bi) {
            loads[bi * 6 + 2] = -0.5 * W_total;
        }
    }
    
    return loads;
}

Vec12 ConductorElement::gravity_loads(Real g) const {
    Vec12 loads{};
    if (!material_ || !section_) return loads;
    
    Real L_cat = catenary_length();
    Real w = material_->rho * section_->A * g;
    if (section_->A_ice > 0.0) {
        w += 917.0 * section_->A_ice * g;
    }
    Real W_total = w * L_cat;
    
    for (Index bi = 0; bi < 2; ++bi) {
        loads[bi * 6 + 2] = -0.5 * W_total;
    }
    
    return loads;
}

Real ConductorElement::current_tension(const Vec12& displacement) const {
    Vec3 x1 = nodes_[0]->displaced_coords({displacement[0], displacement[1], displacement[2]});
    Vec3 x2 = nodes_[1]->displaced_coords({displacement[6], displacement[7], displacement[8]});
    Vec3 dx = math::vec3_sub(x2, x1);
    Real L_chord = math::vec3_norm(dx);
    
    if (L_chord < EPS || L0_ < EPS || !section_ || !material_) return T0_;
    
    Real L_cat = catenary_length();
    Real epsilon = (L_cat - L0_) / L0_;
    Real delta_eps = (L_chord - L0_) / L0_;
    epsilon += std::max(0.0, delta_eps * 0.5);
    
    if (delta_T_ != 0.0) {
        epsilon -= material_->alpha * delta_T_;
    }
    
    Real T = T0_;
    if (section_->A > EPS) {
        T = T0_ + material_->E * section_->A * epsilon;
    }
    return std::max(T, 0.0);
}

Real ConductorElement::catenary_length() const {
    if (sag_ <= 0.0 || L0_ <= 0.0) return L0_;
    
    Real S = L0_;
    Real h = sag_;
    Real L_cat = S * (1.0 + 8.0 / 3.0 * (h * h) / (S * S)
                        - 32.0 / 5.0 * std::pow(h / S, 4)
                        + 256.0 / 7.0 * std::pow(h / S, 6));
    return L_cat;
}

void ConductorElement::catenary_shape(Index n_points,
                                       std::vector<Vec3>& points,
                                       std::vector<Real>& tensions) const
{
    n_points = std::max<Index>(2, n_points);
    points.resize(n_points);
    tensions.resize(n_points);
    
    Vec3 x1 = nodes_[0]->coords();
    Vec3 x2 = nodes_[1]->coords();
    Vec3 dx = math::vec3_sub(x2, x1);
    
    Vec3 e_horiz = math::vec3_normalize(Vec3{dx[0], dx[1], 0.0});
    if (math::vec3_norm2(e_horiz) < EPS) {
        e_horiz = {1.0, 0.0, 0.0};
    }
    Vec3 e_vert = {0.0, 0.0, 1.0};
    
    Real S = L0_;
    Real h = sag_;
    if (h <= 0.0 && T0_ > 0.0 && section_ && material_) {
        Real w = material_->rho * section_->A * GRAVITY;
        if (section_->A_ice > 0.0) {
            w += 917.0 * section_->A_ice * GRAVITY;
        }
        h = w * S * S / (8.0 * T0_);
    }
    
    for (Index i = 0; i < n_points; ++i) {
        Real t = static_cast<Real>(i) / static_cast<Real>(n_points - 1);
        Real x_t = t * S;
        Real y_t = 4.0 * h * t * (1.0 - t);
        
        Vec3 horiz_comp = math::vec3_scale(e_horiz, x_t);
        Vec3 vert_comp = math::vec3_scale(e_vert, -y_t);
        points[i] = math::vec3_add(x1, math::vec3_add(horiz_comp, vert_comp));
        
        Real slope = 4.0 * h / S * (1.0 - 2.0 * t);
        Real angle = std::atan(slope);
        tensions[i] = T0_ / std::cos(angle);
    }
}

void ConductorElement::compute_truss_equivalent_stiffness(
    Mat12x12& K, const Vec3& e,
    Real EA_L_eff, Real T_eff_over_L, bool use_geo) const
{
    Mat3x3 K3{};
    for (Index i = 0; i < 3; ++i) {
        for (Index j = 0; j < 3; ++j) {
            K3[i][j] = EA_L_eff * e[i] * e[j];
            if (use_geo) {
                K3[i][j] += T_eff_over_L * ((i == j ? 1.0 : 0.0) - e[i] * e[j]);
            }
        }
    }
    
    for (Index i = 0; i < 3; ++i) {
        for (Index j = 0; j < 3; ++j) {
            K[i][j]             =  K3[i][j];
            K[i][j + 6]         = -K3[i][j];
            K[i + 6][j]         = -K3[i][j];
            K[i + 6][j + 6]     =  K3[i][j];
        }
    }
}

Real ConductorElement::compute_effective_modulus(Real current_tension) const {
    if (!material_ || !section_) return 0.0;
    
    Real E0 = material_->E;
    if (T0_ <= 0.0 || section_->A <= 0.0) return E0;
    
    Real sigma0 = T0_ / section_->A;
    Real sigma = current_tension / section_->A;
    Real sigma_avg = 0.5 * (sigma0 + sigma);
    Real E_eff = E0;
    
    if (sigma_avg > material_->sigma_y * 0.5) {
        Real ratio = std::min(1.0, (sigma_avg - material_->sigma_y * 0.5) / 
                              (material_->sigma_y * 0.5));
        E_eff = E0 * (1.0 - 0.2 * ratio);
    }
    
    return E_eff;
}

Real ConductorElement::compute_Lame_parameter(Real stress_axial) const {
    if (!material_) return 0.0;
    Real E = material_->E;
    Real nu = material_->nu;
    return E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
}

} // namespace fem
} // namespace hvdc
