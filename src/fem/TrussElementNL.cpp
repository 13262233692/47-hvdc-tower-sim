#include "fem/TrussElementNL.hpp"
#include "common/MathUtils.hpp"
#include <cmath>

namespace hvdc {
namespace fem {

TrussElementNL::TrussElementNL() {
    type_ = ElementType::TrussNonlinear;
}

TrussElementNL::TrussElementNL(Index id, Node* n1, Node* n2,
                                const Material* mat, const TrussSection* sec)
    : material_(mat), section_(sec)
{
    type_ = ElementType::TrussNonlinear;
    id_ = id;
    nodes_[0] = n1;
    nodes_[1] = n2;
    initialize();
}

void TrussElementNL::initialize() {
    if (!nodes_[0] || !nodes_[1]) return;
    Vec3 x1 = nodes_[0]->coords();
    Vec3 x2 = nodes_[1]->coords();
    Vec3 dx = math::vec3_sub(x2, x1);
    L0_ = math::vec3_norm(dx);
    if (material_) material_id_ = material_->id;
    if (section_) section_id_ = section_->id;
}

void TrussElementNL::stiffness_matrix(Mat12x12& K_out) const {
    Vec12 zero{};
    tangent_stiffness_matrix(K_out, zero, false);
}

void TrussElementNL::tangent_stiffness_matrix(
    Mat12x12& K_T_out, const Vec12& displacement,
    bool include_geometric) const
{
    math::mat_set_zero(K_T_out);
    if (!material_ || !section_) return;
    
    Mat3x3 T;
    Vec3 dir;
    Real L_cur, EA_L, N_force;
    compute_stiffness_components(displacement, T, dir, L_cur, EA_L, N_force);
    
    if (L_cur < EPS) return;
    
    Real N_L = include_geometric ? (N_force / L_cur) : 0.0;
    build_stiffness_from_dir(K_T_out, dir, L_cur, EA_L, N_L, include_geometric);
}

void TrussElementNL::mass_matrix(Mat12x12& M_out) const {
    math::mat_set_zero(M_out);
    if (!material_ || !section_) return;
    
    Real m = material_->rho * section_->total_A() * L0_;
    
    for (Index bi = 0; bi < 2; ++bi) {
        for (Index i = 0; i < 3; ++i) {
            M_out[bi * 6 + i][bi * 6 + i] = 0.5 * m;
        }
    }
}

Vec12 TrussElementNL::internal_forces(const Vec12& displacement) const {
    Vec12 f_int{};
    if (!material_ || !section_) return f_int;
    
    Mat3x3 T;
    Vec3 dir;
    Real L_cur, EA_L, N;
    compute_stiffness_components(displacement, T, dir, L_cur, EA_L, N);
    
    if (L_cur < EPS) return f_int;
    
    Vec3 f1 = math::vec3_scale(dir, -N);
    Vec3 f2 = math::vec3_scale(dir, N);
    
    for (Index i = 0; i < 3; ++i) {
        f_int[i]     = f1[i];
        f_int[6 + i] = f2[i];
    }
    
    return f_int;
}

Vec12 TrussElementNL::equivalent_nodal_loads() const {
    return Vec12{};
}

Vec12 TrussElementNL::gravity_loads(Real g) const {
    Vec12 loads{};
    if (!material_ || !section_) return loads;
    
    Real total_A = section_->total_A();
    Real m_total = material_->rho * total_A * L0_;
    if (section_->thickness_ice > 0.0 && section_->A_ice > 0.0) {
        m_total += 917.0 * section_->A_ice * L0_;
    }
    
    for (Index bi = 0; bi < 2; ++bi) {
        loads[bi * 6 + 2] = -0.5 * m_total * g;
    }
    
    return loads;
}

Real TrussElementNL::axial_force(const Vec12& displacement) const {
    Mat3x3 T;
    Vec3 dir;
    Real L_cur, EA_L, N;
    compute_stiffness_components(displacement, T, dir, L_cur, EA_L, N);
    return N;
}

Real TrussElementNL::axial_stress(const Vec12& displacement) const {
    Real N = axial_force(displacement);
    if (section_) {
        return N / section_->total_A();
    }
    return 0.0;
}

Real TrussElementNL::axial_strain(const Vec12& displacement) const {
    Vec3 x1 = nodes_[0]->displaced_coords({displacement[0], displacement[1], displacement[2]});
    Vec3 x2 = nodes_[1]->displaced_coords({displacement[6], displacement[7], displacement[8]});
    Real L_cur = math::vec3_norm(math::vec3_sub(x2, x1));
    return (L_cur - L0_) / L0_;
}

void TrussElementNL::compute_axial_response(const Vec12& disp,
                                             Real& strain, Real& stress, Real& force) const {
    strain = axial_strain(disp);
    stress = material_->E * strain;
    if (delta_T_ != 0.0) {
        stress -= material_->E * material_->alpha * delta_T_;
    }
    force = stress * section_->total_A();
}

void TrussElementNL::compute_stiffness_components(
    const Vec12& displacement,
    Mat3x3& T, Vec3& dir, Real& L_cur,
    Real& EA_L, Real& N_force) const
{
    Vec3 x1 = nodes_[0]->displaced_coords({displacement[0], displacement[1], displacement[2]});
    Vec3 x2 = nodes_[1]->displaced_coords({displacement[6], displacement[7], displacement[8]});
    Vec3 dx = math::vec3_sub(x2, x1);
    L_cur = math::vec3_norm(dx);
    
    if (L_cur < EPS) {
        dir = {1.0, 0.0, 0.0};
        T = math::mat3_identity();
        EA_L = 0.0;
        N_force = 0.0;
        return;
    }
    
    dir = math::vec3_scale(dx, 1.0 / L_cur);
    
    Vec3 e1_L = {1.0, 0.0, 0.0};
    T = math::rotation_matrix_from_vectors(e1_L, dir);
    
    Real epsilon = (L_cur - L0_) / L0_;
    if (delta_T_ != 0.0) {
        epsilon -= material_->alpha * delta_T_;
    }
    Real sigma = material_->E * epsilon;
    N_force = sigma * section_->total_A();
    EA_L = material_->E * section_->total_A() / L_cur;
}

void TrussElementNL::build_stiffness_from_dir(
    Mat12x12& K, const Vec3& e, Real L,
    Real EA_L, Real N_L, bool use_geo) const
{
    Mat3x3 K_mat_block{};
    for (Index i = 0; i < 3; ++i) {
        for (Index j = 0; j < 3; ++j) {
            K_mat_block[i][j] = EA_L * e[i] * e[j];
        }
    }
    
    Mat3x3 K_geo_block{};
    if (use_geo) {
        for (Index i = 0; i < 3; ++i) {
            for (Index j = 0; j < 3; ++j) {
                K_geo_block[i][j] = N_L * (i == j ? 1.0 : 0.0) - N_L * e[i] * e[j];
            }
        }
    }
    
    Mat3x3 K3 = K_mat_block;
    if (use_geo) {
        math::mat3_add_inplace(K3, K_geo_block);
    }
    
    for (Index i = 0; i < 3; ++i) {
        for (Index j = 0; j < 3; ++j) {
            K[i][j]             =  K3[i][j];
            K[i][j + 6]         = -K3[i][j];
            K[i + 6][j]         = -K3[i][j];
            K[i + 6][j + 6]     =  K3[i][j];
        }
    }
    (void)L;
}

} // namespace fem
} // namespace hvdc
