#include "fem/BeamElementNL.hpp"
#include "common/MathUtils.hpp"
#include <cmath>

namespace hvdc {
namespace fem {

BeamElementNL::BeamElementNL() {
    type_ = ElementType::BeamNonlinear;
    math::mat_set_zero(T0_);
    math::vec_set_zero(dx0_);
}

BeamElementNL::BeamElementNL(Index id, Node* n1, Node* n2,
                              const Material* mat, const BeamSection* sec)
    : material_(mat), section_(sec)
{
    type_ = ElementType::BeamNonlinear;
    id_ = id;
    nodes_[0] = n1;
    nodes_[1] = n2;
    math::mat_set_zero(T0_);
    math::vec_set_zero(dx0_);
    initialize();
}

void BeamElementNL::initialize() {
    if (!nodes_[0] || !nodes_[1]) return;
    Vec3 x1 = nodes_[0]->coords();
    Vec3 x2 = nodes_[1]->coords();
    dx0_ = math::vec3_sub(x2, x1);
    L0_ = math::vec3_norm(dx0_);
    
    Vec3 e1_L = {1.0, 0.0, 0.0};
    Vec3 e1_G = math::vec3_normalize(dx0_);
    T0_ = math::rotation_matrix_from_vectors(e1_L, e1_G);
    
    if (material_) material_id_ = material_->id;
    if (section_) section_id_ = section_->id;
}

void BeamElementNL::stiffness_matrix(Mat12x12& K_out) const {
    Vec12 zero_disp{};
    tangent_stiffness_matrix(K_out, zero_disp, false);
}

void BeamElementNL::tangent_stiffness_matrix(
    Mat12x12& K_T_out, const Vec12& displacement,
    bool include_geometric) const
{
    math::mat_set_zero(K_T_out);
    if (!material_ || !section_) return;
    
    Mat3x3 T_LG;
    Vec3 dx_current;
    compute_transformation(T_LG, dx_current, displacement);
    
    Real L_current = math::vec3_norm(dx_current);
    if (L_current < EPS) return;
    
    Mat12x12 K_mat_local;
    compute_material_stiffness_local(K_mat_local);
    
    Vec3 disp1_global = {displacement[0], displacement[1], displacement[2]};
    Vec3 disp2_global = {displacement[6], displacement[7], displacement[8]};
    Vec3 rot1_global = {displacement[3], displacement[4], displacement[5]};
    Vec3 rot2_global = {displacement[9], displacement[10], displacement[11]};
    
    Vec3 disp1_local = math::mat3_mul_vec3(math::mat3_transpose(T_LG), disp1_global);
    Vec3 disp2_local = math::mat3_mul_vec3(math::mat3_transpose(T_LG), disp2_global);
    Vec3 rot1_local = math::mat3_mul_vec3(math::mat3_transpose(T_LG), rot1_global);
    Vec3 rot2_local = math::mat3_mul_vec3(math::mat3_transpose(T_LG), rot2_global);
    
    Vec12 disp_local = {
        disp1_local[0], disp1_local[1], disp1_local[2],
        rot1_local[0],  rot1_local[1],  rot1_local[2],
        disp2_local[0], disp2_local[1], disp2_local[2],
        rot2_local[0],  rot2_local[1],  rot2_local[2]
    };
    
    Real N, Vy, Vz, T, My, Mz;
    compute_axial_force_and_moments(disp_local, N, Vy, Vz, T, My, Mz);
    
    Mat12x12 K_total_local = K_mat_local;
    
    if (include_geometric) {
        Mat12x12 K_geo_local;
        compute_geometric_stiffness_local(K_geo_local, N, Vy, Vz, T, My, Mz);
        math::mat_add_inplace(K_total_local, K_geo_local);
    }
    
    transform_12x12(K_T_out, K_total_local, T_LG);
}

void BeamElementNL::mass_matrix(Mat12x12& M_out) const {
    math::mat_set_zero(M_out);
    if (!material_ || !section_) return;
    
    Real m = section_->total_mass_per_length(material_->rho);
    Real L = L0_;
    Real A = section_->A;
    
    Real rho = material_->rho;
    Real Ixx_r = (section_->Iz + section_->Iy) / A;
    Real Iyy_r = section_->Iy / A;
    Real Izz_r = section_->Iz / A;
    Real Jxx_r = section_->J / A;
    
    Mat12x12 M_lumped_local{};
    
    M_lumped_local[0][0]   = 0.5 * rho * A * L;
    M_lumped_local[1][1]   = 0.5 * rho * A * L;
    M_lumped_local[2][2]   = 0.5 * rho * A * L;
    M_lumped_local[3][3]   = 0.5 * rho * A * Jxx_r * L;
    M_lumped_local[4][4]   = 0.5 * rho * A * Iyy_r * L;
    M_lumped_local[5][5]   = 0.5 * rho * A * Izz_r * L;
    M_lumped_local[6][6]   = 0.5 * rho * A * L;
    M_lumped_local[7][7]   = 0.5 * rho * A * L;
    M_lumped_local[8][8]   = 0.5 * rho * A * L;
    M_lumped_local[9][9]   = 0.5 * rho * A * Jxx_r * L;
    M_lumped_local[10][10] = 0.5 * rho * A * Iyy_r * L;
    M_lumped_local[11][11] = 0.5 * rho * A * Izz_r * L;
    
    (void)m; (void)Ixx_r;
    
    Mat3x3 T = T0_;
    if (T[0][0] == 0.0 && T[1][1] == 0.0 && T[2][2] == 0.0) {
        T = math::mat3_identity();
    }
    transform_12x12(M_out, M_lumped_local, T);
}

Vec12 BeamElementNL::internal_forces(const Vec12& displacement) const {
    Vec12 f_int{};
    if (!material_ || !section_) return f_int;
    
    Mat3x3 T_LG;
    Vec3 dx_current;
    compute_transformation(T_LG, dx_current, displacement);
    
    Vec3 disp1_global = {displacement[0], displacement[1], displacement[2]};
    Vec3 disp2_global = {displacement[6], displacement[7], displacement[8]};
    Vec3 rot1_global = {displacement[3], displacement[4], displacement[5]};
    Vec3 rot2_global = {displacement[9], displacement[10], displacement[11]};
    
    Vec3 T_transposed = {};
    Mat3x3 Tt = math::mat3_transpose(T_LG);
    
    Vec3 disp1_local = math::mat3_mul_vec3(Tt, disp1_global);
    Vec3 disp2_local = math::mat3_mul_vec3(Tt, disp2_global);
    Vec3 rot1_local = math::mat3_mul_vec3(Tt, rot1_global);
    Vec3 rot2_local = math::mat3_mul_vec3(Tt, rot2_global);
    
    Vec12 disp_local = {
        disp1_local[0], disp1_local[1], disp1_local[2],
        rot1_local[0],  rot1_local[1],  rot1_local[2],
        disp2_local[0], disp2_local[1], disp2_local[2],
        rot2_local[0],  rot2_local[1],  rot2_local[2]
    };
    
    Real N, Vy, Vz, T_t, My, Mz;
    compute_axial_force_and_moments(disp_local, N, Vy, Vz, T_t, My, Mz);
    
    Vec12 f_local = {
        -N, -Vy, -Vz,
        -T_t, -My, -Mz,
        N, Vy, Vz,
        T_t, My, Mz
    };
    
    transform_12_vec(f_int, f_local, T_LG);
    
    return f_int;
}

Vec12 BeamElementNL::equivalent_nodal_loads() const {
    Vec12 loads{};
    if (!section_) return loads;
    
    Mat3x3 T = T0_;
    if (T[0][0] == 0.0 && T[1][1] == 0.0 && T[2][2] == 0.0) {
        T = math::mat3_identity();
    }
    
    Real L = L0_;
    Vec3 q = distributed_load_;
    Vec12 q_local = {
        q[0] * L / 2.0, q[1] * L / 2.0, q[2] * L / 2.0,
        0.0, -q[2] * L * L / 12.0, q[1] * L * L / 12.0,
        q[0] * L / 2.0, q[1] * L / 2.0, q[2] * L / 2.0,
        0.0, q[2] * L * L / 12.0, -q[1] * L * L / 12.0
    };
    
    transform_12_vec(loads, q_local, T);
    return loads;
}

Vec12 BeamElementNL::gravity_loads(Real g) const {
    Vec12 loads{};
    if (!material_ || !section_) return loads;
    
    Real mass_per_len = section_->total_mass_per_length(material_->rho);
    Vec3 g_global = {0.0, 0.0, -g * mass_per_len};
    
    Mat3x3 T = T0_;
    if (T[0][0] == 0.0 && T[1][1] == 0.0 && T[2][2] == 0.0) {
        T = math::mat3_identity();
    }
    
    Vec3 g_local = math::mat3_mul_vec3(math::mat3_transpose(T), g_global);
    
    Real L = L0_;
    Vec12 q_local = {
        g_local[0] * L / 2.0, g_local[1] * L / 2.0, g_local[2] * L / 2.0,
        0.0, -g_local[2] * L * L / 12.0, g_local[1] * L * L / 12.0,
        g_local[0] * L / 2.0, g_local[1] * L / 2.0, g_local[2] * L / 2.0,
        0.0, g_local[2] * L * L / 12.0, -g_local[1] * L * L / 12.0
    };
    
    transform_12_vec(loads, q_local, T);
    return loads;
}

void BeamElementNL::compute_material_stiffness_local(Mat12x12& K_mat_local) const {
    math::mat_set_zero(K_mat_local);
    
    Real E = material_->E;
    Real G = material_->G;
    Real L = L0_;
    Real A = section_->A;
    Real Iy = section_->Iy;
    Real Iz = section_->Iz;
    Real J = section_->J;
    
    Real EA_L = E * A / L;
    Real GJ_L = G * J / L;
    Real EIy_L3 = E * Iy / (L * L * L);
    Real EIz_L3 = E * Iz / (L * L * L);
    Real EIy_L2 = E * Iy / (L * L);
    Real EIz_L2 = E * Iz / (L * L);
    Real EIy_L = E * Iy / L;
    Real EIz_L = E * Iz / L;
    
    K_mat_local[0][0] = EA_L;
    K_mat_local[0][6] = -EA_L;
    K_mat_local[6][0] = -EA_L;
    K_mat_local[6][6] = EA_L;
    
    K_mat_local[1][1] = 12.0 * EIz_L3;
    K_mat_local[1][5] = 6.0 * EIz_L2;
    K_mat_local[1][7] = -12.0 * EIz_L3;
    K_mat_local[1][11] = 6.0 * EIz_L2;
    K_mat_local[5][1] = 6.0 * EIz_L2;
    K_mat_local[5][5] = 4.0 * EIz_L;
    K_mat_local[5][7] = -6.0 * EIz_L2;
    K_mat_local[5][11] = 2.0 * EIz_L;
    K_mat_local[7][1] = -12.0 * EIz_L3;
    K_mat_local[7][5] = -6.0 * EIz_L2;
    K_mat_local[7][7] = 12.0 * EIz_L3;
    K_mat_local[7][11] = -6.0 * EIz_L2;
    K_mat_local[11][1] = 6.0 * EIz_L2;
    K_mat_local[11][5] = 2.0 * EIz_L;
    K_mat_local[11][7] = -6.0 * EIz_L2;
    K_mat_local[11][11] = 4.0 * EIz_L;
    
    K_mat_local[2][2] = 12.0 * EIy_L3;
    K_mat_local[2][4] = -6.0 * EIy_L2;
    K_mat_local[2][8] = -12.0 * EIy_L3;
    K_mat_local[2][10] = -6.0 * EIy_L2;
    K_mat_local[4][2] = -6.0 * EIy_L2;
    K_mat_local[4][4] = 4.0 * EIy_L;
    K_mat_local[4][8] = 6.0 * EIy_L2;
    K_mat_local[4][10] = 2.0 * EIy_L;
    K_mat_local[8][2] = -12.0 * EIy_L3;
    K_mat_local[8][4] = 6.0 * EIy_L2;
    K_mat_local[8][8] = 12.0 * EIy_L3;
    K_mat_local[8][10] = 6.0 * EIy_L2;
    K_mat_local[10][2] = -6.0 * EIy_L2;
    K_mat_local[10][4] = 2.0 * EIy_L;
    K_mat_local[10][8] = 6.0 * EIy_L2;
    K_mat_local[10][10] = 4.0 * EIy_L;
    
    K_mat_local[3][3] = GJ_L;
    K_mat_local[3][9] = -GJ_L;
    K_mat_local[9][3] = -GJ_L;
    K_mat_local[9][9] = GJ_L;
}

void BeamElementNL::compute_geometric_stiffness_local(
    Mat12x12& K_geo_local, Real N, Real Vy, Real Vz,
    Real T, Real My, Real Mz) const
{
    math::mat_set_zero(K_geo_local);
    Real L = L0_;
    Real L_inv = 1.0 / L;
    
    Real N_L = N * L_inv;
    Real Vy_L = Vy * L_inv;
    Real Vz_L = Vz * L_inv;
    
    K_geo_local[1][1] += 6.0 / 5.0 * N_L;
    K_geo_local[1][5] += N / 10.0;
    K_geo_local[1][7] += -6.0 / 5.0 * N_L;
    K_geo_local[1][11] += N / 10.0;
    K_geo_local[5][1] += N / 10.0;
    K_geo_local[5][5] += 2.0 * L / 15.0 * N;
    K_geo_local[5][7] += -N / 10.0;
    K_geo_local[5][11] += -L / 30.0 * N;
    K_geo_local[7][1] += -6.0 / 5.0 * N_L;
    K_geo_local[7][5] += -N / 10.0;
    K_geo_local[7][7] += 6.0 / 5.0 * N_L;
    K_geo_local[7][11] += -N / 10.0;
    K_geo_local[11][1] += N / 10.0;
    K_geo_local[11][5] += -L / 30.0 * N;
    K_geo_local[11][7] += -N / 10.0;
    K_geo_local[11][11] += 2.0 * L / 15.0 * N;
    
    K_geo_local[2][2] += 6.0 / 5.0 * N_L;
    K_geo_local[2][4] += -N / 10.0;
    K_geo_local[2][8] += -6.0 / 5.0 * N_L;
    K_geo_local[2][10] += -N / 10.0;
    K_geo_local[4][2] += -N / 10.0;
    K_geo_local[4][4] += 2.0 * L / 15.0 * N;
    K_geo_local[4][8] += N / 10.0;
    K_geo_local[4][10] += -L / 30.0 * N;
    K_geo_local[8][2] += -6.0 / 5.0 * N_L;
    K_geo_local[8][4] += N / 10.0;
    K_geo_local[8][8] += 6.0 / 5.0 * N_L;
    K_geo_local[8][10] += N / 10.0;
    K_geo_local[10][2] += -N / 10.0;
    K_geo_local[10][4] += -L / 30.0 * N;
    K_geo_local[10][8] += N / 10.0;
    K_geo_local[10][10] += 2.0 * L / 15.0 * N;
    
    K_geo_local[0][1] += Vy_L / 2.0;
    K_geo_local[0][2] += Vz_L / 2.0;
    K_geo_local[0][7] += -Vy_L / 2.0;
    K_geo_local[0][8] += -Vz_L / 2.0;
    K_geo_local[6][1] += -Vy_L / 2.0;
    K_geo_local[6][2] += -Vz_L / 2.0;
    K_geo_local[6][7] += Vy_L / 2.0;
    K_geo_local[6][8] += Vz_L / 2.0;
    
    (void)T; (void)My; (void)Mz;
}

void BeamElementNL::compute_local_forces(
    const Vec12& disp_local, Vec12& f_int_local,
    Real& N_axial, Real& Vy_shear, Real& Vz_shear,
    Real& T_torsion, Real& My_bend, Real& Mz_bend) const
{
    compute_axial_force_and_moments(disp_local, N_axial, Vy_shear, Vz_shear,
                                     T_torsion, My_bend, Mz_bend);
    
    f_int_local = {
        -N_axial, -Vy_shear, -Vz_shear,
        -T_torsion, -My_bend, -Mz_bend,
        N_axial, Vy_shear, Vz_shear,
        T_torsion, My_bend, Mz_bend
    };
}

void BeamElementNL::compute_axial_force_and_moments(
    const Vec12& d_local,
    Real& N, Real& Vy, Real& Vz, Real& T, Real& My, Real& Mz) const
{
    Real E = material_->E;
    Real G = material_->G;
    Real L = L0_;
    Real A = section_->A;
    Real Iy = section_->Iy;
    Real Iz = section_->Iz;
    Real J = section_->J;
    
    Real du = d_local[6] - d_local[0];
    N = E * A / L * du;
    
    if (delta_T_ != 0.0 && material_->alpha != 0.0) {
        N -= E * A * material_->alpha * delta_T_;
    }
    
    Real theta_y1 = d_local[4];
    Real theta_y2 = d_local[10];
    Real theta_z1 = d_local[5];
    Real theta_z2 = d_local[11];
    
    Real L2 = L * L;
    Real L3 = L2 * L;
    
    Vy = 12.0 * E * Iz / L3 * (d_local[1] - d_local[7])
       + 6.0 * E * Iz / L2 * (theta_z1 + theta_z2);
    
    Vz = 12.0 * E * Iy / L3 * (d_local[8] - d_local[2])
       + 6.0 * E * Iy / L2 * (theta_y1 + theta_y2);
    
    T = G * J / L * (d_local[9] - d_local[3]);
    
    My = 6.0 * E * Iy / L2 * (d_local[2] - d_local[8])
       + 2.0 * E * Iy / L * (2.0 * theta_y1 + theta_y2);
    
    Mz = 6.0 * E * Iz / L2 * (d_local[7] - d_local[1])
       + 2.0 * E * Iz / L * (2.0 * theta_z2 + theta_z1);
}

void BeamElementNL::extract_node_rotations(const Vec12& disp_global,
                                            Mat3x3& R1, Mat3x3& R2) const {
    Vec3 rot1 = {disp_global[3], disp_global[4], disp_global[5]};
    Vec3 rot2 = {disp_global[9], disp_global[10], disp_global[11]};
    
    compute_current_rotation(T0_, rot1, R1);
    compute_current_rotation(T0_, rot2, R2);
}

void BeamElementNL::compute_current_rotation(
    const Mat3x3& R0, const Vec3& rot_incr, Mat3x3& R_current) const
{
    Real angle = math::vec3_norm(rot_incr);
    if (angle < EPS) {
        R_current = R0;
        return;
    }
    
    Mat3x3 R_incr = math::rotation_matrix_axis_angle(rot_incr, angle);
    R_current = math::mat3_mul(R0, R_incr);
}

} // namespace fem
} // namespace hvdc
