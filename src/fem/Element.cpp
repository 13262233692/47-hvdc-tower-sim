#include "fem/Element.hpp"
#include "common/MathUtils.hpp"

namespace hvdc {
namespace fem {

Real Element2N::length_current(const Vec12& disp) const {
    Vec3 x1 = nodes_[0]->displaced_coords({disp[0], disp[1], disp[2]});
    Vec3 x2 = nodes_[1]->displaced_coords({disp[6], disp[7], disp[8]});
    return math::vec3_norm(math::vec3_sub(x2, x1));
}

void Element2N::compute_transformation(Mat3x3& T_L_to_G, Vec3& dx_global,
                                        const Vec12& disp) const {
    Vec3 x1 = nodes_[0]->displaced_coords({disp[0], disp[1], disp[2]});
    Vec3 x2 = nodes_[1]->displaced_coords({disp[6], disp[7], disp[8]});
    dx_global = math::vec3_sub(x2, x1);
    Vec3 e1_L = {1.0, 0.0, 0.0};
    Vec3 e1_G = math::vec3_normalize(dx_global);
    T_L_to_G = math::rotation_matrix_from_vectors(e1_L, e1_G);
}

void Element2N::compute_transformation_initial(Mat3x3& T_L_to_G, Vec3& dx_global) const {
    Vec3 x1 = nodes_[0]->coords();
    Vec3 x2 = nodes_[1]->coords();
    dx_global = math::vec3_sub(x2, x1);
    Vec3 e1_L = {1.0, 0.0, 0.0};
    Vec3 e1_G = math::vec3_normalize(dx_global);
    T_L_to_G = math::rotation_matrix_from_vectors(e1_L, e1_G);
}

void Element2N::assemble_12x12_from_6x6_blocks(
    Mat12x12& K,
    const Mat6x6& K11, const Mat6x6& K12,
    const Mat6x6& K21, const Mat6x6& K22) const
{
    for (Index i = 0; i < 6; ++i) {
        for (Index j = 0; j < 6; ++j) {
            K[i][j]         = K11[i][j];
            K[i][j + 6]     = K12[i][j];
            K[i + 6][j]     = K21[i][j];
            K[i + 6][j + 6] = K22[i][j];
        }
    }
}

void Element2N::transform_12x12(
    Mat12x12& K_global, const Mat12x12& K_local,
    const Mat3x3& T) const
{
    std::array<Mat3x3, 16> block_T;
    for (Index bi = 0; bi < 4; ++bi) {
        for (Index bj = 0; bj < 4; ++bj) {
            Mat3x3 tmp{};
            if (bi < 2 && bj < 2) {
                for (Index i = 0; i < 3; ++i)
                    for (Index j = 0; j < 3; ++j)
                        for (Index k = 0; k < 3; ++k)
                            tmp[i][j] += T[i][k] * K_local[bi*3+i][bj*3+j];
                for (Index i = 0; i < 3; ++i)
                    for (Index j = 0; j < 3; ++j) {
                        Real sum = 0.0;
                        for (Index k = 0; k < 3; ++k)
                            sum += tmp[i][k] * T[j][k];
                        K_global[bi*3+i][bj*3+j] = sum;
                    }
            } else {
                for (Index i = 0; i < 3; ++i)
                    for (Index j = 0; j < 3; ++j)
                        K_global[bi*3+i][bj*3+j] = K_local[bi*3+i][bj*3+j];
            }
        }
    }
    
    Mat12x12 full_T{};
    for (Index bi = 0; bi < 4; ++bi) {
        for (Index i = 0; i < 3; ++i) {
            for (Index j = 0; j < 3; ++j) {
                full_T[bi*3+i][bi*3+j] = T[i][j];
            }
        }
    }
    
    Mat12x12 tmp1{};
    for (Index i = 0; i < 12; ++i)
        for (Index j = 0; j < 12; ++j) {
            Real sum = 0.0;
            for (Index k = 0; k < 12; ++k)
                sum += full_T[i][k] * K_local[k][j];
            tmp1[i][j] = sum;
        }
    
    for (Index i = 0; i < 12; ++i)
        for (Index j = 0; j < 12; ++j) {
            Real sum = 0.0;
            for (Index k = 0; k < 12; ++k)
                sum += tmp1[i][k] * full_T[j][k];
            K_global[i][j] = sum;
        }
}

void Element2N::transform_12_vec(
    Vec12& v_global, const Vec12& v_local,
    const Mat3x3& T) const
{
    for (Index bi = 0; bi < 4; ++bi) {
        Vec3 local_block = {v_local[bi*3], v_local[bi*3+1], v_local[bi*3+2]};
        Vec3 global_block = math::mat3_mul_vec3(T, local_block);
        v_global[bi*3]   = global_block[0];
        v_global[bi*3+1] = global_block[1];
        v_global[bi*3+2] = global_block[2];
    }
}

Vec3 Element::midpoint() const {
    Vec3 center = {0, 0, 0};
    for (Index i = 0; i < nnodes_; ++i) {
        const Node* n = node(i);
        if (n) {
            for (Index d = 0; d < 3; ++d) center[d] += n->coords()[d];
        }
    }
    if (nnodes_ > 0) {
        for (Index d = 0; d < 3; ++d) center[d] /= static_cast<Real>(nnodes_);
    }
    return center;
}

Vec3 Element::centroid(const Vec12& disp) const {
    Vec3 center = {0, 0, 0};
    for (Index i = 0; i < nnodes_; ++i) {
        const Node* n = node(i);
        if (n) {
            Vec3 coords = n->coords();
            for (Index d = 0; d < 3; ++d) coords[d] += disp[i * 6 + d];
            for (Index d = 0; d < 3; ++d) center[d] += coords[d];
        }
    }
    if (nnodes_ > 0) {
        for (Index d = 0; d < 3; ++d) center[d] /= static_cast<Real>(nnodes_);
    }
    return center;
}

} // namespace fem
} // namespace hvdc
