#include "common/MathUtils.hpp"
#include <cmath>
#include <limits>

namespace hvdc {
namespace math {

Mat3x3 rotation_matrix_from_vectors(const Vec3& a, const Vec3& b) {
    Vec3 a_unit = vec3_normalize(a);
    Vec3 b_unit = vec3_normalize(b);
    
    Real dot = vec3_dot(a_unit, b_unit);
    
    if (dot > 1.0 - EPS) {
        return mat3_identity();
    }
    if (dot < -1.0 + EPS) {
        Vec3 perp;
        if (std::fabs(a_unit[0]) < 0.9) {
            perp = vec3_normalize(vec3_cross(a_unit, {1.0, 0.0, 0.0}));
        } else {
            perp = vec3_normalize(vec3_cross(a_unit, {0.0, 1.0, 0.0}));
        }
        Mat3x3 R = mat3_identity();
        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < 3; ++j)
                R[i][j] = 2.0 * perp[i] * perp[j] - (i == j ? 1.0 : 0.0);
        return R;
    }
    
    Vec3 v = vec3_cross(a_unit, b_unit);
    Real s = vec3_norm(v);
    Real c = dot;
    
    Mat3x3 vx{};
    vx[0][1] = -v[2]; vx[0][2] = v[1];
    vx[1][0] = v[2];  vx[1][2] = -v[0];
    vx[2][0] = -v[1]; vx[2][1] = v[0];
    
    Mat3x3 vxvx{};
    for (Index i = 0; i < 3; ++i)
        for (Index j = 0; j < 3; ++j)
            for (Index k = 0; k < 3; ++k)
                vxvx[i][j] += vx[i][k] * vx[k][j];
    
    Real factor = (1.0 - c) / (s * s);
    
    Mat3x3 R = mat3_identity();
    for (Index i = 0; i < 3; ++i) {
        for (Index j = 0; j < 3; ++j) {
            R[i][j] += vx[i][j] + factor * vxvx[i][j];
        }
    }
    
    return R;
}

Mat3x3 rotation_matrix_axis_angle(const Vec3& axis, Real angle) {
    Vec3 n = vec3_normalize(axis);
    Real c = std::cos(angle);
    Real s = std::sin(angle);
    Real t = 1.0 - c;
    
    Mat3x3 R{};
    R[0][0] = t*n[0]*n[0] + c;
    R[0][1] = t*n[0]*n[1] - s*n[2];
    R[0][2] = t*n[0]*n[2] + s*n[1];
    R[1][0] = t*n[1]*n[0] + s*n[2];
    R[1][1] = t*n[1]*n[1] + c;
    R[1][2] = t*n[1]*n[2] - s*n[0];
    R[2][0] = t*n[2]*n[0] - s*n[1];
    R[2][1] = t*n[2]*n[1] + s*n[0];
    R[2][2] = t*n[2]*n[2] + c;
    
    return R;
}

void mat3_to_euler_zyx(const Mat3x3& R, Real& rx, Real& ry, Real& rz) {
    ry = std::asin(-R[2][0]);
    Real cy = std::cos(ry);
    
    if (std::fabs(cy) > EPS) {
        rz = std::atan2(R[1][0], R[0][0]);
        rx = std::atan2(R[2][1], R[2][2]);
    } else {
        rz = 0.0;
        rx = std::atan2(-R[1][2], R[1][1]);
    }
}

Mat3x3 euler_zyx_to_mat3(Real rx, Real ry, Real rz) {
    Real cx = std::cos(rx), sx = std::sin(rx);
    Real cy = std::cos(ry), sy = std::sin(ry);
    Real cz = std::cos(rz), sz = std::sin(rz);
    
    Mat3x3 R{};
    R[0][0] = cy*cz;
    R[0][1] = sx*sy*cz - cx*sz;
    R[0][2] = cx*sy*cz + sx*sz;
    R[1][0] = cy*sz;
    R[1][1] = sx*sy*sz + cx*cz;
    R[1][2] = cx*sy*sz - sx*cz;
    R[2][0] = -sy;
    R[2][1] = sx*cy;
    R[2][2] = cx*cy;
    
    return R;
}

void newton_raphson_solve(
    const std::function<void(const Vec&, Vec&, std::vector<std::vector<Real>>&)>& assemble,
    Vec& x, Index max_iter, Real tol, bool& converged, Index& iters
) {
    converged = false;
    const Index n = static_cast<Index>(x.size());
    
    Vec residual(n);
    std::vector<std::vector<Real>> jacobian(n, std::vector<Real>(n, 0.0));
    Vec dx(n);
    
    for (iters = 0; iters < max_iter; ++iters) {
        assemble(x, residual, jacobian);
        
        Real res_norm = 0.0;
        for (Index i = 0; i < n; ++i) {
            res_norm += residual[i] * residual[i];
        }
        res_norm = std::sqrt(res_norm);
        
        if (res_norm < tol) {
            converged = true;
            return;
        }
        
        for (Index i = 0; i < n; ++i) {
            dx[i] = -residual[i];
        }
        
        for (Index i = 0; i < n; ++i) {
            Real pivot = jacobian[i][i];
            if (std::fabs(pivot) < EPS_MEDIUM) continue;
            for (Index j = i + 1; j < n; ++j) {
                Real factor = jacobian[j][i] / pivot;
                for (Index k = i + 1; k < n; ++k) {
                    jacobian[j][k] -= factor * jacobian[i][k];
                }
                dx[j] -= factor * dx[i];
            }
        }
        
        for (Index i = n - 1; i >= 0; --i) {
            Real sum = dx[i];
            for (Index j = i + 1; j < n; ++j) {
                sum -= jacobian[i][j] * dx[j];
            }
            if (std::fabs(jacobian[i][i]) > EPS_MEDIUM) {
                dx[i] = sum / jacobian[i][i];
            }
        }
        
        for (Index i = 0; i < n; ++i) {
            x[i] += dx[i];
        }
    }
}

} // namespace math
} // namespace hvdc
