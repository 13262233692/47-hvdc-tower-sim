#ifndef HVDC_COMMON_MATH_UTILS_HPP
#define HVDC_COMMON_MATH_UTILS_HPP

#include "common/Types.hpp"
#include <cmath>
#include <algorithm>

namespace hvdc {
namespace math {

inline Real sign(Real x) {
    return (Real(0) < x) - (x < Real(0));
}

inline Real clamp(Real x, Real lo, Real hi) {
    return std::max(lo, std::min(hi, x));
}

inline Real lerp(Real a, Real b, Real t) {
    return a + t * (b - a);
}

inline Real smoothstep(Real edge0, Real edge1, Real x) {
    Real t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

inline bool approx_equal(Real a, Real b, Real tol = EPS_LARGE) {
    return std::fabs(a - b) <= tol * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

inline bool is_finite(Real x) {
    return std::isfinite(x);
}

inline Vec3 vec3_zero() {
    return {0.0, 0.0, 0.0};
}

inline Vec3 vec3(Real x, Real y, Real z) {
    return {x, y, z};
}

inline Real vec3_dot(const Vec3& a, const Vec3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline Vec3 vec3_cross(const Vec3& a, const Vec3& b) {
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

inline Real vec3_norm2(const Vec3& v) {
    return vec3_dot(v, v);
}

inline Real vec3_norm(const Vec3& v) {
    return std::sqrt(vec3_norm2(v));
}

inline Vec3 vec3_normalize(const Vec3& v) {
    Real n = vec3_norm(v);
    if (n < EPS) return {1.0, 0.0, 0.0};
    Real inv = 1.0 / n;
    return {v[0]*inv, v[1]*inv, v[2]*inv};
}

inline Vec3 vec3_add(const Vec3& a, const Vec3& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

inline Vec3 vec3_sub(const Vec3& a, const Vec3& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

inline Vec3 vec3_scale(const Vec3& v, Real s) {
    return {v[0]*s, v[1]*s, v[2]*s};
}

inline void vec3_add_inplace(Vec3& a, const Vec3& b) {
    a[0] += b[0]; a[1] += b[1]; a[2] += b[2];
}

inline void vec3_sub_inplace(Vec3& a, const Vec3& b) {
    a[0] -= b[0]; a[1] -= b[1]; a[2] -= b[2];
}

inline void vec3_scale_inplace(Vec3& v, Real s) {
    v[0] *= s; v[1] *= s; v[2] *= s;
}

inline Vec3 vec3_neg(const Vec3& v) {
    return {-v[0], -v[1], -v[2]};
}

inline Mat3x3 mat3_identity() {
    Mat3x3 m{};
    m[0][0] = m[1][1] = m[2][2] = 1.0;
    return m;
}

inline Mat3x3 mat3_zero() {
    Mat3x3 m{};
    return m;
}

inline Mat3x3 mat3_transpose(const Mat3x3& m) {
    return {{
        {m[0][0], m[1][0], m[2][0]},
        {m[0][1], m[1][1], m[2][1]},
        {m[0][2], m[1][2], m[2][2]}
    }};
}

inline Vec3 mat3_mul_vec3(const Mat3x3& m, const Vec3& v) {
    return {
        m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
        m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
        m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2]
    };
}

inline Mat3x3 mat3_mul(const Mat3x3& a, const Mat3x3& b) {
    Mat3x3 c{};
    for (Index i = 0; i < 3; ++i)
        for (Index j = 0; j < 3; ++j)
            for (Index k = 0; k < 3; ++k)
                c[i][j] += a[i][k] * b[k][j];
    return c;
}

inline void mat3_add_inplace(Mat3x3& a, const Mat3x3& b) {
    for (Index i = 0; i < 3; ++i)
        for (Index j = 0; j < 3; ++j)
            a[i][j] += b[i][j];
}

inline void mat3_scale_inplace(Mat3x3& m, Real s) {
    for (Index i = 0; i < 3; ++i)
        for (Index j = 0; j < 3; ++j)
            m[i][j] *= s;
}

Mat3x3 rotation_matrix_from_vectors(const Vec3& a, const Vec3& b);

Mat3x3 rotation_matrix_axis_angle(const Vec3& axis, Real angle);

void mat3_to_euler_zyx(const Mat3x3& R, Real& rx, Real& ry, Real& rz);

Mat3x3 euler_zyx_to_mat3(Real rx, Real ry, Real rz);

template <Index N>
void vec_set_zero(VecN<N>& v) {
    v.fill(0.0);
}

template <Index N>
Real vec_dot(const VecN<N>& a, const VecN<N>& b) {
    Real sum = 0.0;
    for (Index i = 0; i < N; ++i) sum += a[i] * b[i];
    return sum;
}

template <Index N>
void vec_add(VecN<N>& c, const VecN<N>& a, const VecN<N>& b) {
    for (Index i = 0; i < N; ++i) c[i] = a[i] + b[i];
}

template <Index N>
void vec_sub(VecN<N>& c, const VecN<N>& a, const VecN<N>& b) {
    for (Index i = 0; i < N; ++i) c[i] = a[i] - b[i];
}

template <Index N>
void vec_scale(VecN<N>& c, const VecN<N>& a, Real s) {
    for (Index i = 0; i < N; ++i) c[i] = a[i] * s;
}

template <Index N>
void vec_axpy(VecN<N>& y, Real alpha, const VecN<N>& x) {
    for (Index i = 0; i < N; ++i) y[i] += alpha * x[i];
}

template <Index M, Index N>
void mat_set_zero(MatMN<M,N>& m) {
    for (auto& row : m) row.fill(0.0);
}

template <Index M, Index N>
void mat_add_inplace(MatMN<M,N>& a, const MatMN<M,N>& b) {
    for (Index i = 0; i < M; ++i)
        for (Index j = 0; j < N; ++j)
            a[i][j] += b[i][j];
}

template <Index M, Index N>
void mat_scale_inplace(MatMN<M,N>& m, Real s) {
    for (Index i = 0; i < M; ++i)
        for (Index j = 0; j < N; ++j)
            m[i][j] *= s;
}

template <Index M, Index N, Index K>
void mat_mul(MatMN<M,K>& C, const MatMN<M,N>& A, const MatMN<N,K>& B) {
    mat_set_zero(C);
    for (Index i = 0; i < M; ++i)
        for (Index j = 0; j < K; ++j)
            for (Index k = 0; k < N; ++k)
                C[i][j] += A[i][k] * B[k][j];
}

template <Index M, Index N>
void mat_mul_vec(VecN<M>& y, const MatMN<M,N>& A, const VecN<N>& x) {
    vec_set_zero(y);
    for (Index i = 0; i < M; ++i)
        for (Index j = 0; j < N; ++j)
            y[i] += A[i][j] * x[j];
}

template <Index N>
void mat_transpose(MatMN<N,N>& At, const MatMN<N,N>& A) {
    for (Index i = 0; i < N; ++i)
        for (Index j = 0; j < N; ++j)
            At[j][i] = A[i][j];
}

void newton_raphson_solve(
    const std::function<void(const Vec&, Vec&, std::vector<std::vector<Real>>&)>& assemble,
    Vec& x, Index max_iter, Real tol, bool& converged, Index& iters
);

} // namespace math
} // namespace hvdc

#endif // HVDC_COMMON_MATH_UTILS_HPP
