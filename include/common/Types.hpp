#ifndef HVDC_COMMON_TYPES_HPP
#define HVDC_COMMON_TYPES_HPP

#include <cstdint>
#include <cstddef>
#include <complex>
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <functional>

namespace hvdc {

using Real = double;
using Complex = std::complex<Real>;
using Index = std::ptrdiff_t;
using UInt = std::uint64_t;
using UInt8 = std::uint8_t;
using UInt16 = std::uint16_t;
using UInt32 = std::uint32_t;
using UInt64 = std::uint64_t;

constexpr Real PI = 3.14159265358979323846;
constexpr Real TWO_PI = 6.28318530717958647692;
constexpr Real HALF_PI = 1.57079632679489661923;
constexpr Real RAD_TO_DEG = 57.2957795130823208768;
constexpr Real DEG_TO_RAD = 0.01745329251994329577;
constexpr Real EPS = 1.0e-15;
constexpr Real EPS_MEDIUM = 1.0e-10;
constexpr Real EPS_LARGE = 1.0e-6;
constexpr Real INF = 1.0e30;
constexpr Real GRAVITY = 9.80665;

template <Index N>
using VecN = std::array<Real, N>;

using Vec2 = VecN<2>;
using Vec3 = VecN<3>;
using Vec6 = VecN<6>;
using Vec12 = VecN<12>;

template <Index M, Index N>
using MatMN = std::array<std::array<Real, N>, M>;

using Mat2x2 = MatMN<2, 2>;
using Mat3x3 = MatMN<3, 3>;
using Mat6x6 = MatMN<6, 6>;
using Mat12x12 = MatMN<12, 12>;

using Vec = std::vector<Real>;
using IndexVec = std::vector<Index>;

enum class ElementType : UInt8 {
    Unknown = 0,
    TrussLinear = 1,
    TrussNonlinear = 2,
    BeamLinear = 3,
    BeamNonlinear = 4,
    ConductorCable = 5,
    ShellLinear = 6,
    SolidHex = 7
};

enum class DOFType : UInt8 {
    UX = 0, UY = 1, UZ = 2,
    ROTX = 3, ROTY = 4, ROTZ = 5,
    TEMP = 6, PRESSURE = 7
};

enum class BoundaryType : UInt8 {
    Free = 0,
    Fixed = 1,
    Displacement = 2,
    Coupled = 3
};

enum class AnalysisType : UInt8 {
    StaticLinear = 0,
    StaticNonlinear = 1,
    Modal = 2,
    TransientLinear = 3,
    TransientNonlinear = 4,
    Buckling = 5,
    FSI_Coupled = 6
};

enum class CouplingScheme : UInt8 {
    Explicit = 0,
    Implicit_Newton = 1,
    Implicit_Aitken = 2,
    Implicit_IQN_ILS = 3
};

enum class MaterialModel : UInt8 {
    LinearElastic = 0,
    ElasticPlastic = 1,
    Hyperelastic = 2,
    Orthotropic = 3
};

struct SimulationConfig {
    Real time_start = 0.0;
    Real time_end = 10.0;
    Real dt = 0.01;
    Real dt_min = 1.0e-6;
    Real dt_max = 0.1;
    Index max_newton_iters = 50;
    Real newton_tol = 1.0e-8;
    Real newton_rel_tol = 1.0e-6;
    Index max_coupling_iters = 20;
    Real coupling_tol = 1.0e-6;
    Real relaxation_param = 0.5;
    Real newmark_beta = 0.25;
    Real newmark_gamma = 0.5;
    Real beta = 0.25;
    Real gamma = 0.5;
    AnalysisType analysis = AnalysisType::TransientNonlinear;
    AnalysisType analysis_type = AnalysisType::TransientNonlinear;
    CouplingScheme coupling = CouplingScheme::Implicit_Aitken;
    bool use_geometric_nonlinearity = true;
    bool use_material_nonlinearity = false;
    bool enable_fsi = true;
    Real t_start = 0.0;
    Real t_end = 10.0;
};

} // namespace hvdc

#endif // HVDC_COMMON_TYPES_HPP
