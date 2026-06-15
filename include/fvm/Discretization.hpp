#ifndef HVDC_FVM_DISCRETIZATION_HPP
#define HVDC_FVM_DISCRETIZATION_HPP

#include "common/Types.hpp"
#include "common/Vector.hpp"
#include "common/SparseMatrix.hpp"
#include "fvm/Grid.hpp"
#include "fvm/Field.hpp"

namespace hvdc {
namespace fvm {

enum class ConvectionScheme : UInt8 {
    FirstOrderUpwind = 0,
    SecondOrderUpwind = 1,
    QUICK = 2,
    CentralDifferencing = 3,
    MUSCL = 4,
    MinMod = 5
};

enum class InterpolationScheme : UInt8 {
    Linear = 0,
    MidPoint = 1,
    Weighted = 2
};

enum class GradientScheme : UInt8 {
    GreenGauss = 0,
    LeastSquares = 1
};

struct DiscretizationConfig {
    ConvectionScheme conv_scheme = ConvectionScheme::SecondOrderUpwind;
    GradientScheme grad_scheme = GradientScheme::GreenGauss;
    InterpolationScheme interp_scheme = InterpolationScheme::Linear;
    Real upwind_blending = 0.0;
    Real relaxation_U = 0.7;
    Real relaxation_p = 0.3;
    Real relaxation_k = 0.8;
    Real relaxation_eps = 0.8;
    bool use_rhie_chow = true;
    Real underrelax_factor = 1.0;
};

class Discretization {
public:
    Discretization() = default;
    explicit Discretization(Grid* grid,
                            const DiscretizationConfig& cfg = DiscretizationConfig());
    
    void set_grid(Grid* grid) { grid_ = grid; }
    void set_config(const DiscretizationConfig& cfg) { config_ = cfg; }
    
    const DiscretizationConfig& config() const { return config_; }
    DiscretizationConfig& mutable_config() { return config_; }
    
    void interpolate_face_scalar(
        const ScalarField& cell_field,
        const std::vector<Vec3>& cell_grad,
        Index face_id,
        Real& face_value,
        Real& face_upwind_value) const;
    
    void interpolate_face_vector(
        const VectorField& cell_field,
        const std::vector<std::array<Vec3, 3>>& cell_grad,
        Index face_id,
        Vec3& face_value,
        Vec3& face_upwind_value) const;
    
    void compute_scalar_gradient(
        const ScalarField& phi,
        const std::vector<Real>& face_phi,
        std::vector<Vec3>& grad) const;
    
    void compute_vector_gradient(
        const VectorField& U,
        const std::vector<Vec3>& face_U,
        std::vector<std::array<Vec3, 3>>& grad) const;
    
    Real face_mass_flux(
        Index face_id,
        const Vec3& face_U,
        Real face_rho) const;
    
    void discretize_convection_diffusion_scalar(
        const ScalarField& phi,
        const ScalarField& rho,
        const ScalarField& Gamma,
        const std::vector<Vec3>& grad_phi,
        SparseMatrix& A,
        Vector& b,
        Real transient_rho = 0.0,
        const ScalarField* phi_prev = nullptr,
        Real dt = 1.0,
        const VectorField* U = nullptr) const;
    
    void discretize_momentum_equation(
        const VelocityField& U,
        const ScalarField& p,
        const ScalarField& rho,
        const ScalarField& nu_eff,
        const std::vector<std::array<Vec3, 3>>& grad_U,
        const std::vector<Vec3>& grad_p,
        Index component,
        SparseMatrix& A,
        Vector& b,
        Real transient_rho = 0.0,
        const VelocityField* U_prev = nullptr,
        Real dt = 1.0) const;
    
    void discretize_pressure_equation(
        const VelocityField& U,
        const ScalarField& rho,
        const SparseMatrix& Ap,
        const Vector& diag_inv_Ap,
        SparseMatrix& A_p,
        Vector& b_p,
        const VectorField* HbyA = nullptr) const;
    
    void compute_face_velocity_rhie_chow(
        const VelocityField& U,
        const ScalarField& p,
        const SparseMatrix& Ap,
        const Vector& diag_inv_Ap,
        std::vector<Vec3>& face_velocity,
        std::vector<Real>& face_mass_flux,
        const std::vector<Vec3>& grad_p) const;
    
    void apply_relaxation(
        SparseMatrix& A,
        Vector& b,
        Vector& field,
        Real relax_factor) const;
    
    void correct_velocity(
        VelocityField& U,
        const SparseMatrix& Ap,
        const Vector& diag_inv_Ap,
        const ScalarField& rho,
        const VectorField& grad_p,
        const VectorField& dp_correction) const;
    
    void under_relax_field(
        ScalarField& field,
        const ScalarField& prev_field,
        Real relax_factor) const;
    
    void under_relax_vector_field(
        VectorField& field,
        const VectorField& prev_field,
        Real relax_factor) const;

protected:
    Grid* grid_ = nullptr;
    DiscretizationConfig config_;
    
    Real compute_convective_weight(
        Index face_id,
        Real m_dot,
        Real Pe) const;
    
    void compute_TVD_limiter(
        Real phi_r,
        ConvectionScheme scheme,
        Real& psi) const;
    
    Real interpolate_central(
        Index face_id,
        const ScalarField& phi) const;
    
    Real interpolate_second_order_upwind(
        Index face_id,
        const ScalarField& phi,
        const std::vector<Vec3>& grad,
        Real m_dot_sign) const;
    
    Vec3 interpolate_central_vector(
        Index face_id,
        const VectorField& U) const;
};

} // namespace fvm
} // namespace hvdc

#endif // HVDC_FVM_DISCRETIZATION_HPP
