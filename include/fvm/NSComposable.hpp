#ifndef HVDC_FVM_NSCOMPOSABLE_HPP
#define HVDC_FVM_NSCOMPOSABLE_HPP

#include "common/Types.hpp"
#include "common/Vector.hpp"
#include "common/SparseMatrix.hpp"
#include "fvm/Grid.hpp"
#include "fvm/Field.hpp"
#include "fvm/Discretization.hpp"

namespace hvdc {
namespace fvm {

struct FlowConfig {
    Real rho_ref = 1.225;
    Real mu_ref = 1.789e-5;
    Real p_ref = 101325.0;
    Real T_ref = 288.15;
    Real gravity[3] = {0.0, 0.0, -9.80665};
    bool use_boussinesq = false;
    Real beta_thermal = 0.0034;
    
    bool enable_turbulence = false;
    std::string turbulence_model = "kEpsilon";
    
    Real C_mu = 0.09;
    Real C1_eps = 1.44;
    Real C2_eps = 1.92;
    Real sigma_k = 1.0;
    Real sigma_eps = 1.3;
    Real sigma_omega = 2.0;
    Real beta_k_omega = 0.075;
    
    bool is_compressible = false;
    bool include_body_force = true;
    bool use_simplec = false;
    
    Index pressure_outer_iters = 4;
    Index pressure_inner_iters = 2;
    Index momentum_predictor_iters = 1;
    
    Real transient_term(Real field, Real field_prev, Real volume,
                        Real rho, Real dt) const {
        return rho * volume * (field - field_prev) / dt;
    }
};

class NSComposable {
public:
    NSComposable() = default;
    NSComposable(Grid* grid, const FlowConfig& flow_cfg,
                 const DiscretizationConfig& disc_cfg = DiscretizationConfig());
    
    void set_grid(Grid* grid) { grid_ = grid; disc_.set_grid(grid); }
    void set_flow_config(const FlowConfig& cfg) { flow_cfg_ = cfg; }
    void set_disc_config(const DiscretizationConfig& cfg) { disc_.set_config(cfg); }
    
    const FlowConfig& flow_config() const { return flow_cfg_; }
    DiscretizationConfig& disc_config() { return disc_.mutable_config(); }
    
    void initialize_state(FluidState& state) const;
    
    void compute_effective_viscosity(
        const FluidState& state,
        ScalarField& nu_eff,
        ScalarField& Gamma_U,
        ScalarField& Gamma_k,
        ScalarField& Gamma_eps) const;
    
    void assemble_momentum_predictor(
        FluidState& state,
        const FluidState* state_prev,
        Real dt,
        Index component,
        SparseMatrix& A_U,
        Vector& b_U,
        VectorField& HbyA,
        Vector& diag_inv_Ap,
        VectorField& grad_p) const;
    
    void assemble_pressure_correction(
        FluidState& state,
        const SparseMatrix& Ap,
        const Vector& diag_inv_Ap,
        const VectorField& HbyA,
        SparseMatrix& A_p,
        Vector& b_p,
        VectorField& face_velocity,
        std::vector<Real>& face_mass_flux,
        const VectorField& grad_p) const;
    
    void correct_velocity_and_pressure(
        FluidState& state,
        const SparseMatrix& Ap,
        const Vector& diag_inv_Ap,
        const VectorField& HbyA,
        const ScalarField& dp,
        VectorField& face_velocity,
        std::vector<Real>& face_mass_flux,
        const VectorField& grad_p) const;
    
    ScalarField compute_kinetic_energy_dissipation(
        const FluidState& state,
        const std::vector<std::array<Vec3, 3>>& grad_U) const;
    
    void assemble_turbulence_equations(
        FluidState& state,
        const FluidState* state_prev,
        Real dt,
        const std::vector<std::array<Vec3, 3>>& grad_U,
        SparseMatrix& A_k, Vector& b_k,
        SparseMatrix& A_eps, Vector& b_eps) const;
    
    void update_turbulence_from_solution(
        FluidState& state,
        const ScalarField& k_new,
        const ScalarField& eps_new) const;
    
    Real compute_cfl_number(
        const FluidState& state,
        const std::vector<Real>& face_flux,
        Real dt) const;
    
    Vec3 compute_total_force_on_boundary(
        Index patch_id,
        const FluidState& state,
        const VectorField& grad_p,
        const std::vector<std::array<Vec3, 3>>& grad_U) const;
    
    Real compute_mass_imbalance(
        const std::vector<Real>& face_flux) const;
    
    std::array<Real, 6> compute_flow_stats(const FluidState& state) const;

protected:
    Grid* grid_ = nullptr;
    FlowConfig flow_cfg_;
    Discretization disc_;
    
    ScalarField compute_production_k(
        const FluidState& state,
        const std::vector<std::array<Vec3, 3>>& grad_U) const;
    
    void compute_strain_rate_tensor(
        const std::array<Vec3, 3>& grad,
        Mat3x3& S) const;
    
    Real compute_strain_rate_magnitude(const Mat3x3& S) const;
};

} // namespace fvm
} // namespace hvdc

#endif // HVDC_FVM_NSCOMPOSABLE_HPP
