#ifndef HVDC_FSI_COUPLINGMANAGER_HPP
#define HVDC_FSI_COUPLINGMANAGER_HPP

#include "common/Types.hpp"
#include "common/Timer.hpp"
#include "fem/FEModel.hpp"
#include "fem/Assembler.hpp"
#include "fem/Solver.hpp"
#include "fvm/FluidSolver.hpp"
#include "fsi/InterfaceMapper.hpp"
#include "fsi/AerodynamicLoads.hpp"
#include "fsi/LoadTransfer.hpp"
#include <vector>
#include <functional>
#include <memory>

namespace hvdc {
namespace fsi {

enum class CouplingIterationMethod : UInt8 {
    Explicit = 0,
    Implicit_Aitken = 1,
    Implicit_Anderson = 2,
    Implicit_IQN_ILS = 3,
    Implicit_Broyden = 4
};

enum class FSIConvergenceCriterion : UInt8 {
    Displacement = 0,
    Force = 1,
    DisplacementAndForce = 2,
    ResidualNorm = 3
};

struct CouplingConvergenceData {
    Real disp_norm = 0.0;
    Real disp_norm_prev = 0.0;
    Real force_norm = 0.0;
    Real force_norm_prev = 0.0;
    Real residual_norm = 0.0;
    Real energy_balance_error = 0.0;
    Real mass_conservation_error = 0.0;
    Real total_force_x = 0.0;
    Real total_force_y = 0.0;
    Real total_force_z = 0.0;
    Index iterations = 0;
    bool converged = false;
};

struct CouplingHistory {
    std::vector<Vector> displacement_history;
    std::vector<Vector> force_history;
    std::vector<Real> time_history;
};

class CouplingManager {
public:
    CouplingManager() = default;
    CouplingManager(
        fem::FEModel* struct_model,
        fem::Assembler* assembler,
        fem::Solver* struct_solver,
        fvm::FluidSolver* fluid_solver,
        InterfaceMapper* mapper,
        AerodynamicLoads* aero,
        LoadTransfer* transfer,
        const SimulationConfig& cfg = SimulationConfig());
    
    void set_config(const SimulationConfig& cfg) { config_ = cfg; }
    void set_coupling_scheme(CouplingIterationMethod m) { method_ = m; }
    void set_convergence_criterion(FSIConvergenceCriterion c) { conv_criterion_ = c; }
    
    void set_timestep(Real dt) { dt_ = dt; }
    void set_max_coupling_iterations(Index n) { max_iter_ = n; }
    void set_convergence_tolerance(Real tol) { tol_ = tol; }
    
    const SimulationConfig& config() const { return config_; }
    
    CouplingConvergenceData perform_fsi_step(
        Vector& struct_displacement,
        Vector& struct_velocity,
        Vector& struct_acceleration,
        fvm::FluidState& fluid_state,
        Real time);
    
    void run_coupled_simulation(
        Vector& struct_displacement,
        Vector& struct_velocity,
        Vector& struct_acceleration,
        fvm::FluidState& fluid_state,
        std::function<void(Index, Real, const Vector&, const fvm::FluidState&)> 
            callback = nullptr);
    
    CouplingConvergenceData perform_explicit_step(
        Vector& struct_disp,
        Vector& struct_vel,
        Vector& struct_acc,
        fvm::FluidState& fluid_state,
        Real time);
    
    CouplingConvergenceData perform_implicit_aitken_step(
        Vector& struct_disp,
        Vector& struct_vel,
        Vector& struct_acc,
        fvm::FluidState& fluid_state,
        Real time);
    
    CouplingConvergenceData perform_implicit_anderson_step(
        Vector& struct_disp,
        Vector& struct_vel,
        Vector& struct_acc,
        fvm::FluidState& fluid_state,
        Real time);
    
    CouplingConvergenceData perform_implicit_iqn_ils_step(
        Vector& struct_disp,
        Vector& struct_vel,
        Vector& struct_acc,
        fvm::FluidState& fluid_state,
        Real time);
    
    CouplingConvergenceData perform_implicit_broyden_step(
        Vector& struct_disp,
        Vector& struct_vel,
        Vector& struct_acc,
        fvm::FluidState& fluid_state,
        Real time);
    
    void compute_interface_residual(
        const Vector& disp_old, const Vector& disp_new,
        const Vector& force_old, const Vector& force_new,
        CouplingConvergenceData& data) const;
    
    bool check_convergence(const CouplingConvergenceData& data) const;
    
    void apply_aitken_relaxation(
        const Vector& x_old, const Vector& x_new,
        Vector& x_relaxed, Real& omega,
        const std::vector<Vector>& x_history) const;
    
    void apply_anderson_mixing(
        const std::vector<Vector>& x_history,
        const std::vector<Vector>& r_history,
        Vector& x_new,
        Index m = 5) const;
    
    void apply_iqn_ils(
        const std::vector<Vector>& displacement_history,
        const std::vector<Vector>& force_history,
        Vector& next_displacement) const;
    
    void update_broyden_jacobian(
        SparseMatrix& J_approx,
        const Vector& delta_x, const Vector& delta_F,
        Real gamma = 1.0) const;
    
    const CouplingHistory& history() const { return history_; }
    void clear_history() {
        history_.displacement_history.clear();
        history_.force_history.clear();
        history_.time_history.clear();
    }
    
    void compute_coupling_diagnostics(
        const Vector& struct_disp,
        const fvm::FluidState& fluid_state,
        Real time) const;
    
    Real estimate_interface_added_mass(
        const Vector& trial_disp,
        const fvm::FluidState& fluid_state,
        Real dt) const;

protected:
    fem::FEModel* struct_model_ = nullptr;
    fem::Assembler* assembler_ = nullptr;
    fem::Solver* struct_solver_ = nullptr;
    fvm::FluidSolver* fluid_solver_ = nullptr;
    InterfaceMapper* mapper_ = nullptr;
    AerodynamicLoads* aero_ = nullptr;
    LoadTransfer* transfer_ = nullptr;
    
    SimulationConfig config_;
    CouplingIterationMethod method_ = CouplingIterationMethod::Implicit_Aitken;
    FSIConvergenceCriterion conv_criterion_ = FSIConvergenceCriterion::DisplacementAndForce;
    
    Real dt_ = 0.001;
    Real tol_ = 1.0e-6;
    Index max_iter_ = 50;
    Real aitken_omega_ = 0.1;
    
    CouplingHistory history_;
    Timer timer_;
    
    void do_one_sided_coupling(
        Vector& struct_disp, Vector& struct_vel, Vector& struct_acc,
        fvm::FluidState& fluid_state,
        Real time,
        Vector& interface_forces,
        Vector& interface_displacements,
        bool compute_geometry_updates = true);
    
    void compute_updated_interface_displacements(
        const Vector& struct_disp,
        Vector& interface_disp_vec) const;
    
    void reassemble_structure_with_aero_loads(
        Vector& struct_disp, Vector& struct_vel, Vector& struct_acc,
        const Vector& aero_forces,
        Real time);
    
    void update_fluid_mesh_boundary_motion(
        const Vector& struct_disp);
    
    Real compute_normalized_residual(
        const Vector& current, const Vector& previous) const;
};

} // namespace fsi
} // namespace hvdc

#endif // HVDC_FSI_COUPLINGMANAGER_HPP
