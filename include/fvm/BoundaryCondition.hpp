#ifndef HVDC_FVM_BOUNDARYCONDITION_HPP
#define HVDC_FVM_BOUNDARYCONDITION_HPP

#include "common/Types.hpp"
#include "common/SparseMatrix.hpp"
#include "fvm/Field.hpp"
#include <vector>
#include <string>
#include <functional>
#include <memory>
#include <cmath>

namespace hvdc {
namespace fvm {

class Grid;

enum class BCType : UInt8 {
    InletVelocity = 0,
    InletPressure = 1,
    InletTotalPressure = 2,
    OutletPressure = 3,
    OutletFlow = 4,
    WallNoSlip = 5,
    WallSlip = 6,
    WallMoving = 7,
    Symmetry = 8,
    Cyclic = 9,
    FarField = 10,
    InterfaceFSI = 11,
    Empty = 12,
    ZeroGradient = 13
};

class BoundaryCondition {
public:
    BoundaryCondition() = default;
    BoundaryCondition(const std::string& name, Index patch_id, BCType type);
    
    const std::string& name() const { return name_; }
    Index patch_id() const { return patch_id_; }
    BCType type() const { return type_; }
    
    void set_velocity(const Vec3& U) { U_inlet_ = U; }
    void set_pressure(Real p) { p_value_ = p; }
    void set_total_pressure(Real p0) { p0_ = p0; }
    void set_temperature(Real T) { T_value_ = T; }
    void set_turbulence(Real k, Real eps, Real omega = 0.0) {
        k_inlet_ = k; eps_inlet_ = eps; omega_inlet_ = omega;
    }
    
    Vec3 velocity() const { return U_inlet_; }
    Real pressure() const { return p_value_; }
    Real total_pressure() const { return p0_; }
    
    virtual void apply_velocity(
        const Grid& grid,
        const FluidState& state,
        VelocityField& U,
        std::vector<Vec3>& face_velocity) const;
    
    virtual void apply_pressure(
        const Grid& grid,
        const FluidState& state,
        PressureField& p,
        std::vector<Real>& face_pressure) const;
    
    virtual void apply_turbulence(
        const Grid& grid,
        FluidState& state) const;
    
    virtual void apply_temperature(
        const Grid& grid,
        ScalarField& T) const;
    
    void update_matrix_rhs_for_bc(
        const Grid& grid,
        Index component,
        SparseMatrix& A,
        Vector& b,
        const VectorField& field) const;

protected:
    std::string name_;
    Index patch_id_ = -1;
    BCType type_ = BCType::ZeroGradient;
    
    Vec3 U_inlet_ = {0.0, 0.0, 0.0};
    Real p_value_ = 0.0;
    Real p0_ = 101325.0;
    Real T_value_ = 288.15;
    Real k_inlet_ = 0.01;
    Real eps_inlet_ = 0.1;
    Real omega_inlet_ = 10.0;
};

class BCPressureOutlet : public BoundaryCondition {
public:
    BCPressureOutlet(const std::string& name, Index patch_id, Real p_outlet)
        : BoundaryCondition(name, patch_id, BCType::OutletPressure)
    {
        p_value_ = p_outlet;
    }
    
    void apply_pressure(
        const Grid& grid,
        const FluidState& state,
        PressureField& p,
        std::vector<Real>& face_pressure) const override;
};

class BCInletVelocity : public BoundaryCondition {
public:
    BCInletVelocity(const std::string& name, Index patch_id, const Vec3& U_in)
        : BoundaryCondition(name, patch_id, BCType::InletVelocity)
    {
        U_inlet_ = U_in;
        compute_turbulence_from_velocity();
    }
    
    void apply_velocity(
        const Grid& grid,
        const FluidState& state,
        VelocityField& U,
        std::vector<Vec3>& face_velocity) const override;
    
    void apply_turbulence(
        const Grid& grid,
        FluidState& state) const override;

protected:
    void compute_turbulence_from_velocity();
};

class BCWallNoSlip : public BoundaryCondition {
public:
    BCWallNoSlip(const std::string& name, Index patch_id)
        : BoundaryCondition(name, patch_id, BCType::WallNoSlip)
    {
        wall_height_ = 0.0;
        use_wall_functions_ = true;
    }
    
    void apply_velocity(
        const Grid& grid,
        const FluidState& state,
        VelocityField& U,
        std::vector<Vec3>& face_velocity) const override;
    
    void set_wall_roughness(Real height) { wall_height_ = height; }
    void set_wall_motion(const Vec3& velocity) { wall_velocity_ = velocity; }

protected:
    Real wall_height_ = 0.0;
    Vec3 wall_velocity_ = {0, 0, 0};
    bool use_wall_functions_ = true;
    
    Real compute_y_plus(Real u_tau, Real y, Real nu) const;
    Real compute_wall_shear_stress(
        Real U_tangential, Real y, Real nu,
        Real& u_tau, Real& y_plus) const;
};

class BCSymmetry : public BoundaryCondition {
public:
    BCSymmetry(const std::string& name, Index patch_id)
        : BoundaryCondition(name, patch_id, BCType::Symmetry) {}
    
    void apply_velocity(
        const Grid& grid,
        const FluidState& state,
        VelocityField& U,
        std::vector<Vec3>& face_velocity) const override;
};

class BCFarField : public BoundaryCondition {
public:
    BCFarField(const std::string& name, Index patch_id,
               Real p_inf, const Vec3& U_inf, Real Ma = 0.1)
        : BoundaryCondition(name, patch_id, BCType::FarField)
    {
        p_value_ = p_inf;
        U_inlet_ = U_inf;
        mach_number_ = Ma;
    }
    
    void apply_velocity(
        const Grid& grid,
        const FluidState& state,
        VelocityField& U,
        std::vector<Vec3>& face_velocity) const override;
    
    void apply_pressure(
        const Grid& grid,
        const FluidState& state,
        PressureField& p,
        std::vector<Real>& face_pressure) const override;

protected:
    Real mach_number_ = 0.1;
};

class BCInterfaceFSI : public BoundaryCondition {
public:
    BCInterfaceFSI(const std::string& name, Index patch_id)
        : BoundaryCondition(name, patch_id, BCType::InterfaceFSI)
    {
        is_moving_wall_ = true;
    }
    
    void apply_velocity(
        const Grid& grid,
        const FluidState& state,
        VelocityField& U,
        std::vector<Vec3>& face_velocity) const override;
    
    void set_interface_velocity(const std::vector<Vec3>& vel) {
        interface_velocity_ = vel;
    }
    
    void set_interface_displacement(const std::vector<Vec3>& disp) {
        interface_displacement_ = disp;
    }
    
    const std::vector<Vec3>& interface_velocity() const {
        return interface_velocity_;
    }
    
    void set_moving_wall(bool b) { is_moving_wall_ = b; }

protected:
    bool is_moving_wall_ = true;
    std::vector<Vec3> interface_velocity_;
    std::vector<Vec3> interface_displacement_;
};

class BoundaryConditionManager {
public:
    BoundaryConditionManager() = default;
    
    Index add_bc(std::unique_ptr<BoundaryCondition> bc);
    
    BoundaryCondition* bc(Index i) { return bcs_[i].get(); }
    const BoundaryCondition* bc(Index i) const { return bcs_[i].get(); }
    
    BoundaryCondition* bc_by_name(const std::string& name);
    BoundaryCondition* bc_by_patch(Index patch_id);
    
    Index num_bcs() const { return static_cast<Index>(bcs_.size()); }
    
    void apply_all_velocity_bc(
        const Grid& grid,
        const FluidState& state,
        VelocityField& U,
        std::vector<Vec3>& face_velocity) const;
    
    void apply_all_pressure_bc(
        const Grid& grid,
        const FluidState& state,
        PressureField& p,
        std::vector<Real>& face_pressure) const;
    
    void apply_all_turbulence_bc(
        const Grid& grid,
        FluidState& state) const;
    
    void apply_all_temperature_bc(
        const Grid& grid,
        ScalarField& T) const;
    
    void update_linear_system_bc(
        const Grid& grid,
        const std::string& field_name,
        Index component,
        SparseMatrix& A,
        Vector& b,
        const VectorField& field) const;

protected:
    std::vector<std::unique_ptr<BoundaryCondition>> bcs_;
    std::unordered_map<std::string, Index> name_map_;
    std::unordered_map<Index, Index> patch_map_;
};

} // namespace fvm
} // namespace hvdc

#endif // HVDC_FVM_BOUNDARYCONDITION_HPP
