#ifndef HVDC_FSI_LOADTRANSFER_HPP
#define HVDC_FSI_LOADTRANSFER_HPP

#include "common/Types.hpp"
#include "fem/FEModel.hpp"
#include "fvm/Grid.hpp"
#include "fvm/Field.hpp"
#include "fsi/AerodynamicLoads.hpp"
#include "fsi/InterfaceMapper.hpp"
#include <vector>
#include <array>
#include <memory>

namespace hvdc {
namespace fsi {

enum class LoadTransferScheme : UInt8 {
    DirectNearest = 0,
    ConservativeInterpolation = 1,
    LumpedNodal = 2,
    ConsistentNodal = 3,
    Mortar = 4
};

enum class DisplacementTransferScheme : UInt8 {
    DirectNearest = 0,
    InterpolatingSpline = 1,
    RadialBasis = 2,
    Transfinite = 3
};

class LoadTransfer {
public:
    LoadTransfer() = default;
    LoadTransfer(fem::FEModel* struct_model,
                 fvm::Grid* fluid_grid,
                 InterfaceMapper* mapper,
                 AerodynamicLoads* aero = nullptr);
    
    void set_struct_model(fem::FEModel* m) { struct_model_ = m; }
    void set_fluid_grid(fvm::Grid* g) { fluid_grid_ = g; }
    void set_mapper(InterfaceMapper* m) { mapper_ = m; }
    void set_aerodynamic_loads(AerodynamicLoads* a) { aero_ = a; }
    
    void set_load_scheme(LoadTransferScheme s) { load_scheme_ = s; }
    void set_disp_scheme(DisplacementTransferScheme s) { disp_scheme_ = s; }
    void set_conservation_tol(Real t) { conservation_tol_ = t; }
    
    void transfer_fluid_to_structural(
        const std::vector<Vec3>& fluid_face_forces,
        Vector& structural_node_forces) const;
    
    void transfer_fluid_pressure_to_structural(
        const Vector& fluid_pressure_cell,
        Vector& structural_node_forces) const;
    
    void transfer_structural_to_fluid(
        const Vector& structural_displacement,
        std::vector<Vec3>& fluid_node_displacement) const;
    
    void transfer_structural_velocity_to_fluid(
        const Vector& structural_velocity,
        std::vector<Vec3>& fluid_face_velocity) const;
    
    void transfer_face_aero_to_nodal(
        const std::vector<FaceAeroLoad>& face_loads,
        Vector& structural_node_forces) const;
    
    void transfer_element_aero_to_nodal(
        const std::vector<ElementAeroLoad>& elem_loads,
        Vector& structural_node_forces) const;
    
    void compute_consistent_nodal_forces(
        const std::vector<Vec3>& distributed_load_per_length,
        Vector& structural_node_forces) const;
    
    void compute_lumped_nodal_forces(
        const std::vector<Vec3>& distributed_load_per_length,
        Vector& structural_node_forces) const;
    
    void compute_mortar_integration_masses(
        const std::vector<Index>& interface_faces,
        std::vector<Real>& mass_source,
        std::vector<Real>& mass_target) const;
    
    void redistribute_forces_for_conservation(
        const Vector& input_forces,
        Vector& output_forces,
        Real desired_total_force_per_direction[3]) const;
    
    void transfer_temperature_from_structural_to_fluid(
        const Vector& structural_temp,
        Vector& fluid_wall_temperature) const;
    
    void transfer_heat_flux_from_fluid_to_structural(
        const Vector& fluid_heat_flux,
        Vector& structural_heat_load) const;
    
    Real check_force_conservation(
        const std::vector<Vec3>& fluid_forces,
        const Vector& struct_forces) const;
    
    Real check_energy_conservation(
        const Vector& struct_vel,
        const Vector& struct_force,
        const std::vector<Vec3>& fluid_wall_vel,
        const std::vector<Vec3>& fluid_force) const;

protected:
    fem::FEModel* struct_model_ = nullptr;
    fvm::Grid* fluid_grid_ = nullptr;
    InterfaceMapper* mapper_ = nullptr;
    AerodynamicLoads* aero_ = nullptr;
    
    LoadTransferScheme load_scheme_ = LoadTransferScheme::ConsistentNodal;
    DisplacementTransferScheme disp_scheme_ = DisplacementTransferScheme::RadialBasis;
    Real conservation_tol_ = 1.0e-4;
    
    void transfer_direct_nearest(
        const std::vector<Vec3>& source,
        Vector& target,
        bool force_to_disp) const;
    
    void transfer_conservative(
        const std::vector<Vec3>& source,
        Vector& target,
        const std::vector<Real>& source_areas,
        const std::vector<Real>& target_areas) const;
    
    void transfer_lumped(
        const std::vector<Vec3>& distributed_per_len,
        Vector& target,
        bool consistent) const;
};

} // namespace fsi
} // namespace hvdc

#endif // HVDC_FSI_LOADTRANSFER_HPP
