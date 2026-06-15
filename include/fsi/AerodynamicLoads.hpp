#ifndef HVDC_FSI_AERODYNAMICLOADS_HPP
#define HVDC_FSI_AERODYNAMICLOADS_HPP

#include "common/Types.hpp"
#include "fem/FEModel.hpp"
#include "fvm/Grid.hpp"
#include "fvm/Field.hpp"
#include "fvm/FluidSolver.hpp"
#include <vector>
#include <array>
#include <string>

namespace hvdc {
namespace fsi {

struct AerodynamicConfig {
    Real rho_air = 1.225;
    Real mu_air = 1.789e-5;
    Real gust_amplitude = 0.0;
    Real gust_frequency = 1.0;
    Real wind_direction[3] = {1.0, 0.0, 0.0};
    Real reference_height = 100.0;
    Real terrain_roughness = 0.03;
    bool use_power_law_wind_profile = true;
    Real power_law_exponent = 0.14;
    bool compute_ice_effects = true;
    Real ice_drag_coefficient_increase = 1.5;
    bool use_rain_effects = false;
    Real rain_intensity = 0.0;
    Real drag_coefficient_cylinder = 1.2;
    Real lift_coefficient_cylinder = 0.0;
    Real strouhal_number = 0.2;
    bool compute_vortex_induced_vibration = false;
    bool compute_buffeting = true;
    Real turbulence_intensity = 0.1;
    Real turbulence_length_scale = 100.0;
};

struct FaceAeroLoad {
    Vec3 pressure_force = {0, 0, 0};
    Vec3 viscous_force = {0, 0, 0};
    Vec3 total_force = {0, 0, 0};
    Real pressure = 0.0;
    Real shear_stress_mag = 0.0;
    Real c_p = 0.0;
    Real c_f = 0.0;
    Real reynolds = 0.0;
    Real mach = 0.0;
};

struct ElementAeroLoad {
    Index element_id = -1;
    Vec3 force_per_unit_length = {0, 0, 0};
    Vec3 moment_per_unit_length = {0, 0, 0};
    Real drag_coeff = 0.0;
    Real lift_coeff = 0.0;
    Real reference_diameter = 0.0;
    Vec3 wind_velocity = {0, 0, 0};
    Real relative_speed = 0.0;
};

class AerodynamicLoads {
public:
    AerodynamicLoads() = default;
    AerodynamicLoads(fem::FEModel* struct_model,
                     fvm::FluidSolver* fluid_solver,
                     const AerodynamicConfig& cfg = AerodynamicConfig());
    
    void set_struct_model(fem::FEModel* m) { struct_model_ = m; }
    void set_fluid_solver(fvm::FluidSolver* s) { fluid_solver_ = s; }
    void set_config(const AerodynamicConfig& cfg) { config_ = cfg; }
    
    const AerodynamicConfig& config() const { return config_; }
    
    void compute_face_loads(
        const fvm::FluidState& state,
        const std::vector<std::array<Vec3, 3>>& grad_U,
        const std::vector<Vec3>& grad_p,
        std::vector<FaceAeroLoad>& face_loads,
        const std::vector<Index>& interface_faces) const;
    
    void compute_element_loads(
        const Vector& structural_displacement,
        const Vector& structural_velocity,
        std::vector<ElementAeroLoad>& elem_loads,
        Real time = 0.0) const;
    
    void compute_conductor_aerodynamic_loads(
        const std::vector<Index>& conductor_element_ids,
        const Vector& structural_displacement,
        const Vector& structural_velocity,
        std::vector<Vec3>& nodal_loads,
        Real time = 0.0) const;
    
    void compute_tower_aerodynamic_loads(
        const std::vector<Index>& tower_element_ids,
        const Vector& structural_displacement,
        const Vector& structural_velocity,
        std::vector<Vec3>& nodal_loads,
        Real time = 0.0) const;
    
    void compute_equivalent_nodal_loads(
        const fvm::FluidState& state,
        const std::vector<std::array<Vec3, 3>>& grad_U,
        const std::vector<Vec3>& grad_p,
        Vector& struct_nodal_loads,
        Real time = 0.0) const;
    
    Vec3 compute_wind_velocity_at_point(const Vec3& point, Real time) const;
    
    Real compute_drag_coefficient(
        Real Reynolds,
        Real relative_roughness,
        Real ice_thickness_ratio) const;
    
    Real compute_lift_coefficient(
        Real angle_of_attack,
        Real Reynolds,
        bool galloping_possible = false) const;
    
    void compute_crosswind_forces(
        Real speed, Real diameter, Real length,
        Real angle_of_attack, Real Re,
        Vec3& lift, Vec3& moment) const;
    
    void compute_vortex_shedding_force(
        Real flow_speed, Real diameter, Real length,
        Real time, Vec3& force, Real& vortex_freq) const;
    
    Vec3 compute_buffeting_force(
        const Vec3& mean_velocity, Real time,
        Real length_scale, Real turbulence_intensity,
        Real element_length) const;
    
    void compute_total_forces_and_moments(
        const std::vector<ElementAeroLoad>& loads,
        Vec3& total_force,
        Vec3& total_moment,
        const Vec3& reference_point = {0, 0, 0}) const;
    
    Real compute_RANS_pressure_coefficient(
        Real pressure,
        const Vec3& velocity,
        const Vec3& free_stream_velocity) const;
    
    Real compute_friction_coefficient(
        Real wall_shear_stress,
        const Vec3& free_stream_velocity) const;
    
    void apply_ice_shape_effects(
        Real ice_thickness, Real diameter,
        Real& Cd_multiplier, Real& Cl_offset,
        Real& effective_diameter) const;

protected:
    fem::FEModel* struct_model_ = nullptr;
    fvm::FluidSolver* fluid_solver_ = nullptr;
    AerodynamicConfig config_;
    
    std::vector<Index> conductor_ids_;
    std::vector<Index> tower_ids_;
    
    Real compute_power_law_wind(Real height) const;
    
    Vec3 rotate_2D(const Vec3& vec, Real angle) const;
    
    Real find_angle_of_attack(
        const Vec3& flow_velocity,
        const Vec3& element_direction,
        const Vec3& element_normal) const;
    
    void resolve_force_components(
        Real Cd, Real Cl, Real dynamic_pressure,
        Real diameter, Real length,
        const Vec3& flow_dir, const Vec3& lift_dir,
        Vec3& force) const;
};

} // namespace fsi
} // namespace hvdc

#endif // HVDC_FSI_AERODYNAMICLOADS_HPP
