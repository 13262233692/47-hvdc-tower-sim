#ifndef HVDC_FEM_SECTION_HPP
#define HVDC_FEM_SECTION_HPP

#include "common/Types.hpp"
#include <string>
#include <cmath>
#include <memory>

namespace hvdc {
namespace fem {

class Material;

struct BeamSection {
    Index id = -1;
    std::string name = "unnamed";
    
    Real A = 0.0;
    Real Ix = 0.0;
    Real Iy = 0.0;
    Real Iz = 0.0;
    Real J = 0.0;
    Real Iw = 0.0;
    Real h = 0.0;
    Real b = 0.0;
    Real t_w = 0.0;
    Real t_f = 0.0;
    
    Real A_ice = 0.0;
    Real D_ice = 0.0;
    Real thickness_ice = 0.0;
    Real rho_ice = 917.0;
    
    static BeamSection CircularTube(Real D, Real t);
    static BeamSection circular_pipe(Real D, Real t, std::shared_ptr<Material> mat_ptr = nullptr);
    static BeamSection RectangularTube(Real h, Real b, Real t);
    static BeamSection I_Section(Real h, Real b, Real t_w, Real t_f);
    static BeamSection i_beam(Real h, Real b, Real t_w, Real t_f, std::shared_ptr<Material> mat_ptr = nullptr);
    static BeamSection L_Angle(Real leg_a, Real leg_b, Real t);
    static BeamSection Conductor(Real D, Real A_strand);
    static BeamSection EqualAngle(Real side, Real thickness);
    
    void add_ice_coating(Real ice_thickness);
    Real total_mass_per_length(Real rho_material) const;
    Real wind_drag_diameter() const;
    Real aerodynamic_area_per_length() const;
    Real ice_area() const { return A_ice; }
};

struct TrussSection {
    Index id = -1;
    std::string name = "unnamed";
    
    Real A = 0.0;
    Real D = 0.0;
    
    Real A_ice = 0.0;
    Real D_ice = 0.0;
    Real thickness_ice = 0.0;
    
    static TrussSection CircularBar(Real D) {
        TrussSection s;
        s.D = D;
        s.A = PI * D * D / 4.0;
        s.name = "CircularBar_D" + std::to_string(D);
        return s;
    }
    
    static TrussSection Angle(Real leg, Real t) {
        TrussSection s;
        s.A = t * (2.0 * leg - t);
        s.D = leg * std::sqrt(2.0);
        s.name = "Angle_" + std::to_string(leg) + "x" + std::to_string(t);
        return s;
    }
    
    void add_ice_coating(Real ice_thickness) {
        thickness_ice = ice_thickness;
        if (D > 0.0) {
            D_ice = D + 2.0 * ice_thickness;
            A_ice = PI * (D_ice * D_ice - D * D) / 4.0;
        } else {
            A_ice = A * 0.3;
            D_ice = 0.05 + 2.0 * ice_thickness;
        }
    }
    
    Real total_A() const { return A + A_ice; }
    Real ice_area() const { return A_ice; }
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_SECTION_HPP
