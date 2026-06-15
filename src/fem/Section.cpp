#include "fem/Section.hpp"
#include "fem/Material.hpp"

namespace hvdc {
namespace fem {

BeamSection BeamSection::CircularTube(Real D, Real t) {
    BeamSection s;
    s.name = "CircularTube_D" + std::to_string(D) + "_t" + std::to_string(t);
    s.h = s.b = D;
    s.t_w = s.t_f = t;
    Real r = D / 2.0;
    Real r_i = r - t;
    s.A = PI * (r * r - r_i * r_i);
    s.Iy = s.Iz = PI * (r * r * r * r - r_i * r_i * r_i * r_i) / 4.0;
    s.Ix = s.Iy + s.Iz;
    s.J = PI * (r * r * r * r - r_i * r_i * r_i * r_i) / 2.0;
    return s;
}

BeamSection BeamSection::RectangularTube(Real h, Real b, Real t) {
    BeamSection s;
    s.name = "RectTube_h" + std::to_string(h) + "_b" + std::to_string(b) + "_t" + std::to_string(t);
    s.h = h; s.b = b;
    s.t_w = s.t_f = t;
    Real hi = h - 2.0 * t;
    Real bi = b - 2.0 * t;
    s.A = h * b - hi * bi;
    s.Iy = (b * h * h * h - bi * hi * hi * hi) / 12.0;
    s.Iz = (h * b * b * b - hi * bi * bi * bi) / 12.0;
    s.Ix = s.Iy + s.Iz;
    Real Am = (h - t) * (b - t);
    s.J = 4.0 * Am * Am * t / (2.0 * (h + b - 2.0 * t));
    return s;
}

BeamSection BeamSection::I_Section(Real h, Real b, Real t_w, Real t_f) {
    BeamSection s;
    s.name = "I_h" + std::to_string(h) + "_b" + std::to_string(b) + 
             "_tw" + std::to_string(t_w) + "_tf" + std::to_string(t_f);
    s.h = h; s.b = b; s.t_w = t_w; s.t_f = t_f;
    s.A = b * t_f * 2.0 + t_w * (h - 2.0 * t_f);
    Real h_w = h - 2.0 * t_f;
    s.Iy = t_w * h_w * h_w * h_w / 12.0 + 2.0 * (b * t_f * t_f * t_f / 12.0 + 
           b * t_f * (h / 2.0 - t_f / 2.0) * (h / 2.0 - t_f / 2.0));
    s.Iz = 2.0 * t_f * b * b * b / 12.0 + h_w * t_w * t_w * t_w / 12.0;
    s.Ix = s.Iy + s.Iz;
    s.J = 1.0 / 3.0 * (2.0 * b * t_f * t_f * t_f + h_w * t_w * t_w * t_w);
    return s;
}

BeamSection BeamSection::L_Angle(Real leg_a, Real leg_b, Real t) {
    BeamSection s;
    s.name = "Angle_" + std::to_string(leg_a) + "x" + std::to_string(leg_b) + "x" + std::to_string(t);
    s.h = leg_a; s.b = leg_b; s.t_w = s.t_f = t;
    s.A = t * (leg_a + leg_b - t);
    Real cx = t * (2.0 * leg_b - t) / (2.0 * (leg_a + leg_b - t));
    Real cy = t * (2.0 * leg_a - t) / (2.0 * (leg_a + leg_b - t));
    s.Iy = leg_a * t * t * t / 12.0 + leg_a * t * (cy - t / 2.0) * (cy - t / 2.0)
         + t * (leg_b - t) * (leg_b - t) * (leg_b - t) / 12.0
         + t * (leg_b - t) * (t / 2.0 + (leg_b - t) / 2.0 - cx) * 
           (t / 2.0 + (leg_b - t) / 2.0 - cx);
    s.Iz = t * leg_b * t * t / 12.0 + t * leg_b * (cx - t / 2.0) * (cx - t / 2.0)
         + (leg_a - t) * t * t * t / 12.0
         + (leg_a - t) * t * (t + (leg_a - t) / 2.0 - cy) * 
           (t + (leg_a - t) / 2.0 - cy);
    s.Ix = s.Iy + s.Iz;
    s.J = 1.0 / 3.0 * t * t * t * (leg_a + leg_b - t);
    return s;
}

BeamSection BeamSection::Conductor(Real D, Real A_strand) {
    BeamSection s;
    s.name = "Conductor_D" + std::to_string(D);
    s.h = s.b = D;
    s.A = A_strand;
    s.Iy = s.Iz = PI * std::pow(D / 2.0, 4) / 4.0 * 0.4;
    s.Ix = s.Iy + s.Iz;
    s.J = PI * std::pow(D / 2.0, 4) / 2.0 * 0.3;
    return s;
}

BeamSection BeamSection::EqualAngle(Real side, Real thickness) {
    return L_Angle(side, side, thickness);
}

void BeamSection::add_ice_coating(Real ice_thickness) {
    thickness_ice = ice_thickness;
    Real D_outer = std::max(h, b);
    D_ice = D_outer + 2.0 * ice_thickness;
    A_ice = PI * (std::pow(D_ice / 2.0, 2) - std::pow(D_outer / 2.0, 2));
}

Real BeamSection::total_mass_per_length(Real rho_material) const {
    return rho_material * A + rho_ice * A_ice;
}

Real BeamSection::wind_drag_diameter() const {
    return D_ice > 0.0 ? D_ice : std::max(h, b);
}

Real BeamSection::aerodynamic_area_per_length() const {
    return wind_drag_diameter();
}

BeamSection BeamSection::circular_pipe(Real D, Real t, std::shared_ptr<Material> mat_ptr) {
    auto s = CircularTube(D, t);
    if (mat_ptr) s.name = mat_ptr->name + "_Pipe_" + s.name;
    return s;
}

BeamSection BeamSection::i_beam(Real h, Real b, Real t_w, Real t_f, std::shared_ptr<Material> mat_ptr) {
    auto s = I_Section(h, b, t_w, t_f);
    if (mat_ptr) s.name = mat_ptr->name + "_IBeam_" + s.name;
    return s;
}

} // namespace fem
} // namespace hvdc
