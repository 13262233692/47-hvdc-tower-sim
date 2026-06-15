#ifndef HVDC_FEM_CONDUCTOR_ELEMENT_HPP
#define HVDC_FEM_CONDUCTOR_ELEMENT_HPP

#include "fem/Element.hpp"
#include "fem/Material.hpp"
#include "fem/Section.hpp"

namespace hvdc {
namespace fem {

class ConductorElement : public Element2N {
public:
    ConductorElement();
    ConductorElement(Index id, Node* n1, Node* n2,
                     const Material* mat, const BeamSection* sec,
                     Real initial_tension = 0.0, Real span_sag = 0.0);
    
    void set_material(const Material* mat) { material_ = mat; }
    void set_section(const BeamSection* sec) { section_ = sec; }
    const Material* material() const { return material_; }
    const BeamSection* section() const { return section_; }
    
    void set_initial_tension(Real T0) { T0_ = T0; }
    Real initial_tension() const { return T0_; }
    
    void set_catenary_sag(Real sag) { sag_ = sag; }
    Real catenary_sag() const { return sag_; }
    
    void initialize();
    
    void stiffness_matrix(Mat12x12& K_out) const override;
    void tangent_stiffness_matrix(
        Mat12x12& K_T_out, const Vec12& displacement,
        bool include_geometric = true) const override;
    void mass_matrix(Mat12x12& M_out) const override;
    
    Vec12 internal_forces(const Vec12& displacement) const override;
    Vec12 equivalent_nodal_loads() const override;
    Vec12 gravity_loads(Real g = GRAVITY) const override;
    
    Real current_tension(const Vec12& displacement) const;
    Real axial_strain(const Vec12& displacement) const;
    Real axial_stress(const Vec12& displacement) const;
    Real catenary_length() const;
    void catenary_shape(Index n_points,
                        std::vector<Vec3>& points,
                        std::vector<Real>& tensions) const;
    
    Real dynamic_amplification_factor() const { return DAF_; }
    void set_DAF(Real daf) { DAF_ = daf; }

protected:
    const Material* material_ = nullptr;
    const BeamSection* section_ = nullptr;
    
    Real T0_ = 0.0;
    Real sag_ = 0.0;
    Real DAF_ = 1.0;
    
    void compute_truss_equivalent_stiffness(
        Mat12x12& K, const Vec3& e,
        Real EA_L_eff, Real T_eff_over_L, bool use_geo) const;
    
    Real compute_effective_modulus(Real current_tension) const;
    Real compute_Lame_parameter(Real stress_axial) const;
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_CONDUCTOR_ELEMENT_HPP
