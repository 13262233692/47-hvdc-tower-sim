#ifndef HVDC_FEM_TRUSS_ELEMENT_NL_HPP
#define HVDC_FEM_TRUSS_ELEMENT_NL_HPP

#include "fem/Element.hpp"
#include "fem/Material.hpp"
#include "fem/Section.hpp"

namespace hvdc {
namespace fem {

class TrussElementNL : public Element2N {
public:
    TrussElementNL();
    TrussElementNL(Index id, Node* n1, Node* n2,
                   const Material* mat, const TrussSection* sec);
    
    void set_material(const Material* mat) { material_ = mat; }
    void set_section(const TrussSection* sec) { section_ = sec; }
    const Material* material() const { return material_; }
    const TrussSection* section() const { return section_; }
    
    void initialize();
    
    void stiffness_matrix(Mat12x12& K_out) const override;
    void tangent_stiffness_matrix(
        Mat12x12& K_T_out, const Vec12& displacement,
        bool include_geometric = true) const override;
    void mass_matrix(Mat12x12& M_out) const override;
    
    Vec12 internal_forces(const Vec12& displacement) const override;
    Vec12 equivalent_nodal_loads() const override;
    Vec12 gravity_loads(Real g = GRAVITY) const override;
    
    Real axial_force(const Vec12& displacement) const;
    Real axial_stress(const Vec12& displacement) const;
    Real axial_strain(const Vec12& displacement) const;
    
    void compute_axial_response(const Vec12& disp,
                                 Real& strain, Real& stress, Real& force) const;

protected:
    const Material* material_ = nullptr;
    const TrussSection* section_ = nullptr;
    
    void compute_stiffness_components(
        const Vec12& displacement,
        Mat3x3& T, Vec3& dir, Real& L_cur,
        Real& E_A_L, Real& N_force) const;
    
    void build_stiffness_from_dir(
        Mat12x12& K, const Vec3& e, Real L,
        Real EA_L, Real N_L, bool use_geo) const;
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_TRUSS_ELEMENT_NL_HPP
