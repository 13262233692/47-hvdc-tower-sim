#ifndef HVDC_FEM_BEAM_ELEMENT_NL_HPP
#define HVDC_FEM_BEAM_ELEMENT_NL_HPP

#include "fem/Element.hpp"
#include "fem/Material.hpp"
#include "fem/Section.hpp"

namespace hvdc {
namespace fem {

class BeamElementNL : public Element2N {
public:
    BeamElementNL();
    BeamElementNL(Index id, Node* n1, Node* n2,
                  const Material* mat, const BeamSection* sec);
    
    void set_material(const Material* mat) { material_ = mat; }
    void set_section(const BeamSection* sec) { section_ = sec; }
    const Material* material() const { return material_; }
    const BeamSection* section() const { return section_; }
    
    void initialize();
    
    void stiffness_matrix(Mat12x12& K_out) const override;
    void tangent_stiffness_matrix(
        Mat12x12& K_T_out, const Vec12& displacement,
        bool include_geometric = true) const override;
    void mass_matrix(Mat12x12& M_out) const override;
    
    Vec12 internal_forces(const Vec12& displacement) const override;
    Vec12 equivalent_nodal_loads() const override;
    Vec12 gravity_loads(Real g = GRAVITY) const override;
    
    void set_distributed_load(const Vec3& q_local) {
        distributed_load_ = q_local;
    }
    
    void compute_local_forces(const Vec12& disp_local, Vec12& f_int_local,
                               Real& N_axial, Real& Vy_shear, Real& Vz_shear,
                               Real& T_torsion, Real& My_bend, Real& Mz_bend) const;
    
    void compute_material_stiffness_local(Mat12x12& K_mat_local) const;
    void compute_geometric_stiffness_local(
        Mat12x12& K_geo_local, Real N, Real Vy, Real Vz,
        Real T, Real My, Real Mz) const;

    void extract_node_rotations(const Vec12& disp_global,
                                 Mat3x3& R1, Mat3x3& R2) const;
    
    void compute_current_rotation(const Mat3x3& R0, const Vec3& rot_incr,
                                   Mat3x3& R_current) const;

protected:
    const Material* material_ = nullptr;
    const BeamSection* section_ = nullptr;
    Vec3 distributed_load_ = {0.0, 0.0, 0.0};
    
    Mat3x3 T0_;
    Vec3 dx0_;
    
    void compute_axial_force_and_moments(
        const Vec12& d_local,
        Real& N, Real& Vy, Real& Vz, Real& T, Real& My, Real& Mz) const;
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_BEAM_ELEMENT_NL_HPP
