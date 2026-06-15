#ifndef HVDC_FEM_MATERIAL_HPP
#define HVDC_FEM_MATERIAL_HPP

#include "common/Types.hpp"
#include <string>

namespace hvdc {
namespace fem {

struct Material {
    Index id = -1;
    std::string name = "unnamed";
    
    Real E = 2.1e11;
    Real nu = 0.3;
    Real rho = 7850.0;
    Real alpha = 1.2e-5;
    
    Real G = 8.077e10;
    Real lambda = 1.212e11;
    Real mu = 8.077e10;
    
    Real sigma_y = 2.5e8;
    Real E_tangent = 1.0e9;
    
    Real yield_stress = 3.55e8;
    Real ultimate_stress = 5.10e8;
    Real elongation = 0.20;
    
    MaterialModel model = MaterialModel::LinearElastic;
    
    void compute_lame_from_E_nu();
    
    Real bulk_modulus() const {
        return E / (3.0 * (1.0 - 2.0 * nu));
    }
    
    Real shear_modulus() const {
        return E / (2.0 * (1.0 + nu));
    }
    
    static Material Q345_Steel();
    static Material Q420_Steel();
    static Material ACSR_Conductor();
};

struct SteelQ345 : public Material {
    SteelQ345() {
        name = "Q345 Steel";
        E = 2.06e11;
        nu = 0.3;
        rho = 7850.0;
        alpha = 1.2e-5;
        sigma_y = 3.45e8;
        yield_stress = 3.45e8;
        ultimate_stress = 4.70e8;
        elongation = 0.21;
        compute_lame_from_E_nu();
    }
};

struct SteelQ420 : public Material {
    SteelQ420() {
        name = "Q420 Steel";
        E = 2.06e11;
        nu = 0.3;
        rho = 7850.0;
        alpha = 1.2e-5;
        sigma_y = 4.20e8;
        yield_stress = 4.20e8;
        ultimate_stress = 5.40e8;
        elongation = 0.19;
        compute_lame_from_E_nu();
    }
};

struct AluminumConductor : public Material {
    AluminumConductor() {
        name = "Aluminum Conductor (ACSR)";
        E = 6.30e10;
        nu = 0.33;
        rho = 2710.0;
        alpha = 2.3e-5;
        sigma_y = 1.65e8;
        yield_stress = 1.65e8;
        ultimate_stress = 2.75e8;
        elongation = 0.015;
        compute_lame_from_E_nu();
    }
};

inline Material Material::Q345_Steel() { return SteelQ345(); }
inline Material Material::Q420_Steel() { return SteelQ420(); }
inline Material Material::ACSR_Conductor() { return AluminumConductor(); }

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_MATERIAL_HPP
