#ifndef HVDC_FSI_INTERFACEMAPPER_HPP
#define HVDC_FSI_INTERFACEMAPPER_HPP

#include "common/Types.hpp"
#include "fem/FEModel.hpp"
#include "fvm/Grid.hpp"
#include <vector>
#include <array>
#include <memory>

namespace hvdc {
namespace fsi {

enum class MappingMethod : UInt8 {
    NearestNeighbor = 0,
    RBF_Gaussian = 1,
    RBF_MultiQuadric = 2,
    Conservative = 3,
    Projection = 4,
    Interpolating = 5
};

struct InterfacePair {
    Index struct_node_id = -1;
    Index fluid_face_id = -1;
    Index fluid_cell_id = -1;
    Real weight = 1.0;
    Real distance = 0.0;
    Vec3 barycentric = {0, 0, 0};
};

class InterfaceMapper {
public:
    InterfaceMapper() = default;
    InterfaceMapper(fem::FEModel* struct_model,
                    fvm::Grid* fluid_grid,
                    MappingMethod method = MappingMethod::RBF_Gaussian);
    
    void set_structural_model(fem::FEModel* model) { struct_model_ = model; }
    void set_fluid_grid(fvm::Grid* grid) { fluid_grid_ = grid; }
    void set_mapping_method(MappingMethod method) { method_ = method; }
    
    void set_structural_interface_nodes(const std::vector<Index>& node_ids);
    void set_fluid_interface_faces(const std::vector<Index>& face_ids);
    void set_fluid_interface_patch(Index boundary_patch_id);
    
    void build_mapping(Real search_radius_factor = 3.0);
    void build_conservative_mapping();
    void build_rbf_mapping(Real support_radius, Real shape_param = 1.0);
    
    Index num_struct_interface_nodes() const { return static_cast<Index>(struct_interface_nodes_.size()); }
    Index num_fluid_interface_faces() const { return static_cast<Index>(fluid_interface_faces_.size()); }
    
    const std::vector<Index>& struct_interface_nodes() const { return struct_interface_nodes_; }
    const std::vector<Index>& fluid_interface_faces() const { return fluid_interface_faces_; }
    
    const std::vector<std::vector<InterfacePair>>& struct_to_fluid_map() const { return struct_to_fluid_; }
    const std::vector<std::vector<InterfacePair>>& fluid_to_struct_map() const { return fluid_to_struct_; }
    
    void map_struct_displacement_to_fluid(
        const Vector& struct_disp,
        std::vector<Vec3>& fluid_face_displacement) const;
    
    void map_struct_velocity_to_fluid(
        const Vector& struct_vel,
        std::vector<Vec3>& fluid_face_velocity) const;
    
    void map_fluid_force_to_struct(
        const std::vector<Vec3>& fluid_face_force,
        Vector& struct_node_force) const;
    
    void map_fluid_pressure_to_struct(
        const std::vector<Real>& fluid_pressure,
        Vector& struct_node_force) const;
    
    Real compute_conservation_error(
        const std::vector<Vec3>& face_forces,
        const Vector& struct_forces) const;
    
    void deform_fluid_mesh(
        const Vector& struct_disp,
        std::vector<Vec3>& new_node_coords,
        Real spring_stiffness = 1.0) const;
    
    void set_rbf_shape_parameter(Real eps) { rbf_eps_ = eps; }
    void set_mapping_tol(Real tol) { mapping_tol_ = tol; }

protected:
    fem::FEModel* struct_model_ = nullptr;
    fvm::Grid* fluid_grid_ = nullptr;
    MappingMethod method_ = MappingMethod::RBF_Gaussian;
    
    std::vector<Index> struct_interface_nodes_;
    std::vector<Index> fluid_interface_faces_;
    
    std::vector<std::vector<InterfacePair>> struct_to_fluid_;
    std::vector<std::vector<InterfacePair>> fluid_to_struct_;
    
    std::vector<Vec3> struct_positions_;
    std::vector<Vec3> fluid_face_positions_;
    std::vector<Vec3> fluid_face_normals_;
    std::vector<Real> fluid_face_areas_;
    
    Real rbf_eps_ = 1.0;
    Real mapping_tol_ = 1.0e-6;
    
    Real rbf_gaussian(Real r, Real eps) const;
    Real rbf_multiquadric(Real r, Real eps) const;
    Real rbf_tps(Real r) const;
    
    Real compute_influence_radius() const;
    
    void normalize_weights(std::vector<InterfacePair>& pairs) const;
    
    void compute_nearest_pairs(
        const std::vector<Vec3>& source,
        const std::vector<Vec3>& target,
        std::vector<std::vector<InterfacePair>>& mapping,
        Real max_radius,
        Index max_neighbors = 20) const;
};

} // namespace fsi
} // namespace hvdc

#endif // HVDC_FSI_INTERFACEMAPPER_HPP
