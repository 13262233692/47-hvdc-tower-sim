#ifndef HVDC_FEM_FEMODEL_HPP
#define HVDC_FEM_FEMODEL_HPP

#include "common/Types.hpp"
#include "common/Vector.hpp"
#include "common/SparseMatrix.hpp"
#include "fem/Node.hpp"
#include "fem/Element.hpp"
#include "fem/BeamElementNL.hpp"
#include "fem/TrussElementNL.hpp"
#include "fem/ConductorElement.hpp"
#include "fem/Material.hpp"
#include "fem/Section.hpp"

#include <vector>
#include <memory>
#include <unordered_map>
#include <functional>

namespace hvdc {
namespace fem {

class FEModel {
public:
    FEModel() = default;
    ~FEModel() = default;
    
    Node* add_node(Index id, Real x, Real y, Real z);
    Node* add_node(const Vec3& coords);
    Node* get_node(Index id);
    const Node* get_node(Index id) const;
    Node* node(Index i) { return (i < static_cast<Index>(nodes_.size())) ? nodes_[i] : nullptr; }
    const Node* node(Index i) const { return (i < static_cast<Index>(nodes_.size())) ? nodes_[i] : nullptr; }
    Index num_nodes() const { return static_cast<Index>(nodes_.size()); }
    
    Index add_material(const Material& mat);
    Index add_material(std::shared_ptr<Material> mat_ptr);
    Material* get_material(Index id);
    const Material* get_material(Index id) const;
    
    Index add_beam_section(const BeamSection& sec);
    Index add_section(std::shared_ptr<BeamSection> sec_ptr);
    Index add_truss_section(std::shared_ptr<TrussSection> sec_ptr);
    BeamSection* get_beam_section(Index id);
    const BeamSection* get_beam_section(Index id) const;
    
    Index add_truss_section(const TrussSection& sec);
    TrussSection* get_truss_section(Index id);
    const TrussSection* get_truss_section(Index id) const;
    
    Element* add_element(std::unique_ptr<Element> elem);
    
    BeamElementNL* add_beam_element(
        Index node1_id, Index node2_id,
        Index material_id, Index section_id);
    
    TrussElementNL* add_truss_element(
        Index node1_id, Index node2_id,
        Index material_id, Index section_id);
    
    ConductorElement* add_conductor_element(
        Index node1_id, Index node2_id,
        Index material_id, Index section_id,
        Real initial_tension = 0.0);
    
    Index num_elements() const { return static_cast<Index>(elements_.size()); }
    Element* element(Index i) { return elements_[i].get(); }
    const Element* element(Index i) const { return elements_[i].get(); }
    
    void setup_dofs(Index dofs_per_node = 6);
    Index total_dofs() const { return total_dofs_; }
    Index num_constrained_dofs() const { return num_constrained_dofs_; }
    Index num_free_dofs() const { return total_dofs_ - num_constrained_dofs_; }
    
    const IndexVec& free_dofs() const { return free_dofs_; }
    const IndexVec& constrained_dofs() const { return constrained_dofs_; }
    
    void apply_fixed_bc(Index node_id, const std::vector<Index>& dofs = {0,1,2,3,4,5});
    void apply_displacement_bc(Index node_id, Index dof, Real value);
    void apply_pin_bc(Index node_id) { apply_fixed_bc(node_id, {0,1,2}); }
    void set_node_bc_all(Index node_id, BoundaryType type);
    void set_node_bc(Index node_id, BoundaryType type,
                     bool tx, bool ty, bool tz, bool rx, bool ry, bool rz);
    
    void gather_displacement_vec(Vec& u_global) const;
    void scatter_displacement_vec(const Vec& u_global);
    
    void apply_body_force(Vector& F_ext, Real g = GRAVITY) const;
    void apply_ice_load(Vector& F_ext, Real rho_ice = 917.0) const;
    
    void apply_ice_loads(Real ice_thickness, Real rho_ice = 917.0);
    void apply_temperature_loads(Real delta_T_conductor, Real delta_T_tower = 0.0);
    
    void mark_fsi_interface_nodes(const std::vector<Index>& node_ids);
    const std::vector<Node*>& interface_nodes() const { return interface_nodes_; }
    std::vector<Index> interface_node_dofs() const;
    
    void partition_domain(int num_procs);
    int element_partition(Index elem_id) const;
    
    void compute_stats(Real& max_sag, Real& total_mass,
                       Vec3& center_of_mass) const;
    
    std::vector<Index> conductor_element_ids() const;
    std::vector<Index> tower_element_ids() const;

private:
    std::vector<std::unique_ptr<Node>> nodes_storage_;
    std::unordered_map<Index, Node*> node_id_map_;
    std::vector<Node*> nodes_;
    
    std::vector<Material> materials_;
    std::vector<BeamSection> beam_sections_;
    std::vector<TrussSection> truss_sections_;
    std::unordered_map<Index, Index> material_id_map_;
    std::unordered_map<Index, Index> beam_sec_id_map_;
    std::unordered_map<Index, Index> truss_sec_id_map_;
    
    std::vector<std::unique_ptr<Element>> elements_;
    
    Index total_dofs_ = 0;
    Index num_constrained_dofs_ = 0;
    IndexVec free_dofs_;
    IndexVec constrained_dofs_;
    
    std::vector<Node*> interface_nodes_;
    
    Index next_node_id_ = 0;
    Index next_material_id_ = 0;
    Index next_beam_sec_id_ = 0;
    Index next_truss_sec_id_ = 0;
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_FEMODEL_HPP
