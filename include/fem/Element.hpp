#ifndef HVDC_FEM_ELEMENT_HPP
#define HVDC_FEM_ELEMENT_HPP

#include "common/Types.hpp"
#include "fem/Node.hpp"
#include "fem/Material.hpp"
#include "fem/Section.hpp"
#include <array>
#include <memory>

namespace hvdc {
namespace fem {

class Element {
public:
    Element() = default;
    virtual ~Element() = default;
    
    Index id() const { return id_; }
    void set_id(Index id) { id_ = id; }
    
    ElementType type() const { return type_; }
    
    Index num_nodes() const { return nnodes_; }
    Index num_dofs() const { return ndofs_; }
    
    virtual Node* node(Index i) = 0;
    virtual const Node* node(Index i) const = 0;
    virtual void set_node(Index i, Node* n) = 0;
    
    virtual Index local_dof_id(Index node_idx, Index dof_idx) const = 0;
    virtual Index global_dof_id(Index node_idx, Index dof_idx) const = 0;
    
    virtual void stiffness_matrix(Mat12x12& K_out) const = 0;
    virtual void tangent_stiffness_matrix(
        Mat12x12& K_T_out, const Vec12& displacement, 
        bool include_geometric = true) const = 0;
    virtual void mass_matrix(Mat12x12& M_out) const = 0;
    
    virtual Vec12 internal_forces(const Vec12& displacement) const = 0;
    virtual Vec12 equivalent_nodal_loads() const = 0;
    virtual Vec12 gravity_loads(Real g = GRAVITY) const = 0;
    
    virtual void apply_temperature(Real delta_T) { delta_T_ = delta_T; }
    Real temperature_change() const { return delta_T_; }
    
    int partition() const { return partition_; }
    void set_partition(int p) { partition_ = p; }
    
    Real length_initial() const { return L0_; }
    virtual Real length_current(const Vec12& disp) const = 0;
    
    virtual Vec3 midpoint() const;
    virtual Vec3 centroid(const Vec12& disp) const;
    
    Index material_id() const { return material_id_; }
    Index section_id() const { return section_id_; }
    void set_material_id(Index id) { material_id_ = id; }
    void set_section_id(Index id) { section_id_ = id; }

protected:
    Index id_ = -1;
    ElementType type_ = ElementType::Unknown;
    Index nnodes_ = 0;
    Index ndofs_ = 0;
    int partition_ = 0;
    Real L0_ = 0.0;
    Real delta_T_ = 0.0;
    Index material_id_ = 0;
    Index section_id_ = 0;
};

class Element2N : public Element {
public:
    Element2N() { nnodes_ = 2; ndofs_ = 12; }
    virtual ~Element2N() = default;
    
    Node* node(Index i) override { return nodes_[i]; }
    const Node* node(Index i) const override { return nodes_[i]; }
    void set_node(Index i, Node* n) override { nodes_[i] = n; }
    
    Index local_dof_id(Index node_idx, Index dof_idx) const override {
        return node_idx * 6 + dof_idx;
    }
    Index global_dof_id(Index node_idx, Index dof_idx) const override {
        return nodes_[node_idx]->dof_start() + dof_idx;
    }
    
    Real length_current(const Vec12& disp) const override;
    
    void compute_transformation(Mat3x3& T_L_to_G, Vec3& dx_global,
                                 const Vec12& disp) const;
    void compute_transformation_initial(Mat3x3& T_L_to_G, Vec3& dx_global) const;
    
    void assemble_12x12_from_6x6_blocks(
        Mat12x12& K,
        const Mat6x6& K11, const Mat6x6& K12,
        const Mat6x6& K21, const Mat6x6& K22) const;
    
    void transform_12x12(
        Mat12x12& K_global, const Mat12x12& K_local,
        const Mat3x3& T) const;
    
    void transform_12_vec(
        Vec12& v_global, const Vec12& v_local,
        const Mat3x3& T) const;

protected:
    std::array<Node*, 2> nodes_ = {};
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_ELEMENT_HPP
