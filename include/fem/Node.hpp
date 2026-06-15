#ifndef HVDC_FEM_NODE_HPP
#define HVDC_FEM_NODE_HPP

#include "common/Types.hpp"
#include <array>

namespace hvdc {
namespace fem {

class Node {
public:
    Node() = default;
    Node(Index id, Real x, Real y, Real z);
    Node(Index id, const Vec3& coords);
    
    Index id() const { return id_; }
    
    const Vec3& coords() const { return coords_; }
    Real x() const { return coords_[0]; }
    Real y() const { return coords_[1]; }
    Real z() const { return coords_[2]; }
    
    void set_coords(const Vec3& c) { coords_ = c; }
    void set_coords(Real x, Real y, Real z) { coords_ = {x, y, z}; }
    
    Vec3 displaced_coords(const Vec& displacement) const;
    Vec3 displaced_coords(const Vec3& disp) const;
    Vec3 displaced_coords(const class Vector& displacement) const;
    
    Index dof_start() const { return dof_start_; }
    void set_dof_start(Index start) { dof_start_ = start; }
    
    Index num_dofs() const { return ndofs_; }
    void set_num_dofs(Index n) { ndofs_ = n; }
    
    bool is_dof_active(Index local_dof) const {
        return bc_types_[local_dof] == BoundaryType::Free;
    }
    
    void set_bc_type(Index local_dof, BoundaryType type) {
        bc_types_[local_dof] = type;
    }
    
    void set_bc_type_all(BoundaryType type) {
        for (auto& t : bc_types_) t = type;
    }
    
    BoundaryType bc_type() const { return bc_types_[0]; }
    BoundaryType bc_type(Index local_dof) const {
        return bc_types_[local_dof];
    }
    
    Real bc_value(Index local_dof) const {
        return bc_values_[local_dof];
    }
    
    void set_bc_value(Index local_dof, Real val) {
        bc_values_[local_dof] = val;
    }
    
    void fix_all() {
        for (auto& t : bc_types_) t = BoundaryType::Fixed;
    }
    
    void free_all() {
        for (auto& t : bc_types_) t = BoundaryType::Free;
    }
    
    bool is_on_interface() const { return on_interface_; }
    void set_on_interface(bool b) { on_interface_ = b; }
    
    Index interface_id() const { return interface_id_; }
    void set_interface_id(Index id) { interface_id_ = id; }
    
    int partition() const { return partition_; }
    void set_partition(int p) { partition_ = p; }

private:
    Index id_ = -1;
    Vec3 coords_ = {0.0, 0.0, 0.0};
    Index dof_start_ = 0;
    Index ndofs_ = 6;
    std::array<BoundaryType, 7> bc_types_ = {};
    std::array<Real, 7> bc_values_ = {};
    bool on_interface_ = false;
    Index interface_id_ = -1;
    int partition_ = 0;
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_NODE_HPP
