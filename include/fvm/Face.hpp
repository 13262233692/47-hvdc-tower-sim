#ifndef HVDC_FVM_FACE_HPP
#define HVDC_FVM_FACE_HPP

#include "common/Types.hpp"
#include <array>
#include <vector>

namespace hvdc {
namespace fvm {

class Face {
public:
    Face() = default;
    Face(Index id, Index num_nodes = 4);
    
    Index id() const { return id_; }
    
    Index num_nodes() const { return static_cast<Index>(node_ids_.size()); }
    const IndexVec& node_ids() const { return node_ids_; }
    IndexVec& node_ids() { return node_ids_; }
    Index node_id(Index i) const { return node_ids_[i]; }
    
    Index owner() const { return owner_; }
    Index& owner() { return owner_; }
    void set_owner(Index id) { owner_ = id; }
    
    Index neighbor() const { return neighbor_; }
    Index& neighbor() { return neighbor_; }
    void set_neighbor(Index id) { neighbor_ = id; }
    
    bool is_on_boundary() const { return neighbor_ < 0; }
    
    const Vec3& normal() const { return normal_; }
    Vec3& normal() { return normal_; }
    void set_normal(const Vec3& n) { normal_ = n; }
    
    Real area() const { return area_; }
    Real& area() { return area_; }
    void set_area(Real a) { area_ = a; }
    
    const Vec3& centroid() const { return centroid_; }
    Vec3& centroid() { return centroid_; }
    void set_centroid(const Vec3& c) { centroid_ = c; }
    
    Vec3 delta() const { return delta_; }
    Vec3& delta() { return delta_; }
    void set_delta(const Vec3& d) { delta_ = d; }
    
    Real delta_mag() const { return delta_mag_; }
    Real& delta_mag() { return delta_mag_; }
    void set_delta_mag(Real m) { delta_mag_ = m; }
    
    Real weight() const { return weight_; }
    Real& weight() { return weight_; }
    void set_weight(Real w) { weight_ = w; }
    
    Real cf() const { return cf_; }
    Real& cf() { return cf_; }
    void set_cf(Real c) { cf_ = c; }
    
    Index boundary_id() const { return boundary_id_; }
    void set_boundary_id(Index id) { boundary_id_ = id; }
    
    bool is_interface() const { return is_interface_; }
    void set_is_interface(bool b) { is_interface_ = b; }
    
    Index interface_group() const { return interface_group_; }
    void set_interface_group(Index g) { interface_group_ = g; }

private:
    Index id_ = -1;
    IndexVec node_ids_;
    Index owner_ = -1;
    Index neighbor_ = -1;
    Vec3 normal_ = {0.0, 0.0, 0.0};
    Real area_ = 0.0;
    Vec3 centroid_ = {0.0, 0.0, 0.0};
    Vec3 delta_ = {0.0, 0.0, 0.0};
    Real delta_mag_ = 0.0;
    Real weight_ = 0.5;
    Real cf_ = 0.0;
    Index boundary_id_ = -1;
    bool is_interface_ = false;
    Index interface_group_ = -1;
};

} // namespace fvm
} // namespace hvdc

#endif // HVDC_FVM_FACE_HPP
