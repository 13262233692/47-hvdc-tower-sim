#ifndef HVDC_FVM_CELL_HPP
#define HVDC_FVM_CELL_HPP

#include "common/Types.hpp"
#include <array>
#include <vector>

namespace hvdc {
namespace fvm {

enum class CellType : UInt8 {
    Hexa = 0,
    Tetra = 1,
    Prism = 2,
    Pyramid = 3,
    Poly = 4
};

class Cell {
public:
    Cell() = default;
    Cell(Index id, CellType type);
    
    Index id() const { return id_; }
    CellType type() const { return type_; }
    
    Index num_nodes() const { return static_cast<Index>(node_ids_.size()); }
    Index num_faces() const { return static_cast<Index>(face_ids_.size()); }
    Index num_neighbors() const { return static_cast<Index>(neighbor_ids_.size()); }
    
    const IndexVec& node_ids() const { return node_ids_; }
    IndexVec& node_ids() { return node_ids_; }
    Index node_id(Index i) const { return node_ids_[i]; }
    
    const IndexVec& face_ids() const { return face_ids_; }
    IndexVec& face_ids() { return face_ids_; }
    Index face_id(Index i) const { return face_ids_[i]; }
    
    const IndexVec& neighbor_ids() const { return neighbor_ids_; }
    IndexVec& neighbor_ids() { return neighbor_ids_; }
    Index neighbor_id(Index i) const { return neighbor_ids_[i]; }
    
    const Vec3& centroid() const { return centroid_; }
    Vec3& centroid() { return centroid_; }
    void set_centroid(const Vec3& c) { centroid_ = c; }
    
    Real volume() const { return volume_; }
    Real& volume() { return volume_; }
    void set_volume(Real v) { volume_ = v; }
    
    bool is_on_boundary() const {
        for (Index nid : neighbor_ids_) {
            if (nid < 0) return true;
        }
        return false;
    }
    
    Index boundary_faces_count() const {
        Index count = 0;
        for (Index nid : neighbor_ids_) {
            if (nid < 0) count++;
        }
        return count;
    }
    
    bool is_in_interface() const { return in_interface_; }
    void set_in_interface(bool b) { in_interface_ = b; }
    
    Index interface_id() const { return interface_id_; }
    void set_interface_id(Index id) { interface_id_ = id; }
    
    int partition() const { return partition_; }
    void set_partition(int p) { partition_ = p; }

private:
    Index id_ = -1;
    CellType type_ = CellType::Hexa;
    IndexVec node_ids_;
    IndexVec face_ids_;
    IndexVec neighbor_ids_;
    Vec3 centroid_ = {0.0, 0.0, 0.0};
    Real volume_ = 0.0;
    bool in_interface_ = false;
    Index interface_id_ = -1;
    int partition_ = 0;
};

} // namespace fvm
} // namespace hvdc

#endif // HVDC_FVM_CELL_HPP
