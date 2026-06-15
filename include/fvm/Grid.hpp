#ifndef HVDC_FVM_GRID_HPP
#define HVDC_FVM_GRID_HPP

#include "common/Types.hpp"
#include "fvm/Cell.hpp"
#include "fvm/Face.hpp"
#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <memory>

namespace hvdc {
namespace fvm {

class BoundaryCondition;

struct GridStats {
    Index num_cells = 0;
    Index num_faces = 0;
    Index num_nodes = 0;
    Index num_boundary_faces = 0;
    Index num_interface_cells = 0;
    Real min_volume = INF;
    Real max_volume = 0.0;
    Real avg_volume = 0.0;
    Real min_face_area = INF;
    Real max_face_area = 0.0;
    Real min_orthogonality = 1.0;
    Real max_aspect_ratio = 0.0;
    Vec3 bounding_min = {INF, INF, INF};
    Vec3 bounding_max = {-INF, -INF, -INF};
};

class Grid {
public:
    Grid() = default;
    ~Grid() = default;
    
    static Grid create_structured_hex(
        const Vec3& origin,
        const Vec3& lengths,
        const std::array<Index, 3>& divisions,
        bool include_outer_boundary = true);
    
    static Grid create_channel_flow(
        Real length, Real height, Real width,
        Index n_x, Index n_y, Index n_z,
        const Vec3& origin = {0, 0, 0});
    
    static Grid create_cylinder_wake(
        Real cylinder_diameter,
        Real domain_length, Real domain_height, Real domain_width,
        Index refinement_level = 2);
    
    bool load_from_vtk(const std::string& filename);
    bool save_to_vtk(const std::string& filename) const;
    bool save_boundary_to_vtk(const std::string& filename) const;
    
    Index add_node(const Vec3& coords);
    Index add_cell(CellType type, const IndexVec& node_ids);
    Face* add_face(Index owner, Index neighbor, const IndexVec& node_ids);
    
    void build_connectivity();
    void compute_geometric_quantities();
    void compute_face_geometry(Face& face);
    void compute_cell_geometry(Cell& cell);
    void compute_interpolation_weights();
    void compute_nonorthogonality_correction();
    
    Index num_cells() const { return static_cast<Index>(cells_.size()); }
    Index num_faces() const { return static_cast<Index>(faces_.size()); }
    Index num_nodes() const { return static_cast<Index>(nodes_.size()); }
    Index num_internal_faces() const { return num_internal_faces_; }
    Index num_boundary_faces() const { return static_cast<Index>(faces_.size()) - num_internal_faces_; }
    
    Cell& cell(Index i) { return cells_[i]; }
    const Cell& cell(Index i) const { return cells_[i]; }
    
    Face& face(Index i) { return faces_[i]; }
    const Face& face(Index i) const { return faces_[i]; }
    
    Vec3& node(Index i) { return nodes_[i]; }
    const Vec3& node(Index i) const { return nodes_[i]; }
    
    const std::vector<Cell>& cells() const { return cells_; }
    std::vector<Cell>& cells() { return cells_; }
    
    const std::vector<Face>& faces() const { return faces_; }
    std::vector<Face>& faces() { return faces_; }
    
    const std::vector<Vec3>& nodes() const { return nodes_; }
    std::vector<Vec3>& nodes() { return nodes_; }
    
    Index add_boundary_patch(const std::string& name, Index start_face, Index n_faces);
    Index boundary_patch_id(const std::string& name) const;
    std::pair<Index, Index> boundary_patch_range(Index patch_id) const;
    const std::vector<std::string>& boundary_names() const { return boundary_names_; }
    
    void set_boundary_patch(Index patch_id, Index start, Index count);
    std::shared_ptr<fvm::BoundaryCondition>& mutable_face_bc(Index face_id);
    
    void mark_fsi_interface_cells(const std::vector<Index>& cell_ids);
    const std::vector<Index>& interface_cell_ids() const { return interface_cells_; }
    
    std::vector<Index> cells_in_box(const Vec3& box_min, const Vec3& box_max) const;
    std::vector<Index> cells_near_point(const Vec3& point, Real radius) const;
    std::vector<Index> cells_on_boundary_patch(Index patch_id) const;
    
    void compute_stats(GridStats& stats) const;
    
    void partition(int num_procs);
    int cell_partition(Index cell_id) const { return cells_[cell_id].partition(); }
    
    void refine_uniform(Index levels = 1);

protected:
    std::vector<Vec3> nodes_;
    std::vector<Cell> cells_;
    std::vector<Face> faces_;
    std::vector<std::shared_ptr<BoundaryCondition>> face_bcs_;
    
    Index num_internal_faces_ = 0;
    
    std::vector<std::string> boundary_names_;
    std::vector<std::pair<Index, Index>> boundary_ranges_;
    std::unordered_map<std::string, Index> boundary_name_map_;
    
    std::vector<Index> interface_cells_;
    
    Index next_node_id_ = 0;
    
    void add_face_to_boundary(Index face_id, Index boundary_id);
    
    void compute_face_normal_and_area(Face& face);
    void compute_face_centroid(Face& face);
    void compute_cell_centroid_and_volume(Cell& cell);
};

} // namespace fvm
} // namespace hvdc

#endif // HVDC_FVM_GRID_HPP
