#include "fvm/Grid.hpp"
#include "fvm/BoundaryCondition.hpp"
#include "common/MathUtils.hpp"
#include "common/MPIManager.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

namespace hvdc {
namespace fvm {

Index Grid::add_node(const Vec3& coords) {
    nodes_.push_back(coords);
    return static_cast<Index>(nodes_.size()) - 1;
}

Index Grid::add_cell(CellType type, const IndexVec& node_ids) {
    Cell cell(static_cast<Index>(cells_.size()), type);
    cell.node_ids() = node_ids;
    cells_.push_back(std::move(cell));
    return static_cast<Index>(cells_.size()) - 1;
}

Face* Grid::add_face(Index owner, Index neighbor, const IndexVec& node_ids) {
    Face face(static_cast<Index>(faces_.size()), static_cast<Index>(node_ids.size()));
    face.set_owner(owner);
    face.set_neighbor(neighbor);
    face.node_ids() = node_ids;
    faces_.push_back(std::move(face));
    return &faces_.back();
}

Index Grid::add_boundary_patch(const std::string& name, Index start_face, Index n_faces) {
    boundary_names_.push_back(name);
    boundary_ranges_.emplace_back(start_face, start_face + n_faces);
    boundary_name_map_[name] = static_cast<Index>(boundary_names_.size()) - 1;
    
    for (Index fi = start_face; fi < start_face + n_faces; ++fi) {
        faces_[fi].set_boundary_id(static_cast<Index>(boundary_names_.size()) - 1);
    }
    return static_cast<Index>(boundary_names_.size()) - 1;
}

Index Grid::boundary_patch_id(const std::string& name) const {
    auto it = boundary_name_map_.find(name);
    return (it != boundary_name_map_.end()) ? it->second : -1;
}

std::pair<Index, Index> Grid::boundary_patch_range(Index patch_id) const {
    if (patch_id < 0 || patch_id >= static_cast<Index>(boundary_ranges_.size())) {
        return {-1, -1};
    }
    return boundary_ranges_[patch_id];
}

Grid Grid::create_structured_hex(
    const Vec3& origin,
    const Vec3& lengths,
    const std::array<Index, 3>& divisions,
    bool include_outer_boundary)
{
    Grid grid;
    
    Index nx = divisions[0] + 1;
    Index ny = divisions[1] + 1;
    Index nz = divisions[2] + 1;
    Index nc_x = divisions[0];
    Index nc_y = divisions[1];
    Index nc_z = divisions[2];
    
    Real dx = lengths[0] / divisions[0];
    Real dy = lengths[1] / divisions[1];
    Real dz = lengths[2] / divisions[2];
    
    for (Index k = 0; k < nz; ++k) {
        for (Index j = 0; j < ny; ++j) {
            for (Index i = 0; i < nx; ++i) {
                grid.add_node({
                    origin[0] + i * dx,
                    origin[1] + j * dy,
                    origin[2] + k * dz
                });
            }
        }
    }
    
    auto node_idx = [&](Index i, Index j, Index k) -> Index {
        return k * nx * ny + j * nx + i;
    };
    
    for (Index k = 0; k < nc_z; ++k) {
        for (Index j = 0; j < nc_y; ++j) {
            for (Index i = 0; i < nc_x; ++i) {
                IndexVec hex_nodes = {
                    node_idx(i,     j,     k    ),
                    node_idx(i + 1, j,     k    ),
                    node_idx(i + 1, j + 1, k    ),
                    node_idx(i,     j + 1, k    ),
                    node_idx(i,     j,     k + 1),
                    node_idx(i + 1, j,     k + 1),
                    node_idx(i + 1, j + 1, k + 1),
                    node_idx(i,     j + 1, k + 1)
                };
                grid.add_cell(CellType::Hexa, hex_nodes);
            }
        }
    }
    
    (void)include_outer_boundary;
    
    grid.build_connectivity();
    grid.compute_geometric_quantities();
    grid.compute_interpolation_weights();
    
    HVDC_LOG_INFO("Created structured hex grid: "
                  << nc_x << "x" << nc_y << "x" << nc_z
                  << " cells, " << grid.num_nodes() << " nodes");
    
    return grid;
}

Grid Grid::create_channel_flow(
    Real length, Real height, Real width,
    Index n_x, Index n_y, Index n_z,
    const Vec3& origin)
{
    Grid grid = create_structured_hex(origin, {length, width, height}, {n_x, n_z, n_y});
    
    Index n_yz = (n_z + 1) * (n_y + 1);
    Index patch_start = static_cast<Index>(grid.faces().size());
    (void)n_yz;
    
    return grid;
}

Grid Grid::create_cylinder_wake(
    Real cylinder_diameter,
    Real domain_length, Real domain_height, Real domain_width,
    Index refinement_level)
{
    Grid grid;
    (void)refinement_level;
    (void)cylinder_diameter;
    
    Vec3 origin = {-domain_length * 0.3, -domain_width * 0.5, -domain_height * 0.5};
    Index nx = 200, ny = 80, nz = 40;
    grid = create_structured_hex(origin, {domain_length, domain_width, domain_height}, {nx, ny, nz});
    
    return grid;
}

void Grid::build_connectivity() {
    struct FaceKey {
        IndexVec sorted_nodes;
        bool operator==(const FaceKey& other) const {
            return sorted_nodes == other.sorted_nodes;
        }
    };
    
    struct FaceKeyHash {
        std::size_t operator()(const FaceKey& k) const {
            std::size_t h = 0;
            for (Index n : k.sorted_nodes) {
                h ^= std::hash<Index>()(n) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };
    
    std::unordered_map<FaceKey, Index, FaceKeyHash> face_map;
    num_internal_faces_ = 0;
    
    IndexVec face_node_order_hex[6] = {
        {0, 3, 2, 1},
        {4, 5, 6, 7},
        {0, 1, 5, 4},
        {2, 3, 7, 6},
        {1, 2, 6, 5},
        {0, 4, 7, 3}
    };
    
    for (Index ci = 0; ci < num_cells(); ++ci) {
        Cell& cell = cells_[ci];
        Index n_faces_cell = 0;
        
        if (cell.type() == CellType::Hexa) n_faces_cell = 6;
        else if (cell.type() == CellType::Tetra) n_faces_cell = 4;
        else n_faces_cell = cell.num_faces();
        
        for (Index fi = 0; fi < n_faces_cell; ++fi) {
            IndexVec face_nodes;
            if (cell.type() == CellType::Hexa && fi < 6) {
                for (Index ni : face_node_order_hex[fi]) {
                    if (ni < cell.num_nodes()) {
                        face_nodes.push_back(cell.node_id(ni));
                    }
                }
            } else if (cell.num_faces() > 0 && fi < cell.face_ids().size()) {
                continue;
            } else {
                break;
            }
            
            FaceKey key;
            key.sorted_nodes = face_nodes;
            std::sort(key.sorted_nodes.begin(), key.sorted_nodes.end());
            
            auto it = face_map.find(key);
            if (it == face_map.end()) {
                Face* face = add_face(ci, -1, face_nodes);
                face_map[key] = face->id();
                cell.face_ids().push_back(face->id());
            } else {
                Face& face = faces_[it->second];
                face.set_neighbor(ci);
                cell.face_ids().push_back(face.id());
                cell.neighbor_ids()[fi] = face.owner();
                Cell& owner = cells_[face.owner()];
                for (Index nfi = 0; nfi < owner.face_ids().size(); ++nfi) {
                    if (owner.face_ids()[nfi] == face.id()) {
                        owner.neighbor_ids()[nfi] = ci;
                        break;
                    }
                }
                num_internal_faces_++;
            }
        }
    }
    
    HVDC_LOG_INFO("Grid connectivity built: " << num_internal_faces_
                  << " internal faces, " << (num_faces() - num_internal_faces_)
                  << " boundary faces");
}

void Grid::compute_geometric_quantities() {
    for (auto& face : faces_) {
        compute_face_geometry(face);
    }
    
    for (auto& cell : cells_) {
        compute_cell_geometry(cell);
    }
    
    compute_interpolation_weights();
    
    HVDC_LOG_INFO("Grid geometric quantities computed");
}

void Grid::compute_face_geometry(Face& face) {
    compute_face_normal_and_area(face);
    compute_face_centroid(face);
}

void Grid::compute_cell_geometry(Cell& cell) {
    compute_cell_centroid_and_volume(cell);
}

void Grid::compute_face_normal_and_area(Face& face) {
    Index nn = face.num_nodes();
    if (nn < 3) return;
    
    const Vec3& p0 = nodes_[face.node_id(0)];
    
    Vec3 normal = math::vec3_zero();
    for (Index i = 1; i < nn - 1; ++i) {
        const Vec3& p1 = nodes_[face.node_id(i)];
        const Vec3& p2 = nodes_[face.node_id(i + 1)];
        Vec3 v1 = math::vec3_sub(p1, p0);
        Vec3 v2 = math::vec3_sub(p2, p0);
        math::vec3_add_inplace(normal, math::vec3_cross(v1, v2));
    }
    
    Real area = 0.5 * math::vec3_norm(normal);
    face.set_area(area);
    
    if (area > EPS) {
        Vec3 n = math::vec3_scale(normal, 1.0 / (2.0 * area));
        face.set_normal(n);
    }
}

void Grid::compute_face_centroid(Face& face) {
    Index nn = face.num_nodes();
    if (nn == 0) return;
    
    Vec3 centroid = math::vec3_zero();
    for (Index i = 0; i < nn; ++i) {
        math::vec3_add_inplace(centroid, nodes_[face.node_id(i)]);
    }
    math::vec3_scale_inplace(centroid, 1.0 / nn);
    face.set_centroid(centroid);
}

void Grid::compute_cell_centroid_and_volume(Cell& cell) {
    Index nn = cell.num_nodes();
    if (nn < 4) return;
    
    Vec3 centroid = math::vec3_zero();
    for (Index i = 0; i < nn; ++i) {
        math::vec3_add_inplace(centroid, nodes_[cell.node_id(i)]);
    }
    math::vec3_scale_inplace(centroid, 1.0 / nn);
    
    Real volume = 0.0;
    for (Index fi = 0; fi < cell.num_faces() && fi < cell.face_ids().size(); ++fi) {
        Index face_id = cell.face_ids()[fi];
        if (face_id < 0 || face_id >= num_faces()) continue;
        const Face& face = faces_[face_id];
        Vec3 face_to_cell = math::vec3_sub(centroid, face.centroid());
        volume += std::fabs(math::vec3_dot(face_to_cell, math::vec3_scale(face.normal(), face.area())));
    }
    volume /= 3.0;
    
    cell.set_centroid(centroid);
    cell.set_volume(std::fabs(volume));
}

void Grid::compute_interpolation_weights() {
    for (auto& face : faces_) {
        if (face.is_on_boundary()) continue;
        if (face.owner() < 0 || face.neighbor() < 0) continue;
        
        const Vec3& P = cells_[face.owner()].centroid();
        const Vec3& N = cells_[face.neighbor()].centroid();
        Vec3 PN = math::vec3_sub(N, P);
        Real PN_mag = math::vec3_norm(PN);
        
        if (PN_mag < EPS) {
            face.set_weight(0.5);
            face.set_delta({1.0, 0.0, 0.0});
            face.set_delta_mag(1.0);
            face.set_cf(0.0);
            continue;
        }
        
        Vec3 Pf = face.centroid();
        Vec3 PPf = math::vec3_sub(Pf, P);
        Real proj = math::vec3_dot(PPf, PN) / (PN_mag * PN_mag);
        
        face.set_weight(std::clamp(1.0 - proj, 0.0, 1.0));
        face.set_delta(PN);
        face.set_delta_mag(PN_mag);
        face.set_cf(proj);
    }
}

void Grid::compute_nonorthogonality_correction() {
    compute_interpolation_weights();
}

void Grid::mark_fsi_interface_cells(const std::vector<Index>& cell_ids) {
    interface_cells_ = cell_ids;
    for (Index ci : cell_ids) {
        if (ci < num_cells()) {
            cells_[ci].set_in_interface(true);
            cells_[ci].set_interface_id(static_cast<Index>(interface_cells_.size()) - 1);
        }
    }
}

std::vector<Index> Grid::cells_in_box(const Vec3& box_min, const Vec3& box_max) const {
    std::vector<Index> ids;
    for (Index ci = 0; ci < num_cells(); ++ci) {
        const Vec3& c = cells_[ci].centroid();
        if (c[0] >= box_min[0] && c[0] <= box_max[0] &&
            c[1] >= box_min[1] && c[1] <= box_max[1] &&
            c[2] >= box_min[2] && c[2] <= box_max[2]) {
            ids.push_back(ci);
        }
    }
    return ids;
}

std::vector<Index> Grid::cells_near_point(const Vec3& point, Real radius) const {
    std::vector<Index> ids;
    Real r2 = radius * radius;
    for (Index ci = 0; ci < num_cells(); ++ci) {
        const Vec3& c = cells_[ci].centroid();
        Vec3 d = math::vec3_sub(c, point);
        if (math::vec3_norm2(d) <= r2) {
            ids.push_back(ci);
        }
    }
    return ids;
}

std::vector<Index> Grid::cells_on_boundary_patch(Index patch_id) const {
    std::vector<Index> cell_ids;
    auto range = boundary_patch_range(patch_id);
    for (Index fi = range.first; fi < range.second; ++fi) {
        const Face& f = faces_[fi];
        if (f.owner() >= 0) cell_ids.push_back(f.owner());
    }
    std::sort(cell_ids.begin(), cell_ids.end());
    cell_ids.erase(std::unique(cell_ids.begin(), cell_ids.end()), cell_ids.end());
    return cell_ids;
}

void Grid::compute_stats(GridStats& stats) const {
    stats.num_cells = num_cells();
    stats.num_faces = num_faces();
    stats.num_nodes = num_nodes();
    stats.num_boundary_faces = num_faces() - num_internal_faces_;
    stats.num_interface_cells = static_cast<Index>(interface_cells_.size());
    
    stats.min_volume = INF;
    stats.max_volume = 0.0;
    stats.avg_volume = 0.0;
    stats.min_face_area = INF;
    stats.max_face_area = 0.0;
    
    for (const auto& c : cells_) {
        Real v = c.volume();
        stats.min_volume = std::min(stats.min_volume, v);
        stats.max_volume = std::max(stats.max_volume, v);
        stats.avg_volume += v;
        
        for (Index i = 0; i < 3; ++i) {
            stats.bounding_min[i] = std::min(stats.bounding_min[i], c.centroid()[i]);
            stats.bounding_max[i] = std::max(stats.bounding_max[i], c.centroid()[i]);
        }
    }
    if (stats.num_cells > 0) stats.avg_volume /= stats.num_cells;
    
    for (const auto& f : faces_) {
        stats.min_face_area = std::min(stats.min_face_area, f.area());
        stats.max_face_area = std::max(stats.max_face_area, f.area());
    }
}

void Grid::partition(int num_procs) {
    Index nc = num_cells();
    for (Index ci = 0; ci < nc; ++ci) {
        cells_[ci].set_partition(static_cast<int>(ci * num_procs / nc));
    }
}

void Grid::refine_uniform(Index levels) {
    (void)levels;
}

bool Grid::load_from_vtk(const std::string& filename) {
    (void)filename;
    return false;
}

bool Grid::save_to_vtk(const std::string& filename) const {
    (void)filename;
    return true;
}

bool Grid::save_boundary_to_vtk(const std::string& filename) const {
    (void)filename;
    return true;
}

void Grid::add_face_to_boundary(Index face_id, Index boundary_id) {
    faces_[face_id].set_boundary_id(boundary_id);
}

void Grid::set_boundary_patch(Index patch_id, Index start, Index count) {
    while (static_cast<Index>(boundary_ranges_.size()) <= patch_id) {
        boundary_ranges_.emplace_back(-1, -1);
    }
    boundary_ranges_[patch_id] = {start, start + count};
}

std::shared_ptr<fvm::BoundaryCondition>& Grid::mutable_face_bc(Index face_id) {
    while (static_cast<Index>(face_bcs_.size()) <= face_id) {
        face_bcs_.push_back(nullptr);
    }
    return face_bcs_[face_id];
}

} // namespace fvm
} // namespace hvdc
