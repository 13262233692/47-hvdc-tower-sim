#include "fem/FEModel.hpp"
#include "common/MathUtils.hpp"
#include "common/Logger.hpp"
#include <algorithm>
#include <numeric>

namespace hvdc {
namespace fem {

Node* FEModel::add_node(Index id, Real x, Real y, Real z) {
    auto node = std::make_unique<Node>(id, x, y, z);
    Node* ptr = node.get();
    nodes_storage_.push_back(std::move(node));
    node_id_map_[id] = ptr;
    nodes_.push_back(ptr);
    if (id >= next_node_id_) next_node_id_ = id + 1;
    return ptr;
}

Node* FEModel::add_node(const Vec3& coords) {
    return add_node(next_node_id_, coords[0], coords[1], coords[2]);
}

Node* FEModel::get_node(Index id) {
    auto it = node_id_map_.find(id);
    return (it != node_id_map_.end()) ? it->second : nullptr;
}

const Node* FEModel::get_node(Index id) const {
    auto it = node_id_map_.find(id);
    return (it != node_id_map_.end()) ? it->second : nullptr;
}

Index FEModel::add_material(const Material& mat) {
    Index id = (mat.id >= 0) ? mat.id : next_material_id_++;
    Material stored = mat;
    stored.id = id;
    materials_.push_back(stored);
    material_id_map_[id] = static_cast<Index>(materials_.size()) - 1;
    return id;
}

Material* FEModel::get_material(Index id) {
    auto it = material_id_map_.find(id);
    return (it != material_id_map_.end()) ? &materials_[it->second] : nullptr;
}

const Material* FEModel::get_material(Index id) const {
    auto it = material_id_map_.find(id);
    return (it != material_id_map_.end()) ? &materials_[it->second] : nullptr;
}

Index FEModel::add_beam_section(const BeamSection& sec) {
    Index id = (sec.id >= 0) ? sec.id : next_beam_sec_id_++;
    BeamSection stored = sec;
    stored.id = id;
    beam_sections_.push_back(stored);
    beam_sec_id_map_[id] = static_cast<Index>(beam_sections_.size()) - 1;
    return id;
}

BeamSection* FEModel::get_beam_section(Index id) {
    auto it = beam_sec_id_map_.find(id);
    return (it != beam_sec_id_map_.end()) ? &beam_sections_[it->second] : nullptr;
}

const BeamSection* FEModel::get_beam_section(Index id) const {
    auto it = beam_sec_id_map_.find(id);
    return (it != beam_sec_id_map_.end()) ? &beam_sections_[it->second] : nullptr;
}

Index FEModel::add_truss_section(const TrussSection& sec) {
    Index id = (sec.id >= 0) ? sec.id : next_truss_sec_id_++;
    TrussSection stored = sec;
    stored.id = id;
    truss_sections_.push_back(stored);
    truss_sec_id_map_[id] = static_cast<Index>(truss_sections_.size()) - 1;
    return id;
}

Index FEModel::add_material(std::shared_ptr<Material> mat_ptr) {
    if (!mat_ptr) return -1;
    return add_material(*mat_ptr);
}

Index FEModel::add_section(std::shared_ptr<BeamSection> sec_ptr) {
    if (!sec_ptr) return -1;
    return add_beam_section(*sec_ptr);
}

Index FEModel::add_truss_section(std::shared_ptr<TrussSection> sec_ptr) {
    if (!sec_ptr) return -1;
    return add_truss_section(*sec_ptr);
}

TrussSection* FEModel::get_truss_section(Index id) {
    auto it = truss_sec_id_map_.find(id);
    return (it != truss_sec_id_map_.end()) ? &truss_sections_[it->second] : nullptr;
}

const TrussSection* FEModel::get_truss_section(Index id) const {
    auto it = truss_sec_id_map_.find(id);
    return (it != truss_sec_id_map_.end()) ? &truss_sections_[it->second] : nullptr;
}

Element* FEModel::add_element(std::unique_ptr<Element> elem) {
    Element* ptr = elem.get();
    elements_.push_back(std::move(elem));
    return ptr;
}

BeamElementNL* FEModel::add_beam_element(
    Index node1_id, Index node2_id,
    Index material_id, Index section_id)
{
    Node* n1 = get_node(node1_id);
    Node* n2 = get_node(node2_id);
    Material* mat = get_material(material_id);
    BeamSection* sec = get_beam_section(section_id);
    if (!n1 || !n2 || !mat || !sec) return nullptr;
    
    auto elem = std::make_unique<BeamElementNL>(
        static_cast<Index>(elements_.size()), n1, n2, mat, sec);
    BeamElementNL* ptr = elem.get();
    elements_.push_back(std::move(elem));
    return ptr;
}

TrussElementNL* FEModel::add_truss_element(
    Index node1_id, Index node2_id,
    Index material_id, Index section_id)
{
    Node* n1 = get_node(node1_id);
    Node* n2 = get_node(node2_id);
    Material* mat = get_material(material_id);
    TrussSection* sec = get_truss_section(section_id);
    if (!n1 || !n2 || !mat || !sec) return nullptr;
    
    auto elem = std::make_unique<TrussElementNL>(
        static_cast<Index>(elements_.size()), n1, n2, mat, sec);
    TrussElementNL* ptr = elem.get();
    elements_.push_back(std::move(elem));
    return ptr;
}

ConductorElement* FEModel::add_conductor_element(
    Index node1_id, Index node2_id,
    Index material_id, Index section_id,
    Real initial_tension)
{
    Node* n1 = get_node(node1_id);
    Node* n2 = get_node(node2_id);
    Material* mat = get_material(material_id);
    BeamSection* sec = get_beam_section(section_id);
    if (!n1 || !n2 || !mat || !sec) return nullptr;
    
    auto elem = std::make_unique<ConductorElement>(
        static_cast<Index>(elements_.size()), n1, n2, mat, sec, initial_tension);
    ConductorElement* ptr = elem.get();
    elements_.push_back(std::move(elem));
    return ptr;
}

void FEModel::setup_dofs(Index dofs_per_node) {
    total_dofs_ = 0;
    num_constrained_dofs_ = 0;
    free_dofs_.clear();
    constrained_dofs_.clear();
    
    for (auto* node : nodes_) {
        node->set_dof_start(total_dofs_);
        node->set_num_dofs(dofs_per_node);
        for (Index d = 0; d < dofs_per_node; ++d) {
            Index gdof = total_dofs_ + d;
            if (node->is_dof_active(d)) {
                free_dofs_.push_back(gdof);
            } else {
                constrained_dofs_.push_back(gdof);
                num_constrained_dofs_++;
            }
        }
        total_dofs_ += dofs_per_node;
    }
    
    HVDC_LOG_INFO("FEM DOF setup: total=" << total_dofs_ 
                  << ", free=" << free_dofs_.size()
                  << ", constrained=" << constrained_dofs_.size());
}

void FEModel::apply_fixed_bc(Index node_id, const std::vector<Index>& dofs) {
    Node* n = get_node(node_id);
    if (!n) {
        HVDC_LOG_WARNING("apply_fixed_bc: node " << node_id << " not found");
        return;
    }
    for (Index d : dofs) {
        if (d < 7) {
            n->set_bc_type(d, BoundaryType::Fixed);
            n->set_bc_value(d, 0.0);
        }
    }
}

void FEModel::set_node_bc_all(Index node_id, BoundaryType type) {
    Node* n = get_node(node_id);
    if (!n) {
        HVDC_LOG_WARNING("set_node_bc_all: node " << node_id << " not found");
        return;
    }
    n->set_bc_type_all(type);
}

void FEModel::set_node_bc(Index node_id, BoundaryType type,
                          bool tx, bool ty, bool tz,
                          bool rx, bool ry, bool rz) {
    Node* n = get_node(node_id);
    if (!n) {
        HVDC_LOG_WARNING("set_node_bc: node " << node_id << " not found");
        return;
    }
    bool active[6] = {tx, ty, tz, rx, ry, rz};
    for (Index d = 0; d < 6; ++d) {
        if (active[d]) n->set_bc_type(d, type);
    }
}

void FEModel::apply_body_force(Vector& F_ext, Real g) const {
    Vec3 gravity = {0.0, 0.0, -g};
    for (const auto& elem : elements_) {
        if (!elem) continue;
        Vec12 loads = elem->gravity_loads(g);
        for (Index ni = 0; ni < elem->num_nodes(); ++ni) {
            const Node* node = elem->node(ni);
            if (!node) continue;
            Index start = node->dof_start();
            for (Index d = 0; d < 3; ++d) {
                Index gdof = start + d;
                if (gdof < F_ext.size() && node->is_dof_active(d)) {
                    F_ext.add(gdof, loads[ni * 6 + d]);
                }
            }
        }
    }
}

void FEModel::apply_ice_load(Vector& F_ext, Real rho_ice) const {
    Vec3 gravity = {0.0, 0.0, -GRAVITY};
    for (const auto& elem : elements_) {
        if (!elem) continue;
        Real A_ice = 0.0, D_ice = 0.0;
        if (auto* beam = dynamic_cast<const BeamElementNL*>(elem.get())) {
            if (beam->section()) {
                A_ice = beam->section()->ice_area();
                D_ice = beam->section()->wind_drag_diameter();
            }
        } else if (auto* truss = dynamic_cast<const TrussElementNL*>(elem.get())) {
            if (truss->section()) {
                A_ice = truss->section()->ice_area();
            }
        } else if (auto* cond = dynamic_cast<const ConductorElement*>(elem.get())) {
            if (cond->section()) {
                A_ice = cond->section()->ice_area();
                D_ice = cond->section()->wind_drag_diameter();
            }
        }
        if (A_ice < EPS) continue;
        Real L = elem->length_initial();
        Real w_ice = rho_ice * GRAVITY * A_ice * L;
        for (Index ni = 0; ni < elem->num_nodes(); ++ni) {
            const Node* node = elem->node(ni);
            if (!node) continue;
            Index start = node->dof_start();
            Index gdof = start + 2;
            if (gdof < F_ext.size() && node->is_dof_active(2)) {
                F_ext.add(gdof, -0.5 * w_ice);
            }
        }
    }
}

void FEModel::apply_displacement_bc(Index node_id, Index dof, Real value) {
    Node* n = get_node(node_id);
    if (!n) {
        HVDC_LOG_WARNING("apply_displacement_bc: node " << node_id << " not found");
        return;
    }
    if (dof < 7) {
        n->set_bc_type(dof, BoundaryType::Displacement);
        n->set_bc_value(dof, value);
    }
}

void FEModel::gather_displacement_vec(Vec& u_global) const {
    u_global.assign(total_dofs_, 0.0);
    for (const auto* node : nodes_) {
        Index start = node->dof_start();
        Index ndofs = node->num_dofs();
        for (Index d = 0; d < ndofs; ++d) {
            if (!node->is_dof_active(d)) {
                u_global[start + d] = node->bc_value(d);
            }
        }
    }
}

void FEModel::scatter_displacement_vec(const Vec& u_global) {
    if (static_cast<Index>(u_global.size()) != total_dofs_) return;
    (void)u_global;
}

void FEModel::apply_ice_loads(Real ice_thickness, Real rho_ice) {
    for (auto& elem : elements_) {
        if (auto* beam = dynamic_cast<BeamElementNL*>(elem.get())) {
            if (beam->section()) {
                const_cast<BeamSection*>(beam->section())->add_ice_coating(ice_thickness);
                const_cast<BeamSection*>(beam->section())->rho_ice = rho_ice;
            }
        } else if (auto* truss = dynamic_cast<TrussElementNL*>(elem.get())) {
            if (truss->section()) {
                const_cast<TrussSection*>(truss->section())->add_ice_coating(ice_thickness);
            }
        } else if (auto* cond = dynamic_cast<ConductorElement*>(elem.get())) {
            if (cond->section()) {
                const_cast<BeamSection*>(cond->section())->add_ice_coating(ice_thickness);
            }
        }
    }
    
    HVDC_LOG_INFO("Applied ice coating: thickness=" << ice_thickness 
                  << "m, rho=" << rho_ice << "kg/m^3");
}

void FEModel::apply_temperature_loads(Real delta_T_conductor, Real delta_T_tower) {
    for (auto& elem : elements_) {
        if (auto* cond = dynamic_cast<ConductorElement*>(elem.get())) {
            cond->apply_temperature(delta_T_conductor);
        } else {
            elem->apply_temperature(delta_T_tower);
        }
    }
}

void FEModel::mark_fsi_interface_nodes(const std::vector<Index>& node_ids) {
    interface_nodes_.clear();
    for (Index id : node_ids) {
        Node* n = get_node(id);
        if (n) {
            n->set_on_interface(true);
            n->set_interface_id(static_cast<Index>(interface_nodes_.size()));
            interface_nodes_.push_back(n);
        }
    }
    HVDC_LOG_INFO("Marked " << interface_nodes_.size() << " FSI interface nodes");
}

std::vector<Index> FEModel::interface_node_dofs() const {
    std::vector<Index> dofs;
    for (const auto* node : interface_nodes_) {
        Index start = node->dof_start();
        for (Index d = 0; d < 3; ++d) {
            if (node->is_dof_active(d)) {
                dofs.push_back(start + d);
            }
        }
    }
    return dofs;
}

void FEModel::partition_domain(int num_procs) {
    Index n_elems = static_cast<Index>(elements_.size());
    for (Index i = 0; i < n_elems; ++i) {
        elements_[i]->set_partition(static_cast<int>(i * num_procs / n_elems));
    }
}

int FEModel::element_partition(Index elem_id) const {
    if (elem_id < 0 || elem_id >= static_cast<Index>(elements_.size())) return 0;
    return elements_[elem_id]->partition();
}

void FEModel::compute_stats(Real& max_sag, Real& total_mass,
                              Vec3& center_of_mass) const {
    max_sag = 0.0;
    total_mass = 0.0;
    center_of_mass = {0.0, 0.0, 0.0};
    Real z_max = -INF;
    Real z_min = INF;
    
    for (const auto* n : nodes_) {
        z_max = std::max(z_max, n->z());
        z_min = std::min(z_min, n->z());
    }
    max_sag = z_max - z_min;
    
    for (const auto& elem : elements_) {
        Real mass = 0.0;
        if (auto* beam = dynamic_cast<const BeamElementNL*>(elem.get())) {
            if (beam->material() && beam->section()) {
                mass = beam->section()->total_mass_per_length(beam->material()->rho) 
                       * beam->length_initial();
            }
        } else if (auto* truss = dynamic_cast<const TrussElementNL*>(elem.get())) {
            if (truss->material() && truss->section()) {
                mass = truss->material()->rho * truss->section()->total_A() 
                       * truss->length_initial();
                if (truss->section()->thickness_ice > 0.0) {
                    mass += 917.0 * truss->section()->A_ice * truss->length_initial();
                }
            }
        } else if (auto* cond = dynamic_cast<const ConductorElement*>(elem.get())) {
            if (cond->material() && cond->section()) {
                mass = cond->material()->rho * cond->section()->A * cond->catenary_length();
                if (cond->section()->A_ice > 0.0) {
                    mass += 917.0 * cond->section()->A_ice * cond->catenary_length();
                }
            }
        }
        
        total_mass += mass;
        
        const Node* n1 = elem->node(0);
        const Node* n2 = elem->node(1);
        if (n1 && n2) {
            Vec3 c = {0.5 * (n1->x() + n2->x()),
                      0.5 * (n1->y() + n2->y()),
                      0.5 * (n1->z() + n2->z())};
            math::vec3_add_inplace(center_of_mass, math::vec3_scale(c, mass));
        }
    }
    
    if (total_mass > EPS) {
        math::vec3_scale_inplace(center_of_mass, 1.0 / total_mass);
    }
}

std::vector<Index> FEModel::conductor_element_ids() const {
    std::vector<Index> ids;
    for (Index i = 0; i < static_cast<Index>(elements_.size()); ++i) {
        if (dynamic_cast<const ConductorElement*>(elements_[i].get())) {
            ids.push_back(i);
        }
    }
    return ids;
}

std::vector<Index> FEModel::tower_element_ids() const {
    std::vector<Index> ids;
    for (Index i = 0; i < static_cast<Index>(elements_.size()); ++i) {
        if (!dynamic_cast<const ConductorElement*>(elements_[i].get())) {
            ids.push_back(i);
        }
    }
    return ids;
}

} // namespace fem
} // namespace hvdc
