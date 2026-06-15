#include "common/Types.hpp"
#include "common/MPIManager.hpp"
#include "common/PETScManager.hpp"
#include "common/Logger.hpp"
#include "fem/FEModel.hpp"
#include "fem/Assembler.hpp"
#include "fem/Solver.hpp"
#include <memory>
#include <cmath>

using namespace hvdc;
using namespace hvdc::fem;

int main(int argc, char** argv) {
    MPIManager mpi(&argc, &argv);
    PETScManager petsc(&argc, &argv);
    Logger::instance().set_level(LogLevel::INFO);
    Logger::instance().set_mpi_rank(mpi.rank());

    HVDC_LOG_INFO("=== Example 1: UHV Transmission Tower Static Analysis ===");

    auto model = std::make_shared<FEModel>();
    auto steel = std::make_shared<Material>(Material::Q420_Steel());
    model->add_material(steel);

    auto leg = std::make_shared<BeamSection>(BeamSection::circular_pipe(0.60, 0.020, steel));
    auto brace = std::make_shared<BeamSection>(BeamSection::circular_pipe(0.25, 0.010, steel));
    auto crossarm = std::make_shared<BeamSection>(BeamSection::i_beam(0.30, 0.40, 0.012, 0.018, steel));
    model->add_section(leg);
    model->add_section(brace);
    model->add_section(crossarm);

    const Real H = 108.0;
    const int NLEV = 18;
    std::vector<Index> legs[4];

    for (int li = 0; li <= NLEV; ++li) {
        Real z = H * li / NLEV;
        Real w = 10.0 - (10.0 - 2.0) * li / NLEV;
        Vec3 corners[4] = {
            { w/2,  w/2, z}, {-w/2,  w/2, z},
            {-w/2, -w/2, z}, { w/2, -w/2, z}
        };
        for (int c = 0; c < 4; ++c) {
            Index nid = model->add_node(corners[c]);
            legs[c].push_back(nid);
            if (li == 0) {
                model->set_node_bc_all(nid, BoundaryType::Fixed);
            }
        }
    }

    for (int c = 0; c < 4; ++c) {
        for (size_t i = 0; i + 1 < legs[c].size(); ++i) {
            model->add_element(std::make_shared<BeamElementNL>(
                model->get_node(legs[c][i]), model->get_node(legs[c][i+1]), leg));
        }
    }
    for (int li = 0; li < NLEV; ++li) {
        for (int c = 0; c < 4; ++c) {
            int c2 = (c + 1) % 4;
            model->add_element(std::make_shared<BeamElementNL>(
                model->get_node(legs[c][li]), model->get_node(legs[c2][li]), brace));
            model->add_element(std::make_shared<BeamElementNL>(
                model->get_node(legs[c][li]), model->get_node(legs[c2][li+1]), brace));
        }
    }

    for (int al : {10, 14, 17}) {
        Real z = H * al / NLEV;
        Real w = 10.0 - 8.0 * al / NLEV;
        for (Real off : {-15.0, -10.0, 10.0, 15.0}) {
            Index a1 = model->add_node({w/2 + off,  w/2, z});
            Index a2 = model->add_node({w/2 + off, -w/2, z});
            model->add_element(std::make_shared<BeamElementNL>(
                model->get_node(legs[0][al]), model->get_node(a1), crossarm));
            model->add_element(std::make_shared<BeamElementNL>(
                model->get_node(legs[3][al]), model->get_node(a2), crossarm));
        }
    }

    model->setup_dofs();
    HVDC_LOG_INFO("Tower model: " << model->num_nodes() << " nodes, "
                  << model->num_elements() << " elements, DOFs="
                  << model->total_dofs());

    Assembler assembler;
    Solver solver(model.get());
    SparseMatrix K;
    Vector F(model->total_dofs()), u(model->total_dofs());
    u.zero(); F.zero();

    assembler.assemble_tangent_stiffness(u, K, model.get());
    model->apply_body_force(F);
    assembler.apply_dirichlet_bc(K, F, u, *model);
    Vector du(model->total_dofs());
    solver.solve_linear_system(K, du, F);

    for (Index i = 0; i < u.size(); ++i) u.add(i, du[i]);

    Real max_d = 0.0, top_d = 0.0;
    for (Index i = 0; i < model->num_nodes(); ++i) {
        const Node* n = model->node(i);
        Index s = n->dof_start();
        Real dx = (s+0 < u.size()) ? u[s+0] : 0.0;
        Real dy = (s+1 < u.size()) ? u[s+1] : 0.0;
        Real dz = (s+2 < u.size()) ? u[s+2] : 0.0;
        Real dmag = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (dmag > max_d) max_d = dmag;
        if (std::fabs(n->coords()[2] - H) < 0.01) top_d = dmag;
    }
    HVDC_LOG_INFO("Gravity-only: max_displacement=" << max_d << " m, "
                  "top_tower_disp=" << top_d << " m");

    for (Index i = 0; i < model->num_nodes(); ++i) {
        const Node* n = model->node(i);
        if (n->bc_type() == BoundaryType::Fixed) continue;
        Index s = n->dof_start();
        Real q_wind = 0.5 * 1.225 * 40.0 * 40.0;
        Real Cd = 1.2;
        Real A = 0.6 * (H / NLEV);
        Real Fx = Cd * q_wind * A / 4.0;
        if (s < F.size()) F.add(s, Fx);
    }

    solver.solve_linear_system(K, du, F);
    for (Index i = 0; i < u.size(); ++i) u.add(i, du[i]);

    max_d = 0.0; top_d = 0.0;
    for (Index i = 0; i < model->num_nodes(); ++i) {
        const Node* n = model->node(i);
        Index s = n->dof_start();
        Vec3 d = {0,0,0};
        for (int k = 0; k < 3 && s + k < u.size(); ++k) d[k] = u[s + k];
        Real dmag = math::vec3_norm(d);
        if (dmag > max_d) max_d = dmag;
        if (std::fabs(n->coords()[2] - H) < 0.01) top_d = dmag;
    }
    HVDC_LOG_INFO("Gravity + 40m/s wind: max_disp=" << max_d
                  << " m, top_tower_disp=" << top_d << " m");

    HVDC_LOG_INFO("=== Example 1 Complete ===");
    return 0;
}
