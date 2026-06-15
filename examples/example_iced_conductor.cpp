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

    HVDC_LOG_INFO("=== Example 2: Iced Conductor Aerodynamic Analysis ===");

    auto model = std::make_shared<FEModel>();
    auto acsr = std::make_shared<Material>(Material::ACSR_Conductor());
    model->add_material(acsr);

    auto cond_sec = std::make_shared<BeamSection>(BeamSection::circular_pipe(0.0336, 0.0112, acsr));
    cond_sec->thickness_ice = 0.010;
    model->add_section(cond_sec);

    const Real span = 400.0;
    const Real sag = 12.0;
    const int NSEG = 80;
    const Real z_left = 70.0;

    std::vector<Index> nodes;
    nodes.push_back(model->add_node({0.0, 0.0, z_left}));
    for (int i = 1; i < NSEG; ++i) {
        Real xi = static_cast<Real>(i) / NSEG;
        Real xc = span * xi;
        Real yc = 0.0;
        Real zc = z_left - 4.0 * sag * xi * (1.0 - xi);
        nodes.push_back(model->add_node({xc, yc, zc}));
    }
    nodes.push_back(model->add_node({span, 0.0, z_left - sag}));

    model->set_node_bc_all(nodes.front(), BoundaryType::Fixed);
    model->set_node_bc_all(nodes.back(), BoundaryType::Fixed);

    for (size_t i = 0; i + 1 < nodes.size(); ++i) {
        model->add_element(std::make_shared<BeamElementNL>(
            model->get_node(nodes[i]), model->get_node(nodes[i+1]), cond_sec));
    }

    model->setup_dofs();
    HVDC_LOG_INFO("Conductor model: " << model->num_nodes() << " nodes, "
                  << model->num_elements() << " elements, DOFs="
                  << model->total_dofs());

    Assembler assembler;
    Solver solver(model.get());
    SparseMatrix K;
    Vector F(model->total_dofs()), u(model->total_dofs());
    u.zero(); F.zero();

    assembler.assemble_tangent_stiffness(u, K, model.get());
    HVDC_LOG_INFO("Tangent stiffness assembled, nonzeros=" << K.num_nonzeros());

    model->apply_body_force(F);
    assembler.apply_dirichlet_bc(K, F, u, *model);
    Vector du(model->total_dofs());
    solver.solve_linear_system(K, du, F);
    for (Index i = 0; i < u.size(); ++i) u.add(i, du[i]);

    Index mid_n = nodes[NSEG / 2];
    const Node* mnode = model->get_node(mid_n);
    Index s = mnode->dof_start();
    HVDC_LOG_INFO("Mid-span deflection (gravity only): "
                  << ((s+2 < u.size()) ? u[s+2] : 0.0) << " m (expected ~" << sag << " m)");

    F.zero();
    Real rho = 1.225;
    for (Index ei = 0; ei < model->num_elements(); ++ei) {
        const Element* elem = model->element(ei);
        if (!elem) continue;
        Vec3 midp = elem->midpoint();
        Real h = midp[2];
        Real U_10m = 40.0;
        Real U_h = U_10m * std::pow(h / 10.0, 0.14);
        Real q = 0.5 * rho * U_h * U_h;
        Real Cd = 1.6;
        Real D = cond_sec->h + 2.0 * cond_sec->thickness_ice;
        Real L = elem->length_initial();
        Real Fy = Cd * q * D * L;
        for (Index ni = 0; ni < 2; ++ni) {
            const Node* n = elem->node(ni);
            Index ns = n->dof_start();
            if (ns + 1 < F.size()) F.add(ns + 1, Fy * 0.5);
        }
    }

    assembler.apply_dirichlet_bc(K, F, u, *model);
    solver.solve_linear_system(K, du, F);
    for (Index i = 0; i < u.size(); ++i) u.add(i, du[i]);

    Real max_sway = 0.0, max_sway_z = 0.0;
    for (Index ni : nodes) {
        const Node* n = model->get_node(ni);
        Index ns = n->dof_start();
        if (ns + 1 < u.size() && std::fabs(u[ns+1]) > max_sway) {
            max_sway = std::fabs(u[ns+1]);
            max_sway_z = n->coords()[2];
        }
    }
    HVDC_LOG_INFO("40m/s wind + gravity: max_sway=" << max_sway
                  << " m at z=" << max_sway_z << " m");

    HVDC_LOG_INFO("=== Example 2 Complete ===");
    return 0;
}
