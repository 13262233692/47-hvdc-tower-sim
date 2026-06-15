#include "common/Types.hpp"
#include "common/MPIManager.hpp"
#include "common/PETScManager.hpp"
#include "common/Logger.hpp"
#include "fem/FEModel.hpp"
#include "fem/Assembler.hpp"
#include "fem/Solver.hpp"
#include "fvm/Grid.hpp"
#include "fvm/FluidSolver.hpp"
#include "fsi/InterfaceMapper.hpp"
#include "fsi/AerodynamicLoads.hpp"
#include "fsi/LoadTransfer.hpp"
#include "fsi/CouplingManager.hpp"
#include <memory>
#include <cmath>

using namespace hvdc;
using namespace hvdc::fem;
using namespace hvdc::fvm;
using namespace hvdc::fsi;

int main(int argc, char** argv) {
    MPIManager mpi(&argc, &argv);
    PETScManager petsc(&argc, &argv);
    Logger::instance().set_level(LogLevel::INFO);
    Logger::instance().set_mpi_rank(mpi.rank());

    HVDC_LOG_INFO("=== Example 3: Typhoon-Ice Coupled FSI Simulation ===");

    const Real WIND_SPEED = 50.0;
    const Real ICE_THICKNESS = 0.015;
    const Real SIM_TIME = 5.0;
    const Real DT = 0.005;

    auto struct_model = std::make_shared<FEModel>();
    auto steel = std::make_shared<Material>(Material::Q420_Steel());
    auto acsr = std::make_shared<Material>(Material::ACSR_Conductor());
    struct_model->add_material(steel);
    struct_model->add_material(acsr);

    auto leg_sec = std::make_shared<BeamSection>(BeamSection::circular_pipe(0.60, 0.020, steel));
    auto brace_sec = std::make_shared<BeamSection>(BeamSection::circular_pipe(0.25, 0.010, steel));
    auto cond_sec = std::make_shared<BeamSection>(BeamSection::circular_pipe(0.0336, 0.0112, acsr));
    leg_sec->thickness_ice = ICE_THICKNESS;
    brace_sec->thickness_ice = ICE_THICKNESS;
    cond_sec->thickness_ice = ICE_THICKNESS;
    struct_model->add_section(leg_sec);
    struct_model->add_section(brace_sec);
    struct_model->add_section(cond_sec);

    const Real H = 108.0;
    const int NLEV = 18;
    std::vector<Index> legs[4];
    for (int li = 0; li <= NLEV; ++li) {
        Real z = H * li / NLEV;
        Real w = 10.0 - 8.0 * li / NLEV;
        Vec3 corners[4] = {{w/2,w/2,z},{-w/2,w/2,z},{-w/2,-w/2,z},{w/2,-w/2,z}};
        for (int c = 0; c < 4; ++c) {
            Index nid = struct_model->add_node(corners[c]);
            legs[c].push_back(nid);
            if (li == 0) struct_model->set_node_bc_all(nid, BoundaryType::Fixed);
        }
    }
    for (int c = 0; c < 4; ++c)
        for (size_t i = 0; i+1 < legs[c].size(); ++i)
            struct_model->add_element(std::make_shared<BeamElementNL>(
                struct_model->get_node(legs[c][i]), struct_model->get_node(legs[c][i+1]), leg_sec));
    for (int li = 0; li < NLEV; ++li)
        for (int c = 0; c < 4; ++c) {
            int c2 = (c+1)%4;
            struct_model->add_element(std::make_shared<BeamElementNL>(
                struct_model->get_node(legs[c][li]), struct_model->get_node(legs[c2][li]), brace_sec));
            struct_model->add_element(std::make_shared<BeamElementNL>(
                struct_model->get_node(legs[c][li]), struct_model->get_node(legs[c2][li+1]), brace_sec));
        }

    Real span = 300.0, sag = 9.0;
    std::vector<Index> cnodes;
    int NCS = 30;
    for (int i = 0; i <= NCS; ++i) {
        Real xi = static_cast<Real>(i) / NCS;
        cnodes.push_back(struct_model->add_node(
            {-50.0 + span * xi, 0.0, 80.0 - 4.0 * sag * xi * (1.0 - xi)}));
    }
    struct_model->set_node_bc_all(cnodes.front(), BoundaryType::Fixed);
    struct_model->set_node_bc_all(cnodes.back(), BoundaryType::Fixed);
    for (size_t i = 0; i+1 < cnodes.size(); ++i)
        struct_model->add_element(std::make_shared<BeamElementNL>(
            struct_model->get_node(cnodes[i]), struct_model->get_node(cnodes[i+1]), cond_sec));

    struct_model->setup_dofs();
    Index ndofs = struct_model->total_dofs();
    HVDC_LOG_INFO("Structural model: " << struct_model->num_nodes() << " nodes, "
                  << struct_model->num_elements() << " elements, DOFs=" << ndofs);

    AerodynamicConfig aero_cfg;
    aero_cfg.wind_direction[0] = WIND_SPEED;
    aero_cfg.wind_direction[1] = 0.0;
    aero_cfg.wind_direction[2] = 0.0;
    aero_cfg.gust_amplitude = 0.3 * WIND_SPEED;
    aero_cfg.gust_frequency = 0.5;
    aero_cfg.reference_height = 80.0;
    aero_cfg.power_law_exponent = 0.14;
    aero_cfg.compute_buffeting = true;
    aero_cfg.turbulence_intensity = 0.15;
    aero_cfg.compute_ice_effects = true;
    aero_cfg.ice_drag_coefficient_increase = 1.5;
    auto aero = std::make_shared<AerodynamicLoads>(struct_model.get(), nullptr, aero_cfg);

    auto fluid_grid = std::make_shared<Grid>();
    Index nx = 30, ny = 25, nz = 30;
    fluid_grid->create_structured_hex(nx, ny, nz,
        {-150, -150, 0}, {350, 150, 180});
    fluid_grid->set_boundary_patch(0, 0, nx * ny * nz);

    auto fluid_solver = std::make_shared<FluidSolver>(fluid_grid.get());
    fluid_solver->set_free_stream_velocity({WIND_SPEED, 0.0, 0.0});
    fluid_solver->initialize_fields();
    HVDC_LOG_INFO("Fluid grid: " << fluid_grid->num_cells() << " cells, "
                  << fluid_grid->num_nodes() << " nodes");

    auto mapper = std::make_shared<InterfaceMapper>(struct_model.get(), fluid_grid.get());
    mapper->build_mapping(3.0);
    HVDC_LOG_INFO("FSI interface: struct_nodes=" << mapper->num_struct_interface_nodes()
                  << ", fluid_faces=" << mapper->num_fluid_interface_faces());

    auto transfer = std::make_shared<LoadTransfer>(struct_model.get(), fluid_grid.get(),
                                                    mapper.get(), aero.get());
    auto assembler = std::make_shared<Assembler>();
    auto solver = std::make_shared<Solver>(struct_model.get());

    SimulationConfig sim_cfg;
    sim_cfg.dt = DT;
    sim_cfg.t_start = 0.0;
    sim_cfg.t_end = SIM_TIME;
    sim_cfg.coupling_scheme = CouplingScheme::StronglyImplicit;
    sim_cfg.max_coupling_iter = 20;
    sim_cfg.coupling_tol = 1e-5;
    sim_cfg.max_fluid_iter = 15;
    sim_cfg.fluid_res_tol = 1e-4;
    sim_cfg.analysis_type = AnalysisType::Transient;
    sim_cfg.beta = 0.25;
    sim_cfg.gamma = 0.5;

    CouplingManager coupling(
        struct_model.get(), assembler.get(), solver.get(),
        fluid_solver.get(), mapper.get(), aero.get(), transfer.get(), sim_cfg);
    coupling.set_coupling_scheme(CouplingIterationMethod::Implicit_Aitken);
    coupling.set_max_coupling_iterations(sim_cfg.max_coupling_iter);
    coupling.set_convergence_tolerance(sim_cfg.coupling_tol);

    Vector disp(ndofs), vel(ndofs), acc(ndofs);
    disp.zero(); vel.zero(); acc.zero();

    coupling.run_coupled_simulation(disp, vel, acc, fluid_solver->current_state(),
        [&](Index step, Real t, const Vector& d, const fvm::FluidState& s) {
            if (step % 50 == 0) {
                Real max_tower = 0.0, max_cond = 0.0;
                for (Index i = 0; i+2 < d.size() && i < (legs[0].size() + legs[1].size() + legs[2].size() + legs[3].size()) * 8; i += 3) {
                    Real m = std::sqrt(d[i]*d[i] + d[i+1]*d[i+1] + d[i+2]*d[i+2]);
                    if (m > max_tower) max_tower = m;
                }
                Index cn0 = cnodes[0] * 6;
                Index cnn = cnodes.back() * 6;
                for (Index gi = cn0; gi + 2 < d.size() && gi < cnn + 6; gi += 6) {
                    Real m = std::sqrt(d[gi]*d[gi] + d[gi+1]*d[gi+1] + d[gi+2]*d[gi+2]);
                    if (m > max_cond) max_cond = m;
                }
                HVDC_LOG_INFO("  step=" << step << " t=" << t
                              << "s, max_tower_disp=" << max_tower
                              << "m, max_cond_disp=" << max_cond << "m");
            }
            (void)s;
        });

    HVDC_LOG_INFO("=== Example 3 Complete ===");
    return 0;
}
