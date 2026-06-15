#include "common/Types.hpp"
#include "common/MPIManager.hpp"
#include "common/PETScManager.hpp"
#include "common/Logger.hpp"
#include "common/Timer.hpp"
#include "fem/FEModel.hpp"
#include "fem/Assembler.hpp"
#include "fem/Solver.hpp"
#include "fvm/Grid.hpp"
#include "fvm/FluidSolver.hpp"
#include "fsi/InterfaceMapper.hpp"
#include "fsi/AerodynamicLoads.hpp"
#include "fsi/LoadTransfer.hpp"
#include "fsi/CouplingManager.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <cmath>
#include <cstring>

using namespace hvdc;
using namespace hvdc::fem;
using namespace hvdc::fvm;
using namespace hvdc::fsi;

struct CLIOptions {
    std::string mode = "help";
    std::string config_file;
    std::string output_file = "hvdc_results.vtk";
    Real wind_speed = 40.0;
    Real ice_thickness = 0.01;
    Real sim_time = 10.0;
    Real dt = 0.01;
    Index nx = 40;
    Index ny = 30;
    Index nz = 40;
    bool use_coupled = true;
    CouplingIterationMethod coupling_scheme = CouplingIterationMethod::Implicit_Aitken;
};

void print_help(const char* prog_name) {
    std::cout << "================================================================\n";
    std::cout << "HVDC Tower-Line System Nonlinear FSI Simulator v1.0.0\n";
    std::cout << "特高压输电塔线体系非线性流固耦合数值模拟内核\n";
    std::cout << "================================================================\n";
    std::cout << "用法:\n";
    std::cout << "  " << prog_name << " <command> [options]\n\n";
    std::cout << "命令:\n";
    std::cout << "  structural     纯结构有限元静力/模态分析\n";
    std::cout << "  fluid          纯流场稳态/瞬态求解\n";
    std::cout << "  coupled        非线性流固耦合仿真\n";
    std::cout << "  typhoon        典型台风-覆冰联合作用仿真\n";
    std::cout << "  tower          特高压塔独立静力分析\n";
    std::cout << "  conductor      覆冰导线气动弹性分析\n";
    std::cout << "  example        内置示例程序\n";
    std::cout << "  help           显示此帮助信息\n\n";
    std::cout << "选项:\n";
    std::cout << "  -w, --wind <speed>        风速 m/s (默认 40.0)\n";
    std::cout << "  -i, --ice <thickness>     覆冰厚度 m (默认 0.01)\n";
    std::cout << "  -t, --time <seconds>      仿真总时间 s (默认 10.0)\n";
    std::cout << "  -dt <step>                时间步长 s (默认 0.01)\n";
    std::cout << "  -o, --output <file>       输出文件 (默认 hvdc_results.vtk)\n";
    std::cout << "  -c, --config <file>       配置文件路径\n";
    std::cout << "  -nx, -ny, -nz <n>         流体网格尺寸 (默认 40x30x40)\n";
    std::cout << "  --scheme <name>           耦合方案: explicit|aitken|anderson|iqn-ils\n";
    std::cout << "  --uncoupled               关闭流固耦合\n";
    std::cout << "================================================================\n";
}

CLIOptions parse_args(int argc, char** argv) {
    CLIOptions opts;
    if (argc < 2) return opts;
    opts.mode = argv[1];

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-w" || arg == "--wind") && i + 1 < argc) {
            opts.wind_speed = std::atof(argv[++i]);
        } else if ((arg == "-i" || arg == "--ice") && i + 1 < argc) {
            opts.ice_thickness = std::atof(argv[++i]);
        } else if ((arg == "-t" || arg == "--time") && i + 1 < argc) {
            opts.sim_time = std::atof(argv[++i]);
        } else if (arg == "-dt" && i + 1 < argc) {
            opts.dt = std::atof(argv[++i]);
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            opts.output_file = argv[++i];
        } else if ((arg == "-c" || arg == "--config") && i + 1 < argc) {
            opts.config_file = argv[++i];
        } else if (arg == "-nx" && i + 1 < argc) {
            opts.nx = std::atoi(argv[++i]);
        } else if (arg == "-ny" && i + 1 < argc) {
            opts.ny = std::atoi(argv[++i]);
        } else if (arg == "-nz" && i + 1 < argc) {
            opts.nz = std::atoi(argv[++i]);
        } else if (arg == "--scheme" && i + 1 < argc) {
            std::string scheme = argv[++i];
            if (scheme == "explicit")      opts.coupling_scheme = CouplingIterationMethod::Explicit;
            else if (scheme == "aitken")   opts.coupling_scheme = CouplingIterationMethod::Implicit_Aitken;
            else if (scheme == "anderson") opts.coupling_scheme = CouplingIterationMethod::Implicit_Anderson;
            else if (scheme == "iqn-ils")  opts.coupling_scheme = CouplingIterationMethod::Implicit_IQN_ILS;
        } else if (arg == "--uncoupled") {
            opts.use_coupled = false;
        }
    }
    return opts;
}

void build_uhv_tower_model(FEModel& model, Real ice_thickness) {
    HVDC_LOG_INFO("构建特高压输电塔模型...");

    auto steel = std::make_shared<Material>(Material::Q420_Steel());
    model.add_material(steel);

    auto tower_leg = std::make_shared<BeamSection>(BeamSection::circular_pipe(0.60, 0.020, steel));
    auto tower_brace = std::make_shared<BeamSection>(BeamSection::circular_pipe(0.25, 0.010, steel));
    auto crossarm_main = std::make_shared<BeamSection>(BeamSection::i_beam(0.30, 0.40, 0.012, 0.018, steel));
    tower_leg->thickness_ice = ice_thickness;
    tower_brace->thickness_ice = ice_thickness;
    crossarm_main->thickness_ice = ice_thickness;
    model.add_section(tower_leg);
    model.add_section(tower_brace);
    model.add_section(crossarm_main);

    const Real H = 108.0;
    const int NLEV = 18;
    std::vector<Index> leg_nodes[4];
    for (int li = 0; li <= NLEV; ++li) {
        Real z = H * li / NLEV;
        Real w = 10.0 - (10.0 - 2.0) * li / NLEV;
        Vec3 corners[4] = {
            { w / 2,  w / 2, z},
            {-w / 2,  w / 2, z},
            {-w / 2, -w / 2, z},
            { w / 2, -w / 2, z}
        };
        for (int c = 0; c < 4; ++c) {
            Index nid = model.add_node(corners[c]);
            leg_nodes[c].push_back(nid);
            if (li == 0) {
                model.set_node_bc_all(nid, BoundaryType::Fixed);
            }
        }
    }
    for (int c = 0; c < 4; ++c) {
        for (size_t li = 0; li + 1 < leg_nodes[c].size(); ++li) {
            model.add_element(std::make_shared<BeamElementNL>(
                model.get_node(leg_nodes[c][li]),
                model.get_node(leg_nodes[c][li + 1]),
                tower_leg));
        }
    }
    for (int li = 0; li < NLEV; ++li) {
        for (int c = 0; c < 4; ++c) {
            int c2 = (c + 1) % 4;
            model.add_element(std::make_shared<BeamElementNL>(
                model.get_node(leg_nodes[c][li]),
                model.get_node(leg_nodes[c2][li]),
                tower_brace));
            if (li + 1 <= NLEV) {
                model.add_element(std::make_shared<BeamElementNL>(
                    model.get_node(leg_nodes[c][li]),
                    model.get_node(leg_nodes[c2][li + 1]),
                    tower_brace));
            }
        }
    }
    std::vector<Index> arm_nodes;
    for (int arm_lev : {10, 14, 17}) {
        Real z = H * arm_lev / NLEV;
        Real w_at_lev = 10.0 - (10.0 - 2.0) * arm_lev / NLEV;
        for (Real off : {-15.0, -10.0, 10.0, 15.0}) {
            Vec3 p1 = {w_at_lev / 2 + off,  w_at_lev / 2, z};
            Vec3 p2 = {w_at_lev / 2 + off, -w_at_lev / 2, z};
            Index a1 = model.add_node(p1);
            Index a2 = model.add_node(p2);
            arm_nodes.push_back(a1);
            arm_nodes.push_back(a2);
            Index cn1 = leg_nodes[0][arm_lev];
            Index cn2 = leg_nodes[3][arm_lev];
            model.add_element(std::make_shared<BeamElementNL>(
                model.get_node(cn1), model.get_node(a1), crossarm_main));
            model.add_element(std::make_shared<BeamElementNL>(
                model.get_node(cn2), model.get_node(a2), crossarm_main));
        }
    }
    HVDC_LOG_INFO("特高压塔模型: " << model.num_nodes() << " 节点, "
                  << model.num_elements() << " 单元");
}

void build_conductor_line_model(FEModel& model, Real ice_thickness,
                                Real span = 400.0, Real sag = 12.0) {
    HVDC_LOG_INFO("构建覆冰架空导线模型 (跨距 " << span << "m, 垂度 " << sag << "m)...");
    auto acsr = std::make_shared<Material>(Material::ACSR_Conductor());
    model.add_material(acsr);
    auto cond_sec = std::make_shared<BeamSection>(BeamSection::circular_pipe(0.0336, 0.0112, acsr));
    cond_sec->thickness_ice = ice_thickness;
    model.add_section(cond_sec);

    Index nL = model.add_node({0.0, 0.0, 70.0});
    Index nR = model.add_node({span, 0.0, 70.0 - sag});
    model.set_node_bc_all(nL, BoundaryType::Fixed);
    model.set_node_bc_all(nR, BoundaryType::Fixed);

    const int NSEG = 40;
    std::vector<Index> seg_nodes;
    seg_nodes.push_back(nL);
    for (int i = 1; i < NSEG; ++i) {
        Real xi = static_cast<Real>(i) / NSEG;
        Real xc = span * xi;
        Real yc = 0.0;
        Real zc = 70.0 - 4.0 * sag * xi * (1.0 - xi);
        seg_nodes.push_back(model.add_node({xc, yc, zc}));
    }
    seg_nodes.push_back(nR);

    for (size_t i = 0; i + 1 < seg_nodes.size(); ++i) {
        model.add_element(std::make_shared<BeamElementNL>(
            model.get_node(seg_nodes[i]),
            model.get_node(seg_nodes[i + 1]),
            cond_sec));
    }
    HVDC_LOG_INFO("架空导线模型: " << model.num_nodes() << " 节点, "
                  << model.num_elements() << " 单元");
}

void build_full_tower_line_system(FEModel& model, Real ice_thickness) {
    build_uhv_tower_model(model, ice_thickness);
    build_conductor_line_model(model, ice_thickness, 400.0, 12.0);
}

int run_structural_analysis(const CLIOptions& opts) {
    HVDC_LOG_INFO("=== 执行纯结构有限元静力分析 ===");
    auto model = std::make_shared<FEModel>();
    build_full_tower_line_system(*model, opts.ice_thickness);
    model->setup_dofs();
    HVDC_LOG_INFO("总自由度: " << model->total_dofs());

    Assembler assembler;
    Solver solver(model.get());
    SparseMatrix K;
    Vector F_ext(model->total_dofs());
    Vector disp(model->total_dofs());
    disp.zero(); F_ext.zero();

    assembler.assemble_tangent_stiffness(disp, K, model.get());
    HVDC_LOG_INFO("切线刚度矩阵非零元: " << K.num_nonzeros());
    model->apply_body_force(F_ext);
    model->apply_ice_load(F_ext);
    assembler.apply_dirichlet_bc(K, F_ext, disp, *model);
    Vector du(model->total_dofs());
    solver.solve_linear_system(K, du, F_ext);
    for (Index i = 0; i < disp.size(); ++i) disp.add(i, du[i]);

    Real max_d = 0.0;
    Index max_node = -1;
    for (Index i = 0; i < model->num_nodes(); ++i) {
        const auto* n = model->node(i);
        if (!n) continue;
        Index s = n->dof_start();
        Vec3 d = {0, 0, 0};
        for (int k = 0; k < 3; ++k)
            if (s + k < disp.size()) d[k] = disp[s + k];
        Real norm = math::vec3_norm(d);
        if (norm > max_d) { max_d = norm; max_node = i; }
    }
    HVDC_LOG_INFO("最大节点位移: " << max_d << " m (节点 " << max_node << ")");
    return 0;
}

int run_coupled_simulation(const CLIOptions& opts) {
    HVDC_LOG_INFO("=== 执行非线性流固耦合仿真 ===");
    HVDC_LOG_INFO("  风速: " << opts.wind_speed << " m/s, 覆冰: "
                  << opts.ice_thickness * 1000 << " mm");

    auto struct_model = std::make_shared<FEModel>();
    build_full_tower_line_system(*struct_model, opts.ice_thickness);

    AerodynamicConfig aero_cfg;
    aero_cfg.wind_direction[0] = opts.wind_speed;
    aero_cfg.rho_air = 1.225;
    aero_cfg.mu_air = 1.789e-5;
    aero_cfg.power_law_exponent = 0.14;
    aero_cfg.reference_height = 80.0;
    aero_cfg.compute_buffeting = true;
    aero_cfg.turbulence_intensity = 0.12;
    aero_cfg.compute_ice_effects = (opts.ice_thickness > EPS);
    aero_cfg.ice_drag_coefficient_increase = 1.4;

    auto aero = std::make_shared<AerodynamicLoads>(struct_model.get(), nullptr, aero_cfg);
    struct_model->setup_dofs();
    Index ndofs = struct_model->total_dofs();
    HVDC_LOG_INFO("结构自由度: " << ndofs);

    auto fluid_grid = std::make_shared<Grid>();
    Real xmin = -100.0, xmax = 600.0;
    Real ymin = -200.0, ymax = 200.0;
    Real zmin = 0.0, zmax = 180.0;
    fluid_grid->create_structured_hex(opts.nx, opts.ny, opts.nz,
                                      {xmin, ymin, zmin}, {xmax, ymax, zmax});
    fluid_grid->set_boundary_patch(0, 0, opts.nx * opts.ny * opts.nz);

    auto fluid_solver = std::make_shared<FluidSolver>(fluid_grid.get());
    fluid_solver->set_free_stream_velocity({opts.wind_speed, 0.0, 0.0});
    fluid_solver->initialize_fields();

    auto mapper = std::make_shared<InterfaceMapper>(struct_model.get(), fluid_grid.get());
    mapper->build_mapping(3.0);

    auto transfer = std::make_shared<LoadTransfer>(struct_model.get(), fluid_grid.get(),
                                                    mapper.get(), aero.get());

    auto assembler = std::make_shared<Assembler>();
    auto solver = std::make_shared<Solver>(struct_model.get());

    SimulationConfig sim_cfg;
    sim_cfg.dt = opts.dt;
    sim_cfg.t_start = 0.0;
    sim_cfg.t_end = opts.sim_time;
    sim_cfg.coupling_scheme = CouplingScheme::StronglyImplicit;
    sim_cfg.max_coupling_iter = 30;
    sim_cfg.coupling_tol = 1e-5;
    sim_cfg.max_fluid_iter = 20;
    sim_cfg.fluid_res_tol = 1e-4;
    sim_cfg.analysis_type = AnalysisType::Transient;
    sim_cfg.beta = 0.25;
    sim_cfg.gamma = 0.5;

    CouplingManager coupling(
        struct_model.get(), assembler.get(), solver.get(),
        fluid_solver.get(), mapper.get(), aero.get(), transfer.get(), sim_cfg);
    coupling.set_coupling_scheme(opts.coupling_scheme);
    coupling.set_max_coupling_iterations(sim_cfg.max_coupling_iter);
    coupling.set_convergence_tolerance(sim_cfg.coupling_tol);

    Vector disp(ndofs), vel(ndofs), acc(ndofs);
    disp.zero(); vel.zero(); acc.zero();

    coupling.run_coupled_simulation(disp, vel, acc, fluid_solver->current_state(),
        [&](Index step, Real t, const Vector& d, const fvm::FluidState& s) {
            if (step % 100 == 0) {
                Real max_d = 0.0;
                for (Index i = 0; i < d.size() && i < ndofs; i += 3) {
                    Vec3 v = {d[i], d[i+1], d[i+2]};
                    max_d = std::max(max_d, math::vec3_norm(v));
                }
                HVDC_LOG_INFO("  进度 step=" << step << ", t=" << t
                              << "s, max_disp=" << max_d << "m");
            }
        });

    return 0;
}

int main(int argc, char** argv) {
    MPIManager mpi(&argc, &argv);
    PETScManager petsc(&argc, &argv);

    const char* log_file = (mpi.rank() == 0) ? "hvdc_sim.log" : nullptr;
    Logger::instance().set_level(LogLevel::INFO);
    Logger::instance().set_file_output(log_file);
    Logger::instance().set_mpi_rank(mpi.rank());

    CLIOptions opts = parse_args(argc, argv);

    if (opts.mode == "help") {
        if (mpi.rank() == 0) print_help(argv[0]);
        return 0;
    }

    if (mpi.rank() == 0) {
        HVDC_LOG_INFO("================================================================");
        HVDC_LOG_INFO("特高压输电塔线体系非线性流固耦合数值模拟内核 (HVDC-FSI)");
        HVDC_LOG_INFO("MPI 进程数: " << mpi.size() << ", 本进程: " << mpi.rank());
        HVDC_LOG_INFO("计算模式: " << opts.mode);
        HVDC_LOG_INFO("================================================================");
    }

    int status = 0;

    if (opts.mode == "structural" || opts.mode == "tower" || opts.mode == "conductor") {
        status = run_structural_analysis(opts);
    } else if (opts.mode == "coupled" || opts.mode == "typhoon" || opts.mode == "example") {
        status = run_coupled_simulation(opts);
    } else if (opts.mode == "fluid") {
        auto grid = std::make_shared<Grid>();
        grid->create_structured_hex(opts.nx, opts.ny, opts.nz,
                                     {-100, -100, 0}, {500, 100, 150});
        auto fs = std::make_shared<FluidSolver>(grid.get());
        fs->set_free_stream_velocity({opts.wind_speed, 0, 0});
        fs->initialize_fields();
        HVDC_LOG_INFO("流体网格初始化完成: "
                      << grid->num_cells() << " 单元, "
                      << grid->num_faces() << " 面, "
                      << grid->num_nodes() << " 节点");
        fs->solve_steady(500, 1e-5);
    } else {
        if (mpi.rank() == 0) {
            std::cerr << "未知命令: " << opts.mode << std::endl;
            print_help(argv[0]);
        }
        status = 1;
    }

    if (mpi.rank() == 0) {
        HVDC_LOG_INFO("================================================================");
        HVDC_LOG_INFO("计算完成，退出状态: " << status);
        HVDC_LOG_INFO("================================================================");
    }

    return status;
}
