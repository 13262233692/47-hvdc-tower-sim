cmake_minimum_required(VERSION 3.16)
project(hvdc_tower_sim
    VERSION 1.0.0
    LANGUAGES CXX
    DESCRIPTION "UHV Transmission Tower-Line System Nonlinear FSI Numerical Simulation Kernel"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(ENABLE_MPI "Enable MPI parallel computing" ON)
option(ENABLE_PETSC "Enable PETSc for large-scale sparse linear solver" ON)
option(ENABLE_OPENMP "Enable OpenMP shared-memory parallelism" ON)
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_EXAMPLES "Build example cases" ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")

if(MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

if(ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    message(STATUS "MPI found: ${MPI_CXX_VERSION}")
    add_compile_definitions(HVDC_USE_MPI)
endif()

if(ENABLE_PETSC)
    set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
    find_package(PETSc REQUIRED)
    message(STATUS "PETSc found: ${PETSC_VERSION}")
    add_compile_definitions(HVDC_USE_PETSC)
endif()

if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
        add_compile_definitions(HVDC_USE_OPENMP)
    endif()
endif()

find_package(Eigen3 3.3 QUIET)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 found: ${Eigen3_VERSION}")
    add_compile_definitions(HVDC_USE_EIGEN)
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)

set(HVDC_SOURCES
    src/common/Types.cpp
    src/common/MathUtils.cpp
    src/common/MPIManager.cpp
    src/common/PETScManager.cpp
    src/common/Timer.cpp
    src/common/Logger.cpp
    src/fem/Node.cpp
    src/fem/Material.cpp
    src/fem/Section.cpp
    src/fem/Element.cpp
    src/fem/BeamElementNL.cpp
    src/fem/TrussElementNL.cpp
    src/fem/ConductorElement.cpp
    src/fem/FEModel.cpp
    src/fem/Assembler.cpp
    src/fem/Solver.cpp
    src/fvm/Cell.cpp
    src/fvm/Face.cpp
    src/fvm/Grid.cpp
    src/fvm/Field.cpp
    src/fvm/Discretization.cpp
    src/fvm/NSComposable.cpp
    src/fvm/FluidSolver.cpp
    src/fvm/BoundaryCondition.cpp
    src/fsi/InterfaceMapper.cpp
    src/fsi/AerodynamicLoads.cpp
    src/fsi/CouplingManager.cpp
    src/fsi/LoadTransfer.cpp
)

set(HVDC_HEADERS
    include/common/Types.hpp
    include/common/MathUtils.hpp
    include/common/MPIManager.hpp
    include/common/PETScManager.hpp
    include/common/Timer.hpp
    include/common/Logger.hpp
    include/common/SparseMatrix.hpp
    include/common/Vector.hpp
    include/fem/Node.hpp
    include/fem/Material.hpp
    include/fem/Section.hpp
    include/fem/Element.hpp
    include/fem/BeamElementNL.hpp
    include/fem/TrussElementNL.hpp
    include/fem/ConductorElement.hpp
    include/fem/FEModel.hpp
    include/fem/Assembler.hpp
    include/fem/Solver.hpp
    include/fvm/Cell.hpp
    include/fvm/Face.hpp
    include/fvm/Grid.hpp
    include/fvm/Field.hpp
    include/fvm/Discretization.hpp
    include/fvm/NSComposable.hpp
    include/fvm/FluidSolver.hpp
    include/fvm/BoundaryCondition.hpp
    include/fsi/InterfaceMapper.hpp
    include/fsi/AerodynamicLoads.hpp
    include/fsi/CouplingManager.hpp
    include/fsi/LoadTransfer.hpp
)

add_library(hvdc_core STATIC ${HVDC_SOURCES})

target_include_directories(hvdc_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(MPI_CXX_FOUND)
    target_link_libraries(hvdc_core PUBLIC MPI::MPI_CXX)
endif()

if(PETSC_FOUND)
    target_include_directories(hvdc_core PUBLIC ${PETSC_INCLUDES})
    target_link_libraries(hvdc_core PUBLIC ${PETSC_LIBRARIES})
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(hvdc_core PUBLIC OpenMP::OpenMP_CXX)
endif()

if(Eigen3_FOUND)
    target_link_libraries(hvdc_core PUBLIC Eigen3::Eigen)
endif()

add_executable(hvdc_sim src/main.cpp)
target_link_libraries(hvdc_sim PRIVATE hvdc_core)

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

message(STATUS "")
message(STATUS "========================================")
message(STATUS "HVDC Tower Simulation Configuration:")
message(STATUS "  MPI:      ${ENABLE_MPI}")
message(STATUS "  PETSc:    ${ENABLE_PETSC}")
message(STATUS "  OpenMP:   ${ENABLE_OPENMP}")
message(STATUS "  Tests:    ${BUILD_TESTS}")
message(STATUS "  Examples: ${BUILD_EXAMPLES}")
message(STATUS "========================================")
message(STATUS "")
