find_package(GTest QUIET)

if(NOT GTest_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
    )
    FetchContent_MakeAvailable(googletest)
endif()

add_executable(test_common
    test_main.cpp
    test_math_utils.cpp
    test_vector.cpp
    test_sparse_matrix.cpp
)
target_link_libraries(test_common PRIVATE hvdc_core gtest gtest_main)
add_test(NAME test_common COMMAND test_common)

add_executable(test_fem
    test_main.cpp
    test_node.cpp
    test_material.cpp
    test_beam_element_nl.cpp
    test_assembler.cpp
)
target_link_libraries(test_fem PRIVATE hvdc_core gtest gtest_main)
add_test(NAME test_fem COMMAND test_fem)

add_executable(test_fvm
    test_main.cpp
    test_grid.cpp
    test_field.cpp
    test_discretization.cpp
)
target_link_libraries(test_fvm PRIVATE hvdc_core gtest gtest_main)
add_test(NAME test_fvm COMMAND test_fvm)

add_executable(test_fsi
    test_main.cpp
    test_interface_mapper.cpp
    test_aerodynamic_loads.cpp
    test_coupling_manager.cpp
)
target_link_libraries(test_fsi PRIVATE hvdc_core gtest gtest_main)
add_test(NAME test_fsi COMMAND test_fsi)
