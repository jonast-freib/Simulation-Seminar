cmake_minimum_required(VERSION 3.10)
project(SPHFluidSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The solver is far too slow unoptimised; default to an optimised build.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Find SFML
find_package(SFML 3 COMPONENTS Graphics REQUIRED)

include_directories(include)

add_executable(sph_solver
    src/main.cpp
    src/Scenes.cpp
    src/SPHSolver.cpp
    src/NeighborSearch.cpp
    src/Renderer.cpp
)

target_link_libraries(sph_solver PRIVATE SFML::Graphics)

# NeighborSearch unit test target
add_executable(neighbor_search_test
    tests/NeighborSearchTest.cpp
    src/NeighborSearch.cpp
)

target_link_libraries(neighbor_search_test PRIVATE SFML::Graphics)

add_executable(sph_solver_test
    tests/SPHSolverTest.cpp
    src/SPHSolver.cpp
    src/NeighborSearch.cpp
)

target_link_libraries(sph_solver_test PRIVATE SFML::Graphics)

enable_testing()
add_test(NAME NeighborSearchTest COMMAND neighbor_search_test)
add_test(NAME SPHSolverTest COMMAND sph_solver_test)
add_executable(sph_density_test
    tests/SPHDensityTest.cpp
    src/SPHSolver.cpp
    src/NeighborSearch.cpp
)

target_link_libraries(sph_density_test PRIVATE SFML::Graphics)
add_test(NAME SPHDensityTest COMMAND sph_density_test)

add_executable(hydrostatic_test
    tests/HydrostaticTest.cpp
    src/Scenes.cpp
    src/SPHSolver.cpp
    src/NeighborSearch.cpp
)

target_link_libraries(hydrostatic_test PRIVATE SFML::Graphics)
add_test(NAME HydrostaticTest COMMAND hydrostatic_test)

# Headless experiment runner that writes the CSV data for the paper's figures.
find_package(Threads REQUIRED)
add_executable(sph_experiments
    tools/experiments.cpp
    src/Scenes.cpp
    src/SPHSolver.cpp
    src/NeighborSearch.cpp
)
target_link_libraries(sph_experiments PRIVATE SFML::Graphics Threads::Threads)
