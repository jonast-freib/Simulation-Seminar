cmake_minimum_required(VERSION 3.10)
project(SPHFluidSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find SFML
find_package(SFML 3 COMPONENTS Graphics REQUIRED)

include_directories(include)

add_executable(sph_solver
    src/main.cpp
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

enable_testing()
add_test(NAME NeighborSearchTest COMMAND neighbor_search_test)