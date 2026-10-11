# SPH Fluid Solver

A small 2D smoothed-particle hydrodynamics (SPH) experiment written in C++17. The application uses SFML for rendering and exposes several solver parameters as interactive controls while the simulation is running.

## Features

- 2D particle simulation with density, pressure, viscosity, and gravity forces.
- Radius-based neighbor search for local particle interactions.
- Fixed particles for simple container boundaries.
- Runtime controls for solver parameters.
- A reset action that restores the active scene to its initial positions.
- A standalone neighbor-search test executable registered with CTest.

## Requirements

- CMake 3.10 or newer.
- A C++17 compiler.
- SFML 3 with the Graphics component available to CMake.
- A desktop environment capable of opening an SFML window.

## Build and run

From the repository root:

```sh
cmake -S . -B build
cmake --build build
./build/sph_solver
```

The exact executable location can vary by generator. Use `cmake --build build --target sph_solver` to rebuild only the application.

## Tests

```sh
ctest --test-dir build --output-on-failure
```

See [docs/BUILDING.md](docs/BUILDING.md) for a fuller setup walkthrough and [docs/TESTING.md](docs/TESTING.md) for the current test scope.

## Project layout

| Path | Responsibility |
| --- | --- |
| `src/main.cpp` | Application setup, presets, event loop, and simulation stepping |
| `src/SPHSolver.cpp` | SPH density, pressure, force, and integration calculations |
| `src/NeighborSearch.cpp` | Neighbor list construction |
| `src/Renderer.cpp` | SFML window, particle drawing, sliders, and buttons |
| `include/` | Public headers for the simulation components |
| `tests/` | Executable test programs |
| `docs/` | Developer and maintenance notes |
