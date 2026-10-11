# Architecture Notes

The simulator is split into three small runtime responsibilities.

## Particle state

`Particle` is the shared data record for the simulation. It stores position, velocity, acceleration, mass, density, pressure, fixed-state information, and a render color.

## Neighbor search

`NeighborSearch` receives the current particle collection and builds a list of nearby particle indices for each particle. The search radius is the solver's smoothing radius. The solver reuses these lists while calculating local density and forces.

## SPH solver

`SPHSolver::step` advances the state in three conceptual stages:

1. Compute density and pressure from neighboring particles.
2. Compute pressure, viscosity, gravity, and boundary-related forces.
3. Integrate acceleration into velocity and position for the requested time step.

The solver owns the neighbor-search helper and keeps the physical parameters separate from rendering concerns.

## Renderer

`Renderer` owns the SFML window and translates particle state into drawable vertices. It also owns the interactive sliders and buttons used by the application loop. The renderer reports input state to `main.cpp`; it does not perform simulation calculations.

## Application loop

`main.cpp` creates the renderer and solver, builds the selected particle preset, then repeatedly:

- Handles window and control events.
- Applies changed control values to the solver.
- Handles reset input.
- Advances the solver using configured substeps.
- Draws the current particle state.

This arrangement keeps the simulation model testable without requiring an SFML window in every test.
