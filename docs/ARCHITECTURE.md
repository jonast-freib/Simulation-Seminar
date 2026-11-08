# Architecture Notes

The simulator is split into three small runtime responsibilities.

## Particle state

`Particle` is the shared data record for the simulation. It stores position, velocity, acceleration, mass, density, pressure, fixed-state information, and a render color.

## Neighbor search

`NeighborSearch` receives the current particle collection and builds a list of nearby particle indices for each particle. The search radius is twice the solver's smoothing radius (the support of the cubic spline kernel). The solver reuses these lists while calculating local density and forces.

## SPH solver

`SPHSolver::step` advances the state in three conceptual stages:

1. Compute density and pressure from neighboring particles.
2. Compute pressure, viscosity, gravity, and boundary-related forces.
3. Integrate acceleration into velocity and position for the requested time step.

Wall particles (boundary samples) contribute to the density like fluid at rest density. In the pressure acceleration their pressure and density are mirrored from the fluid particle (p_ib = p_i, rho_ib = rho_i), as in the course notes. Mirroring is not energy-consistent and feeds energy into the fluid near walls, so the default viscosity is 600; `setBoundaryPressureMirroring(false)` switches to p_ib = 0 for comparison. A short-range projection at half a particle spacing remains as a safety net against penetration.

The smoothing length is h = 1.3 dx (`SPHSolver::kernelSpacingRatio`); larger ratios trigger the pairing instability of the cubic spline, which shows up as fluid particles clumping into widely spaced horizontal rows.

Pressure follows the state equation of the course notes, p_i = max(k((rho_i/rho_0)^7 - 1), 0), with the stiffness coefficient k. `SPHSolver::advance` advances a given duration in adaptive substeps limited by the CFL condition dt <= lambda h / (sqrt(7k/rho_0) + |v_max|) with lambda = 0.4, the force condition and the viscous condition.

The solver owns the neighbor-search helper and keeps the physical parameters separate from rendering concerns.

## Renderer

`Renderer` owns the SFML window and translates particle state into drawable vertices. It also owns the interactive sliders and buttons used by the application loop. The renderer reports input state to `main.cpp`; it does not perform simulation calculations.

## Application loop

`Scenes.cpp` builds the test scenes (dam break, hydrostatic column): tank walls made of enough layers of fixed particles to fill the 2h kernel support (three at h = 1.3 dx), the fluid block, and the matching solver parameters. The stiffness coefficient k is chosen so that pressure waves, sqrt(7 k / rho_0), are ten times faster than the expected maximum flow speed.

`main.cpp` creates the renderer and solver, loads the selected scene, then repeatedly:

- Handles window and control events.
- Applies changed control values to the solver.
- Handles reset input.
- Advances the solver by one frame of simulated time (adaptive substeps).
- Draws the current particle state.

This arrangement keeps the simulation model testable without requiring an SFML window in every test.
