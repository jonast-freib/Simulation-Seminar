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

Wall particles get a pressure extrapolated from the surrounding fluid (Adami, Hu & Adams 2012), so the fluid-wall pressure force has the same symmetric form as the fluid-fluid force. A short-range projection at half a particle spacing remains as a safety net against penetration.

The smoothing length is h = 1.3 dx (`SPHSolver::kernelSpacingRatio`); larger ratios trigger the pairing instability of the cubic spline, which shows up as fluid particles clumping into widely spaced horizontal rows.

Pressure follows the Tait equation of state, p = B((rho/rho0)^7 - 1) with B = rho0 c^2 / 7, clamped at zero. `SPHSolver::advance` advances a given duration in adaptive substeps limited by the acoustic CFL condition, the force condition and the viscous condition.

The solver owns the neighbor-search helper and keeps the physical parameters separate from rendering concerns.

## Renderer

`Renderer` owns the SFML window and translates particle state into drawable vertices. It also owns the interactive sliders and buttons used by the application loop. The renderer reports input state to `main.cpp`; it does not perform simulation calculations.

## Application loop

`Scenes.cpp` builds the test scenes (dam break, hydrostatic column): tank walls made of enough layers of fixed particles to fill the 2h kernel support (three at h = 1.3 dx), the fluid block, and the matching solver parameters. The speed of sound is chosen as ten times the expected maximum flow speed.

`main.cpp` creates the renderer and solver, loads the selected scene, then repeatedly:

- Handles window and control events.
- Applies changed control values to the solver.
- Handles reset input.
- Advances the solver by one frame of simulated time (adaptive substeps).
- Draws the current particle state.

This arrangement keeps the simulation model testable without requiring an SFML window in every test.
