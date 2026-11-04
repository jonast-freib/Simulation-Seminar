#pragma once

#include <string>
#include <vector>
#include "Particle.h"
#include "SPHSolver.h"

// Physical setup shared by the test scenes.
struct SceneParameters {
    double smoothingRadius;
    double restDensity;
    double speedOfSound;
    double viscosity;
    double gravity;
};

// Rectangular tank: inner region [left, left + width] x [top, top + height]
// in window coordinates (y points down). Walls are layers of fixed particles
// outside the inner region (left, right, bottom and optionally a lid).
struct Tank {
    double left;
    double top;
    double width;
    double height;
    bool closedTop;
};

struct Scene {
    std::string name;
    SceneParameters parameters;
    Tank tank;
    // Height of the initial fluid body; used to choose the speed of sound.
    double fluidHeight;
    std::vector<Particle> particles;
};

// Dam break: a fluid column in the lower-left corner of a closed tank collapses.
Scene makeDamBreakScene(double windowWidth, double windowHeight);

// Hydrostatic column: the tank bottom is covered with fluid at rest.
// Density should stay close to rho0 and pressure should grow linearly with depth.
Scene makeHydrostaticScene(double windowWidth, double windowHeight);

// Push the scene parameters into the solver.
void applySceneParameters(SPHSolver& solver, const SceneParameters& parameters);
