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
    double kernelSpacingRatio;
};

// Knobs for building a scene (window units: pixels and seconds). The defaults
// are the setup used by the interactive app; experiments vary them.
struct SceneOptions {
    double particleSpacing = 10.0;     // initial spacing dx
    double kernelSpacingRatio = SPHSolver::defaultKernelSpacingRatio; // h / dx
    double speedOfSoundFactor = 10.0;  // c = factor * sqrt(2 g H)
    double restDensity = 300.0;
    double viscosity = 150.0;
    double gravity = 200.0;
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
    // Initial fluid body, placed in the lower-left corner of the tank.
    double fluidWidth;
    double fluidHeight;
    std::vector<Particle> particles;
};

// Dam break: a fluid column in the lower-left corner of a closed tank collapses.
Scene makeDamBreakScene(double windowWidth, double windowHeight, const SceneOptions& options = SceneOptions());

// Hydrostatic column: the tank bottom is covered with fluid at rest.
// Density should stay close to rho0 and pressure should grow linearly with depth.
Scene makeHydrostaticScene(double windowWidth, double windowHeight, const SceneOptions& options = SceneOptions());

// Push the scene parameters into the solver.
void applySceneParameters(SPHSolver& solver, const SceneParameters& parameters);
