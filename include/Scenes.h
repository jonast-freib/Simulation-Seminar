#pragma once

#include <string>
#include <vector>
#include "Particle.h"
#include "SPHSolver.h"

// Physical setup shared by the test scenes.
struct SceneParameters
{
    double smoothingRadius;
    double restDensity;
    double stiffness;
    double viscosity;
    double gravity;
    double kernelSpacingRatio;
    bool mirrorBoundaryPressure;
};

struct SceneOptions
{
    double particleSpacing = 10.0; // initial spacing
    double kernelSpacingRatio = SPHSolver::defaultKernelSpacingRatio;
    double stiffnessScale = 1.0;
    double restDensity = 300.0;
    double viscosity = 600.0;
    double gravity = 200.0;
    bool mirrorBoundaryPressure = true;
};

struct Tank
{
    double left;
    double top;
    double width;
    double height;
    bool closedTop;
};

struct Scene
{
    std::string name;
    SceneParameters parameters;
    Tank tank;
    double fluidWidth;
    double fluidHeight;
    std::vector<Particle> particles;
};

Scene makeDamBreakScene(double windowWidth, double windowHeight, const SceneOptions &options = SceneOptions());

Scene makeHydrostaticScene(double windowWidth, double windowHeight, const SceneOptions &options = SceneOptions());

void applySceneParameters(SPHSolver &solver, const SceneParameters &parameters);
