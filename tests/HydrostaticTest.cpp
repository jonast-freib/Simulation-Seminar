#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "Scenes.h"

static void expect(bool condition, const char *message)
{
    if (!condition)
    {
        std::fprintf(stderr, "Test failed: %s\n", message);
        std::exit(EXIT_FAILURE);
    }
}

// The Tait EOS with B = rho0 c^2 / 7 must give dp/drho = c^2 at rest density.
static void testTaitSpeedOfSound()
{
    SPHSolver solver(20.0, 300.0, 1000.0, 0.0, Vector2D(0.0, 0.0));
    const double rho0 = 300.0;
    const double eps = 1e-4 * rho0;
    const double dpdrho = (solver.pressureFromDensity(rho0 + eps) - solver.pressureFromDensity(rho0)) / eps;
    expect(std::abs(dpdrho - 1000.0 * 1000.0) < 1e-3 * 1000.0 * 1000.0, "Tait EOS slope equals c^2 at rho0");
    expect(solver.pressureFromDensity(0.9 * rho0) == 0.0, "negative pressure is clamped");
}

// Run the hydrostatic column and check density error, containment, the
// adaptive timestep and the hydrostatic pressure profile p = rho0 * g * depth.
static void testHydrostaticColumn()
{
    Scene scene = makeHydrostaticScene(1700.0, 900.0);
    SPHSolver solver;
    applySceneParameters(solver, scene.parameters);
    solver.setBounds(1700.0, 900.0);

    const SceneParameters &params = scene.parameters;
    const double cflLimit = 0.4 * params.smoothingRadius / params.speedOfSound;
    const double frameTime = 1.0 / 60.0;
    const int frames = 180; // 3 s of simulated time

    for (int frame = 0; frame < frames; ++frame)
    {
        const int substeps = solver.advance(scene.particles, frameTime);
        expect(substeps >= 1, "advance takes at least one substep");
        expect(solver.getLastTimeStep() <= cflLimit * (1.0 + 1e-9), "substep respects the acoustic CFL limit");
    }

    const Tank &tank = scene.tank;
    const double bottom = tank.top + tank.height;
    double maxDensityError = 0.0;
    double sumPressure = 0.0;
    double sumExpected = 0.0;
    for (const Particle &p : scene.particles)
    {
        if (p.fixed)
        {
            continue;
        }
        expect(p.position.x > tank.left && p.position.x < tank.left + tank.width &&
                   p.position.y < bottom,
               "fluid stays inside the tank");
        maxDensityError = std::max(maxDensityError, p.density / params.restDensity - 1.0);

        // Compare pressure away from the walls and the free surface.
        const double depth = bottom - p.position.y;
        const double fromSide = std::min(p.position.x - tank.left, tank.left + tank.width - p.position.x);
        if (depth > 2.0 * params.smoothingRadius && fromSide > 2.0 * params.smoothingRadius &&
            depth < scene.fluidHeight - 4.0 * params.smoothingRadius)
        {
            sumPressure += p.pressure;
            const double surfaceY = bottom - scene.fluidHeight;
            sumExpected += params.restDensity * params.gravity * (p.position.y - surfaceY);
        }
    }

    const double pressureRatio = sumPressure / sumExpected;
    std::printf("max density error = %.2f %%, mean p / (rho0 g d) = %.3f\n",
                100.0 * maxDensityError, pressureRatio);
    expect(maxDensityError < 0.03, "density error stays below 3%");
    expect(std::abs(pressureRatio - 1.0) < 0.15, "pressure follows rho0 * g * depth within 15%");
}

int main()
{
    testTaitSpeedOfSound();
    testHydrostaticColumn();
    return 0;
}
