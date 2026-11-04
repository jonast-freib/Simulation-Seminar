#include "Scenes.h"

#include <cmath>

namespace
{
// Default physical setup (window units: pixels and seconds).
constexpr double smoothingRadius = 13.0; // particle spacing 10 px at h / dx = 1.3
constexpr double restDensity = 300.0;
constexpr double gravity = 200.0;
constexpr double viscosity = 150.0;

double particleSpacing(double h)
{
    return SPHSolver::particleSpacingFor(h);
}

double particleMass(double rho0, double h)
{
    const double spacing = particleSpacing(h);
    return rho0 * spacing * spacing;
}

// Number of fixed-particle layers so the walls fill the full kernel support
// (2h) of a fluid particle touching the wall.
int wallLayerCount(double h)
{
    return static_cast<int>(std::ceil(2.0 * h / particleSpacing(h)));
}

// Weakly compressible SPH: c = 10 * v_max keeps density errors around 1%.
// The fastest flow is roughly free fall from the top of the fluid body.
double chooseSpeedOfSound(double g, double fluidHeight)
{
    const double maxFlowSpeed = std::sqrt(2.0 * g * fluidHeight);
    return 10.0 * maxFlowSpeed;
}

Tank centeredTank(double windowWidth, double windowHeight, double width, double height, bool closedTop)
{
    const double bottomMargin = 80.0;
    return Tank{
        0.5 * (windowWidth - width),
        windowHeight - bottomMargin - height,
        width,
        height,
        closedTop};
}

// The wall surface coincides with the tank edge: fluid particles sit half a
// spacing inside it and the first wall layer half a spacing outside it, so
// neighbouring fluid and wall particles are one spacing apart.
void addTankWalls(std::vector<Particle> &particles, const Tank &tank, double h, double mass)
{
    const double spacing = particleSpacing(h);
    const int wallLayers = wallLayerCount(h);
    const double right = tank.left + tank.width;
    const double bottom = tank.top + tank.height;
    const double outer = (wallLayers - 0.5) * spacing;

    const int sideRows = static_cast<int>(std::floor(tank.height / spacing));
    const int horizontalColumns =
        static_cast<int>(std::round((tank.width + 2.0 * outer) / spacing)) + 1;

    for (int layer = 0; layer < wallLayers; ++layer)
    {
        const double offset = (layer + 0.5) * spacing;

        // Bottom wall (and lid), extended over the side walls to close the corners.
        for (int col = 0; col < horizontalColumns; ++col)
        {
            const double x = tank.left - outer + col * spacing;
            particles.emplace_back(Vector2D(x, bottom + offset), mass, true);
            if (tank.closedTop)
            {
                particles.emplace_back(Vector2D(x, tank.top - offset), mass, true);
            }
        }

        // Side walls.
        for (int row = 0; row < sideRows; ++row)
        {
            const double y = bottom - (row + 0.5) * spacing;
            particles.emplace_back(Vector2D(tank.left - offset, y), mass, true);
            particles.emplace_back(Vector2D(right + offset, y), mass, true);
        }
    }
}

// Fill [x0, x0 + width] x [bottom - height, bottom] with fluid on a square lattice
// whose outermost particles sit half a spacing inside the region.
void addFluidBlock(std::vector<Particle> &particles, double x0, double bottom,
                   double width, double height, double h, double mass)
{
    const double spacing = particleSpacing(h);
    const int columns = static_cast<int>(std::floor(width / spacing));
    const int rows = static_cast<int>(std::floor(height / spacing));

    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < columns; ++col)
        {
            particles.emplace_back(
                Vector2D(x0 + (col + 0.5) * spacing, bottom - (row + 0.5) * spacing),
                mass,
                false);
        }
    }
}

Scene makeTankScene(const std::string &name, const Tank &tank,
                    double fluidWidth, double fluidHeight)
{
    Scene scene;
    scene.name = name;
    scene.tank = tank;
    scene.fluidHeight = fluidHeight;
    scene.parameters = SceneParameters{
        smoothingRadius,
        restDensity,
        chooseSpeedOfSound(gravity, fluidHeight),
        viscosity,
        gravity};

    const double mass = particleMass(restDensity, smoothingRadius);
    addTankWalls(scene.particles, tank, smoothingRadius, mass);
    addFluidBlock(scene.particles, tank.left, tank.top + tank.height,
                  fluidWidth, fluidHeight, smoothingRadius, mass);
    return scene;
}
} // namespace

Scene makeDamBreakScene(double windowWidth, double windowHeight)
{
    const Tank tank = centeredTank(windowWidth, windowHeight, 800.0, 520.0, true);
    return makeTankScene("Dam break", tank, 200.0, 400.0);
}

Scene makeHydrostaticScene(double windowWidth, double windowHeight)
{
    const Tank tank = centeredTank(windowWidth, windowHeight, 300.0, 520.0, false);
    return makeTankScene("Hydrostatic column", tank, 300.0, 400.0);
}

void applySceneParameters(SPHSolver &solver, const SceneParameters &parameters)
{
    solver.setRestDensity(parameters.restDensity);
    solver.setSmoothingRadius(parameters.smoothingRadius);
    solver.setSpeedOfSound(parameters.speedOfSound);
    solver.setViscosity(parameters.viscosity);
    solver.setGravity(Vector2D(0.0, parameters.gravity));
}
