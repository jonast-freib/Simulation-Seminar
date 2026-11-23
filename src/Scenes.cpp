#include "Scenes.h"

#include <cmath>

namespace
{

    int wallLayerCount(double h, double spacing)
    {
        return static_cast<int>(std::ceil(2.0 * h / spacing));
    }

    double referenceStiffness(double rho0, double g, double fluidHeight)
    {
        const double maxFlowSpeed = std::sqrt(2.0 * g * fluidHeight);
        const double waveSpeed = 10.0 * maxFlowSpeed;
        return rho0 * waveSpeed * waveSpeed / 7.0;
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

    void addTankWalls(std::vector<Particle> &particles, const Tank &tank,
                      double h, double spacing, double mass)
    {
        const int wallLayers = wallLayerCount(h, spacing);
        const double right = tank.left + tank.width;
        const double bottom = tank.top + tank.height;
        const double outer = (wallLayers - 0.5) * spacing;

        const int sideRows = static_cast<int>(std::floor(tank.height / spacing));
        const int horizontalColumns =
            static_cast<int>(std::round((tank.width + 2.0 * outer) / spacing)) + 1;

        for (int layer = 0; layer < wallLayers; ++layer)
        {
            const double offset = (layer + 0.5) * spacing;

            for (int col = 0; col < horizontalColumns; ++col)
            {
                const double x = tank.left - outer + col * spacing;
                particles.emplace_back(Vector2D(x, bottom + offset), mass, true);
                if (tank.closedTop)
                {
                    particles.emplace_back(Vector2D(x, tank.top - offset), mass, true);
                }
            }

            for (int row = 0; row < sideRows; ++row)
            {
                const double y = bottom - (row + 0.5) * spacing;
                particles.emplace_back(Vector2D(tank.left - offset, y), mass, true);
                particles.emplace_back(Vector2D(right + offset, y), mass, true);
            }
        }
    }

    void addFluidBlock(std::vector<Particle> &particles, double x0, double bottom,
                       double width, double height, double spacing, double mass)
    {
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
                        double fluidWidth, double fluidHeight, const SceneOptions &options)
    {
        const double spacing = options.particleSpacing;
        const double h = options.kernelSpacingRatio * spacing;

        Scene scene;
        scene.name = name;
        scene.tank = tank;
        scene.fluidWidth = fluidWidth;
        scene.fluidHeight = fluidHeight;
        scene.parameters = SceneParameters{
            h,
            options.restDensity,
            options.stiffnessScale * referenceStiffness(options.restDensity, options.gravity, fluidHeight),
            options.viscosity,
            options.gravity,
            options.kernelSpacingRatio,
            options.mirrorBoundaryPressure};

        const double mass = options.restDensity * spacing * spacing;
        addTankWalls(scene.particles, tank, h, spacing, mass);
        addFluidBlock(scene.particles, tank.left, tank.top + tank.height,
                      fluidWidth, fluidHeight, spacing, mass);
        return scene;
    }
} // namespace

Scene makeDamBreakScene(double windowWidth, double windowHeight, const SceneOptions &options)
{

    const Tank tank = centeredTank(windowWidth, windowHeight, 800.0, 520.0, true);
    return makeTankScene("Dam break", tank, 200.0, 400.0, options);
}

Scene makeHydrostaticScene(double windowWidth, double windowHeight, const SceneOptions &options)
{
    const Tank tank = centeredTank(windowWidth, windowHeight, 300.0, 520.0, false);
    return makeTankScene("Hydrostatic column", tank, 300.0, 400.0, options);
}

void applySceneParameters(SPHSolver &solver, const SceneParameters &parameters)
{
    solver.setKernelSpacingRatio(parameters.kernelSpacingRatio);
    solver.setBoundaryPressureMirroring(parameters.mirrorBoundaryPressure);
    solver.setRestDensity(parameters.restDensity);
    solver.setSmoothingRadius(parameters.smoothingRadius);
    solver.setStiffness(parameters.stiffness);
    solver.setViscosity(parameters.viscosity);
    solver.setGravity(Vector2D(0.0, parameters.gravity));
}
