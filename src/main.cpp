#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>
#include "Particle.h"
#include "Renderer.h"
#include "SPHSolver.h"
#include "Scenes.h"

int main()
{
    const float particleRadius = 5.f; // visual radius only, physics uses smoothing radius

    const unsigned int windowWidth = 1700;
    const unsigned int windowHeight = 900;

    // Simulated time advanced per rendered frame. The solver splits it into
    // as many substeps as the stability limits require.
    constexpr double frameTime = 1.0 / 60.0;

    Renderer renderer(windowWidth, windowHeight, "SPH Fluid Solver", particleRadius);
    SPHSolver sphSolver;
    sphSolver.setBounds(windowWidth, windowHeight);

    // Debug sliders for tuning parameters at runtime.
    enum SliderIndex
    {
        RestDensitySlider,
        SpeedOfSoundSlider,
        ViscositySlider,
        SmoothingRadiusSlider,
        GravitySlider
    };
    const float sliderX = 20.f;
    float sliderY = 50.f;
    const float sliderGap = 40.f;
    renderer.addSlider("restDensity", 100.f, 800.f, 300.f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("speedOfSound", 500.f, 10000.f, 4000.f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("viscosity", 0.f, 1000.f, 150.f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("smoothingRadius", 5.f, 40.f, 13.f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("gravity", 0.f, 1000.f, 200.f, sliderX, sliderY);

    enum ButtonIndex
    {
        ResetButton,
        DamBreakButton,
        ColumnButton,
        ColorButton
    };
    const float buttonX = sliderX + 240.f;
    renderer.addButton("Reset", buttonX, 50.f, 80.f, 28.f);
    renderer.addButton("Dam", buttonX, 90.f, 80.f, 28.f);
    renderer.addButton("Column", buttonX, 130.f, 80.f, 28.f);
    renderer.addButton("Color", buttonX, 170.f, 80.f, 28.f);

    enum class SceneKind
    {
        DamBreak,
        HydrostaticColumn
    };
    SceneKind activeScene = SceneKind::DamBreak;
    std::vector<Particle> particles;
    double simulationTime = 0.0;

    auto loadScene = [&](SceneKind kind)
    {
        Scene scene = (kind == SceneKind::DamBreak)
                          ? makeDamBreakScene(windowWidth, windowHeight)
                          : makeHydrostaticScene(windowWidth, windowHeight);

        applySceneParameters(sphSolver, scene.parameters);
        renderer.setSliderValue(RestDensitySlider, static_cast<float>(scene.parameters.restDensity));
        renderer.setSliderValue(SpeedOfSoundSlider, static_cast<float>(scene.parameters.speedOfSound));
        renderer.setSliderValue(ViscositySlider, static_cast<float>(scene.parameters.viscosity));
        renderer.setSliderValue(SmoothingRadiusSlider, static_cast<float>(scene.parameters.smoothingRadius));
        renderer.setSliderValue(GravitySlider, static_cast<float>(scene.parameters.gravity));

        // Colour scale: hydrostatic pressure at the bottom of the initial fluid body.
        renderer.setReferencePressure(scene.parameters.restDensity * scene.parameters.gravity * scene.fluidHeight);
        renderer.resetColorSmoothing();

        particles = std::move(scene.particles);
        activeScene = kind;
        simulationTime = 0.0;
        std::printf("Loaded scene '%s': %zu particles, c = %.0f\n",
                    scene.name.c_str(), particles.size(), scene.parameters.speedOfSound);
    };

    loadScene(SceneKind::DamBreak);

    int frameCounter = 0;
    int substepsSinceReport = 0;
    while (renderer.isOpen())
    {
        renderer.handleEvents();

        if (renderer.consumeButtonPress(ResetButton))
        {
            loadScene(activeScene);
        }
        if (renderer.consumeButtonPress(DamBreakButton))
        {
            loadScene(SceneKind::DamBreak);
        }
        if (renderer.consumeButtonPress(ColumnButton))
        {
            loadScene(SceneKind::HydrostaticColumn);
        }
        if (renderer.consumeButtonPress(ColorButton))
        {
            // Toggle fluid colouring: pressure map <-> plain blue (walls stay grey).
            renderer.setColorByPressure(!renderer.isColorByPressure());
        }

        ++frameCounter;
        if (frameCounter >= 60)
        {
            // Read slider values and apply to solver at a reduced cadence.
            sphSolver.setRestDensity(renderer.getSliderValue(RestDensitySlider));
            sphSolver.setSpeedOfSound(renderer.getSliderValue(SpeedOfSoundSlider));
            sphSolver.setViscosity(renderer.getSliderValue(ViscositySlider));
            sphSolver.setSmoothingRadius(renderer.getSliderValue(SmoothingRadiusSlider));
            sphSolver.setGravity(Vector2D(0.0, renderer.getSliderValue(GravitySlider)));
            const double currentParticleMass = sphSolver.getParticleMass();
            for (auto &particle : particles)
            {
                particle.mass = currentParticleMass;
            }

            // Diagnostics: the density error should stay around 1% for c = 10 * v_max.
            double maxDensityError = 0.0;
            double maxSpeed = 0.0;
            for (const auto &p : particles)
            {
                if (p.fixed)
                {
                    continue;
                }
                maxDensityError = std::max(maxDensityError, p.density / sphSolver.getRestDensity() - 1.0);
                maxSpeed = std::max(maxSpeed, std::hypot(p.velocity.x, p.velocity.y));
            }
            std::printf("t = %6.2f s | substeps/frame = %5.1f | dt = %.2e | max density error = %5.2f %% | "
                        "max speed = %6.1f (Mach %.3f)\n",
                        simulationTime,
                        substepsSinceReport / static_cast<double>(frameCounter),
                        sphSolver.getLastTimeStep(),
                        100.0 * maxDensityError,
                        maxSpeed,
                        maxSpeed / sphSolver.getSpeedOfSound());

            frameCounter = 0;
            substepsSinceReport = 0;
        }

        substepsSinceReport += sphSolver.advance(particles, frameTime);
        simulationTime += frameTime;

        renderer.clear();
        renderer.drawParticles(particles);
        renderer.display();
    }

    return 0;
}
