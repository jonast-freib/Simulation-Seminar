#include <cmath>
#include <iostream>
#include <cstddef>
#include <vector>
#include "Particle.h"
#include "Renderer.h"
#include "SPHSolver.h"

int main()
{
    const float particleRadius = 5.f; // visual radius only, physics uses smoothing radius

    const unsigned int windowWidth = 1700;
    const unsigned int windowHeight = 900;

    Renderer renderer(windowWidth, windowHeight, "SPH Fluid Solver", particleRadius);
    // Use a tuned SPH parameter set for stable motion in 2D.
    // Increase rest density to produce stronger pressure repulsion
    // and use a smaller timestep for better stability.
    SPHSolver sphSolver(32.0, 300.0, 2000.0, 200.0, Vector2D(0.0, 10.0));
    sphSolver.setBounds(windowWidth, windowHeight);

    // Add runtime sliders for tuning parameters
    const float sliderX = 20.f;
    float sliderY = 50.f;
    const float sliderGap = 40.f;
    renderer.addSlider("restDensity", 100.f, 800.f, 300.f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("gasConstant", 500.f, 5000.f, 2000.f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("viscosity", 0.f, 300.f, 200.f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("smoothingRadius", 10.f, 30.f, 13.f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("timestep", 0.0001f, 0.1f, 0.0007f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("gravity", 0.f, 20.f, 10.f, sliderX, sliderY);

    // Add reset and preset buttons
    const float buttonX = sliderX + 240.f;
    renderer.addButton("Reset", buttonX, 50.f, 80.f, 28.f);

    std::vector<Particle> particles;
    particles.reserve(2500);
    std::vector<Vector2D> initialPositions;


    auto saveInitialPositions = [&]()
    {
        initialPositions.clear();
        initialPositions.reserve(particles.size());
        for (const auto &p : particles)
        {
            initialPositions.push_back(p.position);
        }
    };

    struct SolverPreset
    {
        float restDensity;
        float gasConstant;
        float viscosity;
        float smoothingRadius;
        float timestep;
        float gravity;
    };

    auto applyPreset = [&](const SolverPreset &preset)
    {
        renderer.setSliderValue(0, preset.restDensity);
        renderer.setSliderValue(1, preset.gasConstant);
        renderer.setSliderValue(2, preset.viscosity);
        renderer.setSliderValue(3, preset.smoothingRadius);
        renderer.setSliderValue(4, preset.timestep);
        renderer.setSliderValue(5, preset.gravity);

        sphSolver.setRestDensity(preset.restDensity);
        sphSolver.setGasConstant(preset.gasConstant);
        sphSolver.setViscosity(preset.viscosity);
        sphSolver.setSmoothingRadius(preset.smoothingRadius);
        sphSolver.setGravity(Vector2D(0.0, preset.gravity));
    };

    auto getPhysicalSpacing = [](double smoothingRadius)
    {
        return static_cast<float>(smoothingRadius * 0.5);
    };

    auto addFixedBlock = [&](float cx, float cy, int cols, int rows, float spacing)
    {
        const double particleMass = sphSolver.getParticleMass();
        const float left = cx - (cols - 1) * 0.5f * spacing;
        for (int row = 0; row < rows; ++row)
        {
            for (int col = 0; col < cols; ++col)
            {
                particles.emplace_back(
                    Vector2D(left + col * spacing, cy + row * spacing),
                    particleMass,
                    true);
            }
        }
    };

    auto addUShapedBoundary = [&]()
    {
        const float physicalSpacing = getPhysicalSpacing(sphSolver.getSmoothingRadius());
        const float floorY = windowHeight - 200.f;
        const float floorLeft = windowWidth * 0.45f;
        const float floorRight = windowWidth * 0.55f;
        const int floorLength = 36;
        const int wallHeight = 12;

        addFixedBlock(
            floorLeft + floorLength * 0.5f * physicalSpacing,
            floorY,
            floorLength,
            2,
            physicalSpacing);
        addFixedBlock(
            floorLeft - 1.5f * physicalSpacing,
            floorY - 10.0f * physicalSpacing,
            2,
            wallHeight,
            physicalSpacing);
        addFixedBlock(
            floorRight + 1.5f * physicalSpacing,
            floorY - 10.0f * physicalSpacing,
            2,
            wallHeight,
            physicalSpacing);
    };

    auto addFluidCloud = [&](int columns, int rows, float top)
    {
        const double particleMass = sphSolver.getParticleMass();
        const float spacing = getPhysicalSpacing(sphSolver.getSmoothingRadius());
        const float left =
            windowWidth * 0.5f -
            (columns - 1) * 0.5f * spacing;

        for (int row = 0; row < rows; ++row)
        {
            for (int col = 0; col < columns; ++col)
            {
                particles.emplace_back(
                    Vector2D(
                        left + col * spacing +
                            (row % 2 ? spacing * 0.5f : 0.f),
                        top + row * spacing),
                    particleMass,
                    false);
            }
        }
    };

    auto makeDropOnBlock = [&](bool singleParticle)
    {
        applyPreset({300.f, 500.f, 200.f, 10.f, 0.005f, 10.f});
        particles.clear();
        addUShapedBoundary();
        if (singleParticle)
        {
            particles.emplace_back(
                Vector2D(windowWidth * 0.5f, windowHeight - 300.f),
                sphSolver.getParticleMass(),
                false);
        }
        else
        {
            addFluidCloud(5, 5, windowHeight - 250.f);
        }
        saveInitialPositions();
    };

    auto makeParticleCloud = [&]()
    {
        applyPreset({300.f, 500.f, 200.f, 32.f, 0.07f, 10.f});
        particles.clear();
        addFluidCloud(50, 20, 200.f);
        saveInitialPositions();
    };

    makeParticleCloud();

    // Main loop
    int frameCounter = 0;
    while (renderer.isOpen())
    {
        renderer.handleEvents();

        ++frameCounter;
        if (frameCounter >= 60)
        {
            // Read slider values and apply to solver at a reduced cadence.
            sphSolver.setRestDensity(static_cast<double>(renderer.getSliderValue(0)));
            sphSolver.setGasConstant(static_cast<double>(renderer.getSliderValue(1)));
            sphSolver.setViscosity(static_cast<double>(renderer.getSliderValue(2)));
            sphSolver.setSmoothingRadius(static_cast<double>(renderer.getSliderValue(3)));
            sphSolver.setGravity(Vector2D(0.0, static_cast<double>(renderer.getSliderValue(5))));
            const double currentParticleMass = sphSolver.getParticleMass();
            for (auto &particle : particles)
            {
                particle.mass = currentParticleMass;
            }
            frameCounter = 0;
        }

        // iterate over all particles and print all values if not fixed
        /*for (const auto& particle : particles) {
            if (particle.fixed) {
                std::cout << "Particle: pos(" << particle.position.x << ", " << particle.position.y << ") "
                          << "vel(" << particle.velocity.x << ", " << particle.velocity.y << ") "
                          << "acc(" << particle.acceleration.x << ", " << particle.acceleration.y << ") "
                          << "mass(" << particle.mass << ") "
                          << "density(" << particle.density << ") "
                          << "pressure(" << particle.pressure << ") "
                          << "cflViolated(" << (particle.cflViolated ? "true" : "false") << ")"
                          << std::endl;
            }
        }*/

        // Handle Reset and preset buttons
        if (renderer.consumeButtonPress(0))
        {
            for (std::size_t i = 0; i < particles.size(); ++i)
            {
                particles[i].position = initialPositions[i];
                particles[i].velocity = Vector2D(0.0, 0.0);
                particles[i].acceleration = Vector2D(0.0, 0.0);
            }
        }

        const double dt = static_cast<double>(renderer.getSliderValue(4));
        constexpr int substeps = 10;
        const double substepDt = dt / static_cast<double>(substeps);
        for (int i = 0; i < substeps; ++i)
        {
            sphSolver.step(particles, substepDt);
        }

        // Calculate and print max velocity, check CFL condition
        double maxVelocity = 0.0;
        const double smoothingRadius = sphSolver.getSmoothingRadius();
        const double cflThreshold = 0.4; // CFL stability threshold
        int cflViolations = 0;

        for (auto &p : particles)
        {
            const double speed = std::sqrt(p.velocity.x * p.velocity.x + p.velocity.y * p.velocity.y);
            maxVelocity = std::max(maxVelocity, speed);

            // Check CFL: velocity * dt / smoothing_radius
            const double cfl = (speed * substepDt) / smoothingRadius;
            p.cflViolated = (cfl > cflThreshold);
            if (p.cflViolated)
            {
                cflViolations++;
            }
        }

        renderer.clear();
        renderer.drawParticles(particles);
        renderer.display();
    }

    return 0;
}