#include <iostream>
#include <vector>
#include "Particle.h"
#include "Renderer.h"
#include "SPHSolver.h"

int main()
{
    const int particlesX = 10;
    const int particlesY = 10;
    const float particleRadius = 5.f; // visual radius only, physics uses smoothing radius

    const unsigned int windowWidth = 1700;
    const unsigned int windowHeight = 900;

    Renderer renderer(windowWidth, windowHeight, "SPH Fluid Solver", particleRadius);
    // Use a tuned SPH parameter set for stable motion in 2D.
    // Increase rest density to produce stronger pressure repulsion
    // and use a smaller timestep for better stability.
    SPHSolver sphSolver(16.0, 300.0, 2000.0, 200.0, Vector2D(0.0, 10.0));
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
    particles.reserve(particlesX * particlesY + 200);
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

    auto getPhysicalSpacing = [&](double smoothingRadius)
    {
        return static_cast<float>(smoothingRadius * 0.5);
    };

    auto applyPresetSliders = [&](float restDensity, float gasConstant, float viscosity, float smoothingRadius, float timestep, float gravity)
    {
        renderer.setSliderValue(0, restDensity);
        renderer.setSliderValue(1, gasConstant);
        renderer.setSliderValue(2, viscosity);
        renderer.setSliderValue(3, smoothingRadius);
        renderer.setSliderValue(4, timestep);
        renderer.setSliderValue(5, gravity);

        sphSolver.setRestDensity(restDensity);
        sphSolver.setGasConstant(gasConstant);
        sphSolver.setViscosity(viscosity);
        sphSolver.setSmoothingRadius(smoothingRadius);
        sphSolver.setGravity(Vector2D(0.0, gravity));
    };

    const float boundaryTop = windowHeight - 180.f;
    const float boundaryLeft = windowWidth * 0.5f;

    auto makeDropOnBlock = [&]()
    {
        applyPresetSliders(300.f, 500.f, 200.f, 10.f, 0.005f, 10.f);
        particles.clear();
        const double particleMass = sphSolver.getParticleMass();
        const float physicalSpacing = getPhysicalSpacing(sphSolver.getSmoothingRadius());

        const float floorY = windowHeight - 200.f;
        const float floorLeft = windowWidth * 0.45f;
        const float floorRight = windowWidth * 0.55f;
        const int floorLength = 18 * 2;
        const int wallHeight = 12;

        auto addFixedBlock = [&](float cx, float cy, int cols, int rows)
        {
            const float left = cx - (cols - 1) * 0.5f * physicalSpacing;
            for (int row = 0; row < rows; ++row)
            {
                for (int col = 0; col < cols; ++col)
                {
                    const float x = left + col * physicalSpacing;
                    const float y = cy + row * physicalSpacing;
                    particles.emplace_back(Vector2D(x, y), particleMass, true);
                }
            }
        };

        // U-shaped boundary: left wall, right wall, and base floor
        addFixedBlock(floorLeft + (floorLength * 0.5f * physicalSpacing), floorY, floorLength, 2);
        addFixedBlock(floorLeft - 1.5f * physicalSpacing, floorY - 10.0f * physicalSpacing, 2, wallHeight);
        addFixedBlock(floorRight + 1.5f * physicalSpacing, floorY - 10.0f * physicalSpacing, 2, wallHeight);

        const int cloudCols = 5;
        const int cloudRows = 5;
        const float cloudLeft = windowWidth * 0.5f - (cloudCols - 1) * 0.5f * physicalSpacing;
        const float cloudTop = windowHeight - 250.f;
        for (int row = 0; row < cloudRows; ++row)
        {
            for (int col = 0; col < cloudCols; ++col)
            {
                const float x = cloudLeft + col * physicalSpacing + (row % 2 ? physicalSpacing * 0.5f : 0.f);
                const float y = cloudTop + row * physicalSpacing;
                particles.emplace_back(Vector2D(x, y), particleMass, false);
            }
        }
        saveInitialPositions();
    };

    auto makeDropOnBlockSingle = [&]()
    {
        applyPresetSliders(300.f, 500.f, 200.f, 10.f, 0.005f, 10.f);
        particles.clear();
        const double particleMass = sphSolver.getParticleMass();
        const float physicalSpacing = getPhysicalSpacing(sphSolver.getSmoothingRadius());

        const float floorY = windowHeight - 200.f;
        const float floorLeft = windowWidth * 0.45f;
        const float floorRight = windowWidth * 0.55f;
        const int floorLength = 18 * 2;
        const int wallHeight = 12;

        auto addFixedBlock = [&](float cx, float cy, int cols, int rows)
        {
            const float left = cx - (cols - 1) * 0.5f * physicalSpacing;
            for (int row = 0; row < rows; ++row)
            {
                for (int col = 0; col < cols; ++col)
                {
                    const float x = left + col * physicalSpacing;
                    const float y = cy + row * physicalSpacing;
                    particles.emplace_back(Vector2D(x, y), particleMass, true);
                }
            }
        };

        // U-shaped boundary: left wall, right wall, and base floor
        addFixedBlock(floorLeft + (floorLength * 0.5f * physicalSpacing), floorY, floorLength, 2);
        addFixedBlock(floorLeft - 1.5f * physicalSpacing, floorY - 10.0f * physicalSpacing, 2, wallHeight);
        addFixedBlock(floorRight + 1.5f * physicalSpacing, floorY - 10.0f * physicalSpacing, 2, wallHeight);

        // add a single particle above the block
        const float singleParticleX = windowWidth * 0.5f;
        const float singleParticleY = windowHeight - 300.f;
        particles.emplace_back(Vector2D(singleParticleX, singleParticleY), particleMass, false);
        saveInitialPositions();
    };

    auto makeParticleCloud = [&]()
    {
        applyPresetSliders(300.f, 500.f, 200.f, 14.f, 0.007f, 10.f);
        particles.clear();
        const double particleMass = sphSolver.getParticleMass();
        const float physicalSpacing = getPhysicalSpacing(sphSolver.getSmoothingRadius());
        const int cloudCols = 5;
        const int cloudRows = 5;
        const float cloudLeft = windowWidth * 0.5f - (cloudCols - 1) * 0.5f * physicalSpacing;
        const float cloudTop = 800.f;
        for (int row = 0; row < cloudRows; ++row)
        {
            for (int col = 0; col < cloudCols; ++col)
            {
                const float x = cloudLeft + col * physicalSpacing + (row % 2 ? physicalSpacing * 0.5f : 0.f);
                const float y = cloudTop + row * physicalSpacing;
                particles.emplace_back(Vector2D(x, y), particleMass, false);
            }
        }
        saveInitialPositions();
    };

    // makeDropOnBlock();
    // makeParticleCloud();
    makeDropOnBlockSingle();

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
        constexpr int substeps = 40;
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