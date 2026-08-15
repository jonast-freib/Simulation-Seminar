#include <iostream>
#include <vector>
#include "Particle.h"
#include "Renderer.h"
#include "SPHSolver.h"


int main() {
    const int particlesX = 10;
    const int particlesY = 10;
    const float particleRadius = 4.f; // visual radius only, physics uses smoothing radius

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

    auto saveInitialPositions = [&]() {
        initialPositions.clear();
        initialPositions.reserve(particles.size());
        for (const auto& p : particles) {
            initialPositions.push_back(p.position);
        }
    };

    auto getPhysicalSpacing = [&](double smoothingRadius) {
        return static_cast<float>(smoothingRadius * 0.7);
    };

    auto applyPresetSliders = [&](float restDensity, float gasConstant, float viscosity, float smoothingRadius, float timestep, float gravity) {
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

    auto makeDropOnBlock = [&]() {
        applyPresetSliders(300.f, 2000.f, 200.f, 10.f, 0.05f, 10.f);
        particles.clear();
        const double particleMass = sphSolver.getParticleMass();
        const float physicalSpacing = getPhysicalSpacing(sphSolver.getSmoothingRadius());
        const int boundaryRows = 4;
        const int boundaryCols = 28;
        const float boundaryLeftPhys = boundaryLeft - (boundaryCols - 1) * 0.5f * physicalSpacing;

        for (int row = 0; row < boundaryRows; ++row) {
            for (int col = 0; col < boundaryCols; ++col) {
                const float x = boundaryLeftPhys + col * physicalSpacing;
                const float y = boundaryTop + row * physicalSpacing;
                particles.emplace_back(Vector2D(x, y), particleMass, true);
            }
        }
        particles.emplace_back(Vector2D(windowWidth * 0.5f, boundaryTop - 2.5f * physicalSpacing), particleMass, false);
        saveInitialPositions();
    };

    auto makeParticleCloud = [&]() {
        applyPresetSliders(300.f, 2000.f, 200.f, 12.f, 0.0007f, 10.f);
        particles.clear();
        const double particleMass = sphSolver.getParticleMass();
        const float physicalSpacing = getPhysicalSpacing(sphSolver.getSmoothingRadius());
        const int cloudCols = 5;
        const int cloudRows = 5;
        const float cloudLeft = windowWidth * 0.5f - (cloudCols - 1) * 0.5f * physicalSpacing;
        const float cloudTop = 800.f;
        for (int row = 0; row < cloudRows; ++row) {
            for (int col = 0; col < cloudCols; ++col) {
                const float x = cloudLeft + col * physicalSpacing + (row % 2 ? physicalSpacing * 0.5f : 0.f);
                const float y = cloudTop + row * physicalSpacing;
                particles.emplace_back(Vector2D(x, y), particleMass, false);
            }
        }
        saveInitialPositions();
    };

    auto makeWalls = [&]() {
        applyPresetSliders(400.f, 4000.f, 250.f, 16.f, 0.004f, 10.f);
        particles.clear();

        const double particleMass = sphSolver.getParticleMass();
        const float physicalSpacing = getPhysicalSpacing(sphSolver.getSmoothingRadius());
        const int wallDepth = 3;
        const int wallHeight = 30;
        const float wallTop = 800.f;
        const float leftWallX = windowWidth * 0.35f;
        const float rightWallX = windowWidth * 0.65f;

        auto addWall = [&](float wallCenterX) {
            for (int depth = 0; depth < wallDepth; ++depth) {
                for (int i = 0; i < wallHeight; ++i) {
                    const float x = wallCenterX + (depth - (wallDepth - 1) * 0.5f) * physicalSpacing - 2.f;
                    const float y = wallTop + i * physicalSpacing - 2.f;
                    particles.emplace_back(Vector2D(x, y), particleMass, true);
                }
            }
        };

        addWall(leftWallX);
        addWall(rightWallX);

        const int cloudCols = 4;
        const int cloudRows = 4;
        const float cloudLeft = windowWidth * 0.5f - (cloudCols - 1) * 0.5f * physicalSpacing;
        const float cloudTop = 700.f;

        for (int row = 0; row < cloudRows; ++row) {
            for (int col = 0; col < cloudCols; ++col) {
                const float x = cloudLeft + col * physicalSpacing + (row % 2 ? physicalSpacing * 0.5f : 0.f);
                const float y = cloudTop + row * physicalSpacing;
                particles.emplace_back(Vector2D(x, y), particleMass, false);
                
            }
        }

        saveInitialPositions();
    };

    auto makeHighGravityImpact = [&]() {
        applyPresetSliders(300.f, 2000.f, 200.f, 16.f, 0.0007f, 20.f);
        particles.clear();
        const double particleMass = sphSolver.getParticleMass();
        const float physicalSpacing = getPhysicalSpacing(sphSolver.getSmoothingRadius());
        const int boundaryRows = 4;
        const int boundaryCols = 28;
        const float boundaryLeftPhys = boundaryLeft - (boundaryCols - 1) * 0.5f * physicalSpacing;

        for (int row = 0; row < boundaryRows; ++row) {
            for (int col = 0; col < boundaryCols; ++col) {
                const float x = boundaryLeftPhys + col * physicalSpacing;
                const float y = boundaryTop + row * physicalSpacing;
                particles.emplace_back(Vector2D(x, y), particleMass, true);
            }
        }
        particles.emplace_back(Vector2D(windowWidth * 0.5f, boundaryTop - 4.0f * physicalSpacing), particleMass, false);
        saveInitialPositions();
    };

    makeDropOnBlock();
    //makeParticleCloud();
    

    // Main loop
    int frameCounter = 0;
    while (renderer.isOpen()) {
        renderer.handleEvents();

        ++frameCounter;
        if (frameCounter >= 60) {
            // Read slider values and apply to solver at a reduced cadence.
            sphSolver.setRestDensity(static_cast<double>(renderer.getSliderValue(0)));
            sphSolver.setGasConstant(static_cast<double>(renderer.getSliderValue(1)));
            sphSolver.setViscosity(static_cast<double>(renderer.getSliderValue(2)));
            sphSolver.setSmoothingRadius(static_cast<double>(renderer.getSliderValue(3)));
            sphSolver.setGravity(Vector2D(0.0, static_cast<double>(renderer.getSliderValue(5))));
            const double currentParticleMass = sphSolver.getParticleMass();
            for (auto& particle : particles) {
            particle.mass = currentParticleMass;
            }
            frameCounter = 0;
        }

        

        // Handle Reset and preset buttons
        if (renderer.consumeButtonPress(0)) {
            for (std::size_t i = 0; i < particles.size(); ++i) {
                particles[i].position = initialPositions[i];
                particles[i].velocity = Vector2D(0.0, 0.0);
                particles[i].acceleration = Vector2D(0.0, 0.0);
            }
        }

        const double dt = static_cast<double>(renderer.getSliderValue(4));
        constexpr int substeps = 10;
        const double substepDt = dt / static_cast<double>(substeps);
        for (int i = 0; i < substeps; ++i) {
            sphSolver.step(particles, substepDt);
        }

        // Calculate and print max velocity, check CFL condition
        double maxVelocity = 0.0;
        const double smoothingRadius = sphSolver.getSmoothingRadius();
        const double cflThreshold = 0.4;  // CFL stability threshold
        int cflViolations = 0;

        for (auto& p : particles) {
            const double speed = std::sqrt(p.velocity.x * p.velocity.x + p.velocity.y * p.velocity.y);
            maxVelocity = std::max(maxVelocity, speed);

            // Check CFL: velocity * dt / smoothing_radius
            const double cfl = (speed * substepDt) / smoothingRadius;
            p.cflViolated = (cfl > cflThreshold);
            if (p.cflViolated) {
                cflViolations++;
            }
        }
        std::cout << "Frame " << frameCounter << " | Max velocity: " << maxVelocity 
                  << " | CFL violations: " << cflViolations << std::endl;

        renderer.clear();
        renderer.drawParticles(particles);
        renderer.display();
    }

    return 0;
}