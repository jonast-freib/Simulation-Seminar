#include <iostream>
#include <vector>
#include "Particle.h"
#include "Renderer.h"
#include "SPHSolver.h"


int main() {
    const int particlesX = 10;
    const int particlesY = 10;
    const float particleRadius = 6.f; // visual radius only, physics uses smoothing radius

    const unsigned int windowWidth = 1800;
    const unsigned int windowHeight = 1000;

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
    renderer.addSlider("smoothingRadius", 10.f, 40.f, 16.f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("timestep", 0.0001f, 0.01f, 0.0007f, sliderX, sliderY);
    sliderY += sliderGap;
    renderer.addSlider("gravity", 0.f, 20.f, 10.f, sliderX, sliderY);

    // Add reset and preset buttons
    const float buttonX = sliderX + 240.f;
    renderer.addButton("Reset", buttonX, 50.f, 80.f, 28.f);
    renderer.addButton("Drop Block", buttonX, 90.f, 100.f, 28.f);
    renderer.addButton("Particle Cloud", buttonX, 130.f, 120.f, 28.f);
    renderer.addButton("High Gravity", buttonX, 170.f, 100.f, 28.f);

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
        return static_cast<float>(smoothingRadius * 0.75);
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
        applyPresetSliders(300.f, 2000.f, 200.f, 16.f, 0.0007f, 10.f);
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
        applyPresetSliders(300.f, 2000.f, 200.f, 16.f, 0.0007f, 10.f);
        particles.clear();
        const double particleMass = sphSolver.getParticleMass();
        const float physicalSpacing = getPhysicalSpacing(sphSolver.getSmoothingRadius());
        const int cloudCols = 14;
        const int cloudRows = 8;
        const float cloudLeft = windowWidth * 0.5f - (cloudCols - 1) * 0.5f * physicalSpacing;
        const float cloudTop = 120.f;
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

    // Main loop
    while (renderer.isOpen()) {
        renderer.handleEvents();

        // Read slider values and apply to solver
        sphSolver.setRestDensity(static_cast<double>(renderer.getSliderValue(0)));
        sphSolver.setGasConstant(static_cast<double>(renderer.getSliderValue(1)));
        sphSolver.setViscosity(static_cast<double>(renderer.getSliderValue(2)));
        sphSolver.setSmoothingRadius(static_cast<double>(renderer.getSliderValue(3)));
        sphSolver.setGravity(Vector2D(0.0, static_cast<double>(renderer.getSliderValue(5))));

        const double currentParticleMass = sphSolver.getParticleMass();
        for (auto& particle : particles) {
            particle.mass = currentParticleMass;
        }

        // Handle Reset and preset buttons
        if (renderer.consumeButtonPress(0)) {
            for (std::size_t i = 0; i < particles.size(); ++i) {
                particles[i].position = initialPositions[i];
                particles[i].velocity = Vector2D(0.0, 0.0);
                particles[i].acceleration = Vector2D(0.0, 0.0);
            }
        }
        if (renderer.consumeButtonPress(1)) {
            makeDropOnBlock();
        }
        if (renderer.consumeButtonPress(2)) {
            makeParticleCloud();
        }
        if (renderer.consumeButtonPress(3)) {
            makeHighGravityImpact();
        }

        const double dt = static_cast<double>(renderer.getSliderValue(4));
        sphSolver.step(particles, dt);

        renderer.clear();
        renderer.drawParticles(particles);
        renderer.display();
    }

    return 0;
}