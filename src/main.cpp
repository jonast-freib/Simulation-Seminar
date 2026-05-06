#include <iostream>
#include <vector>
#include "Particle.h"
#include "Renderer.h"

int main() {
    const int particlesX = 300;
    const int particlesY = 200;
    const float particleRadius = 4.f;
    const float spacing = particleRadius * 2.4f;
    const float margin = 20.f;

    const unsigned int windowWidth = 1800;
    const unsigned int windowHeight = 1000;

    Renderer renderer(windowWidth, windowHeight, "SPH Fluid Solver", particleRadius);

    std::vector<Particle> particles;
    particles.reserve(particlesX * particlesY);

    for (int row = 0; row < particlesY; ++row) {
        for (int col = 0; col < particlesX; ++col) {
            const float x = margin + col * spacing;
            const float y = margin + row * spacing;
            particles.emplace_back(Vector2D(x, y), 1.0);
        }
    }

    // Main loop
    while (renderer.isOpen()) {
        renderer.handleEvents();

        for (auto& particle : particles) {
            particle.velocity.y += 0.1f;
            particle.velocity.x += (rand() % 3 - 1) * 0.05f;
            particle.position = particle.position + particle.velocity;

            if (particle.position.y > static_cast<float>(windowHeight)) {
                particle.position.y = static_cast<float>(windowHeight);
                particle.velocity.y *= -0.8f;
            }
            if (particle.position.x < 0.f) {
                particle.position.x = 0.f;
                particle.velocity.x *= -0.8f;
            }
            if (particle.position.x > static_cast<float>(windowWidth)) {
                particle.position.x = static_cast<float>(windowWidth);
                particle.velocity.x *= -0.8f;
            }
        }

        renderer.clear();
        renderer.drawParticles(particles);
        renderer.display();
    }

    return 0;
}