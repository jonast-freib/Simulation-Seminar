#include <iostream>
#include <vector>
#include "Particle.h"
#include "Renderer.h"

int main() {
    Renderer renderer(1280, 720, "SPH Fluid Solver");

    // Create some particles
    std::vector<Particle> particles;
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            Vector2D pos(200.0 + i * 4, 100.0 + j * 4);
            particles.emplace_back(pos, 1.0);
        }
    }

    // Main loop
    while (renderer.isOpen()) {
        renderer.handleEvents();

        // Update particle positions (simple gravity for now)
        for (auto& particle : particles) {
            particle.velocity.y += 0.1; // gravity
            particle.position = particle.position + particle.velocity;

            // Bounce off walls
            if (particle.position.y > 700) {
                particle.position.y = 700;
                particle.velocity.y *= -0.8; // damping
            }
            if (particle.position.x < 0) {
                particle.position.x = 0;
                particle.velocity.x *= -0.8;
            }
            if (particle.position.x > 1280) {
                particle.position.x = 1280;
                particle.velocity.x *= -0.8;
            }
        }

        renderer.clear();
        renderer.drawParticles(particles);
        renderer.display();
    }

    return 0;
}