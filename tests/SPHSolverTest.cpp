// Keep assert() active in optimised (NDEBUG) builds.
#undef NDEBUG
#include <cassert>
#include <vector>

#include "SPHSolver.h"

int main()
{
    SPHSolver solver(10.0, 300.0, 500.0, 0.0, Vector2D(0.0, 0.0));
    solver.setBounds(100.0, 100.0);

    std::vector<Particle> particles;
    const double mass = solver.getParticleMass();
    for (int x = -4; x <= 4; ++x)
    {
        particles.emplace_back(Vector2D(x * 5.0, 0.0), mass, true);
    }
    // Start inside the wall safety radius (half a particle spacing) so the
    // projection must push the particle out and cancel its inward velocity.
    const double collisionRadius = 0.5 * solver.getParticleSpacing();
    particles.emplace_back(Vector2D(0.0, 0.5 * collisionRadius), mass, false);
    particles.back().velocity.y = -1.0;

    solver.step(particles, 0.001);

    const Particle &fluid = particles.back();
    assert(fluid.position.y >= collisionRadius - 1e-9);
    assert(fluid.velocity.y >= 0.0);
}
