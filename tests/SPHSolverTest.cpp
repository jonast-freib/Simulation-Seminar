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
    particles.emplace_back(Vector2D(0.0, 4.0), mass, false);
    particles.back().velocity.y = -1.0;

    solver.step(particles, 0.001);

    const Particle &fluid = particles.back();
    assert(fluid.position.y >= 0.5 * solver.getSmoothingRadius());
    assert(fluid.velocity.y >= 0.0);
}
