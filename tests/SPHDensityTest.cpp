// Keep assert() active in optimised (NDEBUG) builds.
#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "SPHSolver.h"

// An isolated particle must still carry its own mass in the density sum:
// rho_i = m_i * W(0, h), with W(0, h) = 4 * 5 / (14 * pi * h^2) for the 2D cubic spline.
static void testIsolatedParticleDensity()
{
    const double h = 10.0;
    SPHSolver solver(h, 300.0, 500.0, 0.0, Vector2D(0.0, 0.0));
    solver.setBounds(1000.0, 1000.0);

    std::vector<Particle> particles;
    particles.emplace_back(Vector2D(500.0, 500.0), solver.getParticleMass(), false);

    solver.step(particles, 0.001);

    const double expected = solver.getParticleMass() * 4.0 * 5.0 / (14.0 * M_PI * h * h);
    assert(std::abs(particles[0].density - expected) < 1e-9 * expected);
}

// The particle mass is calibrated as rho0 * dx^2, so a particle in the
// interior of a regular lattice with spacing dx should sit close to rho0.
static void testLatticeInteriorReachesRestDensity()
{
    const double h = 10.0;
    const double rho0 = 300.0;
    SPHSolver solver(h, rho0, 500.0, 0.0, Vector2D(0.0, 0.0));
    solver.setBounds(1000.0, 1000.0);

    const double spacing = SPHSolver::particleSpacingFor(h);
    const int n = 21;
    std::vector<Particle> particles;
    for (int y = 0; y < n; ++y)
    {
        for (int x = 0; x < n; ++x)
        {
            particles.emplace_back(Vector2D(400.0 + x * spacing, 400.0 + y * spacing),
                                   solver.getParticleMass(), false);
        }
    }

    solver.step(particles, 0.0);

    const double centerDensity = particles[(n / 2) * n + n / 2].density;
    std::printf("lattice centre density = %.3f (rho0 = %.1f)\n", centerDensity, rho0);
    assert(std::abs(centerDensity - rho0) < 0.02 * rho0);
}

int main()
{
    testIsolatedParticleDensity();
    testLatticeInteriorReachesRestDensity();
    return 0;
}
