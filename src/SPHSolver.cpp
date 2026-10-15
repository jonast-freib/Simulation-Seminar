#include "SPHSolver.h"
#include <algorithm>
#include <cmath>

SPHSolver::SPHSolver(double smoothingRadius,
                     double restDensity,
                     double gasConstant,
                     double viscosity,
                     const Vector2D &gravity)
    : smoothingRadius(smoothingRadius), smoothingRadiusSquared(smoothingRadius * smoothingRadius), restDensity(restDensity), gasConstant(gasConstant), viscosity(viscosity), gravity(gravity), boundaryWidth(0.0), boundaryHeight(0.0), boundaryDamping(0.5), particleMass(computeParticleMass()), neighborSearch(smoothingRadius)
{
}

void SPHSolver::setBounds(double width, double height)
{
    boundaryWidth = width;
    boundaryHeight = height;
}

void SPHSolver::setRestDensity(double rho0)
{
    restDensity = rho0;
    particleMass = computeParticleMass();
}

void SPHSolver::setGasConstant(double k)
{
    gasConstant = k;
}

void SPHSolver::setViscosity(double mu)
{
    viscosity = mu;
}

void SPHSolver::setSmoothingRadius(double h)
{
    smoothingRadius = h;
    smoothingRadiusSquared = h * h;

    // Density kernel has support 2h.
    neighborSearch.setRadius(2.0 * h);

    particleMass = computeParticleMass();
}

void SPHSolver::setGravity(const Vector2D &g)
{
    gravity = g;
}

double SPHSolver::getParticleMass() const
{
    return particleMass;
}

double SPHSolver::getSmoothingRadius() const
{
    return smoothingRadius;
}

void SPHSolver::step(std::vector<Particle> &particles, double dt)
{
    // Build the neighbor list once per step.
    // The cubic spline density kernel has support 2h.
    neighborSearch.setRadius(2.0 * smoothingRadius);
    neighborSearch.build(particles);

    // Compute density and pressure for every particle.
    computeDensityPressure(particles);

    // Compute pressure, viscosity and gravity accelerations.
    computeForces(particles);

    // Integrate motion and handle boundaries.
    integrate(particles, dt);
}

void SPHSolver::computeDensityPressure(std::vector<Particle> &particles)
{
    const std::size_t particleCount = particles.size();

    for (std::size_t i = 0; i < particleCount; ++i)
    {
        Particle &pi = particles[i];

        if (pi.fixed)
        {
            pi.density = restDensity;
            pi.pressure = 0.0;
            continue;
        }

        double densitySum = 0.0;

        for (std::size_t j : neighborSearch.getNeighbors(i))
        {
            const Particle &pj = particles[j];

            const Vector2D rij = pj.position - pi.position;

            const double r2 =
                rij.x * rij.x +
                rij.y * rij.y;

            const double kernelValue =
                cubicSplineKernel(r2, smoothingRadius);

            densitySum += pj.mass * kernelValue;

            // Fixed particles represent the mirrored solid half-space.
            if (pj.fixed)
            {
                densitySum += pj.mass * kernelValue;
            }
        }
        pi.density = densitySum;
        pi.pressure =
            gasConstant *
            ((pi.density / restDensity) - 1);

        if (pi.pressure < 0.0)
        {
            pi.pressure = 0.0;
        }
    }
}

void SPHSolver::computeForces(std::vector<Particle> &particles)
{
    const std::size_t particleCount = particles.size();

    std::vector<Vector2D> accelNonPressure(
        particleCount,
        Vector2D(0.0, 0.0));

    for (std::size_t i = 0; i < particleCount; ++i)
    {
        const Particle &pi = particles[i];

        if (pi.fixed)
        {
            accelNonPressure[i] = Vector2D(0.0, 0.0);

            continue;
        }

        Vector2D viscositySum(0.0, 0.0);

        for (std::size_t j : neighborSearch.getNeighbors(i))
        {
            const Particle &pj = particles[j];

            const Vector2D rij = pi.position - pj.position;

            const double r2 = rij.x * rij.x + rij.y * rij.y;

            const double r = std::sqrt(r2);

            if (r <= 0.0)
            {
                continue;
            }

            const Vector2D vRel = pi.velocity - pj.velocity;

            const double vijDotXij = vRel.x * rij.x + vRel.y * rij.y;

            const double epsilon = 0.01 * smoothingRadius * smoothingRadius;

            const Vector2D gradW = cubicSplineGradient(rij, r, smoothingRadius);

            const double neighborDensity = std::max(pj.density, 1e-8);

            viscositySum +=
                (pj.mass / neighborDensity) *
                (vijDotXij / (r2 + epsilon)) *
                gradW;
        }
        const Vector2D viscosityAcceleration =
            2.0 * viscosity * viscositySum;

        accelNonPressure[i] =
            viscosityAcceleration +
            gravity;
    }

    for (std::size_t i = 0; i < particleCount; ++i)
    {
        Particle &pi = particles[i];

        if (pi.fixed)
        {
            pi.acceleration =
                Vector2D(0.0, 0.0);

            continue;
        }

        Vector2D pressureForceSumFluid(0.0, 0.0);
        Vector2D pressureForceSumBoundary(0.0, 0.0);

        for (std::size_t j : neighborSearch.getNeighbors(i))
        {
            const Particle &pj = particles[j];

            const Vector2D rij =
                pi.position - pj.position;

            const double r2 =
                rij.x * rij.x +
                rij.y * rij.y;

            const double r = std::sqrt(r2);

            if (r <= 0.0)
            {
                continue;
            }

            const double particleDensity =
                std::max(pi.density, 1e-8);

            const double neighborDensity =
                std::max(pj.density, 1e-8);

            double pressureTermFluid = 0.0;
            double pressureTermBoundary = 0.0;

            if (!pj.fixed)
            {
                // Fluid-Fluid interaction
                pressureTermFluid =
                    (pi.pressure / (particleDensity * particleDensity) +
                     pj.pressure / (neighborDensity * neighborDensity));
            }
            if (pj.fixed)
            {
                // Fluid-Boundary interaction
                pressureTermBoundary =
                    2.0 * pi.pressure /
                    (particleDensity * particleDensity);
            }

            const Vector2D gradW =
                cubicSplineGradient(
                    rij,
                    r,
                    smoothingRadius);

            if (!pj.fixed)
            {
                pressureForceSumFluid +=
                    -pj.mass *
                    pressureTermFluid *
                    gradW;
            }
            if (pj.fixed)
            {
                pressureForceSumBoundary +=
                    -pj.mass *
                    pressureTermBoundary *
                    gradW;
            }
        }

        const Vector2D pressureAcceleration =
            pressureForceSumFluid + pressureForceSumBoundary;
        pi.acceleration =
            accelNonPressure[i] +
            pressureAcceleration;
    }
}

void SPHSolver::integrate(
    std::vector<Particle> &particles,
    double dt)
{
    const double maxSpeed = 2000.0;

    for (std::size_t particleIndex = 0;
         particleIndex < particles.size();
         ++particleIndex)
    {
        Particle &particle = particles[particleIndex];

        particle.velocity +=
            particle.acceleration * dt;

        /*
         * Safety clamp to prevent numerical explosions from
         * completely destroying the simulation.
         */
        const double speedSquared =
            particle.velocity.x *
                particle.velocity.x +
            particle.velocity.y *
                particle.velocity.y;

        if (speedSquared >
            maxSpeed * maxSpeed)
        {

            const double speed =
                std::sqrt(speedSquared);

            particle.velocity *=
                (maxSpeed / speed);
        }

        particle.position +=
            particle.velocity * dt;

        if (!particle.fixed)
        {
            const double collisionRadius = 0.5 * smoothingRadius;
            const double collisionRadiusSquared =
                collisionRadius * collisionRadius;
            std::size_t closestFixed = particles.size();
            double closestDistanceSquared = collisionRadiusSquared;

            for (std::size_t neighborIndex :
                 neighborSearch.getNeighbors(particleIndex))
            {
                if (!particles[neighborIndex].fixed)
                {
                    continue;
                }

                const Vector2D offset =
                    particle.position -
                    particles[neighborIndex].position;
                const double distanceSquared =
                    offset.x * offset.x +
                    offset.y * offset.y;

                if (distanceSquared < closestDistanceSquared)
                {
                    closestDistanceSquared = distanceSquared;
                    closestFixed = neighborIndex;
                }
            }

            if (closestFixed < particles.size() &&
                closestDistanceSquared > 0.0)
            {
                const Vector2D offset =
                    particle.position - particles[closestFixed].position;
                const double distance =
                    std::sqrt(closestDistanceSquared);
                const Vector2D normal = offset / distance;

                particle.position =
                    particles[closestFixed].position +
                    normal * collisionRadius;

                const double velocityAlongNormal =
                    particle.velocity.x * normal.x +
                    particle.velocity.y * normal.y;
                if (velocityAlongNormal < 0.0)
                {
                    particle.velocity -=
                        normal *
                        ((1.0 + boundaryDamping) *
                         velocityAlongNormal);
                }
            }
        }

        /*
         * Boundary collision handling.
         */
        if (particle.position.x < 0.0)
        {
            particle.position.x = 0.0;

            particle.velocity.x =
                -particle.velocity.x *
                boundaryDamping;
        }
        else if (particle.position.x >
                 boundaryWidth)
        {

            particle.position.x =
                boundaryWidth;

            particle.velocity.x =
                -particle.velocity.x *
                boundaryDamping;
        }

        if (particle.position.y < 0.0)
        {
            particle.position.y = 0.0;

            particle.velocity.y =
                -particle.velocity.y *
                boundaryDamping;
        }
        else if (particle.position.y >
                 boundaryHeight)
        {

            particle.position.y =
                boundaryHeight;

            particle.velocity.y =
                -particle.velocity.y *
                boundaryDamping;
        }
    }
}

double SPHSolver::cubicSplineKernel(
    double r2,
    double h) const
{
    if (h <= 0.0)
    {
        return 0.0;
    }

    const double r =
        std::sqrt(r2);
    const double q = r / h;

    if (q >= 2.0)
    {
        return 0.0;
    }

    const double alpha =
        5.0 /
        (14.0 * M_PI * h * h);

    if (1.0 <= q && q < 2.0)
    {
        return alpha *
               ((2 - q) *
                (2 - q) *
                (2 - q));
    }

    if (0.0 <= q && q < 1.0)
    {
        return alpha *
               ((2.0 - q) *
                    (2.0 - q) *
                    (2.0 - q) -
                4.0 * (1.0 - q) *
                    (1.0 - q) *
                    (1.0 - q));
    }

    const double t =
        2.0 - q;

    return alpha *
           0.25 *
           t * t * t;
}

Vector2D SPHSolver::cubicSplineGradient(
    const Vector2D &rij,
    double r,
    double h) const
{
    if (h <= 0.0 || r <= 0.0 || r >= 2.0 * h)
    {
        return Vector2D(0.0, 0.0);
    }

    const double q = r / h;
    const double alpha =
        5.0 /
        (14.0 * M_PI * h * h);

    double dWdq = 0.0;

    if (q < 1.0)
    {
        dWdq = -3.0 * (2.0 - q) * (2.0 - q) + 12.0 * (1.0 - q) * (1.0 - q);
    }
    else if (q < 2.0)
    {
        const double t = 2.0 - q;
        dWdq = -3.0 * t * t;
    }

    const double factor =
        alpha * dWdq / h / r;

    return rij * factor;
}

double SPHSolver::computeParticleMass() const
{
    const double spacing =
        0.5 * smoothingRadius;

    return restDensity *
           spacing *
           spacing;
}