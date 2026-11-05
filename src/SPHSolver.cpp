#include "SPHSolver.h"
#include <algorithm>
#include <cmath>

SPHSolver::SPHSolver(double smoothingRadius,
                     double restDensity,
                     double speedOfSound,
                     double viscosity,
                     const Vector2D &gravity)
    : smoothingRadius(smoothingRadius), smoothingRadiusSquared(smoothingRadius * smoothingRadius), restDensity(restDensity), speedOfSound(speedOfSound), viscosity(viscosity), gravity(gravity), boundaryWidth(0.0), boundaryHeight(0.0), boundaryDamping(0.5), kernelSpacingRatio(defaultKernelSpacingRatio), particleMass(computeParticleMass()), lastTimeStep(0.0), neighborSearch(2.0 * smoothingRadius)
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

void SPHSolver::setSpeedOfSound(double c)
{
    speedOfSound = c;
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

double SPHSolver::getSpeedOfSound() const
{
    return speedOfSound;
}

double SPHSolver::getRestDensity() const
{
    return restDensity;
}

void SPHSolver::setKernelSpacingRatio(double ratio)
{
    kernelSpacingRatio = ratio;
    particleMass = computeParticleMass();
}

double SPHSolver::getKernelSpacingRatio() const
{
    return kernelSpacingRatio;
}

double SPHSolver::getParticleSpacing() const
{
    return smoothingRadius / kernelSpacingRatio;
}

std::size_t SPHSolver::getNeighborCount(std::size_t particleIndex) const
{
    return neighborSearch.getNeighbors(particleIndex).size();
}

double SPHSolver::getLastTimeStep() const
{
    return lastTimeStep;
}

double SPHSolver::pressureFromDensity(double density) const
{
    // Tait equation (Monaghan 1994, Becker & Teschner 2007) with gamma = 7.
    // The stiffness B makes the speed of sound at rest density equal to c.
    constexpr double gamma = 7.0;
    const double B = restDensity * speedOfSound * speedOfSound / gamma;
    const double ratio = density / restDensity;
    const double r2 = ratio * ratio;
    const double r7 = r2 * r2 * r2 * ratio;
    // Negative pressures are clamped to avoid tensile clumping at the free surface.
    return std::max(0.0, B * (r7 - 1.0));
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
    lastTimeStep = dt;
    integrate(particles, dt);
}

int SPHSolver::advance(std::vector<Particle> &particles, double duration)
{
    // Guard against a stall if the state has blown up (dt -> 0).
    constexpr int maxSubsteps = 2000;
    const double minTimeStep = duration / maxSubsteps;

    double elapsed = 0.0;
    int substeps = 0;
    while (elapsed < duration && substeps < maxSubsteps)
    {
        neighborSearch.setRadius(2.0 * smoothingRadius);
        neighborSearch.build(particles);
        computeDensityPressure(particles);
        computeForces(particles);

        double dt = std::max(computeStableTimeStep(particles), minTimeStep);
        const double remaining = duration - elapsed;
        if (dt >= remaining)
        {
            dt = remaining;
        }
        else if (dt > 0.5 * remaining)
        {
            // Split the rest evenly instead of leaving a tiny final step.
            dt = 0.5 * remaining;
        }

        lastTimeStep = dt;
        integrate(particles, dt);
        elapsed += dt;
        ++substeps;
    }
    return substeps;
}

double SPHSolver::computeStableTimeStep(const std::vector<Particle> &particles) const
{
    double maxSpeedSquared = 0.0;
    double maxAccelSquared = 0.0;
    for (const Particle &p : particles)
    {
        if (p.fixed)
        {
            continue;
        }
        maxSpeedSquared = std::max(maxSpeedSquared,
                                   p.velocity.x * p.velocity.x + p.velocity.y * p.velocity.y);
        maxAccelSquared = std::max(maxAccelSquared,
                                   p.acceleration.x * p.acceleration.x + p.acceleration.y * p.acceleration.y);
    }

    const double h = smoothingRadius;

    // Acoustic CFL condition: information travels at c + |v|.
    double dt = 0.4 * h / (speedOfSound + std::sqrt(maxSpeedSquared));

    // Force condition (Monaghan 1992).
    const double maxAccel = std::sqrt(maxAccelSquared);
    if (maxAccel > 0.0)
    {
        dt = std::min(dt, 0.25 * std::sqrt(h / maxAccel));
    }

    // Viscous diffusion condition.
    if (viscosity > 0.0)
    {
        dt = std::min(dt, 0.125 * h * h / viscosity);
    }

    return dt;
}

void SPHSolver::computeDensityPressure(std::vector<Particle> &particles)
{
    const std::size_t particleCount = particles.size();

    for (std::size_t i = 0; i < particleCount; ++i)
    {
        Particle &pi = particles[i];

        if (pi.fixed)
        {
            // Wall pressure is extrapolated from the fluid below.
            continue;
        }

        // Self-contribution: neighbor lists exclude the particle itself,
        // but the SPH density sum runs over all j including j == i.
        double densitySum = pi.mass * cubicSplineKernel(0.0, smoothingRadius);

        for (std::size_t j : neighborSearch.getNeighbors(i))
        {
            const Particle &pj = particles[j];

            const Vector2D rij = pj.position - pi.position;

            const double r2 =
                rij.x * rij.x +
                rij.y * rij.y;

            const double kernelValue =
                cubicSplineKernel(r2, smoothingRadius);

            // Fixed particles contribute like fluid at rest density. Walls must
            // be thick enough (3 layers at spacing h/2) to fill the 2h support.
            densitySum += pj.mass * kernelValue;
        }
        pi.density = densitySum;
        pi.pressure = pressureFromDensity(pi.density);
    }

    computeWallPressure(particles);
}

void SPHSolver::computeWallPressure(std::vector<Particle> &particles)
{
    // Generalized wall boundary condition (Adami, Hu & Adams 2012):
    // each wall particle gets the Shepard-interpolated pressure of the
    // surrounding fluid plus the hydrostatic correction for its offset,
    //   p_b = (sum_f p_f W_bf + g . sum_f rho_f (x_b - x_f) W_bf) / sum_f W_bf.
    // The wall then pushes back with the same symmetric pressure term as
    // the fluid, so the first fluid layer needs no excess pressure.
    constexpr double gamma = 7.0;
    const double B = restDensity * speedOfSound * speedOfSound / gamma;

    for (std::size_t b = 0; b < particles.size(); ++b)
    {
        Particle &wall = particles[b];
        if (!wall.fixed)
        {
            continue;
        }

        double weightSum = 0.0;
        double pressureSum = 0.0;
        double hydrostaticSum = 0.0;
        for (std::size_t f : neighborSearch.getNeighbors(b))
        {
            const Particle &fluid = particles[f];
            if (fluid.fixed)
            {
                continue;
            }
            const Vector2D rbf = wall.position - fluid.position;
            const double w = cubicSplineKernel(rbf.x * rbf.x + rbf.y * rbf.y, smoothingRadius);
            weightSum += w;
            pressureSum += fluid.pressure * w;
            hydrostaticSum += fluid.density * (gravity.x * rbf.x + gravity.y * rbf.y) * w;
        }

        if (weightSum > 0.0)
        {
            wall.pressure = std::max(0.0, (pressureSum + hydrostaticSum) / weightSum);
        }
        else
        {
            wall.pressure = 0.0;
        }
        // Wall density from the inverse equation of state.
        wall.density = restDensity * std::pow(wall.pressure / B + 1.0, 1.0 / gamma);
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
                // Fluid-Boundary interaction: same symmetric form as
                // fluid-fluid, using the extrapolated wall pressure.
                pressureTermBoundary =
                    (pi.pressure / (particleDensity * particleDensity) +
                     pj.pressure / (neighborDensity * neighborDensity));
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
            // Safety net only: resting fluid sits a full spacing away from
            // the wall particles, so this projection does not engage unless
            // a particle is about to penetrate the wall.
            const double collisionRadius = 0.5 * getParticleSpacing();
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
        getParticleSpacing();

    return restDensity *
           spacing *
           spacing;
}