#include "SPHSolver.h"
#include <algorithm>
#include <cmath>

SPHSolver::SPHSolver(double smoothingRadius,
                     double restDensity,
                     double stiffness,
                     double viscosity,
                     const Vector2D &gravity)
    : smoothingRadius(smoothingRadius), restDensity(restDensity), stiffness(stiffness), viscosity(viscosity), gravity(gravity), boundaryDamping(0.5), kernelSpacingRatio(defaultKernelSpacingRatio), mirrorBoundaryPressure(true), particleMass(computeParticleMass()), lastTimeStep(0.0), neighborSearch(2.0 * smoothingRadius)
{
}

void SPHSolver::setRestDensity(double rho0)
{
    restDensity = rho0;
    particleMass = computeParticleMass();
}

void SPHSolver::setStiffness(double k)
{
    stiffness = k;
}

void SPHSolver::setViscosity(double mu)
{
    viscosity = mu;
}

void SPHSolver::setSmoothingRadius(double h)
{
    smoothingRadius = h;

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

double SPHSolver::getPressureWaveSpeed() const
{
    return std::sqrt(7.0 * stiffness / restDensity);
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

void SPHSolver::setBoundaryPressureMirroring(bool enabled)
{
    mirrorBoundaryPressure = enabled;
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
    const double ratio = density / restDensity;
    const double r2 = ratio * ratio;
    const double r7 = r2 * r2 * r2 * ratio;
    return std::max(0.0, stiffness * (r7 - 1.0));
}

void SPHSolver::step(std::vector<Particle> &particles, double dt)
{
    neighborSearch.build(particles);

    computeDensityPressure(particles);

    computeForces(particles);

    lastTimeStep = dt;
    integrate(particles, dt);
}

int SPHSolver::advance(std::vector<Particle> &particles, double duration)
{
    constexpr int maxSubsteps = 2000;
    const double minTimeStep = duration / maxSubsteps;

    double elapsed = 0.0;
    int substeps = 0;
    while (elapsed < duration && substeps < maxSubsteps)
    {
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

    double dt = 0.4 * h / (getPressureWaveSpeed() + std::sqrt(maxSpeedSquared));

    const double maxAccel = std::sqrt(maxAccelSquared);
    if (maxAccel > 0.0)
    {
        dt = std::min(dt, 0.25 * std::sqrt(h / maxAccel));
    }

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
            pi.density = restDensity;
            pi.pressure = 0.0;
            continue;
        }

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

            densitySum += pj.mass * kernelValue;
        }
        pi.density = densitySum;
        pi.pressure = pressureFromDensity(pi.density);
    }
}

void SPHSolver::computeForces(std::vector<Particle> &particles)
{
    const double epsilon = 0.01 * smoothingRadius * smoothingRadius;
    const double boundaryFactor = mirrorBoundaryPressure ? 2.0 : 1.0;

    for (std::size_t i = 0; i < particles.size(); ++i)
    {
        Particle &pi = particles[i];

        if (pi.fixed)
        {
            pi.acceleration = Vector2D(0.0, 0.0);
            continue;
        }

        const double particleDensity = std::max(pi.density, 1e-8);
        const double pressureTermI = pi.pressure / (particleDensity * particleDensity);

        Vector2D viscositySum(0.0, 0.0);
        Vector2D pressureSumFluid(0.0, 0.0);
        Vector2D pressureSumBoundary(0.0, 0.0);

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

            const Vector2D gradW = cubicSplineGradient(rij, r, smoothingRadius);
            const double neighborDensity = std::max(pj.density, 1e-8);

            const Vector2D vRel = pi.velocity - pj.velocity;
            const double vijDotXij = vRel.x * rij.x + vRel.y * rij.y;
            viscositySum += (pj.mass / neighborDensity) * (vijDotXij / (r2 + epsilon)) * gradW;

            if (pj.fixed)
            {
                pressureSumBoundary += -pj.mass * (boundaryFactor * pressureTermI) * gradW;
            }
            else
            {
                pressureSumFluid += -pj.mass * (pressureTermI + pj.pressure / (neighborDensity * neighborDensity)) * gradW;
            }
        }

        const Vector2D nonPressureAcceleration = 2.0 * viscosity * viscositySum + gravity;
        pi.acceleration = nonPressureAcceleration + (pressureSumFluid + pressureSumBoundary);
    }
}

void SPHSolver::integrate(std::vector<Particle> &particles, double dt)
{
    const double collisionRadius = 0.5 * getParticleSpacing();
    const double collisionRadiusSquared = collisionRadius * collisionRadius;

    for (std::size_t particleIndex = 0; particleIndex < particles.size(); ++particleIndex)
    {
        Particle &particle = particles[particleIndex];

        if (particle.fixed)
        {
            continue;
        }

        particle.velocity += particle.acceleration * dt;
        particle.position += particle.velocity * dt;

        std::size_t closestFixed = particles.size();
        double closestDistanceSquared = collisionRadiusSquared;

        for (std::size_t neighborIndex : neighborSearch.getNeighbors(particleIndex))
        {
            if (!particles[neighborIndex].fixed)
            {
                continue;
            }

            const Vector2D offset = particle.position - particles[neighborIndex].position;
            const double distanceSquared = offset.x * offset.x + offset.y * offset.y;

            if (distanceSquared < closestDistanceSquared)
            {
                closestDistanceSquared = distanceSquared;
                closestFixed = neighborIndex;
            }
        }

        if (closestFixed < particles.size() && closestDistanceSquared > 0.0)
        {
            const Vector2D offset = particle.position - particles[closestFixed].position;
            const Vector2D normal = offset / std::sqrt(closestDistanceSquared);

            particle.position = particles[closestFixed].position + normal * collisionRadius;

            const double velocityAlongNormal =
                particle.velocity.x * normal.x + particle.velocity.y * normal.y;
            if (velocityAlongNormal < 0.0)
            {
                particle.velocity -= normal * ((1.0 + boundaryDamping) * velocityAlongNormal);
            }
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

    return alpha *
           ((2.0 - q) *
                (2.0 - q) *
                (2.0 - q) -
            4.0 * (1.0 - q) *
                (1.0 - q) *
                (1.0 - q));
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