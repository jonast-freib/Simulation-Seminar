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

		double densitySum = 0.0;

		for (std::size_t j : neighborSearch.getNeighbors(i))
		{
			const Particle &pj = particles[j];

			const Vector2D rij = pj.position - pi.position;

			const double r2 =
				rij.x * rij.x +
				rij.y * rij.y;

			densitySum +=
				pj.mass *
				cubicSplineKernel(r2, smoothingRadius);
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

	/*
	 * Phase 1:
	 *
	 * Compute non-pressure accelerations:
	 *   - gravity
	 *   - viscosity
	 *
	 * These accelerations are stored temporarily because pressure
	 * acceleration is calculated separately below.
	 */
	std::vector<Vector2D> accelNonPressure(
		particleCount,
		Vector2D(0.0, 0.0));

	for (std::size_t i = 0; i < particleCount; ++i)
	{
		const Particle &pi = particles[i];

		if (pi.fixed)
		{
			accelNonPressure[i] =
				Vector2D(0.0, 0.0);

			continue;
		}

		Vector2D viscositySum(0.0, 0.0);

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

			const Vector2D vRel =
				pi.velocity - pj.velocity;

			const double numerator =
				vRel.x * rij.x +
				vRel.y * rij.y;

			const double epsilon =
				0.01 * smoothingRadius * smoothingRadius;

			const Vector2D gradW =
				cubicSplineGradient(
					rij,
					r,
					smoothingRadius);

			const double particleDensity =
				std::max(pi.density, 1e-8);

			viscositySum +=
				(2.0 * viscosity * pj.mass / particleDensity) *
				((numerator / (r2 + epsilon))) *
				gradW;
		}

		const double particleDensity =
			std::max(pi.density, 1e-8);

		const Vector2D viscosityAcceleration =
			viscositySum / particleDensity;

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

		Vector2D pressureForceSum(0.0, 0.0);

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

			/*
			 * For fixed boundary particles we still use the
			 * symmetric pressure formulation. This requires
			 * boundary particles to have a meaningful density.
			 */
			const double pressureTerm =
				(pi.pressure + pj.pressure) /
				(2.0 *
				 particleDensity *
				 neighborDensity);

			const Vector2D gradW =
				cubicSplineGradient(
					rij,
					r,
					smoothingRadius);

			pressureForceSum +=
				-pj.mass *
				pressureTerm *
				gradW;
		}

		/*
		 * pressureForceSum already contains 1/rho_i,
		 * so there is no additional division by particle density.
		 */
		const Vector2D pressureAcceleration =
			pressureForceSum;
		if (!pi.fixed)
		{
			printf("Particle %zu: non Pressure acceleration = (%f, %f)\n", i, accelNonPressure[i].x, accelNonPressure[i].y);
		}
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

	for (auto &particle : particles)
	{
		if (particle.fixed)
		{
			particle.velocity =
				Vector2D(0.0, 0.0);

			particle.acceleration =
				Vector2D(0.0, 0.0);

			continue;
		}

		/*
		 * Semi-implicit Euler integration:
		 *
		 * v(t + dt) = v(t) + a(t) dt
		 * x(t + dt) = x(t) + v(t + dt) dt
		 */
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

/*
 * Standard 2D cubic B-spline kernel.
 *
 * Support:
 *
 *     0 <= r < 2h
 *
 * With:
 *
 *     q = r / h
 *
 * The normalized kernel is:
 *
 *     W(q) = alpha *
 *
 *       (1 - 1.5 q^2 + 0.75 q^3)     0 <= q < 1
 *
 *       0.25 (2-q)^3                  1 <= q < 2
 *
 * where:
 *
 *     alpha = 10 / (7 pi h^2)
 */
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
			   (2 - q) );
	}

	if (0.0 <= q && q < 1.0)
	{
		return alpha *
			   (2-q) *
			   (2-q) *
			   (2-q) -
			   4.0 *(1-q) *
			   (1-q) *
			   (1-q);
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
		10.0 / (7.0 * M_PI * h * h);

	double dWdq;

	if (q < 1.0)
	{
		dWdq =
			-3.0 * q +
			2.25 * q * q;
	}
	else
	{
		const double t = 2.0 - q;

		dWdq =
			-0.75 * t * t;
	}

	const double factor =
		alpha * dWdq / h / r;

	return rij * factor;
}

/*
 * Standard 2D SPH viscosity Laplacian.
 *
 * Support:
 *
 *     0 <= r < h
 *
 * One commonly used 2D formulation is:
 *
 *     laplacian(W) =
 *         20 / (3 pi h^5) * (h-r)
 */
double SPHSolver::viscosityLaplacian(
	double r,
	double h) const
{
	if (h <= 0.0 ||
		r <= 0.0 ||
		r >= h)
	{

		return 0.0;
	}

	return 20.0 /
		   (3.0 * M_PI * std::pow(h, 5)) *
		   (h - r);
}

/*
 * Compute particle mass from the desired rest density
 * and particle spacing.
 *
 * In 2D:
 *
 *     mass = density * area_per_particle
 *
 * We use:
 *
 *     spacing = 0.7 h
 *
 * so:
 *
 *     mass = rho0 * spacing^2
 */
double SPHSolver::computeParticleMass() const
{
	const double spacing =
		0.5 * smoothingRadius;

	return restDensity *
		   spacing *
		   spacing;
}