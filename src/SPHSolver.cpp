#include "SPHSolver.h"
#include <cmath>

SPHSolver::SPHSolver(double smoothingRadius,
					 double restDensity,
					 double gasConstant,
					 double viscosity,
					 const Vector2D& gravity)
	: smoothingRadius(smoothingRadius)
	, smoothingRadiusSquared(smoothingRadius * smoothingRadius)
	, restDensity(restDensity)
	, gasConstant(gasConstant)
	, viscosity(viscosity)
	, gravity(gravity)
	, boundaryWidth(0.0)
	, boundaryHeight(0.0)
	, boundaryDamping(0.5)
	, particleMass(computeParticleMass())
	, neighborSearch(smoothingRadius) {
}

void SPHSolver::setBounds(double width, double height) {
	boundaryWidth = width;
	boundaryHeight = height;
}

void SPHSolver::setRestDensity(double rho0) {
	restDensity = rho0;
	particleMass = computeParticleMass();
}

void SPHSolver::setGasConstant(double k) {
	gasConstant = k;
}

void SPHSolver::setViscosity(double mu) {
	viscosity = mu;
}

void SPHSolver::setSmoothingRadius(double h) {
	smoothingRadius = h;
	smoothingRadiusSquared = h * h;
	neighborSearch.setRadius(2.0 * h);
	particleMass = computeParticleMass();
}

void SPHSolver::setGravity(const Vector2D& g) {
	gravity = g;
}

double SPHSolver::getParticleMass() const {
	return particleMass;
}

double SPHSolver::getSmoothingRadius() const {
	return smoothingRadius;
}

void SPHSolver::step(std::vector<Particle>& particles, double dt) {
	// Build the neighbor lists once per step so the solver can reuse them.
	neighborSearch.setRadius(2.0 * smoothingRadius);
	neighborSearch.build(particles);

	// Compute density and pressure for every particle.
	computeDensityPressure(particles);

	// Compute the fluid forces that arise from pressure and viscosity.
	computeForces(particles);

	// Integrate motion and apply simple boundary collisions.
	integrate(particles, dt);
}

void SPHSolver::computeDensityPressure(std::vector<Particle>& particles) {
	const std::size_t particleCount = particles.size();

	for (std::size_t i = 0; i < particleCount; ++i) {
		Particle& pi = particles[i];
		pi.density = 0.0;

		// Add self density contribution before neighbor contributions.
		pi.density = pi.mass * cubicSplineKernel(0.0, smoothingRadius);

		for (std::size_t j : neighborSearch.getNeighbors(i)) {
			const Particle& pj = particles[j];
			const Vector2D rij = pi.position - pj.position;
			const double r2 = rij.x * rij.x + rij.y * rij.y;
			// Sum contributions from neighboring particles using the cubic spline kernel.
			pi.density += pj.mass * cubicSplineKernel(r2, smoothingRadius);
		}

		// Use the equation of state to compute pressure from density.
		pi.pressure = gasConstant * (pi.density - restDensity);
		if (pi.pressure < 0.0) {
			// Prevent tensile instability by clamping negative pressure.
			pi.pressure = 0.0;
		}
	}
}

void SPHSolver::computeForces(std::vector<Particle>& particles) {
	const std::size_t particleCount = particles.size();

	// Phase 1: compute non-pressure accelerations (viscosity + gravity)
	std::vector<Vector2D> accelNonPressure(particleCount, Vector2D(0.0, 0.0));
	for (std::size_t i = 0; i < particleCount; ++i) {
		const Particle& pi = particles[i];
		Vector2D viscForceSum(0.0, 0.0);

		for (std::size_t j : neighborSearch.getNeighbors(i)) {
			const Particle& pj = particles[j];
			const Vector2D rij = pi.position - pj.position;
			const double r2 = rij.x * rij.x + rij.y * rij.y;
			const double r = std::sqrt(r2);

			if (r > 0.0) {
				const double neighborDensity = std::max(pj.density, 1e-8);
				const double laplacian = viscosityLaplacian(r, smoothingRadius);
				// accumulate viscosity force term (will divide by density later)
				viscForceSum += viscosity * pj.mass * ((pj.velocity - pi.velocity) / neighborDensity) * laplacian;
			}
		}

		const double particleDensity = std::max(pi.density, 1e-8);
		// viscosity acceleration: divide accumulated viscous force by particle density
		const Vector2D viscAccel = viscForceSum / particleDensity;
		// gravity is an acceleration term per-particle (not mass-scaled)
		accelNonPressure[i] = viscAccel + gravity;
	}

	// Phase 2: compute pressure accelerations and combine
	for (std::size_t i = 0; i < particleCount; ++i) {
		Particle& pi = particles[i];
		Vector2D pressureForceSum(0.0, 0.0);

		for (std::size_t j : neighborSearch.getNeighbors(i)) {
			const Particle& pj = particles[j];
			const Vector2D rij = pi.position - pj.position;
			const double r2 = rij.x * rij.x + rij.y * rij.y;
			const double r = std::sqrt(r2);

			if (r > 0.0) {
				const double neighborDensity = std::max(pj.density, 1e-8);
				const double pressureTerm = (pi.pressure + pj.pressure) / (2.0 * neighborDensity);
				const Vector2D gradW = spikyGradient(rij, r, smoothingRadius);
				// pressure force contribution (note the negative sign for pressure gradient)
				pressureForceSum += -pj.mass * pressureTerm * gradW;
			}
		}

		const double particleDensity = std::max(pi.density, 1e-8);
		const Vector2D pressureAccel = pressureForceSum / particleDensity;

		// final acceleration is non-pressure (visc + gravity) + pressure acceleration
		pi.acceleration = accelNonPressure[i] + pressureAccel;
	}
}

void SPHSolver::integrate(std::vector<Particle>& particles, double dt) {
	const double maxSpeed = 2000.0;

	for (auto& particle : particles) {
		if (particle.fixed) {
			particle.velocity = Vector2D(0.0, 0.0);
			particle.acceleration = Vector2D(0.0, 0.0);
			continue;
		}

		// Semi-explicit Euler integration: update velocity, then position.
		particle.velocity += particle.acceleration * dt;

		// Clamp the velocity to prevent runaway motion.
		const double speedSquared = particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y;
		if (speedSquared > maxSpeed * maxSpeed) {
			const double speed = std::sqrt(speedSquared);
			particle.velocity *= (maxSpeed / speed);
		}

		particle.position += particle.velocity * dt;

		// Simple collision handling with the window boundaries.
		if (particle.position.x < 0.0) {
			particle.position.x = 0.0;
			particle.velocity.x = -particle.velocity.x * boundaryDamping;
		} else if (particle.position.x > boundaryWidth) {
			particle.position.x = boundaryWidth;
			particle.velocity.x = -particle.velocity.x * boundaryDamping;
		}

		if (particle.position.y < 0.0) {
			particle.position.y = 0.0;
			particle.velocity.y = -particle.velocity.y * boundaryDamping;
		} else if (particle.position.y > boundaryHeight) {
			particle.position.y = boundaryHeight;
			particle.velocity.y = -particle.velocity.y * boundaryDamping;
		}
	}
}

// Cubic spline smoothing kernel used for density estimation in 2D.
// The current smoothingRadius is treated as the SPH smoothing length h.
// The kernel support is 2h, and q = r / h.
double SPHSolver::cubicSplineKernel(double r2, double h) const {
	const double r = std::sqrt(r2);
	const double R = 2.0 * h;
	if (r >= R) {
		return 0.0;
	}

	const double q = r / h;
	const double alpha = 5.0 / (14.0 * M_PI * h * h);

	if (q <= 1.0) {
		return alpha * ((2.0 - q) * (2.0 - q) * (2.0 - q) - 4.0 * (1.0 - q) * (1.0 - q) * (1.0 - q));
	}

	const double t = 2.0 - q;
	return alpha * t * t * t;
}

// Cubic spline gradient for pressure force in 2D.
Vector2D SPHSolver::spikyGradient(const Vector2D& rij, double r, double h) const {
	if (r >= 2.0 * h || r <= 0.0) {
		return Vector2D(0.0, 0.0);
	}

	const double q = r / h;
	const double sigma = 10.0 / (7.0 * M_PI * h * h);
	double dWdq;

	if (q <= 1.0) {
		dWdq = -3.0 * q + 2.25 * q * q;
	} else {
		const double t = 2.0 - q;
		dWdq = -0.75 * t * t;
	}

	const double factor = sigma * dWdq / h;
	return rij * (factor / r);
}

double SPHSolver::spikyMagnitude(double r, double h) const {
	if (r >= h || r <= 0.0) {
		return 0.0;
	}
	const double diff = h - r;
	return 30.0 / (M_PI * std::pow(h, 6)) * diff * diff;
}

// Viscosity kernel laplacian for viscosity force in 2D.
double SPHSolver::viscosityLaplacian(double r, double h) const {
	if (r >= 2.0 * h) {
		return 0.0;
	}

	if (r <= h) {
		const double q = r / h;
		// simple piecewise approximation for viscosity support 2h
		const double factor = 40.0 / (M_PI * std::pow(h, 6));
		return factor * (h - r);
	}

	const double factor = 40.0 / (M_PI * std::pow(h, 6));
	const double t = 2.0 - r / h;
	return factor * 0.25 * t * t;
}

double SPHSolver::computeParticleMass() const {
	const double spacing = 0.75 * smoothingRadius;
	return restDensity * spacing * spacing;
}

