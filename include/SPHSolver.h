#pragma once

#include <vector>
#include "Particle.h"
#include "NeighborSearch.h"

class SPHSolver {
public:

	static constexpr double defaultKernelSpacingRatio = 1.3;

	SPHSolver(double smoothingRadius = 30.0,
			  double restDensity = 300.0,
			  double stiffness = 1.0e6,
			  double viscosity = 200.0,
			  const Vector2D& gravity = Vector2D(0.0, 10.0));

	void setRestDensity(double rho0);
	void setStiffness(double k);
	void setViscosity(double mu);
	void setSmoothingRadius(double h);
	void setGravity(const Vector2D& g);
	void setKernelSpacingRatio(double ratio);
	void setBoundaryPressureMirroring(bool enabled);

	double getParticleMass() const;
	double getPressureWaveSpeed() const;
	double getRestDensity() const;

	double getParticleSpacing() const;

	std::size_t getNeighborCount(std::size_t particleIndex) const;

	void step(std::vector<Particle>& particles, double dt);

	int advance(std::vector<Particle>& particles, double duration);

	double getLastTimeStep() const;

	double pressureFromDensity(double density) const;

private:
	// Core SPH steps
	void computeDensityPressure(std::vector<Particle>& particles);
	void computeForces(std::vector<Particle>& particles);
	void integrate(std::vector<Particle>& particles, double dt);

	// Largest stable dt for the current particle state (requires computed forces).
	double computeStableTimeStep(const std::vector<Particle>& particles) const;

	// Smoothing kernels & helpers
	double cubicSplineKernel(double r2, double h) const;
	Vector2D cubicSplineGradient(const Vector2D& rij, double r, double h) const;
	double computeParticleMass() const;

private:
	double smoothingRadius;
	double restDensity;
	double stiffness;
	double viscosity;
	Vector2D gravity;

	double boundaryDamping;

	double kernelSpacingRatio;
	bool mirrorBoundaryPressure;
	double particleMass;
	double lastTimeStep;

	NeighborSearch neighborSearch;
};
