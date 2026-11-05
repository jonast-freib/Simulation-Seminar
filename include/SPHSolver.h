#pragma once

#include <vector>
#include "Particle.h"
#include "NeighborSearch.h"

class SPHSolver {
public:
	// Default ratio of smoothing length to initial particle spacing, h / dx.
	// The cubic spline suffers from the pairing instability (particles clump
	// into rows) when h / dx is too large; 1.3 gives ~20 neighbours in 2D.
	static constexpr double defaultKernelSpacingRatio = 1.3;

	// speedOfSound: numerical speed of sound c of the weakly compressible fluid.
	// Choose c >= 10 * expected maximum flow speed to keep density errors ~1%.
	SPHSolver(double smoothingRadius = 30.0,
			  double restDensity = 300.0,
			  double speedOfSound = 500.0,
			  double viscosity = 200.0,
			  const Vector2D& gravity = Vector2D(0.0, 10.0));

	void setBounds(double width, double height);

	void setRestDensity(double rho0);
	void setSpeedOfSound(double c);
	void setViscosity(double mu);
	void setSmoothingRadius(double h);
	void setGravity(const Vector2D& g);
	void setKernelSpacingRatio(double ratio);

	double getParticleMass() const;
	double getSmoothingRadius() const;
	double getSpeedOfSound() const;
	double getRestDensity() const;
	double getKernelSpacingRatio() const;

	// Initial particle spacing dx = h / (h/dx) that the particle mass is calibrated for.
	double getParticleSpacing() const;

	// Number of neighbours (within the 2h kernel support) of a particle in the last step.
	std::size_t getNeighborCount(std::size_t particleIndex) const;

	// Advance the simulation state by exactly one step of size dt.
	void step(std::vector<Particle>& particles, double dt);

	// Advance the simulation by `duration`, splitting it into as many
	// substeps as the stability limits (CFL, force, viscosity) require.
	// Returns the number of substeps taken.
	int advance(std::vector<Particle>& particles, double duration);

	// Size of the most recent substep and the stability limit it was taken from.
	double getLastTimeStep() const;

	// Tait equation of state: p = B * ((rho / rho0)^7 - 1), B = rho0 * c^2 / 7.
	double pressureFromDensity(double density) const;

private:
	// Core SPH steps
	void computeDensityPressure(std::vector<Particle>& particles);
	void computeWallPressure(std::vector<Particle>& particles);
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
	double smoothingRadiusSquared;
	double restDensity;
	double speedOfSound;
	double viscosity;
	Vector2D gravity;

	double boundaryWidth;
	double boundaryHeight;
	double boundaryDamping;

	double kernelSpacingRatio;
	double particleMass;
	double lastTimeStep;

	NeighborSearch neighborSearch;
};
