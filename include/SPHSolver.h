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

	// stiffness: stiffness coefficient k of the state equation
	// p_i = k ((rho_i / rho_0)^7 - 1). Larger k -> less compressible fluid,
	// smaller time step. Density errors stay ~1% when the pressure-wave speed
	// sqrt(7 k / rho_0) is >= 10 * the maximum flow speed.
	SPHSolver(double smoothingRadius = 30.0,
			  double restDensity = 300.0,
			  double stiffness = 1.0e6,
			  double viscosity = 200.0,
			  const Vector2D& gravity = Vector2D(0.0, 10.0));

	void setBounds(double width, double height);

	void setRestDensity(double rho0);
	void setStiffness(double k);
	void setViscosity(double mu);
	void setSmoothingRadius(double h);
	void setGravity(const Vector2D& g);
	void setKernelSpacingRatio(double ratio);
	// Boundary pressure in the pressure acceleration: mirrored from the fluid
	// particle (p_ib = p_i, course notes; default) or zero (p_ib = 0).
	void setBoundaryPressureMirroring(bool enabled);

	double getParticleMass() const;
	double getSmoothingRadius() const;
	double getStiffness() const;
	// Speed of pressure waves at rest density, sqrt(dp/drho) = sqrt(7 k / rho_0).
	// Used by the CFL condition and for Mach-number diagnostics.
	double getPressureWaveSpeed() const;
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

	// State equation (weakly compressible SPH): p = max(k ((rho / rho_0)^7 - 1), 0).
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
	double smoothingRadiusSquared;
	double restDensity;
	double stiffness;
	double viscosity;
	Vector2D gravity;

	double boundaryWidth;
	double boundaryHeight;
	double boundaryDamping;

	double kernelSpacingRatio;
	bool mirrorBoundaryPressure;
	double particleMass;
	double lastTimeStep;

	NeighborSearch neighborSearch;
};
