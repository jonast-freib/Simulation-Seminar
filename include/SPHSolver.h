#pragma once

#include <vector>
#include "Particle.h"
#include "NeighborSearch.h"

class SPHSolver {
public:
	SPHSolver(double smoothingRadius = 30.0,
			  double restDensity = 300.0,
			  double gasConstant = 2000.0,
			  double viscosity = 200.0,
			  const Vector2D& gravity = Vector2D(0.0, -10.0));

	void setBounds(double width, double height);

	void setRestDensity(double rho0);
	void setGasConstant(double k);
	void setViscosity(double mu);
	void setSmoothingRadius(double h);
	void setGravity(const Vector2D& g);

	double getParticleMass() const;
	double getSmoothingRadius() const;

	// Advance the simulation state by dt for the provided particles.
	void step(std::vector<Particle>& particles, double dt);

private:
	// Core SPH steps
	void computeDensityPressure(std::vector<Particle>& particles);
	void computeForces(std::vector<Particle>& particles);
	void integrate(std::vector<Particle>& particles, double dt);

	// Smoothing kernels & helpers
	double cubicSplineKernel(double r2, double h) const;
	Vector2D spikyGradient(const Vector2D& rij, double r, double h) const;
	double spikyMagnitude(double r, double h) const;
	double viscosityLaplacian(double r, double h) const;
	double computeParticleMass() const;

private:
	double smoothingRadius;
	double smoothingRadiusSquared;
	double restDensity;
	double gasConstant;
	double viscosity;
	Vector2D gravity;

	double boundaryWidth;
	double boundaryHeight;
	double boundaryDamping;

	double particleMass;

	NeighborSearch neighborSearch;
};

