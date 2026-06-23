#pragma once

#include <vector>
#include "Particle.h"
#include "NeighborSearch.h"

class SPHSolver {
public:
    SPHSolver(double smoothingRadius = 16.0,
              double restDensity = 300.0,
              double gasConstant = 2000.0,
              double viscosity = 200.0,
              const Vector2D& gravity = Vector2D(0.0, 10.0));

    // Step the simulation forward by dt seconds for all particles.
    void step(std::vector<Particle>& particles, double dt);

    // Update the boundary size used for collision handling.
    void setBounds(double width, double height);
    // Runtime parameter setters
    void setRestDensity(double rho0);
    void setGasConstant(double k);
    void setViscosity(double mu);
    void setSmoothingRadius(double h);
    // runtime gravity setter
    void setGravity(const Vector2D& g);
    // Particle mass is derived from smoothing radius and rest density.
    double getParticleMass() const;
    double getSmoothingRadius() const;

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

    void computeDensityPressure(std::vector<Particle>& particles);
    void computeForces(std::vector<Particle>& particles);
    void integrate(std::vector<Particle>& particles, double dt);

    double cubicSplineKernel(double r2, double h) const;
    Vector2D spikyGradient(const Vector2D& rij, double r, double h) const;
    double spikyMagnitude(double r, double h) const;
    double viscosityLaplacian(double r, double h) const;
    double computeParticleMass() const;
};
