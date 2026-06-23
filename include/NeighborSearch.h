#pragma once

#include <vector>
#include "Particle.h"

class NeighborSearch {
public:
    NeighborSearch(double smoothingRadius = 30.0);

    // Set the search radius used for neighbor detection.
    void setRadius(double smoothingRadius);

    // Build a neighbor list for every particle in the scene.
    void build(const std::vector<Particle>& particles);

    // Get the indices of neighbors for a given particle.
    const std::vector<std::size_t>& getNeighbors(std::size_t particleIndex) const;

private:
    double smoothingRadius;
    double smoothingRadiusSquared;
    std::vector<std::vector<std::size_t>> neighborLists;
};
