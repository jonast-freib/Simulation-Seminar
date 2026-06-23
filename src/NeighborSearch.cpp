#include "NeighborSearch.h"
#include <cmath>

NeighborSearch::NeighborSearch(double smoothingRadius)
    : smoothingRadius(smoothingRadius)
    , smoothingRadiusSquared(smoothingRadius * smoothingRadius) {
}

void NeighborSearch::setRadius(double newRadius) {
    smoothingRadius = newRadius;
    smoothingRadiusSquared = newRadius * newRadius;
}

void NeighborSearch::build(const std::vector<Particle>& particles) {
    const std::size_t particleCount = particles.size();
    neighborLists.assign(particleCount, {});

    for (std::size_t i = 0; i < particleCount; ++i) {
        for (std::size_t j = i + 1; j < particleCount; ++j) {
            const Vector2D rij = particles[i].position - particles[j].position;
            const double r2 = rij.x * rij.x + rij.y * rij.y;

            if (r2 < smoothingRadiusSquared) {
                neighborLists[i].push_back(j);
                neighborLists[j].push_back(i);
            }
        }
    }
}

const std::vector<std::size_t>& NeighborSearch::getNeighbors(std::size_t particleIndex) const {
    return neighborLists[particleIndex];
}
