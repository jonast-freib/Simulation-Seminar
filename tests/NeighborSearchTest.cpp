#include "NeighborSearch.h"
#include <cassert>
#include <iostream>

static void expect(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "Test failed: " << message << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

static void printNeighbors(std::size_t index, const std::vector<std::size_t>& neighbors) {
    std::cout << "Particle " << index << " neighbors: [";
    for (std::size_t i = 0; i < neighbors.size(); ++i) {
        std::cout << neighbors[i];
        if (i + 1 < neighbors.size()) std::cout << ", ";
    }
    std::cout << "] (count=" << neighbors.size() << ")" << std::endl;
}

static void testEmptyParticles() {
    std::vector<Particle> particles;
    NeighborSearch search(50.0);
    search.build(particles);
    std::cout << "testEmptyParticles: particleCount=0, radius=50" << std::endl;
    expect(true, "Empty particle list should not crash");
}

static void testSingleParticle() {
    std::vector<Particle> particles;
    particles.emplace_back(Vector2D(0.0, 0.0), 1.0, false);
    NeighborSearch search(50.0);
    search.build(particles);
    const auto& neighbors = search.getNeighbors(0);
    printNeighbors(0, neighbors);
    expect(neighbors.empty(), "Single particle should have no neighbors");
}

static void testTwoNeighborsWithinRadius() {
    std::vector<Particle> particles;
    particles.emplace_back(Vector2D(0.0, 0.0), 1.0, false);
    particles.emplace_back(Vector2D(30.0, 0.0), 1.0, false);

    const double radius = 50.0;
    NeighborSearch search(radius);
    search.build(particles);

    const auto& neighbors0 = search.getNeighbors(0);
    const auto& neighbors1 = search.getNeighbors(1);

    std::cout << "testTwoNeighborsWithinRadius: radius=" << radius << ", distance=30" << std::endl;
    printNeighbors(0, neighbors0);
    printNeighbors(1, neighbors1);

    expect(neighbors0.size() == 1, "Particle 0 should have one neighbor");
    expect(neighbors0[0] == 1, "Particle 0 neighbor should be particle 1");
    expect(neighbors1.size() == 1, "Particle 1 should have one neighbor");
    expect(neighbors1[0] == 0, "Particle 1 neighbor should be particle 0");
}

static void testTwoParticlesOutsideRadius() {
    std::vector<Particle> particles;
    particles.emplace_back(Vector2D(0.0, 0.0), 1.0, false);
    particles.emplace_back(Vector2D(60.0, 0.0), 1.0, false);

    const double radius = 50.0;
    NeighborSearch search(radius);
    search.build(particles);

    const auto& neighbors0 = search.getNeighbors(0);
    const auto& neighbors1 = search.getNeighbors(1);

    std::cout << "testTwoParticlesOutsideRadius: radius=" << radius << ", distance=60" << std::endl;
    printNeighbors(0, neighbors0);
    printNeighbors(1, neighbors1);

    expect(neighbors0.empty(), "Particle 0 should have no neighbors outside radius");
    expect(neighbors1.empty(), "Particle 1 should have no neighbors outside radius");
}

static void testVariableRadius() {
    std::vector<Particle> particles;
    particles.emplace_back(Vector2D(0.0, 0.0), 1.0, false);
    particles.emplace_back(Vector2D(40.0, 0.0), 1.0, false);

    NeighborSearch search(30.0);
    search.build(particles);
    std::cout << "testVariableRadius: start radius=30, distance=40" << std::endl;
    printNeighbors(0, search.getNeighbors(0));
    expect(search.getNeighbors(0).empty(), "Particles should not be neighbors at radius 30");

    search.setRadius(50.0);
    search.build(particles);
    const auto& neighbors0 = search.getNeighbors(0);
    std::cout << "testVariableRadius: updated radius=50" << std::endl;
    printNeighbors(0, neighbors0);

    expect(neighbors0.size() == 1, "Particles should be neighbors at radius 50");
    expect(neighbors0[0] == 1, "Neighbor index should be 1 after increasing radius");
}

static void testLargeParticleSet() {
    const int rows = 10;
    const int cols = 10;
    const double spacing = 20.0;
    const double radius = 30.0;

    std::vector<Particle> particles;
    particles.reserve(rows * cols);
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {
            particles.emplace_back(Vector2D(x * spacing, y * spacing), 1.0, false);
        }
    }

    NeighborSearch search(radius);
    search.build(particles);

    const std::size_t centerIndex = (rows / 2) * cols + (cols / 2);
    const std::size_t cornerIndex = 0;
    const std::size_t edgeIndex = (rows / 2) * cols;

    const auto& centerNeighbors = search.getNeighbors(centerIndex);
    const auto& cornerNeighbors = search.getNeighbors(cornerIndex);
    const auto& edgeNeighbors = search.getNeighbors(edgeIndex);

    std::cout << "testLargeParticleSet: grid=" << rows << "x" << cols << ", radius=" << radius << ", spacing=" << spacing << std::endl;
    printNeighbors(centerIndex, centerNeighbors);
    printNeighbors(cornerIndex, cornerNeighbors);
    printNeighbors(edgeIndex, edgeNeighbors);

    expect(centerNeighbors.size() == 8, "Center particle should have 8 neighbors in the grid");
    expect(cornerNeighbors.size() == 3, "Corner particle should have 3 neighbors in the grid");
    expect(edgeNeighbors.size() == 5, "Edge particle should have 5 neighbors in the grid");
}

int main() {
    testEmptyParticles();
    testSingleParticle();
    testTwoNeighborsWithinRadius();
    testTwoParticlesOutsideRadius();
    testVariableRadius();
    testLargeParticleSet();

    std::cout << "NeighborSearch tests passed." << std::endl;
    return 0;
}
