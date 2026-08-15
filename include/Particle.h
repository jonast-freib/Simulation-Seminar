#pragma once

#include <SFML/System/Vector2.hpp>
#include <SFML/Graphics/Color.hpp>
#include <random>

using Vector2D = sf::Vector2<double>;

struct Particle {
    Vector2D position;
    Vector2D velocity;
    Vector2D acceleration;
    double mass;
    double density;
    double pressure;
    bool fixed;
    bool cflViolated;
    sf::Color color;

    Particle(Vector2D pos = Vector2D(), double m = 2.5, bool fixedParticle = false)
        : position(pos), velocity(0, 0), acceleration(0, 0),
          mass(m), density(0.0), pressure(0.0), fixed(fixedParticle), cflViolated(false) {}
};