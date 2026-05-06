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
    sf::Color color;

    Particle(Vector2D pos = Vector2D(), double m = 1.0)
        : position(pos), velocity(0, 0), acceleration(0, 0),
          mass(m), density(0.0), pressure(0.0) {
        // Generate random color
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, 255);
        color = sf::Color(dis(gen), dis(gen), dis(gen), 255);
    }
};