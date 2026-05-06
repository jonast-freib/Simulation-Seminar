#include "Renderer.h"

Renderer::Renderer(unsigned int width, unsigned int height, const std::string& title)
    : window(sf::VideoMode(sf::Vector2u(width, height)), title) {
    window.setFramerateLimit(60);
}

bool Renderer::isOpen() const {
    return window.isOpen();
}

void Renderer::handleEvents() {
    while (auto event = window.pollEvent()) {
        if (event->is<sf::Event::Closed>()) {
            window.close();
        }
    }
}

void Renderer::clear() {
    window.clear(sf::Color::Black);
}

void Renderer::drawParticles(const std::vector<Particle>& particles) {
    for (const auto& particle : particles) {
        sf::CircleShape circle(2.0); // radius 2
        circle.setPosition(sf::Vector2f(particle.position.x - 2, particle.position.y - 2));
        circle.setFillColor(particle.color);
        window.draw(circle);
    }
}

void Renderer::display() {
    window.display();
}
