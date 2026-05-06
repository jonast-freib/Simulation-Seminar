#include "Renderer.h"
#include <cmath>
#include <string>
#include <iostream>

Renderer::Renderer(unsigned int width, unsigned int height, const std::string& title, float particleRadius)
    : window(sf::VideoMode(sf::Vector2u(width, height)), title)
    , particleVertices(sf::PrimitiveType::Triangles)
    , particleRadius(particleRadius)
    , fpsText(font, "", 24) {
    window.setFramerateLimit(60);

    if (!font.openFromFile("/usr/share/fonts/Adwaita/AdwaitaSans-Regular.ttf")) {
        // Fallback: try another common font path
        if (!font.openFromFile("/usr/share/fonts/truetype/freefont/FreeSans.ttf")) {
            // If no font found, just skip FPS display
            std::cerr << "Warning: Could not load font for FPS display" << std::endl;
            fontLoaded = false;
        } else {
            fontLoaded = true;
        }
    } else {
        fontLoaded = true;
    }

    if (fontLoaded) {
        fpsText.setFillColor(sf::Color::White);
        fpsText.setPosition(sf::Vector2f(10.f, 10.f));
    }
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
    const float radius = particleRadius;
    const unsigned int segments = 4;
    const float twoPi = 6.28318530718f;
    const std::size_t verticesPerCircle = segments * 3;

    particleVertices.clear();
    particleVertices.resize(particles.size() * verticesPerCircle);

    for (std::size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const float cx = static_cast<float>(particle.position.x);
        const float cy = static_cast<float>(particle.position.y);
        const sf::Color color = particle.color;
        const std::size_t baseIndex = i * verticesPerCircle;

        for (unsigned int s = 0; s < segments; ++s) {
            const float angle0 = (static_cast<float>(s) / segments) * twoPi;
            const float angle1 = (static_cast<float>(s + 1) / segments) * twoPi;
            const std::size_t triIndex = baseIndex + s * 3;

            particleVertices[triIndex + 0].position = {cx, cy};
            particleVertices[triIndex + 1].position = {cx + std::cos(angle0) * radius, cy + std::sin(angle0) * radius};
            particleVertices[triIndex + 2].position = {cx + std::cos(angle1) * radius, cy + std::sin(angle1) * radius};

            particleVertices[triIndex + 0].color = color;
            particleVertices[triIndex + 1].color = color;
            particleVertices[triIndex + 2].color = color;
        }
    }

    window.draw(particleVertices);
}

void Renderer::display() {
    // Calculate FPS
    timeElapsed += fpsClock.restart().asSeconds();
    frameCount++;

    if (timeElapsed >= 1.f) {
        float fps = frameCount / timeElapsed;
        fpsText.setString("FPS: " + std::to_string(static_cast<int>(fps)));
        frameCount = 0;
        timeElapsed = 0.f;
    }

    // Draw FPS text if font is loaded
    if (fontLoaded) {
        window.draw(fpsText);
    }

    window.display();
}
