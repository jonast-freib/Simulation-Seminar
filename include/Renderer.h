#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include "Particle.h"

class Renderer {
public:
    Renderer(unsigned int width, unsigned int height, const std::string& title, float particleRadius = 2.f);
    bool isOpen() const;
    void handleEvents();
    void clear();
    void drawParticles(const std::vector<Particle>& particles);
    void display();

private:
    sf::RenderWindow window;
    sf::VertexArray particleVertices;
    float particleRadius;
    sf::Font font;
    sf::Text fpsText;
    sf::Clock fpsClock;
    int frameCount = 0;
    float timeElapsed = 0.f;
    bool fontLoaded = false;
};
