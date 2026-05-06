#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include "Particle.h"

class Renderer {
public:
    Renderer(unsigned int width, unsigned int height, const std::string& title);
    bool isOpen() const;
    void handleEvents();
    void clear();
    void drawParticles(const std::vector<Particle>& particles);
    void display();

private:
    sf::RenderWindow window;
};
