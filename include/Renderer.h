#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include "Particle.h"

class Renderer {
public:
    Renderer(unsigned int width, unsigned int height, const std::string& title, float particleRadius = 1.f);
    bool isOpen() const;
    void handleEvents();
    void clear();
    void drawParticles(const std::vector<Particle>& particles);
    void display();
    // UI sliders
    void addSlider(const std::string& name, float minValue, float maxValue, float initialValue, float x, float y, float width = 200.f);
    float getSliderValue(std::size_t index) const;
    void setSliderValue(std::size_t index, float value);

    // UI buttons
    void addButton(const std::string& name, float x, float y, float width = 80.f, float height = 28.f);
    bool consumeButtonPress(std::size_t index);

private:
    void createParticleTexture();

    sf::RenderWindow window;
    sf::VertexBuffer particleBuffer;
    sf::Texture particleTexture;
    std::vector<sf::Vertex> particleVertices;
    float particleRadius;
    sf::Font font;
    sf::Text fpsText;
    sf::Clock fpsClock;
    int frameCount = 0;
    float timeElapsed = 0.f;
    bool fontLoaded = false;
    struct Slider {
        std::string name;
        float minValue;
        float maxValue;
        float value;
        sf::RectangleShape track;
        sf::CircleShape knob;
        // label text will be constructed at draw time (requires font)
        bool dragging = false;
        float width = 200.f;
        void updateKnob() {
            const float t = (value - minValue) / (maxValue - minValue);
            const float trackLeft = track.getPosition().x;
            const float trackTop = track.getPosition().y;
            const float knobX = trackLeft + t * width;
            const float knobY = trackTop + track.getSize().y * 0.5f;
            knob.setPosition(sf::Vector2f(knobX - knob.getRadius(), knobY - knob.getRadius()));
        }
    };
    std::vector<Slider> sliders;

    struct Button {
        std::string name;
        sf::RectangleShape rect;
        bool pressed = false;
    };
    std::vector<Button> buttons;
};
