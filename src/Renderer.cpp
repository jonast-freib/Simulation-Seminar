#include "Renderer.h"
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdint>

Renderer::Renderer(unsigned int width, unsigned int height, const std::string& title, float particleRadius)
    : window(sf::VideoMode(sf::Vector2u(width, height)), title)
    , particleBuffer(sf::PrimitiveType::Triangles, sf::VertexBuffer::Usage::Stream)
    , particleVertices()
    , particleRadius(particleRadius)
    , fpsText(font, "", 24) {
    //window.setFramerateLimit(60);
    createParticleTexture();

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

void Renderer::addSlider(const std::string& name, float minValue, float maxValue, float initialValue, float x, float y, float width) {
    Slider s;
    s.name = name;
    s.minValue = minValue;
    s.maxValue = maxValue;
    s.value = std::max(minValue, std::min(maxValue, initialValue));
    s.width = width;

    s.track.setSize({width, 6.f});
    s.track.setFillColor(sf::Color(100, 100, 100));
    s.track.setPosition(sf::Vector2f(x, y));

    s.knob.setRadius(8.f);
    s.knob.setFillColor(sf::Color(200, 200, 200));
    s.knob.setOutlineColor(sf::Color::Black);
    s.knob.setOutlineThickness(1.f);

    s.updateKnob();
    sliders.push_back(s);
}

float Renderer::getSliderValue(std::size_t index) const {
    if (index >= sliders.size()) return 0.f;
    return sliders[index].value;
}

void Renderer::setSliderValue(std::size_t index, float value) {
    if (index >= sliders.size()) return;
    auto& s = sliders[index];
    s.value = std::max(s.minValue, std::min(s.maxValue, value));
    s.updateKnob();
}

void Renderer::createParticleTexture() {
    const unsigned int textureSize = 12;
    sf::Image image(sf::Vector2u(textureSize, textureSize), sf::Color::Transparent);

    const float center = static_cast<float>(textureSize) * 0.5f;
    const float radius = center - 1.f;
    const float radiusSquared = radius * radius;

    for (unsigned int y = 0; y < textureSize; ++y) {
        for (unsigned int x = 0; x < textureSize; ++x) {
            const float dx = static_cast<float>(x) + 0.5f - center;
            const float dy = static_cast<float>(y) + 0.5f - center;
            const float distanceSquared = dx * dx + dy * dy;
            const uint8_t alpha = (distanceSquared <= radiusSquared) ? 255 : 0;
            image.setPixel({x, y}, sf::Color(255, 255, 255, alpha));
        }
    }

    if (!particleTexture.loadFromImage(image)) {
        std::cerr << "Failed to create particle texture" << std::endl;
    }
    particleTexture.setSmooth(true);
}

bool Renderer::isOpen() const {
    return window.isOpen();
}

void Renderer::handleEvents() {
    while (auto event = window.pollEvent()) {
        if (event->is<sf::Event::Closed>()) {
            window.close();
        }

        // Handle slider interaction
        if (event->is<sf::Event::MouseButtonPressed>()) {
            const sf::Vector2i mp = sf::Mouse::getPosition(window);
            const sf::Vector2f mpf(static_cast<float>(mp.x), static_cast<float>(mp.y));
            for (auto& s : sliders) {
                const sf::FloatRect trackBounds = s.track.getGlobalBounds();
                if (trackBounds.contains(mpf)) {
                    s.dragging = true;
                }
            }
        } else if (event->is<sf::Event::MouseButtonReleased>()) {
            // On mouse release, stop dragging sliders and detect button presses
            const sf::Vector2i mp = sf::Mouse::getPosition(window);
            const sf::Vector2f mpf(static_cast<float>(mp.x), static_cast<float>(mp.y));
            for (auto& s : sliders) s.dragging = false;
            for (auto& b : buttons) {
                if (b.rect.getGlobalBounds().contains(mpf)) {
                    b.pressed = true;
                }
            }
        } else if (event->is<sf::Event::MouseMoved>()) {
            const sf::Vector2i mp = sf::Mouse::getPosition(window);
            const float mx = static_cast<float>(mp.x);
            for (auto& s : sliders) {
                if (!s.dragging) continue;
                const float trackLeft = s.track.getPosition().x;
                const float t = std::max(0.f, std::min(s.width, mx - trackLeft)) / s.width;
                s.value = s.minValue + t * (s.maxValue - s.minValue);
                s.updateKnob();
            }
        }
    }
}

void Renderer::clear() {
    window.clear(sf::Color::Black);
}

void Renderer::drawParticles(const std::vector<Particle>& particles) {
    const float halfSize = particleRadius;
    const std::size_t requiredVertices = particles.size() * 6;

    if (particleVertices.size() != requiredVertices) {
        particleVertices.resize(requiredVertices);
        if (!particleBuffer.create(requiredVertices)) {
            std::cerr << "Failed to create particle vertex buffer" << std::endl;
            return;
        }
    }

    const float texWidth = static_cast<float>(particleTexture.getSize().x);
    const float texHeight = static_cast<float>(particleTexture.getSize().y);

    // determine pressure bounds for coloring
    double minP = std::numeric_limits<double>::infinity();
    double maxP = -std::numeric_limits<double>::infinity();
    for (const auto& p : particles) {
        minP = std::min(minP, p.pressure);
        maxP = std::max(maxP, p.pressure);
    }
    if (minP == std::numeric_limits<double>::infinity()) {
        minP = 0.0; maxP = 1.0;
    }

    auto pressureToColor = [&](const Particle& particle) -> sf::Color {
        // If CFL violated, show purple
        if (particle.cflViolated) {
            return sf::Color(200, 0, 0, 220);  // Purple
        }

        double n = 0.0;
        if (maxP > minP) n = (particle.pressure - minP) / (maxP - minP);
        n = std::max(0.0, std::min(1.0, n));
        // simple blue -> cyan -> green -> yellow -> red mapping
        const double r = std::min(1.0, std::max(0.0, 1.5 * n));
        const double g = std::max(0.0, 1.0 - std::abs(n - 0.5) * 2.0);
        const double b = std::min(1.0, std::max(0.0, 1.5 * (1.0 - n)));
        const uint8_t R = static_cast<uint8_t>(std::round(r * 255.0));
        const uint8_t Gc = static_cast<uint8_t>(std::round(g * 255.0));
        const uint8_t B = static_cast<uint8_t>(std::round(b * 255.0));
        return sf::Color(R, Gc, B, 220);
    };

    for (std::size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const float cx = static_cast<float>(particle.position.x);
        const float cy = static_cast<float>(particle.position.y);
        const sf::Color color = pressureToColor(particle);
        const std::size_t baseIndex = i * 6;

        particleVertices[baseIndex + 0].position = {cx - halfSize, cy - halfSize};
        particleVertices[baseIndex + 1].position = {cx + halfSize, cy - halfSize};
        particleVertices[baseIndex + 2].position = {cx + halfSize, cy + halfSize};
        particleVertices[baseIndex + 3].position = {cx - halfSize, cy - halfSize};
        particleVertices[baseIndex + 4].position = {cx + halfSize, cy + halfSize};
        particleVertices[baseIndex + 5].position = {cx - halfSize, cy + halfSize};

        particleVertices[baseIndex + 0].texCoords = {0.f, 0.f};
        particleVertices[baseIndex + 1].texCoords = {texWidth, 0.f};
        particleVertices[baseIndex + 2].texCoords = {texWidth, texHeight};
        particleVertices[baseIndex + 3].texCoords = {0.f, 0.f};
        particleVertices[baseIndex + 4].texCoords = {texWidth, texHeight};
        particleVertices[baseIndex + 5].texCoords = {0.f, texHeight};

        for (std::size_t v = 0; v < 6; ++v) {
            particleVertices[baseIndex + v].color = color;
        }
    }

    if (!particleBuffer.update(particleVertices.data())) {
        std::cerr << "Failed to update particle vertex buffer" << std::endl;
        return;
    }

    window.draw(particleBuffer, &particleTexture);

    // Draw UI sliders on top
    for (const auto& s : sliders) {
        window.draw(s.track);
        window.draw(s.knob);
        if (fontLoaded) {
            sf::Text label(font, s.name + ": " + std::to_string(s.value), 14);
            label.setFillColor(sf::Color::White);
            label.setPosition(sf::Vector2f(s.track.getPosition().x, s.track.getPosition().y - 20.f));
            window.draw(label);
        }
    }

    // Draw buttons
    for (const auto& b : buttons) {
        window.draw(b.rect);
        if (fontLoaded) {
            sf::Text label(font, b.name, 16);
            label.setFillColor(sf::Color::White);
            // center label in button
            const sf::Vector2f pos = b.rect.getPosition();
            const sf::Vector2f size = b.rect.getSize();
            label.setPosition(sf::Vector2f(pos.x + 6.f, pos.y + (size.y - 16.f) * 0.5f));
            window.draw(label);
        }
    }
}

void Renderer::display() {
    
    window.display();
}

void Renderer::addButton(const std::string& name, float x, float y, float width, float height) {
    Button b;
    b.name = name;
    b.rect.setSize(sf::Vector2f(width, height));
    b.rect.setFillColor(sf::Color(70, 70, 70));
    b.rect.setOutlineColor(sf::Color::White);
    b.rect.setOutlineThickness(1.f);
    b.rect.setPosition(sf::Vector2f(x, y));
    b.pressed = false;
    buttons.push_back(b);
}

bool Renderer::consumeButtonPress(std::size_t index) {
    if (index >= buttons.size()) return false;
    bool p = buttons[index].pressed;
    buttons[index].pressed = false;
    return p;
}
