#include <cstddef>
#include <cmath>
#include <optional>

#include "SFML/Graphics.hpp"

#include "args.hpp"
#include "plotting.hpp"

std::size_t GraphBase::n() const { return this->vertexes.getVertexCount(); }

sf::Vertex * GraphBase::begin() {
  return &(this->vertexes[0]);
}
sf::Vertex * GraphBase::end() {
  return &(this->vertexes[vertexes.getVertexCount() - 1]) + 1;
}

GraphBase::GraphBase(sf::VertexArray const vertexes) : vertexes{vertexes} {}

GraphBase::GraphBase(sf::PrimitiveType const primitive,
    std::size_t const n_vertexes, sf::Color const &init_color) : vertexes{[&](){
      sf::VertexArray vertexes(primitive, n_vertexes);
      for (std::size_t i{0}; i < vertexes.getVertexCount(); ++i)
        vertexes[i].color = init_color;
      return vertexes;
    }()} {}

void draw(sf::RenderTarget &target, GraphBase const &graph) {
  target.draw(graph.vertexes);
}

sf::Vector2f thickness_px_to_2d(float const thickness,
    sf::RenderTarget const &target) {
  auto size_window{target.getSize()};
  auto size_view{target.getView().getSize()};
  return {
    std::abs(thickness * size_view.x) / size_window.x,
    std::abs(thickness * size_view.y) / size_window.y};
}

std::size_t Curve::n() const {
  return this->vertexes.getVertexCount() / 4 + 1;
}

Curve::Graph(std::size_t const n, sf::Color const &init_color,
    sf::Vector2f const &thickness_2d) :
    GraphBase(sf::PrimitiveType::TriangleStrip, (n - 1) * 4, init_color),
    thickness_2d{thickness_2d} {}

std::size_t Fill::n() const {
  return this->vertexes.getVertexCount() / 2;
}

Fill::Graph(std::size_t const n, sf::Color const &init_color,
    float const to, bool const vertical) :
    GraphBase(sf::PrimitiveType::TriangleStrip, n * 2, init_color),
    to{to}, vertical{vertical} {}

