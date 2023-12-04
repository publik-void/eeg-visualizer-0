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

GraphBase::state_t iterative_set_init(GraphBase &self) {
  return self.begin();
}

GraphBase::state_t iterative_set_point(GraphBase &,
    sf::Vector2f const &p, GraphBase::state_t const &s) {
  s->position = p;
  return s + 1;
}

GraphBase::state_t iterative_set_point(GraphBase &self,
    sf::Vector2f const &p, sf::Color const &c, GraphBase::state_t const &s) {
  s->color = c;
  return iterative_set_point(self, p, s);
}

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

sf::Vertex * get_segment(Curve &self, std::size_t const i) {
  return &(self.vertexes[i * 4]);
}

sf::Vertex * prev_segment(Curve const &, sf::Vertex * const segment) {
  return segment - 4;
}

sf::Vertex * next_segment(Curve const &, sf::Vertex * const segment) {
  return segment + 4;
}

sf::Vertex * set_segment(Curve &self, sf::Vertex * segment,
    sf::Vector2f const &p0, sf::Vector2f const &p1) {
  sf::Vector2f const nv{p0.y - p1.y, p1.x - p0.x};
  //sf::Vector2f const nv{p1.y - p0.y, p0.x - p1.x};
  float const length{std::sqrt((nv.x * nv.x) + (nv.y * nv.y))};
  sf::Vector2f const snv{nv.x * self.thickness_2d.x / length,
                         nv.y * self.thickness_2d.y / length};

  (segment++)->position = p0 + snv;
  (segment++)->position = p0 - snv;
  (segment++)->position = p1 + snv;
  (segment++)->position = p1 - snv;
  return segment;
}

sf::Vertex * set_segment(Curve &self, std::size_t const i,
    sf::Vector2f const &p0, sf::Vector2f const &p1) {
  return set_segment(self, get_segment(self, i), p0, p1);
}

sf::Vertex * set_segment(Curve &, sf::Vertex * segment,
    sf::Color const &c0, sf::Color const &c1) {
  (segment++)->color = c0;
  (segment++)->color = c0;
  (segment++)->color = c1;
  (segment++)->color = c1;
  return segment;
}

sf::Vertex * set_segment(Curve &self, std::size_t const i,
    sf::Color const &c0, sf::Color const &c1) {
  return set_segment(self, get_segment(self, i), c0, c1);
}

Curve::state_t iterative_set_init(Curve &self,
    std::optional<sf::Color> const c_init) {
  return {prev_segment(self, get_segment(self)), {},
    c_init.value_or(get_segment(self)->color)};
}

Curve::state_t iterative_set_point(Curve &self, sf::Vector2f const &p,
    Curve::state_t const &s) {
  sf::Vertex * const next{s.segment >= get_segment(self)
    ? set_segment(self, s.segment, s.p, p)
    : next_segment(self, s.segment)};
  return {next, p, s.c};
}

Curve::state_t iterative_set_point(Curve &self, sf::Vector2f const &p,
    sf::Color const &c, Curve::state_t const &s) {
  if (s.segment >= get_segment(self)) set_segment(self, s.segment, s.c, c);
  return iterative_set_point(self, p, {s.segment, s.p, c});
}

std::size_t Fill::n() const {
  return this->vertexes.getVertexCount() / 2;
}

Fill::Graph(std::size_t const n, sf::Color const &init_color,
    float const to, bool const vertical) :
    GraphBase(sf::PrimitiveType::TriangleStrip, n * 2, init_color),
    to{to}, vertical{vertical} {}

sf::Vertex * get_point(Fill &self, std::size_t const i) {
  return &(self.vertexes[i * 2]);
}

sf::Vertex * prev_point(Fill const &, sf::Vertex * const point) {
  return point - 2;
}

sf::Vertex * next_point(Fill const &, sf::Vertex * const point) {
  return point + 2;
}

sf::Vertex * set_point(Fill &self, sf::Vertex * const point,
    sf::Vector2f const &p) {
  auto const point_fill{point};
  auto const point_value{std::next(point)};
  // NOTE: Definitely potential for optimization here, e.g. not branching on
  // every call and not setting one coordinate to the same `to` value on every
  // call.
  if (not self.vertical) {
    point_fill->position.x = p.x;
    point_fill->position.y = self.to;
  } else {
    point_fill->position.x = self.to;
    point_fill->position.y = p.y;
  }
  point_value->position.x = p.x;
  point_value->position.y = p.y;
  return next_point(self, point);
}

sf::Vertex * set_point(Fill &self, sf::Vertex * const point,
    sf::Vector2f const &p, sf::Color const &c0, sf::Color const &c1) {
  auto const point_fill{point};
  auto const point_value{std::next(point)};
  point_fill->color = c1;
  point_value->color = c0;
  return set_point(self, point, p);
}

sf::Vertex * set_point(Fill &self, sf::Vertex * const point,
    sf::Vector2f const &p, sf::Color const &c) {
  return set_point(self, point, p, c, c);
}

Fill::state_t iterative_set_point(Fill &self, sf::Vector2f const &p,
    Fill::state_t const &s) {
  return set_point(self, s, p);
}

Fill::state_t iterative_set_point(Fill &self,
    sf::Vector2f const &p, sf::Color const &c, Fill::state_t const &s) {
  return set_point(self, s, p, c);
}


