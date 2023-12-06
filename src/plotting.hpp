#pragma once

#include "SFML/Graphics.hpp"

#include "args.hpp"

struct GraphBase {
  sf::VertexArray vertexes;
  using state_t = sf::Vertex *;

  std::size_t n() const;

  sf::Vertex * begin();
  sf::Vertex * end();

  GraphBase(sf::VertexArray const vertexes);
  GraphBase(sf::PrimitiveType const primitive, std::size_t const n_vertexes,
      sf::Color const &init_color);
};

inline GraphBase::state_t iterative_set_init(GraphBase &self) {
  return self.begin();
}

inline GraphBase::state_t iterative_set_point(GraphBase &,
    sf::Vector2f const &p, GraphBase::state_t const &s) {
  s->position = p;
  return s + 1;
}

inline GraphBase::state_t iterative_set_point(GraphBase &self,
    sf::Vector2f const &p, sf::Color const &c, GraphBase::state_t const &s) {
  s->color = c;
  return iterative_set_point(self, p, s);
}

void draw(sf::RenderTarget &, GraphBase const &);

sf::Vector2f thickness_px_to_2d(float const, sf::RenderTarget const &);

args::modes_t<2> constexpr modes_graph{"curve", "fill"};

template<std::size_t mode>
struct Graph : public GraphBase {};

template<>
struct Graph<args::get_mode_index(modes_graph, "curve")> : public GraphBase {
  sf::Vector2f thickness_2d;

  struct state_t {
    sf::Vertex * segment;
    sf::Vector2f p;
    sf::Color c;
  };

  std::size_t n() const;

  Graph(std::size_t const n, sf::Color const &, sf::Vector2f const &);
};

using Curve = Graph<args::get_mode_index(modes_graph, "curve")>;

inline sf::Vertex * get_segment(Curve &self, std::size_t const i = 0) {
  return &(self.vertexes[i * 4]);
}

inline sf::Vertex * prev_segment(Curve const &, sf::Vertex * const segment) {
  return segment - 4;
}

inline sf::Vertex * next_segment(Curve const &, sf::Vertex * const segment) {
  return segment + 4;
}

inline sf::Vertex * set_segment(Curve &self, sf::Vertex * segment,
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

inline sf::Vertex * set_segment(Curve &self, std::size_t const i,
    sf::Vector2f const &p0, sf::Vector2f const &p1) {
  return set_segment(self, get_segment(self, i), p0, p1);
}

inline sf::Vertex * set_segment(Curve &, sf::Vertex * segment,
    sf::Color const &c0, sf::Color const &c1) {
  (segment++)->color = c0;
  (segment++)->color = c0;
  (segment++)->color = c1;
  (segment++)->color = c1;
  return segment;
}

inline sf::Vertex * set_segment(Curve &self, std::size_t const i,
    sf::Color const &c0, sf::Color const &c1) {
  return set_segment(self, get_segment(self, i), c0, c1);
}

inline Curve::state_t iterative_set_init(Curve &self,
    std::optional<sf::Color> const c_init = {}) {
  return {prev_segment(self, get_segment(self)), {},
    c_init.value_or(get_segment(self)->color)};
}

inline Curve::state_t iterative_set_point(Curve &self, sf::Vector2f const &p,
    Curve::state_t const &s) {
  sf::Vertex * const next{s.segment >= get_segment(self)
    ? set_segment(self, s.segment, s.p, p)
    : next_segment(self, s.segment)};
  return {next, p, s.c};
}

inline Curve::state_t iterative_set_point(Curve &self, sf::Vector2f const &p,
    sf::Color const &c, Curve::state_t const &s) {
  if (s.segment >= get_segment(self)) set_segment(self, s.segment, s.c, c);
  return iterative_set_point(self, p, {s.segment, s.p, c});
}

template<>
struct Graph<args::get_mode_index(modes_graph, "fill")> : public GraphBase {
  float to;
  bool vertical;

  std::size_t n() const;

  Graph(std::size_t const, sf::Color const &,
      float const = 0.f, bool const = false);
};

using Fill = Graph<args::get_mode_index(modes_graph, "fill")>;

inline sf::Vertex * get_point(Fill &self, std::size_t const i = 0) {
  return &(self.vertexes[i * 2]);
}

inline sf::Vertex * prev_point(Fill const &, sf::Vertex * const point) {
  return point - 2;
}

inline sf::Vertex * next_point(Fill const &, sf::Vertex * const point) {
  return point + 2;
}

inline sf::Vertex * set_point(Fill &self, sf::Vertex * const point,
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

inline sf::Vertex * set_point(Fill &self, sf::Vertex * const point,
    sf::Vector2f const &p, sf::Color const &c0, sf::Color const &c1) {
  auto const point_fill{point};
  auto const point_value{std::next(point)};
  point_fill->color = c1;
  point_value->color = c0;
  return set_point(self, point, p);
}

inline sf::Vertex * set_point(Fill &self, sf::Vertex * const point,
    sf::Vector2f const &p, sf::Color const &c) {
  return set_point(self, point, p, c, c);
}

inline Fill::state_t iterative_set_point(Fill &self, sf::Vector2f const &p,
    Fill::state_t const &s) {
  return set_point(self, s, p);
}

inline Fill::state_t iterative_set_point(Fill &self,
    sf::Vector2f const &p, sf::Color const &c, Fill::state_t const &s) {
  return set_point(self, s, p, c);
}

//struct PlotAnnotations {

//};
