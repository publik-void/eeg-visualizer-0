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

GraphBase::state_t iterative_set_init(GraphBase &);

GraphBase::state_t iterative_set_point(GraphBase &,
    sf::Vector2f const &, GraphBase::state_t const &);

GraphBase::state_t iterative_set_point(GraphBase &,
    sf::Vector2f const &, sf::Color const &, GraphBase::state_t const &);

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

sf::Vertex * get_segment(Curve &, std::size_t const = 0);

sf::Vertex * prev_segment(Curve const &, sf::Vertex * const);

sf::Vertex * next_segment(Curve const &, sf::Vertex * const);

sf::Vertex * set_segment(Curve &, sf::Vertex *,
    sf::Vector2f const &, sf::Vector2f const &);

sf::Vertex * set_segment(Curve &, std::size_t const,
    sf::Vector2f const &, sf::Vector2f const &);

sf::Vertex * set_segment(Curve &, sf::Vertex *,
    sf::Color const &, sf::Color const &);

sf::Vertex * set_segment(Curve &, std::size_t const,
    sf::Color const &, sf::Color const &);

Curve::state_t iterative_set_init(Curve &,
    std::optional<sf::Color> const = {});

Curve::state_t iterative_set_point(Curve &, sf::Vector2f const &,
    Curve::state_t const &);

Curve::state_t iterative_set_point(Curve &, sf::Vector2f const &,
    sf::Color const &, Curve::state_t const &);

template<>
struct Graph<args::get_mode_index(modes_graph, "fill")> : public GraphBase {
  float to;
  bool vertical;

  std::size_t n() const;

  Graph(std::size_t const, sf::Color const &,
      float const = 0.f, bool const = false);
};

using Fill = Graph<args::get_mode_index(modes_graph, "fill")>;

sf::Vertex * get_point(Fill &, std::size_t const = 0);

sf::Vertex * prev_point(Fill const &, sf::Vertex * const);

sf::Vertex * next_point(Fill const &, sf::Vertex * const);

sf::Vertex * set_point(Fill &, sf::Vertex * const, sf::Vector2f const &);

sf::Vertex * set_point(Fill &, sf::Vertex * const,
    sf::Vector2f const &, sf::Color const &, sf::Color const &);

sf::Vertex * set_point(Fill &, sf::Vertex * const,
    sf::Vector2f const &, sf::Color const &);

Fill::state_t iterative_set_point(Fill &, sf::Vector2f const &,
    Fill::state_t const &);

Fill::state_t iterative_set_point(Fill &,
    sf::Vector2f const &, sf::Color const &, Fill::state_t const &);

