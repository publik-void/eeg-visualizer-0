#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>
#include <string>

#include "SFML/Graphics.hpp"

#include "args.hpp"
#include "lerp.hpp"
#include "affine.hpp"
#include "box.hpp"

using ViewTransform = Affine<float, 2, 2>;

std::pair<ViewTransform, ViewTransform> as_view_transform(
    sf::Vector2f const & = {1.f, 1.f}, sf::Vector2f const & = {0.f, 0.f},
    float const = 0.f);

std::pair<array2d<float, 2, 2>, array2d<float, 2, 2>>
  as_rotation_matrix(float const = 0.f);

Box<float> get_render_target_box(sf::RenderTarget const &);
Box<float> sf_view_to_box(sf::View const &);

template <typename T>
struct RGB {
  T r;
  T g;
  T b;
};

using RGB8 = RGB<std::uint8_t>;

inline sf::Color rgb_to_sf_color(RGB8 const &c) { return {c.r, c.g, c.b}; }
inline RGB8 sf_color_to_rgb(sf::Color const &c) { return {c.r, c.g, c.b}; }

inline sf::Vertex * set_vertex(sf::Vertex * const vertex,
    std::optional<float> const &x_opt = {},
    std::optional<float> const &y_opt = {},
    std::optional<RGB8> const &rgb_opt = {},
    std::optional<std::uint8_t> const &a_opt = {}) {
  if (x_opt.has_value()) vertex->position.x = x_opt.value();
  if (y_opt.has_value()) vertex->position.y = y_opt.value();
  if (rgb_opt.has_value()) {
    vertex->color.r = rgb_opt->r;
    vertex->color.g = rgb_opt->g;
    vertex->color.b = rgb_opt->b;
  }
  if (a_opt.has_value()) vertex->color.a = a_opt.value();
  return std::next(vertex);
}

inline sf::Vertex * set_vertex(sf::Vertex * const vertex,
    std::optional<sf::Vector2f> const &position_opt = {},
    std::optional<RGB8> const &rgb_opt = {},
    std::optional<std::uint8_t> const &a_opt = {}) {
  std::optional<float> x_opt = {}, y_opt = {};
  if (position_opt.has_value()) {
    x_opt = position_opt->x;
    y_opt = position_opt->y;
  }
  return set_vertex(vertex, x_opt, y_opt, rgb_opt, a_opt);
}

void set_n(sf::VertexArray &, std::size_t const);
void set_primitive_type(sf::VertexArray &, sf::PrimitiveType const);

struct GraphBase {
  sf::VertexArray vertexes;
  using state_t = sf::Vertex *;

  std::size_t n() const;

  sf::Vertex * begin();
  sf::Vertex * end();

  GraphBase(sf::VertexArray const vertexes);
  GraphBase(sf::PrimitiveType const primitive = sf::PrimitiveType::Points,
      std::size_t const n_vertexes = 0, sf::Color const &init_color = {});
};

void set_n(GraphBase &, std::size_t const);
void set_primitive_type(GraphBase &, sf::PrimitiveType const);

inline GraphBase::state_t iterative_set_init(GraphBase &self) {
  return self.begin();
}

inline GraphBase::state_t iterative_set_point(GraphBase &,
    GraphBase::state_t const &s, ViewTransform const &vt,
    std::optional<sf::Vector2f> const &p = {},
    std::optional<RGB8> const &rgb = {},
    std::optional<std::uint8_t> const &a = {}) {
  return set_vertex(s, vt(p), rgb, a);
}

void draw(sf::RenderTarget &, GraphBase const &);

sf::Vector2f thickness_px_to_2d(float const, sf::RenderTarget const &);

args::modes_t<3> constexpr modes_graph{"lines", "curve", "fill"};

template<std::size_t mode>
struct Graph : public GraphBase {};

template<>
struct Graph<args::get_mode_index(modes_graph, "lines")> : public GraphBase {
  struct state_t {
    sf::Vertex * segment;
    sf::Vector2f r;
    std::optional<RGB8> c;
    std::optional<float> a;
  };

  std::size_t n() const;

  Graph(std::size_t const n = 0, sf::Color const & = {});
};

using Lines = Graph<args::get_mode_index(modes_graph, "lines")>;

void set_n(Lines &, std::size_t const);

inline sf::Vertex * get_segment(Lines &self, std::size_t const i = 0) {
  return &(self.vertexes[i * 6]);
}

inline sf::Vertex * prev_segment(Lines const &, sf::Vertex * const segment) {
  return segment - 6;
}

inline sf::Vertex * next_segment(Lines const &, sf::Vertex * const segment) {
  return segment + 6;
}

inline sf::Vertex * set_segment(Lines &, sf::Vertex * const segment,
    ViewTransform const &vt,
    sf::Vector2f const &p0, sf::Vector2f const &p1,
    sf::Vector2f const &r, std::optional<RGB8> const &c = {},
    std::optional<float> const &a = {}) {
  auto const _p0{vt(p0)}, _p1{vt(p1)};
  sf::Vector2f const nv{_p0.y - _p1.y, _p1.x - _p0.x};
  //sf::Vector2f const nv{_p1.y - _p0.y, _p0.x - _p1.x};
  float const length{std::sqrt(abs2(sf_vec_to_array(nv)))};
  sf::Vector2f const nnv{nv.x / length, nv.y / length};
  sf::Vector2f const rnv{sf::Vector2f{nnv.x * r.x, nnv.y * r.y}};

  auto const p00{_p0 - rnv}, p01{_p0 + rnv}, p10{_p1 - rnv}, p11{_p1 + rnv};

  auto vertex{segment};
  vertex = set_vertex(vertex, p00, c, a);
  vertex = set_vertex(vertex, p01, c, a);
  vertex = set_vertex(vertex, p10, c, a);
  vertex = set_vertex(vertex, p10, c, a);
  vertex = set_vertex(vertex, p01, c, a);
  vertex = set_vertex(vertex, p11, c, a);
  return vertex;
}

inline Lines::state_t iterative_set_init(Lines &self,
    sf::Vector2f const r_init = {},
    std::optional<RGB8> const c_init = {},
    std::optional<float> const a_init = {}) {
  return {get_segment(self), r_init, c_init, a_init};
}

inline Lines::state_t iterative_set_point(Lines &self, Lines::state_t const &s,
    ViewTransform const &vt, sf::Vector2f const &p0,sf::Vector2f const &p1,
    std::optional<sf::Vector2f> const &r = {},
    std::optional<RGB8> const &c = {}, std::optional<float> const &a = {}) {
  sf::Vertex * const next{set_segment(self, s.segment, vt, p0, p1,
    r.value_or(s.r), c.has_value() ? c : s.c, a.has_value() ? a : s.a)};
  return {next, s.r, s.c, s.a};
}

template<>
struct Graph<args::get_mode_index(modes_graph, "curve")> : public GraphBase {
  struct state_t {
    sf::Vertex * segment;
    sf::Vector2f p;
    sf::Vector2f r;
    std::optional<RGB8> c;
  };

  std::size_t n() const;

  Graph(std::size_t const n = 1, RGB8 const & = {});
};

using Curve = Graph<args::get_mode_index(modes_graph, "curve")>;

void set_n(Curve &, std::size_t const);

inline sf::Vertex * get_segment(Curve &self, std::size_t const i = 0) {
  return &(self.vertexes[i * 4]);
}

inline sf::Vertex * prev_segment(Curve const &, sf::Vertex * const segment) {
  return segment - 4;
}

inline sf::Vertex * next_segment(Curve const &, sf::Vertex * const segment) {
  return segment + 4;
}

inline sf::Vertex * set_segment(Curve &, sf::Vertex * const segment,
    ViewTransform const &vt,
    sf::Vector2f const &p0, sf::Vector2f const &p1,
    sf::Vector2f const &r0, sf::Vector2f const &r1,
    std::optional<RGB8> const &c0 = {},
    std::optional<RGB8> const &c1 = {}) {
  auto const _p0{vt(p0)}, _p1{vt(p1)};
  sf::Vector2f const nv{_p0.y - _p1.y, _p1.x - _p0.x};
  //sf::Vector2f const nv{_p1.y - _p0.y, _p0.x - _p1.x};
  float const length{std::sqrt(abs2(sf_vec_to_array(nv)))};
  sf::Vector2f const nnv{nv.x / length, nv.y / length};
  sf::Vector2f const r0nv{sf::Vector2f{nnv.x * r0.x, nnv.y * r0.y}};
  sf::Vector2f const r1nv{sf::Vector2f{nnv.x * r1.x, nnv.y * r1.y}};

  auto vertex{segment};
  vertex = set_vertex(vertex, _p0 + r0nv, c0);
  vertex = set_vertex(vertex, _p0 - r0nv, c0);
  vertex = set_vertex(vertex, _p1 + r1nv, c1);
  vertex = set_vertex(vertex, _p1 - r1nv, c1);
  return vertex;
}

inline Curve::state_t iterative_set_init(Curve &self,
    sf::Vector2f const r_init = {}, std::optional<RGB8> const c_init = {}) {
  return {prev_segment(self, get_segment(self)), {}, r_init, c_init};
}

inline Curve::state_t iterative_set_point(Curve &self, Curve::state_t const &s,
    ViewTransform const &vt, sf::Vector2f const &p,
    std::optional<sf::Vector2f> const &r = {},
    std::optional<RGB8> const &c = {}) {
  sf::Vertex * const next{s.segment >= get_segment(self)
    ? set_segment(self, s.segment, vt, s.p, p, s.r, r.value_or(s.r), s.c, c)
    : next_segment(self, s.segment)};
  return {next, p, r.value_or(s.r), c};
}

template<>
struct Graph<args::get_mode_index(modes_graph, "fill")> : public GraphBase {
  std::size_t n() const;

  Graph(std::size_t const = 0, sf::Color const & = {});
};

using Fill = Graph<args::get_mode_index(modes_graph, "fill")>;

void set_n(Fill &, std::size_t const);

inline sf::Vertex * get_point(Fill &self, std::size_t const i = 0) {
  return &(self.vertexes[i * 2]);
}

inline sf::Vertex * prev_point(Fill const &, sf::Vertex * const point) {
  return point - 2;
}

inline sf::Vertex * next_point(Fill const &, sf::Vertex * const point) {
  return point + 2;
}

inline sf::Vertex * set_point(Fill &, sf::Vertex * const point,
    ViewTransform const &vt,
    std::optional<sf::Vector2f> const &p_v = {},
    std::optional<sf::Vector2f> const &p_f = {},
    std::optional<RGB8> const &rgb_v = {},
    std::optional<RGB8> const &rgb_f = {},
    std::optional<float> const &a_v = {},
    std::optional<float> const &a_f = {}) {
  auto vertex{point};
  vertex = set_vertex(vertex, vt(p_f), rgb_f, a_f);
  vertex = set_vertex(vertex, vt(p_v), rgb_v, a_v);
  return vertex;
}

inline Fill::state_t iterative_set_point(Fill &self, Fill::state_t const &s,
    ViewTransform const &vt,
    sf::Vector2f const &p_v, std::optional<sf::Vector2f> const &p_f = {},
    std::optional<RGB8> const &rgb_v = {},
    std::optional<RGB8> const &rgb_f = {},
    std::optional<float> const &a_v = {},
    std::optional<float> const &a_f = {}) {
  return set_point(self, s, vt, p_v, p_f, rgb_v, rgb_f, a_v, a_f);
}

struct PlotAnnotations {
  using minor_t = float;
  using major_t = std::pair<minor_t, std::optional<std::string>>;
  template <typename T> using container_t = std::vector<T>;
  using majors_t = container_t<major_t>;
  using minors_t = container_t<minor_t>;

  enum struct TicksMode { off, on, alt };

  template <typename T> static T const empty;

  Box<majors_t> majors;
  Box<minors_t> minors;
  Box<std::optional<std::string>> labels;
  Box<TicksMode> tick_modes_major, tick_modes_minor, grid_modes_major,
    grid_modes_minor;
  RGB8 color_frame; std::uint8_t alpha_frame;
  RGB8 color_grid_major; std::uint8_t alpha_grid_major;
  RGB8 color_grid_minor; std::uint8_t alpha_grid_minor;
  RGB8 color_vignette; std::uint8_t alpha_vignette;
  float tick_length_major, tick_length_minor, thickness_frame, thickness_tick,
    thickness_grid, font_size_label, font_size_tick, vignette_factor;
  sf::Font font;
};

Box<PlotAnnotations::majors_t const &> get_major_ticks_modemapped(
    PlotAnnotations const &);
Box<PlotAnnotations::minors_t const &> get_minor_ticks_modemapped(
    PlotAnnotations const &);
Box<PlotAnnotations::majors_t const &> get_major_grids_modemapped(
    PlotAnnotations const &);
Box<PlotAnnotations::minors_t const &> get_minor_grids_modemapped(
    PlotAnnotations const &);
std::size_t n_major_ticks(PlotAnnotations const &);
std::size_t n_minor_ticks(PlotAnnotations const &);
std::size_t n_major_grids(PlotAnnotations const &);
std::size_t n_minor_grids(PlotAnnotations const &);
Box<float> get_frame_box(PlotAnnotations const &, Box<float> const &);

struct PlotAnnotationsStore {
  Lines frame;
  Lines ticks;
  Lines grid;
  GraphBase vignette;
};

void lower(PlotAnnotationsStore &, PlotAnnotations const &,
  Box<float> const &, ViewTransform const &, sf::Vector2f const &);
void draw_underlay(sf::RenderTarget &, PlotAnnotationsStore const &);
void draw_overlay(sf::RenderTarget &, PlotAnnotationsStore const &);
