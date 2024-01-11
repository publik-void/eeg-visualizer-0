#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <memory>
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
Box<T> sf_rect_to_box(sf::Rect<T> const &rect) {
  return {rect.left, rect.top, rect.left + rect.width, rect.top + rect.height};
}

// For debugging
bool constexpr draw_bounding_boxes{false};
void draw_bounding_box(sf::RenderTarget &, Box<float> const &,
    ViewTransform const &vt = ViewTransform::identity(),
    sf::Color const & = {0x00, 0xff, 0xff, 0x7f},
    sf::Color const & = {0xff, 0x00, 0xff, 0x7f});
void draw_bounding_box(sf::RenderTarget &, sf::Text const &,
    ViewTransform const &vt = ViewTransform::identity(),
    sf::Color const & = {0x00, 0xff, 0xff, 0x7f},
    sf::Color const & = {0xff, 0x00, 0xff, 0x7f});

template <typename T>
struct RGB {
  T r;
  T g;
  T b;
};

using RGB8 = RGB<std::uint8_t>;

inline sf::Color rgb_to_sf_color(RGB8 const &c, std::uint8_t const a = 0xff) {
  return {c.r, c.g, c.b, a};
}
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

// NOTE: Subclasses of this are meant to have a unified interface of free
// functions. However, as long as it is not abstract, the value of using
// inheritance like this is probably questionable.
struct GraphInterface {};

struct GraphBase : public GraphInterface {
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

  Graph(std::size_t const = 1, RGB8 const & = {});
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

struct FillCurve : public GraphInterface {
  struct state_t {
    Fill::state_t s_fill;
    Curve::state_t s_curve;
  };

  Fill fill;
  Curve curve;

  FillCurve(std::size_t const = 1, sf::Color const & = {}, RGB8 const & = {});
};

void set_n(FillCurve &, std::size_t const);

void draw(sf::RenderTarget &, FillCurve const &);

inline FillCurve::state_t iterative_set_init(FillCurve &self,
    sf::Vector2f const r_init = {}, std::optional<RGB8> const c_init = {}) {
  return {iterative_set_init(self.fill),
    iterative_set_init(self.curve, r_init, c_init)};
}

inline FillCurve::state_t iterative_set_point(FillCurve &self,
    FillCurve::state_t const &s, ViewTransform const &vt,
    sf::Vector2f const &p_v, sf::Vector2f const &p_f = {},
    std::optional<RGB8> const &rgb_v = {},
    std::optional<RGB8> const &rgb_f = {},
    std::optional<float> const &a_v = {},
    std::optional<float> const &a_f = {},
    std::optional<sf::Vector2f> const &r = {},
    std::optional<RGB8> const &rgb_c = {}) {
  return {iterative_set_point(self.fill, s.s_fill, vt, p_v, p_f, rgb_v, rgb_f,
      a_v, a_f),
    iterative_set_point(self.curve, s.s_curve, vt, p_v, r, rgb_c)};
}

struct PlotAnnotations {
  using label_t = std::optional<std::string>;
  using minor_t = float;
  using major_t = std::pair<minor_t, label_t>;
  template <typename T> using container_t = std::vector<T>;
  using majors_t = container_t<major_t>;
  using minors_t = container_t<minor_t>;

  using mode_t = std::uint8_t;
  static mode_t constexpr off        {0b00000000};
  static mode_t constexpr on         {0b00000001};
  static mode_t constexpr alt        {0b00000010};
  static mode_t constexpr nolabel    {0b00000100};
  static mode_t constexpr bold       {0b00001000};
  static mode_t constexpr italic     {0b00010000};
  static mode_t constexpr underlined {0b00100000};
  static mode_t constexpr vertical   {0b01000000};

  static sf::Uint32 mode_to_sf_text_style(mode_t const);

  template <typename T> static T const empty;

  Box<majors_t> majors;
  Box<minors_t> minors;
  Box<label_t> frame_labels;
  Box<mode_t> modes_ticks_major, modes_ticks_minor, modes_grid_major,
    modes_grid_minor, modes_frame_labels;
  RGB8 color_frame; std::uint8_t alpha_frame;
  RGB8 color_grid_major; std::uint8_t alpha_grid_major;
  RGB8 color_grid_minor; std::uint8_t alpha_grid_minor;
  RGB8 color_frame_labels; std::uint8_t alpha_frame_labels;
  RGB8 color_ticks_labels; std::uint8_t alpha_ticks_labels;
  RGB8 color_vignette; std::uint8_t alpha_vignette;
  float tick_length_major, tick_length_minor, thickness_frame, thickness_ticks,
    thickness_grid, font_size_frame_labels, font_size_ticks_labels,
    gap_ticks_labels, gap_frame_labels, gap_parent_box, vignette_factor;
  std::shared_ptr<sf::Font> font;
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
std::size_t n_frame_labels(PlotAnnotations const &);
std::size_t n_ticks_labels(PlotAnnotations const &);
Box<float> get_frame_box(PlotAnnotations const &, Box<float> const &,
    sf::Vector2f const &);

struct PlotAnnotationsStore {
  Lines frame;
  Lines ticks;
  Lines grid;
  GraphBase vignette;
  std::vector<sf::Text> labels_frame;
  std::vector<sf::Text> labels_ticks;
};

void lower(PlotAnnotationsStore &, PlotAnnotations const &,
  Box<float> const &, ViewTransform const &, sf::Vector2f const &,
  float const);
void draw_underlay(sf::RenderTarget &, PlotAnnotationsStore const &);
void draw_overlay(sf::RenderTarget &, PlotAnnotationsStore const &);
