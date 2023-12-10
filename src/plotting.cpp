#include <cstddef>
#include <cmath>
#include <optional>
#include <utility>
#include <type_traits>

//#include <iostream>

#include "SFML/Graphics.hpp"

#include "args.hpp"
#include "affine.hpp"
#include "plotting.hpp"

namespace {
  template <typename T>
  T const flip_sign(T const &x, bool const flip) {
    return flip ? -x : x;
  }

  template <typename T>
  ViewTransform as_view_transform(T const &scale_x, T const &scale_y,
      T const &offset_x, T const &offset_y, T const &rotation,
      bool const invert_scale, bool const invert_offset,
      bool const invert_rotation) {
    T const cos_theta{std::cos(rotation)}, sin_theta{std::sin(rotation)};
    auto const apply_scale{[](T const &x, T const &scale, bool const invert){
        return invert ? x / scale : x * scale;
      }};

    array2d<T, 2, 2> const ra{{
      {{flip_sign(apply_scale(cos_theta, scale_x, invert_scale),
          false),
        flip_sign(apply_scale(sin_theta, scale_y, invert_scale), 
          not invert_rotation)}},
      {{flip_sign(apply_scale(sin_theta, scale_x, invert_scale),
          invert_rotation),
        flip_sign(apply_scale(cos_theta, scale_y, invert_scale),
          false)}}}};
    array2d<T, 2, 1> const rb{{
      {{flip_sign(offset_x * cos_theta,
          invert_offset) +
        flip_sign(offset_y * sin_theta,
          invert_offset != (not invert_rotation))}},
      {{flip_sign(offset_x * sin_theta,
          invert_offset != invert_rotation) +
        flip_sign(offset_y * cos_theta,
          invert_offset)}}}};

    return {ra, rb};
  }

  template <typename T>
  array2d<T, 2, 2> as_rotation_matrix(T const &rotation,
      bool const invert_rotation) {
    T const cos_theta{std::cos(rotation)}, sin_theta{std::sin(rotation)};
    return {{{{cos_theta, flip_sign(sin_theta, not invert_rotation)}},
             {{flip_sign(sin_theta, invert_rotation), cos_theta}}}};
  }

}

std::pair<ViewTransform, ViewTransform> as_view_transform(
    sf::Vector2f const &scale, sf::Vector2f const &offset,
    float const rotation) {
  return {as_view_transform<float>(scale.x, scale.y, offset.x, offset.y,
      rotation, false, false, false),
    as_view_transform<float>(scale.x, scale.y, offset.x, offset.y,
      rotation, true, true, true)};
}

std::pair<array2d<float, 2, 2>, array2d<float, 2, 2>> as_rotation_matrix(
    float const rotation) {
  return {as_rotation_matrix(rotation, false),
          as_rotation_matrix(rotation, true)};
}

Box<float> get_render_target_box(sf::RenderTarget const &target) {
  sf::Vector2u const target_size{target.getSize()};
  return {0.f, 0.f, target_size.x, target_size.y};
}

Box<float> sf_view_to_box(sf::View const &view) {
  auto const view_center{view.getCenter()};
  auto const view_size_half{view.getSize() * .5f};
  return {sf_vec_to_array(view_center - view_size_half),
          sf_vec_to_array(view_center + view_size_half)};
}

void set_n(sf::VertexArray &self, std::size_t const n) {
  if (self.getVertexCount() != n) self.resize(n);
}

void set_primitive_type(sf::VertexArray &self, sf::PrimitiveType const t) {
  if (self.getPrimitiveType() != t) self.setPrimitiveType(t);
}

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

void set_n(GraphBase &self, std::size_t const n) {
  return set_n(self.vertexes, n);
}

void set_primitive_type(GraphBase &self, sf::PrimitiveType const t) {
  return set_primitive_type(self.vertexes, t);
}

void draw(sf::RenderTarget &target, GraphBase const &graph) {
  target.draw(graph.vertexes);
}

std::size_t Lines::n() const {
  return this->vertexes.getVertexCount() / 6;
}

Lines::Graph(std::size_t const n, sf::Color const &init_color) :
    GraphBase(sf::PrimitiveType::Triangles, n * 6, init_color) {}

void set_n(Lines &self, std::size_t const n) {
  return set_n(self.vertexes, n * 6);
}

std::size_t Curve::n() const {
  return this->vertexes.getVertexCount() / 4 + 1;
}

Curve::Graph(std::size_t const n, RGB8 const &init_color) :
    GraphBase(sf::PrimitiveType::TriangleStrip, (n - 1) * 4,
      rgb_to_sf_color(init_color)) {}

void set_n(Curve &self, std::size_t const n) {
  return set_n(self.vertexes, (n - 1) * 4);
}

std::size_t Fill::n() const {
  return this->vertexes.getVertexCount() / 2;
}

Fill::Graph(std::size_t const n, sf::Color const &init_color) :
    GraphBase(sf::PrimitiveType::TriangleStrip, n * 2, init_color) {}

void set_n(Fill &self, std::size_t const n) {
  return set_n(self.vertexes, n * 2);
}

template <typename T> T const PlotAnnotations::empty{};

template <typename T>
T const & get_modemapped(PlotAnnotations::TicksMode const mode,
    T const &data, T const &alternative) {
  return mode == PlotAnnotations::TicksMode::off ? PlotAnnotations::empty<T> :
    mode == PlotAnnotations::TicksMode::on ? data : alternative;
}

template <typename T>
Box<PlotAnnotations::container_t<T> const &> get_modemapped(
    Box<PlotAnnotations::TicksMode> const &modes,
    Box<PlotAnnotations::container_t<T>> const &data) {
  return {
    get_modemapped(modes.x0, data.x0, data.x1),
    get_modemapped(modes.y0, data.y0, data.y1),
    get_modemapped(modes.x1, data.x1, data.x0),
    get_modemapped(modes.y1, data.y1, data.y0)};
}

Box<PlotAnnotations::majors_t const &> get_major_ticks_modemapped(
    PlotAnnotations const &self) {
  return get_modemapped(self.tick_modes_major, self.majors);
}

Box<PlotAnnotations::minors_t const &> get_minor_ticks_modemapped(
    PlotAnnotations const &self) {
  return get_modemapped(self.tick_modes_minor, self.minors);
}

Box<PlotAnnotations::majors_t const &> get_major_grids_modemapped(
    PlotAnnotations const &self) {
  return get_modemapped(self.grid_modes_major, self.majors);
}

Box<PlotAnnotations::minors_t const &> get_minor_grids_modemapped(
    PlotAnnotations const &self) {
  return get_modemapped(self.grid_modes_minor, self.minors);
}

template <typename T>
std::size_t size_sum(Box<PlotAnnotations::container_t<T> const &> const &x) {
  return x.x0.size() + x.y0.size() + x.x1.size() + x.y1.size();
}

std::size_t n_major_ticks(PlotAnnotations const &self) {
  return size_sum(get_major_ticks_modemapped(self));
}

std::size_t n_minor_ticks(PlotAnnotations const &self) {
  return size_sum(get_minor_ticks_modemapped(self));
}

std::size_t n_major_grids(PlotAnnotations const &self) {
  return size_sum(get_major_grids_modemapped(self));
}

std::size_t n_minor_grids(PlotAnnotations const &self) {
  return size_sum(get_minor_grids_modemapped(self));
}

Box<float> get_frame_box(PlotAnnotations const &,
    Box<float> const &parent_box) {
  return lerp(Box<float>{{.1f, .1f}, {.9f, .9f}}, parent_box);
}

void lower(PlotAnnotationsStore &store, PlotAnnotations const &pa,
    Box<float> const &parent_box, ViewTransform const &vt,
    sf::Vector2f const &view_scale) {
  // TODO: Implement some of these things as methods / free functions on
  // `PlotAnnotations`?
  auto const frame_box{get_frame_box(pa, parent_box)};
  float const r_frame{pa.thickness_frame * .5f};
  std::array<float, 2> const r_frame_scaled{{r_frame / view_scale.x,
                                             r_frame / view_scale.y}};
  Box<float> const frame_corner_offsets{
    -std::get<0>(r_frame_scaled), -std::get<1>(r_frame_scaled),
     std::get<0>(r_frame_scaled),  std::get<1>(r_frame_scaled)};
  auto const frame_edges{frame_box.edges()};
  auto const frame_corner_offset_edges{frame_corner_offsets.edges()};

  // TODO: Split this up into sub-functions?

  { // Frame edges
    set_n(store.frame, 4);
    auto s{iterative_set_init(store.frame, sf::Vector2f{r_frame, r_frame},
      pa.color_frame, pa.alpha_frame)};
    auto const frame_edges_corrected{frame_box.edges(r_frame_scaled)};
    for (auto const &[p0, p1] : frame_edges_corrected)
      s = iterative_set_point(
        store.frame, s, vt, array_to_sf_vec(p0), array_to_sf_vec(p1));
  }

  { // Vignette
    set_n(store.vignette, 4 * 6);
    set_primitive_type(store.vignette, sf::PrimitiveType::Triangles);
    auto const vignette_outer_edges{add(frame_edges, Box<float>{
      -pa.vignette_factor / view_scale.x,
      -pa.vignette_factor / view_scale.y,
       pa.vignette_factor / view_scale.x,
       pa.vignette_factor / view_scale.y}.edges())};
    auto const frame_outer_edges{add(frame_edges, frame_corner_offset_edges)};
    auto s{iterative_set_init(store.vignette)};
    auto const set_next_point{[&](auto const &p, auto const &alpha){
        s = iterative_set_point(store.vignette, s, vt,
          {{std::get<0>(p), std::get<1>(p)}}, pa.color_vignette, alpha);
      }};
    for (std::size_t i{0}; i < 4; ++i) {
      auto const vignette_outer_edge{vignette_outer_edges[i]};
      auto const frame_outer_edge{frame_outer_edges[i]};
      set_next_point(std::get<0>(vignette_outer_edge), 0.f);
      set_next_point(std::get<1>(vignette_outer_edge), 0.f);
      set_next_point(std::get<0>(frame_outer_edge),    pa.alpha_vignette);
      set_next_point(std::get<0>(frame_outer_edge),    pa.alpha_vignette);
      set_next_point(std::get<1>(vignette_outer_edge), 0.f);
      set_next_point(std::get<1>(frame_outer_edge),    pa.alpha_vignette);
    }
  }

  { // Ticks, tick labels
    auto const set_next_ticks_for_edge{[&](auto &s, auto const &ticks,
          std::array<float, 2> const &range, float const r_scaled,
          float const r_frame_scaled, float const r_frame_scaled_normal,
          std::optional<float> const tick_length, float const frame_edge,
          float const frame_edge_opposite, sf::Vector2f const &normal){
        sf::Vector2f const normal_upright{std::abs(normal.x),
                                          std::abs(normal.y)};
        sf::Vector2f const axis{1.f - normal_upright.x,
                                1.f - normal_upright.y};
        std::optional<std::string> const &no_label{};
        std::optional<std::uint8_t> const alpha_enabled{};
        std::optional<std::uint8_t> const alpha_disabled{0};
        for (auto const &tick : ticks) {
          auto const &[pos, label_opt]{[&](){
              if constexpr (std::is_same_v<decltype(tick),
                  PlotAnnotations::major_t const &>) return tick;
              else return std::make_pair(tick, no_label);
            }()};
          bool const enabled{
            pos >= std::get<0>(range) + r_frame_scaled + r_scaled and
            pos <= std::get<1>(range) - r_frame_scaled - r_scaled};
          auto const &alpha{enabled ? alpha_enabled : alpha_disabled};
          sf::Vector2f p0{axis * pos}, p1{p0};
          if (tick_length.has_value()) { // we are drawing ticks
            sf::Vector2f const tick_length_scaled{
              tick_length.value() / view_scale.x,
              tick_length.value() / view_scale.y};
            p0 += normal * r_frame_scaled_normal + normal_upright * frame_edge;
            p1 = p0 + sf::Vector2f{tick_length_scaled.x * normal.x,
                                   tick_length_scaled.y * normal.y};
          } else { // we are drawing grid lines
            p0 += normal_upright * frame_edge - normal * r_frame_scaled_normal;
            p1 += normal_upright * frame_edge_opposite +
              normal * r_frame_scaled_normal;
          }
          s = iterative_set_point(store.ticks, s, vt, p0, p1, {}, {}, alpha);
        }
      }};
    Box<sf::Vector2f> const normals{
      {0.f, 1.f}, {1.f, 0.f}, {0.f, -1.f}, {-1.f, 0.f}};
    sf::Vector2f const r_grid{pa.thickness_grid * .5f, pa.thickness_grid * .5f};
    sf::Vector2f const r_grid_scaled{r_grid.x / view_scale.x,
                                     r_grid.y / view_scale.y};
    auto const major_grids{get_major_grids_modemapped(pa)};
    auto const minor_grids{get_minor_grids_modemapped(pa)};
    auto const n_grid_lines{size_sum(major_grids) + size_sum(minor_grids)};
    set_n(store.grid, n_grid_lines);

    sf::Vector2f const r_tick{pa.thickness_tick * .5f, pa.thickness_tick * .5f};
    sf::Vector2f const r_tick_scaled{r_tick.x / view_scale.x,
                                     r_tick.y / view_scale.y};
    auto const major_ticks{get_major_ticks_modemapped(pa)};
    auto const minor_ticks{get_minor_ticks_modemapped(pa)};
    auto const n_ticks{size_sum(major_ticks) + size_sum(minor_ticks)};
    set_n(store.ticks, n_ticks);

    auto s_grid{iterative_set_init(store.grid, r_grid, pa.color_grid_minor,
        pa.alpha_grid_minor)};
    set_next_ticks_for_edge(s_grid, minor_grids.x0, frame_box.xs(),
      r_grid_scaled.x, std::get<0>(r_frame_scaled), std::get<1>(r_frame_scaled),
      {}, frame_box.y0, frame_box.y1, normals.x0);
    set_next_ticks_for_edge(s_grid, minor_grids.y0, frame_box.ys(),
      r_grid_scaled.y, std::get<1>(r_frame_scaled), std::get<0>(r_frame_scaled),
      {}, frame_box.x0, frame_box.x1, normals.y0);
    set_next_ticks_for_edge(s_grid, minor_grids.x1, frame_box.xs(),
      r_grid_scaled.x, std::get<0>(r_frame_scaled), std::get<1>(r_frame_scaled),
      {}, frame_box.y1, frame_box.y0, normals.x1);
    set_next_ticks_for_edge(s_grid, minor_grids.y1, frame_box.ys(),
      r_grid_scaled.y, std::get<1>(r_frame_scaled), std::get<0>(r_frame_scaled),
      {}, frame_box.x1, frame_box.x0, normals.y1);
    s_grid.c = pa.color_grid_major; s_grid.a = pa.alpha_grid_major;
    set_next_ticks_for_edge(s_grid, major_grids.x0, frame_box.xs(),
      r_grid_scaled.x, std::get<0>(r_frame_scaled), std::get<1>(r_frame_scaled),
      {}, frame_box.y0, frame_box.y1, normals.x0);
    set_next_ticks_for_edge(s_grid, major_grids.y0, frame_box.ys(),
      r_grid_scaled.y, std::get<1>(r_frame_scaled), std::get<0>(r_frame_scaled),
      {}, frame_box.x0, frame_box.x1, normals.y0);
    set_next_ticks_for_edge(s_grid, major_grids.x1, frame_box.xs(),
      r_grid_scaled.x, std::get<0>(r_frame_scaled), std::get<1>(r_frame_scaled),
      {}, frame_box.y1, frame_box.y0, normals.x1);
    set_next_ticks_for_edge(s_grid, major_grids.y1, frame_box.ys(),
      r_grid_scaled.y, std::get<1>(r_frame_scaled), std::get<0>(r_frame_scaled),
      {}, frame_box.x1, frame_box.x0, normals.y1);

    auto s_ticks{iterative_set_init(store.ticks, r_tick, pa.color_frame,
        pa.alpha_frame)};
    set_next_ticks_for_edge(s_ticks, minor_ticks.x0, frame_box.xs(),
      r_tick_scaled.x, std::get<0>(r_frame_scaled), std::get<1>(r_frame_scaled),
      pa.tick_length_minor, frame_box.y0, frame_box.y1, normals.x0);
    set_next_ticks_for_edge(s_ticks, minor_ticks.y0, frame_box.ys(),
      r_tick_scaled.y, std::get<1>(r_frame_scaled), std::get<0>(r_frame_scaled),
      pa.tick_length_minor, frame_box.x0, frame_box.x1, normals.y0);
    set_next_ticks_for_edge(s_ticks, minor_ticks.x1, frame_box.xs(),
      r_tick_scaled.x, std::get<0>(r_frame_scaled), std::get<1>(r_frame_scaled),
      pa.tick_length_minor, frame_box.y1, frame_box.y0, normals.x1);
    set_next_ticks_for_edge(s_ticks, minor_ticks.y1, frame_box.ys(),
      r_tick_scaled.y, std::get<1>(r_frame_scaled), std::get<0>(r_frame_scaled),
      pa.tick_length_minor, frame_box.x1, frame_box.x0, normals.y1);
    set_next_ticks_for_edge(s_ticks, major_ticks.x0, frame_box.xs(),
      r_tick_scaled.x, std::get<0>(r_frame_scaled), std::get<1>(r_frame_scaled),
      pa.tick_length_major, frame_box.y0, frame_box.y1, normals.x0);
    set_next_ticks_for_edge(s_ticks, major_ticks.y0, frame_box.ys(),
      r_tick_scaled.y, std::get<1>(r_frame_scaled), std::get<0>(r_frame_scaled),
      pa.tick_length_major, frame_box.x0, frame_box.x1, normals.y0);
    set_next_ticks_for_edge(s_ticks, major_ticks.x1, frame_box.xs(),
      r_tick_scaled.x, std::get<0>(r_frame_scaled), std::get<1>(r_frame_scaled),
      pa.tick_length_major, frame_box.y1, frame_box.y0, normals.x1);
    set_next_ticks_for_edge(s_ticks, major_ticks.y1, frame_box.ys(),
      r_tick_scaled.y, std::get<1>(r_frame_scaled), std::get<0>(r_frame_scaled),
      pa.tick_length_major, frame_box.x1, frame_box.x0, normals.y1);
  }
}

void draw_underlay(sf::RenderTarget &target,
    PlotAnnotationsStore const &store) {
  draw(target, store.grid);
}

void draw_overlay(sf::RenderTarget &target, PlotAnnotationsStore const &store) {
  draw(target, store.vignette);
  draw(target, store.ticks);
  draw(target, store.frame);
}
