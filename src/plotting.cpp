#include <cstddef>
#include <cmath>
#include <optional>
#include <utility>
#include <numeric>
#include <type_traits>

//#include <iostream>

#include "SFML/Graphics.hpp"

#include "args.hpp"
#include "numbers.hpp"
#include "affine.hpp"
#include "plotting.hpp"

namespace {
  template <typename T>
  T constexpr flip_sign(T const &x, bool const flip) {
    return flip ? -x : x;
  }

  template <typename T>
  ViewTransform constexpr as_view_transform(T const &scale_x, T const &scale_y,
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
  array2d<T, 2, 2> constexpr as_rotation_matrix(T const &rotation,
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

std::pair<array2d<float, 2, 2>, array2d<float, 2, 2>>
    as_rotation_matrix(float const rotation) {
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

// For debugging
void draw_bounding_box(sf::RenderTarget &target, Box<float> const &box,
    ViewTransform const &vt,
    sf::Color const &color_outline, sf::Color const &color_cross) {
  sf::VertexArray outline(sf::PrimitiveType::LineStrip, 5u);
  sf::VertexArray cross(sf::PrimitiveType::Lines, 4u);
  auto p00{array_to_sf_vec(vt(box.p00()))},
       p01{array_to_sf_vec(vt(box.p01()))},
       p10{array_to_sf_vec(vt(box.p10()))},
       p11{array_to_sf_vec(vt(box.p11()))};
  outline[0] = {p00, color_outline};
  outline[1] = {p01, color_outline};
  outline[2] = {p11, color_outline};
  outline[3] = {p10, color_outline};
  outline[4] = {p00, color_outline};
  cross[0] = {(p00 + p01) * .5f, color_cross};
  cross[1] = {(p10 + p11) * .5f, color_cross};
  cross[2] = {(p00 + p10) * .5f, color_cross};
  cross[3] = {(p01 + p11) * .5f, color_cross};
  target.draw(cross);
  target.draw(outline);
}
void draw_bounding_box(sf::RenderTarget &target, sf::Text const &text,
    ViewTransform const &vt,
    sf::Color const &color_outline, sf::Color const &color_cross) {
  return draw_bounding_box(target, sf_rect_to_box(text.getGlobalBounds()), vt,
    color_outline, color_cross);
}

template <typename T>
void set_n(std::vector<T> &self, std::size_t const n) {
  if (self.size() != n) self.resize(n);
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

void draw(sf::RenderTarget &target, sf::Drawable const &x) {
  target.draw(x);
}

void draw(sf::RenderTarget &target, GraphBase const &graph) {
  draw(target, graph.vertexes);
}

template <typename T>
void draw(sf::RenderTarget &target, std::vector<T> const &xs) {
  for (T const &x : xs) draw(target, x);
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

sf::Uint32 PlotAnnotations::mode_to_sf_text_style(mode_t const mode) {
  sf::Uint32 style{sf::Text::Style::Regular};
  if (mode & PlotAnnotations::bold) style |= sf::Text::Style::Bold;
  if (mode & PlotAnnotations::italic) style |= sf::Text::Style::Italic;
  if (mode & PlotAnnotations::underlined) style |= sf::Text::Style::Underlined;
  return style;
}

template <typename T> T const PlotAnnotations::empty{};

template <typename T>
T const & get_modemapped(PlotAnnotations::mode_t const mode,
    T const &data, T const &alternative) {
  return (mode & PlotAnnotations::on) ? (mode & PlotAnnotations::alt) ?
    alternative : data : PlotAnnotations::empty<T>;
}

template <typename T>
Box<T const &> get_modemapped(Box<PlotAnnotations::mode_t> const &modes,
    Box<T> const &data) {
  return {
    get_modemapped(modes.x0, data.x0, data.x1),
    get_modemapped(modes.y0, data.y0, data.y1),
    get_modemapped(modes.x1, data.x1, data.x0),
    get_modemapped(modes.y1, data.y1, data.y0)};
}

Box<PlotAnnotations::majors_t const &> get_major_ticks_modemapped(
    PlotAnnotations const &self) {
  return get_modemapped(self.modes_ticks_major, self.majors);
}

Box<PlotAnnotations::minors_t const &> get_minor_ticks_modemapped(
    PlotAnnotations const &self) {
  return get_modemapped(self.modes_ticks_minor, self.minors);
}

Box<PlotAnnotations::majors_t const &> get_major_grids_modemapped(
    PlotAnnotations const &self) {
  return get_modemapped(self.modes_grid_major, self.majors);
}

Box<PlotAnnotations::minors_t const &> get_minor_grids_modemapped(
    PlotAnnotations const &self) {
  return get_modemapped(self.modes_grid_minor, self.minors);
}

Box<PlotAnnotations::label_t const &> get_frame_labels_modemapped(
    PlotAnnotations const &self) {
  return get_modemapped(self.modes_frame_labels, self.frame_labels);
}

template <typename T>
std::size_t size_sum(Box<PlotAnnotations::container_t<T> const &> const &x) {
  return x.x0.size() + x.y0.size() + x.x1.size() + x.y1.size();
}

template <typename T>
bool has_value(std::optional<T> const &x) {
  return x.has_value();
}

bool has_label(PlotAnnotations::major_t const &x) {
  return has_value(x.second);
}

template <typename T, typename F>
std::size_t size_sum(T const &x, F const &&predicate) {
  return std::size_t{predicate(x)};
}

template <typename T, typename F>
std::size_t size_sum(PlotAnnotations::container_t<T> const &xs,
    F const &&predicate) {
  return std::accumulate(xs.cbegin(), xs.cend(), std::size_t{0u},
    [&](std::size_t const acc, T const &x){
        return acc + size_sum(x, predicate);
      });
}

template <typename T, typename F0>
std::size_t size_sum(Box<T const &> const &x, F0 const &&predicate,
    Box<bool> const &edge_toggles = {true, true, true, true}) {
  std::size_t constexpr zero{0u};
  return (edge_toggles.x0 ? size_sum(x.x0, predicate) : zero) +
         (edge_toggles.y0 ? size_sum(x.y0, predicate) : zero) +
         (edge_toggles.x1 ? size_sum(x.x1, predicate) : zero) +
         (edge_toggles.y1 ? size_sum(x.y1, predicate) : zero);
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

std::size_t n_labels_frame(PlotAnnotations const &self) {
  return size_sum(get_frame_labels_modemapped(self),
    has_value<PlotAnnotations::label_t::value_type>);
}

std::size_t n_labels_ticks(PlotAnnotations const &self) {
  return size_sum(get_major_ticks_modemapped(self), has_label, {
    not (self.modes_ticks_major.x0 & PlotAnnotations::nolabel),
    not (self.modes_ticks_major.y0 & PlotAnnotations::nolabel),
    not (self.modes_ticks_major.x1 & PlotAnnotations::nolabel),
    not (self.modes_ticks_major.y1 & PlotAnnotations::nolabel)});
}

Box<sf::Vector2f> const normals{
  {0.f, 1.f}, {1.f, 0.f}, {0.f, -1.f}, {-1.f, 0.f}};

float get_frame_edge(PlotAnnotations::mode_t const mode_ticks_major,
    PlotAnnotations::mode_t const mode_frame_label,
    float const font_size_frame_labels,
    float const font_size_ticks_labels, float const gap_ticks_labels,
    float const gap_frame_labels, float const gap_parent_edge,
    float const normal, float const parent_edge) {
  return parent_edge + normal * (gap_parent_edge +
      (gap_frame_labels + font_size_frame_labels) *
      (mode_frame_label & PlotAnnotations::on) +
      (gap_ticks_labels + font_size_ticks_labels) *
      (not (mode_ticks_major & PlotAnnotations::nolabel)));
}

float get_frame_edge(PlotAnnotations const &self,
    PlotAnnotations::mode_t const mode_ticks_major,
    PlotAnnotations::mode_t const mode_frame_label,
    float const normal, float const parent_edge) {
  return get_frame_edge(mode_ticks_major, mode_frame_label,
    self.font_size_frame_labels, self.font_size_ticks_labels,
    self.gap_ticks_labels, self.gap_frame_labels, self.gap_parent_box,
    normal, parent_edge);
}

Box<float> get_frame_box(PlotAnnotations const &self,
    Box<float> const &parent_box, sf::Vector2f const &view_scale) {
  return Box<float>{
    get_frame_edge(self, self.modes_ticks_major.y0, self.modes_frame_labels.y0,
      normals.y0.x / view_scale.x, parent_box.x0),
    get_frame_edge(self, self.modes_ticks_major.x0, self.modes_frame_labels.x0,
      normals.x0.y / view_scale.y, parent_box.y0),
    get_frame_edge(self, self.modes_ticks_major.y1, self.modes_frame_labels.y1,
      normals.y1.x / view_scale.x, parent_box.x1),
    get_frame_edge(self, self.modes_ticks_major.x1, self.modes_frame_labels.x1,
      normals.x1.y / view_scale.y, parent_box.y1)};
}

void lower(PlotAnnotationsStore &store, PlotAnnotations const &pa,
    Box<float> const &parent_box, ViewTransform const &vt,
    sf::Vector2f const &view_scale, float const view_rotation) {
  // TODO: Implement some of these things as methods / free functions on
  // `PlotAnnotations`?
  auto const frame_box{get_frame_box(pa, parent_box, view_scale)};
  float const r_frame{pa.thickness_frame * .5f};
  std::array<float, 2> const r_frame_scaled{{r_frame / view_scale.x,
                                             r_frame / view_scale.y}};
  Box<float> const frame_corner_offsets{
    -std::get<0>(r_frame_scaled), -std::get<1>(r_frame_scaled),
     std::get<0>(r_frame_scaled),  std::get<1>(r_frame_scaled)};
  auto const frame_edges{frame_box.edges()};
  auto const frame_corner_offset_edges{frame_corner_offsets.edges()};

  auto const set_next_label{[&](std::vector<sf::Text>::iterator const &itr,
        PlotAnnotations::label_t const &label, unsigned int const size,
        PlotAnnotations::mode_t const mode,
        RGB8 const &color, std::uint8_t const alpha,
        sf::Vector2f const &p, sf::Vector2f const &normal){
      if (not (label.has_value() and (mode & PlotAnnotations::on)))
        return itr;

      float constexpr vertical_rotation{pi_halves<float>};
      array2d<float, 2, 2> constexpr vertical_rotation_matrix{
        {{{0.f, 1.f}}, {{-1.f, 0.f}}}};

      bool const vertical(mode & PlotAnnotations::vertical);
      auto const _normal{sf_vec_to_array2d(normal)};
      auto const _normal_rotated{vertical
        ? muladd(vertical_rotation_matrix, _normal) : _normal};

      sf::Text &text{*itr};
      text.setString(label.value());
      text.setFont(*(pa.font));
      text.setCharacterSize(size);
      text.setFillColor(rgb_to_sf_color(color, alpha));
      text.setStyle(PlotAnnotations::mode_to_sf_text_style(mode));

      float rotation{view_rotation};
      if (vertical) rotation += vertical_rotation;
      bool const flip{rotation  > pi_halves<float> and
                      rotation <= three_pi_halves<float>};
      if (flip) rotation += pi<float>;
      Box<float> const from{-1.f, flip ? -1.f : 1.f, 1.f, flip ? 1.f : -1.f};
      auto const bounds{sf_rect_to_box(text.getLocalBounds())};
      text.setOrigin(
        array2d_to_sf_vec(lerp(_normal_rotated, from, bounds)));

      text.setPosition(vt(p));
      text.setScale(1.f, -1.f);
      text.setRotation(rotation * three_hundred_and_sixty_over_two_pi<float>);

      return std::next(itr);
    }};

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

  { // Frame labels
    auto const frame_labels{get_frame_labels_modemapped(pa)};
    set_n(store.labels_frame, n_labels_frame(pa));

    auto itr{store.labels_frame.begin()};
    auto const font_size{static_cast<unsigned int>(pa.font_size_frame_labels)};

    auto const set_next_frame_label{[&](auto const &itr,
          PlotAnnotations::label_t const &label,
          PlotAnnotations::mode_t const mode_frame_label,
          PlotAnnotations::mode_t const mode_ticks_major,
          std::array<float, 2> const &p0, std::array<float, 2> const &p1,
          sf::Vector2f const &normal, float const view_scale){
        sf::Vector2f const p{(array_to_sf_vec(p0) + array_to_sf_vec(p1)) * .5f};
        float const offset{(pa.gap_frame_labels +
          (pa.gap_ticks_labels + pa.font_size_ticks_labels) *
          (not (mode_ticks_major & PlotAnnotations::nolabel))) / view_scale};
        return set_next_label(itr, label, font_size, mode_frame_label,
          pa.color_frame_labels, pa.alpha_frame_labels, p - normal * offset,
          normal);
      }};

    itr = set_next_frame_label(itr, frame_labels.x0, pa.modes_frame_labels.x0,
      pa.modes_ticks_major.x0, frame_box.p00(), frame_box.p01(), normals.x0,
      view_scale.y);
    itr = set_next_frame_label(itr, frame_labels.y0, pa.modes_frame_labels.y0,
      pa.modes_ticks_major.y0, frame_box.p00(), frame_box.p10(), normals.y0,
      view_scale.x);
    itr = set_next_frame_label(itr, frame_labels.x1, pa.modes_frame_labels.x1,
      pa.modes_ticks_major.x1, frame_box.p10(), frame_box.p11(), normals.x1,
      view_scale.y);
    itr = set_next_frame_label(itr, frame_labels.y1, pa.modes_frame_labels.y1,
      pa.modes_ticks_major.y1, frame_box.p01(), frame_box.p11(), normals.y1,
      view_scale.x);
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
    set_n(store.labels_ticks, n_labels_ticks(pa));
    auto itr_labels{store.labels_ticks.begin()};

    auto const font_size{static_cast<unsigned int>(pa.font_size_ticks_labels)};

    auto const set_next_ticks_for_edge{[&](auto &s, auto const &ticks,
          std::array<float, 2> const &range, float const r_scaled,
          float const r_frame_scaled, float const r_frame_scaled_normal,
          std::optional<float> const tick_length, float const frame_edge,
          float const frame_edge_opposite, sf::Vector2f const &normal,
          PlotAnnotations::mode_t const mode_ticks_major = 0){
        sf::Vector2f const normal_upright{std::abs(normal.x),
                                          std::abs(normal.y)};
        sf::Vector2f normal_scaled{normal.x / view_scale.x,
                                   normal.y / view_scale.y};
        sf::Vector2f const axis{1.f - normal_upright.x,
                                1.f - normal_upright.y};
        std::optional<std::string> const &no_label{};
        std::optional<std::uint8_t> constexpr alpha_enabled{};
        std::optional<std::uint8_t> constexpr alpha_disabled{0};
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
            sf::Vector2f const p{p0 + normal_upright * frame_edge};
            p0 = p + normal * r_frame_scaled_normal;
            // TODO: I think I can use the scaled normal here
            p1 = p0 + sf::Vector2f{tick_length_scaled.x * normal.x,
                                   tick_length_scaled.y * normal.y};

            if (label_opt.has_value() and
                (not (mode_ticks_major & PlotAnnotations::nolabel))) {
              itr_labels = set_next_label(itr_labels, label_opt, font_size,
                mode_ticks_major, pa.color_ticks_labels,
                enabled ? pa.alpha_ticks_labels : 0x00,
                p - normal_scaled * pa.gap_ticks_labels, normal);
            }
          } else { // we are drawing grid lines
            p0 += normal_upright * frame_edge - normal * r_frame_scaled_normal;
            p1 += normal_upright * frame_edge_opposite +
              normal * r_frame_scaled_normal;
          }
          s = iterative_set_point(store.ticks, s, vt, p0, p1, {}, {}, alpha);
        }
      }};

    sf::Vector2f const r_grid{pa.thickness_grid * .5f, pa.thickness_grid * .5f};
    sf::Vector2f const r_grid_scaled{r_grid.x / view_scale.x,
                                     r_grid.y / view_scale.y};
    auto const major_grids{get_major_grids_modemapped(pa)};
    auto const minor_grids{get_minor_grids_modemapped(pa)};
    auto const n_grid_lines{size_sum(major_grids) + size_sum(minor_grids)};
    set_n(store.grid, n_grid_lines);

    sf::Vector2f const r_tick{pa.thickness_ticks * .5f,
                              pa.thickness_ticks * .5f};
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
      pa.tick_length_major, frame_box.y0, frame_box.y1, normals.x0,
      pa.modes_ticks_major.x0);
    set_next_ticks_for_edge(s_ticks, major_ticks.y0, frame_box.ys(),
      r_tick_scaled.y, std::get<1>(r_frame_scaled), std::get<0>(r_frame_scaled),
      pa.tick_length_major, frame_box.x0, frame_box.x1, normals.y0,
      pa.modes_ticks_major.y0);
    set_next_ticks_for_edge(s_ticks, major_ticks.x1, frame_box.xs(),
      r_tick_scaled.x, std::get<0>(r_frame_scaled), std::get<1>(r_frame_scaled),
      pa.tick_length_major, frame_box.y1, frame_box.y0, normals.x1,
      pa.modes_ticks_major.x1);
    set_next_ticks_for_edge(s_ticks, major_ticks.y1, frame_box.ys(),
      r_tick_scaled.y, std::get<1>(r_frame_scaled), std::get<0>(r_frame_scaled),
      pa.tick_length_major, frame_box.x1, frame_box.x0, normals.y1,
      pa.modes_ticks_major.y1);
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
  draw(target, store.labels_frame);
  draw(target, store.labels_ticks);

  if constexpr (draw_bounding_boxes) {
    for (auto const &text : store.labels_frame)
      draw_bounding_box(target, text);
    for (auto const &text : store.labels_ticks)
      draw_bounding_box(target, text);
  }
}
