#include <cstddef>
#include <cmath>
#include <complex>
#include <utility>
#include <memory>
#include <optional>
#include <array>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <atomic>
#include <chrono>

//#include <iostream>

#include "SFML/Graphics.hpp"
#include "SFML/Window.hpp"
#include "SFML/Window/Mouse.hpp"
#include "SFML/Window/Keyboard.hpp"
#include "fftw3.h"

#include "signal-handling.hpp"
#include "args.hpp"
#include "dsp.hpp"
#include "lsl-customized.hpp"
#include "numbers-and-math.hpp"

args::desc_t constexpr desc_visualizer{
  "Visualize incoming EEG data in real time."};

args::desc_t constexpr desc_time_series{
  "Show the latest data up to a certain duration ago."};

args::desc_t constexpr desc_short_time_spectrum{
  "Compute and show the short-time Fourier transform of incoming EEG data."};

args::modes_t<4> constexpr modes_visualizer{"error", "help", "time-series",
  "short-time-spectrum"};

std::array<args::desc_t, 4> constexpr modes_descs_visualizer{nullptr, nullptr,
  desc_time_series, desc_short_time_spectrum};


struct Buffer;
inline std::size_t get_scalar_index(Buffer const &,
    std::array<std::size_t, 2> const);

struct Buffer {
  lsl::stream_info stream;
  lsl::stream_inlet inlet;
  std::size_t n_timepoints;
  std::vector<float> data;
  std::size_t n_pulled;
  std::atomic_size_t i_in;
  std::atomic_size_t i_out;

  using scalar_iterator = decltype(data)::iterator;
  using sample = std::pair<scalar_iterator, scalar_iterator>;

  double sample_rate() const {
    return this->stream.nominal_srate();
  }

  std::size_t n_channels() const {
    return static_cast<size_t>(this->inlet.get_channel_count());
  }

  Buffer(lsl::attachment_t & attachment, double const buffer_duration) :
      stream{std::move(lsl::get_stream(attachment))},
      inlet{std::move(lsl::get_inlet(attachment))},
      n_timepoints{static_cast<std::size_t>(
        std::round(this->sample_rate() * buffer_duration))},
      data{std::vector<float>(this->n_channels() * this->n_timepoints, 0.f)},
      n_pulled{0}, i_in{0}, i_out{this->n_timepoints - 1} {}

  inline float & operator[](std::array<std::size_t, 2> const ji) {
    return this->data[get_scalar_index(*this, ji)];
  }

  inline float const & operator[](std::array<std::size_t, 2> const ji) const {
    return this->data[get_scalar_index(*this, ji)];
  }
};

std::size_t pull(Buffer &self) {
  float * const buf{self.data.data()};
  std::size_t const size{self.data.size()};
  std::size_t const k{self.n_channels()};
  std::size_t i_in_k{self.i_in.load() * k};

  std::size_t n_pulled_k{self.inlet.pull_chunk_multiplexed(
    buf + i_in_k, nullptr, size - i_in_k, 0)};
  i_in_k += n_pulled_k;
  if (i_in_k >= size) {
    i_in_k = self.inlet.pull_chunk_multiplexed(buf, nullptr, size, 0);
    n_pulled_k += i_in_k;
  }

  std::size_t const i_in{i_in_k / k};
  self.n_pulled = n_pulled_k / k;
  self.i_in = i_in;
  self.i_out = i_in == 0 ? self.n_timepoints - 1 : i_in - 1;
  return self.n_pulled;
}

inline std::size_t get_scalar_index(Buffer const &buffer, std::size_t const i) {
  // NOTE: `std::size_t` is guaranteed to be unsigned
  auto const t{buffer.i_out - i};
  auto const t_wrapped{t >= (std::size_t{0} - buffer.n_timepoints)
    ? t + buffer.n_timepoints : t};
  return t_wrapped * buffer.n_channels();
}

inline std::size_t get_scalar_index(Buffer const &buffer,
    std::array<std::size_t, 2> const ji) {
  auto const [j, i]{ji};
  return get_scalar_index(buffer, i) + j;
}

inline Buffer::sample get_sample(Buffer &buffer, std::size_t const i) {
  auto const first{std::begin(buffer.data) + get_scalar_index(buffer, i)};
  auto const last{first + buffer.n_channels()};
  return {first, last};
}

inline Buffer::sample get_sample(Buffer &buffer) {
  return get_sample(buffer, buffer.n_pulled - 1);
}

inline Buffer::sample next_sample(Buffer &buffer,
    Buffer::scalar_iterator const &last) {
  return last < std::end(buffer.data)
    ? std::make_pair(last, last + buffer.n_channels())
    : std::make_pair(std::begin(buffer.data),
      std::begin(buffer.data) + buffer.n_channels());
}

inline Buffer::sample next_sample(Buffer &buffer,
    Buffer::sample const &sample) {
  return next_sample(buffer, sample.second);
}

inline void advance_sample(Buffer &buffer, Buffer::sample &sample) {
  sample = next_sample(buffer, sample);
}

inline bool past_newest_sample(Buffer &buffer,
    Buffer::scalar_iterator const &first) {
  auto const newest_sample{get_sample(buffer, 0)};
  return newest_sample.second == first or
    (newest_sample.second == std::end(buffer.data) and
      first == std::begin(buffer.data));
}

inline bool past_newest_sample(Buffer &buffer,
    Buffer::sample const &sample) {
  return past_newest_sample(buffer, sample.first);
}

struct GraphBase {
  sf::VertexArray vertexes;
  std::size_t n() const { return this->vertexes.getVertexCount(); }

  sf::Vertex * begin() {
    return &(this->vertexes[0]);
  }
  sf::Vertex * end() {
    return &(this->vertexes[vertexes.getVertexCount() - 1]) + 1;
  }

  using state_t = sf::Vertex *;

  GraphBase(sf::VertexArray const vertexes) : vertexes{vertexes} {}
  GraphBase(sf::PrimitiveType const primitive, std::size_t const n_vertexes,
      sf::Color const &init_color) : vertexes{[&](){
        sf::VertexArray vertexes(primitive, n_vertexes);
        for (std::size_t i{0}; i < vertexes.getVertexCount(); ++i)
          vertexes[i].color = init_color;
        return vertexes;
      }()} {}
  virtual ~GraphBase() {}
};

inline GraphBase::state_t iterative_set_init(GraphBase &self) {
  return self.begin();
}

inline GraphBase::state_t iterative_set_point(GraphBase &, sf::Vector2f const p,
    GraphBase::state_t const &s) {
  s->position = p;
  return s + 1;
}

inline GraphBase::state_t iterative_set_point(GraphBase &self,
    sf::Vector2f const p, sf::Color const c, GraphBase::state_t const &s) {
  s->color = c;
  return iterative_set_point(self, p, s);
}

void draw(sf::RenderTarget &target, GraphBase const &graph) {
  target.draw(graph.vertexes);
}

sf::Vector2f thickness_px_to_2d(float const thickness,
    auto/*sf::Window or sf::RenderTarget*/ const &target) {
  auto size{target.getSize()};
  return {thickness / size.x, thickness / size.y};
}

args::modes_t<1> constexpr modes_graph{"curve"};

template<std::size_t mode>
struct Graph : public GraphBase {};

template<>
struct Graph<args::get_mode_index(modes_graph, "curve")> : public GraphBase {
  sf::Vector2f thickness_2d;

  std::size_t n() const {
    return this->vertexes.getVertexCount() / 4 + 1;
  }

  struct state_t {
    sf::Vertex * segment;
    sf::Vector2f p;
    sf::Color c;
  };

  Graph(std::size_t const n, sf::Color const &init_color,
      sf::Vector2f const thickness_2d) :
      GraphBase(sf::PrimitiveType::TriangleStrip, (n - 1) * 4, init_color),
      thickness_2d{thickness_2d} {}
};

using Curve = Graph<args::get_mode_index(modes_graph, "curve")>;

inline sf::Vertex * get_segment(Curve &self, std::size_t const i = 0) {
  return &(self.vertexes[i * 4]);
}

inline sf::Vertex * prev_segment(Curve &, sf::Vertex * const segment) {
  return segment - 4;
}

inline sf::Vertex * next_segment(Curve &, sf::Vertex * const segment) {
  return segment + 4;
}

inline sf::Vertex * set_segment(Curve &self, sf::Vertex * segment,
    sf::Vector2f const p0, sf::Vector2f const p1) {
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
    sf::Vector2f const p0, sf::Vector2f const p1) {
  return set_segment(self, get_segment(self, i), p0, p1);
}

inline sf::Vertex * set_segment(Curve &, sf::Vertex * segment,
    sf::Color const c0, sf::Color const c1) {
  (segment++)->color = c0;
  (segment++)->color = c0;
  (segment++)->color = c1;
  (segment++)->color = c1;
  return segment;
}

inline sf::Vertex * set_segment(Curve &self, std::size_t const i,
    sf::Color const c0, sf::Color const c1) {
  return set_segment(self, get_segment(self, i), c0, c1);
}

inline Curve::state_t iterative_set_init(Curve &self,
    std::optional<sf::Color> const c_init = {}) {
  return {prev_segment(self, get_segment(self)), {},
    c_init.value_or(get_segment(self)->color)};
}

inline Curve::state_t iterative_set_point(Curve &self, sf::Vector2f const p,
    Curve::state_t const &s) {
  sf::Vertex * const next{s.segment >= get_segment(self)
    ? set_segment(self, s.segment, s.p, p)
    : next_segment(self, s.segment)};
  return {next, p, s.c};
}

inline Curve::state_t iterative_set_point(Curve &self, sf::Vector2f const p,
    sf::Color const c, Curve::state_t const &s) {
  if (s.segment >= get_segment(self)) set_segment(self, s.segment, s.c, c);
  return iterative_set_point(self, p, {s.segment, s.p, c});
}

template<std::size_t mode>
struct Visualizer {};

template<std::size_t mode>
void compute(Visualizer<mode> &, Buffer &) {}

template<std::size_t mode>
void draw(Visualizer<mode> &, Buffer const &, sf::RenderTarget &) {}

template<>
struct Visualizer<args::get_mode_index(modes_visualizer, "time-series")> {
  // Parameters
  float y_scale;
  bool scrolling;

  // State
  std::vector<sf::VertexArray> vertexess;

  std::size_t n_timepoints() const {
    return this->vertexess.front().getVertexCount();
  }

  std::size_t n_channels() const {
    return this->vertexess.size();
  }

  Visualizer(std::size_t const n_channels, std::size_t const n_timepoints,
      sf::Color const color, float const y_scale, bool const scrolling = true) :
      y_scale{y_scale}, scrolling{scrolling},
      vertexess{std::vector<sf::VertexArray>(n_channels,
        sf::VertexArray(sf::LineStrip, n_timepoints))} {
    for (auto &vertexes : this->vertexess)
      for (std::size_t i{0}; i < n_timepoints; ++i)
        vertexes[i].color = color;
  }
};

using TimeSeries =
  Visualizer<args::get_mode_index(modes_visualizer, "time-series")>;

void draw(TimeSeries &self, Buffer const &buffer, sf::RenderTarget &target) {
  std::size_t const n_timepoints{
    std::min(self.n_timepoints(), buffer.n_timepoints)};
  std::size_t const n_channels{
    std::min(self.n_channels(), buffer.n_channels())};
  for (std::size_t i{0}; i < n_timepoints; ++i)
    for (std::size_t j{0}; j < n_channels; ++j) {
      float const s{self.scrolling
        ? buffer[{j, n_timepoints - 1 - i}]
        : buffer.data[i * n_channels + j]};
      float const x{i / static_cast<float>(n_timepoints - 1)};
      float const y{(s * (.5f * self.y_scale) + j + .5f) / n_channels};
      self.vertexess[j][i].position = sf::Vector2f{x, y};
    }
  for (auto const &vertexes : self.vertexess) target.draw(vertexes);
}

template<>
struct Visualizer<
      args::get_mode_index(modes_visualizer, "short-time-spectrum")>  {
  // Parameters
  std::size_t n_in;
  float x_min, x_max;
  float y_min_db, y_max_db;
  float thickness;

  // State
  float * window, * in, * aggregate;
  fftwf_complex * out;
  fftwf_plan plan;
  //sf::VertexArray vertexes;
  Curve graph;

  std::size_t n_out() const {
    return this->n_in / 2u + 1u;
  }

  std::complex<float> * out_std() const {
    return reinterpret_cast<std::complex<float> *>(this->out);
  }

  float f_min(Buffer const &buffer) const {
    return buffer.sample_rate() * (1.f / this->n_in);
  }
  float f_max(Buffer const &buffer) const {
    return buffer.sample_rate() *
      (static_cast<float>(this->n_out() - 1) / this->n_in);
  }
  float f_range(Buffer const &buffer) const {
    return this->f_max(buffer) - this->f_min(buffer);
  }

  float x_min_log2() const { return std::log2(this->x_min); }
  float x_max_log2() const { return std::log2(this->x_max); }
  float x_range_log2() const { return this->x_max_log2() - this->x_min_log2(); }

  float y_min() const { return std::exp2(this->y_min_db / 6.f); }
  float y_max() const { return std::exp2(this->y_max_db / 6.f); }
  float y_range_db() const { return this->y_max_db - this->y_min_db; }

  Visualizer(std::size_t const n_in, auto &&window_function,
      sf::Color const color, float const x_min, float const x_max,
      float const y_min_db, float const y_max_db, float const thickness) :
      n_in{n_in}, x_min{x_min}, x_max{x_max},
      y_min_db{y_min_db}, y_max_db{y_max_db},
      thickness{thickness},
      window{[&](){
        float * const window{
          static_cast<float *>(fftwf_malloc(sizeof(float) * n_in))};
        for (std::size_t i{0}; i < n_in; ++i)
          window[i] = window_function(i, n_in);
        return window; }()},
      in{static_cast<float *>(fftwf_malloc(sizeof(float) * n_in))},
      aggregate{static_cast<float *>(
        fftwf_malloc(sizeof(float) * this->n_out()))},
      out{static_cast<fftwf_complex *>(
        fftwf_malloc(sizeof(fftwf_complex) * this->n_out()))},
      plan{fftwf_plan_dft_r2c_1d(n_in, in, out, FFTW_ESTIMATE)},
      //vertexes{[&](){
      //  sf::VertexArray vertexes(sf::LineStrip, this->n_out() - 1);
      //  for (std::size_t i{0}; i < vertexes.getVertexCount(); ++i)
      //    vertexes[i].color = color;
      //  return vertexes; }()}
      graph(this->n_out() - 1, color, {}) {}

  ~Visualizer() {
    fftwf_destroy_plan(this->plan);
    fftwf_free(this->in);
    fftwf_free(this->out);
    fftwf_free(this->aggregate);
  }
};

using ShortTimeSpectrum =
  Visualizer<args::get_mode_index(modes_visualizer, "short-time-spectrum")>;

void compute(ShortTimeSpectrum &self, Buffer &buffer) {
  std::fill(self.aggregate, self.aggregate + self.n_out(), 0.f);
  for (std::size_t j{0}; j < buffer.n_channels(); ++j) {
    for (std::size_t i{0}; i < self.n_in; ++i)
      self.in[i] = self.window[i] * buffer[{j, self.n_in - 1 - i}];
    fftwf_execute(self.plan);
    std::transform(self.out_std(), self.out_std() + self.n_out(),
      self.aggregate, self.aggregate,
      [](auto const &out, auto const &agg){
        auto const abs2{[](auto const &x){ return x * x; }};
        return abs2(std::real(out)) + abs2(std::imag(out)) + agg; });
  }
  std::transform(self.aggregate, self.aggregate + self.n_out(),
    self.aggregate, [&](auto const x){
      return 6.f * std::log2(std::sqrt(x * (1.f / buffer.n_channels()))); });
}

void draw(ShortTimeSpectrum &self, Buffer const &buffer,
    sf::RenderTarget &target) {
  self.graph.thickness_2d = thickness_px_to_2d(self.thickness, target);
  auto s{iterative_set_init(self.graph)};
  for (std::size_t i{1}; i < self.n_out(); ++i) {
    float const x{(std::log2(self.f_min(buffer) + self.f_range(buffer) *
        static_cast<float>(i - 1) / (self.n_out() - 2)) -
      self.x_min_log2()) / self.x_range_log2()};
    float const y{(self.aggregate[i] - self.y_min_db) / self.y_range_db()};
    //self.vertexes[i - 1].position = sf::Vector2f{x, y};
    s = iterative_set_point(self.graph, {x, y}, s);
  }
  //target.draw(self.vertexes);
  draw(target, self.graph);
}

void pan(sf::RenderTarget &target, int const &delta_x, int const &delta_y,
    bool const fine){
  sf::View view{target.getView()};
  sf::Vector2f size_view{view.getSize()};
  sf::Vector2u const size_target{target.getSize()};
  if (fine) size_view *= .1f;
  view.move(delta_x * (size_view.x / size_target.x),
            delta_y * (size_view.y / size_target.y));
  target.setView(view);
}

void zoom_anchored(sf::RenderTarget &target, float const &delta_x,
    float const &delta_y, sf::Vector2i const &anchor_target, bool const fine){
  float const factor{fine ? .005f : .05f};
  sf::Vector2f const alpha{std::exp2(factor * delta_x),
                           std::exp2(factor * delta_y)};

  sf::View view{target.getView()};
  sf::Vector2u const size_target{target.getSize()};
  sf::Vector2f const size_view{view.getSize()};
  sf::Vector2f const center_view{view.getCenter()};
  sf::Vector2f const anchor_view{
    anchor_target.x * (size_view.x / size_target.x) + center_view.x,
    anchor_target.y * (size_view.y / size_target.y) + center_view.y};
  sf::Vector2f const size_view_new{size_view.x * alpha.x,
                                   size_view.y * alpha.y};

  sf::Vector2f const center_view_new{
    (center_view.x - anchor_view.x) * alpha.x +
      anchor_view.x - .5f * (size_view.x - size_view_new.x),
    (center_view.y - anchor_view.y) * alpha.y +
      anchor_view.y - .5f * (size_view.y - size_view_new.y)};
  target.setView({center_view_new, size_view_new});
}

args::flag_descs_t const flag_descs_visualizer{
  {"use-compute-thread", "Separate the data processing from the graphical "
    "processing thread."},
  {"no-vertical-sync", nullptr},
  {"fullscreen", nullptr},
  {"no-highpass", "Disable preprocessing highpass filter."}};

args::opt_descs_t const opt_descs_visualizer{
  {"window-width", {"n", "1024", "Window width in pixels.", {}}},
  {"window-height", {"n", "768", "Window height in pixels.", {}}},
  {"bits-per-pixel", {"n", "32", nullptr, {}}},
  {"antialiasing-level", {"n", "8", nullptr, {}}},
  {"frame-rate-limit", {
    "f", "0", "Limit frame rate to <f> Hz.", {
      {"0", "Disable frame rate limit."}}}},
  {"ref-index", {
    "i", "-1", "Use channel <i> as reference.", {
      {"-1", "No (re-)referencing."}}}},
  {"buffer-duration", {
    "t", "10", "Use a data processing buffer length of <t> seconds.", {}}},
  {"highpass-frequency", {
    "f", ".25", "Cutoff frequency for the preprocessing highpass filter.", {}}},
  {"highpass-q", {
    "q", ".70710678118654752",
    "Quality factor for the preprocessing highpass filter.", {}}},
  {"timeout-resolve", {
    "t", "3", "Wait <t> seconds for input stream to be resolved before "
      "aborting.", {}}},
  {"timeout-open", {
    "t", "3", "Wait <t> seconds for input stream to be opened before aborting.",
    {}}},
  {"lsl-property", {
    "str", "type", "Search LSL streams that match on this metadata property.", {
      {"type", "Stream type."}}}},
  {"lsl-value", {
    "str", "EEG", "Match the LSL metadata property to this value.", {
      {"EEG", "With --lsl-property=type, search for EEG streams."}}}}};

args::flag_descs_t const flag_descs_time_series{{"help", nullptr}};

args::opt_descs_t const opt_descs_time_series{};

args::flag_descs_t const flag_descs_short_time_spectrum{{"help", nullptr}};

args::opt_descs_t const opt_descs_short_time_spectrum{
  {"fft-duration",
    {"t", "1", "Compute FFT on a time-domain window of <t> seconds.", {}}}};

int print_usage_visualizer(args::arg_t const &program_name,
    args::UsageMessageMode const mode) {
  return args::print_usage(program_name, {"client", "visualizer"},
    desc_visualizer, std::next(std::cbegin(modes_visualizer)),
    std::cend(modes_visualizer), std::next(std::cbegin(modes_descs_visualizer)),
    std::cend(modes_descs_visualizer), flag_descs_visualizer, 
    opt_descs_visualizer, mode);
}

int print_usage_time_series(args::arg_t const &program_name,
    args::UsageMessageMode const mode) {
  return args::print_usage(program_name,
    {"client", "visualizer", "time-series"}, desc_time_series, nullptr, nullptr,
    nullptr, nullptr, flag_descs_time_series, opt_descs_time_series, mode);
}

int print_usage_short_time_spectrum(args::arg_t const &program_name,
    args::UsageMessageMode const mode) {
  return args::print_usage(program_name,
    {"client", "visualizer", "time-series"}, desc_short_time_spectrum, nullptr,
    nullptr, nullptr, nullptr, flag_descs_short_time_spectrum,
    opt_descs_short_time_spectrum, mode);
}

template<std::size_t mode_visualizer>
int main_client_visualizer(args::args_t::const_iterator &pos,
    args::args_t::const_iterator const &end, args::arg_t const &program_name,
    args::flags_t &flags_visualizer, args::opts_t &opts_visualizer) {
  args::flags_t flags;
  args::opts_t opts;
  int (*print_usage)(args::arg_t const &, args::UsageMessageMode const);
  if constexpr (mode_visualizer ==
      args::get_mode_index(modes_visualizer, "help"))
    return print_usage_visualizer(program_name, args::UsageMessageMode::help);
  else if constexpr (mode_visualizer ==
      args::get_mode_index(modes_visualizer, "time-series")) {
    flags = args::init_flags(flag_descs_time_series);
    opts = args::init_opts(opt_descs_time_series);
    print_usage = &print_usage_time_series;
  } else if constexpr (mode_visualizer ==
      args::get_mode_index(modes_visualizer, "short-time-spectrum")) {
    flags = args::init_flags(flag_descs_short_time_spectrum);
    opts = args::init_opts(opt_descs_short_time_spectrum);
    print_usage = &print_usage_short_time_spectrum;
  } else // catch-all, includes `mode_visualizer == error`
    return print_usage_visualizer(program_name, args::UsageMessageMode::error);

  args::parse_opts_and_flags(flags, opts, pos, end);
  if (flags["help"] or pos != end) return (*print_usage)(
    program_name, flags["help"] ? args::UsageMessageMode::help :
      args::UsageMessageMode::error);


  // Parameters
  std::string const property{opts_visualizer["lsl-property"].value()};
  std::string const value{opts_visualizer["lsl-value"].value()};
  double const timeout_resolve{args::parse_opt(
    args::parse_double, opts_visualizer["timeout-resolve"])};
  double const timeout_open{args::parse_opt(
    args::parse_double, opts_visualizer["timeout-open"])};

  double const buffer_duration{args::parse_opt(
    args::parse_double, opts_visualizer["buffer-duration"])};
  bool const use_compute_thread{flags_visualizer["use-compute-thread"]};
  double const update_rate{60.};

  long const ref_index{args::parse_opt(
    args::parse_long, opts_visualizer["ref-index"])};
  bool const use_highpass{not flags_visualizer["no-highpass"]};
  double const highpass_frequency{args::parse_opt(
    args::parse_double, opts_visualizer["highpass-frequency"])};
  double const highpass_quality_factor{args::parse_opt(
    args::parse_double, opts_visualizer["highpass-q"])};

  std::string const window_title{std::string{"Visualizer ("} +
    modes_visualizer[mode_visualizer] + ")"};
  auto const window_style{flags_visualizer["fullscreen"] ?
    sf::Style::Fullscreen : sf::Style::Resize | sf::Style::Close};
  sf::VideoMode const window_video_mode{
    args::parse_opt(args::parse_unsigned_int, opts_visualizer["window-width"]),
    args::parse_opt(args::parse_unsigned_int, opts_visualizer["window-height"]),
    args::parse_opt(args::parse_unsigned_int, opts_visualizer["bits-per-pixel"]),
  };
  sf::ContextSettings const window_context_settings{
    0, // depth bits
    0, // stencil bits
    args::parse_opt(args::parse_unsigned_int,
      opts_visualizer["antialiasing-level"]),
  };
  unsigned int const frame_rate_limit{args::parse_opt(
    args::parse_unsigned_int, opts_visualizer["frame-rate-limit"])};
  bool const vertical_sync{not flags_visualizer["no-vertical-sync"]};
  sf::View const default_view{{0.f, 1.f, 1.f, -1.f}};
  auto const mouse_button_move{sf::Mouse::Left};
  auto const mouse_button_scale{sf::Mouse::Right};
  auto const key_reset{sf::Keyboard::Enter};
  auto const key_fine{sf::Keyboard::LShift};

  // Set up LSL input and buffer
  auto lsl_attachment_var{lsl::attach(property, value,
      timeout_resolve, timeout_open)};
  if (lsl::has_failed(lsl_attachment_var))
    return lsl::exit_code(lsl_attachment_var);
  Buffer buffer(lsl::get_attachment(lsl_attachment_var), buffer_duration);

  // Set up preprocessing (TODO: Abstract this)
  auto const [highpass_bs, highpass_as]{
    dsp::biquad_rbj<double, dsp::BiquadShape::highpass>(
      2 * pi<double> * highpass_frequency / buffer.sample_rate(), 
      highpass_quality_factor)};
  dsp::IIRSpec<float, 2, 2> highpass_spec(highpass_bs, highpass_as);
  auto highpass_state{dsp::make_state(highpass_spec, buffer.n_channels())};

  // Set up visualizer
  Visualizer<mode_visualizer> visualizer{[&](){
  if constexpr (mode_visualizer ==
      args::get_mode_index(modes_visualizer, "time-series")) {
    std::size_t const n_timepoints{buffer.n_timepoints};
    sf::Color const color{0xff, 0xff, 0xff, 0xff};
    float const y_scale{1.f};
    bool const scrolling{true};
    return Visualizer<mode_visualizer>(
      buffer.n_channels(), n_timepoints, color, y_scale, scrolling);
  } else if constexpr (mode_visualizer ==
      args::get_mode_index(modes_visualizer, "short-time-spectrum")) {
    double const fft_duration{args::parse_opt(
      args::parse_double, opts["fft-duration"])};
    std::size_t const n_in{static_cast<std::size_t>(
      std::round(buffer.sample_rate() * fft_duration))};
    auto const window_function{[](std::size_t const i, std::size_t const n){
      return dsp::window<float, dsp::WindowShape::blackman>(i, n, .16f); }};
    sf::Color const color{0xff, 0xff, 0xff, 0xff};
    float const x_min{.5f}, x_max{80.f}, y_min_db{-24.f}, y_max_db{48.f};

    if (n_in > buffer.n_timepoints) throw(std::out_of_range{"Requested FFT "
      "size (" + std::to_string(n_in) + ") is bigger than available buffer "
      "size (" + std::to_string(buffer.n_timepoints) + ")"});

    return Visualizer<mode_visualizer>(n_in, window_function, color,
      x_min, x_max, y_min_db, y_max_db, 5.f);
  } else {
    return Visualizer<mode_visualizer>{};
  }}()};

  // Set up computation thread
  auto const compute_step{[&](auto const /*update_rate*/){
      pull(buffer);

      { // Preprocessing (TODO: Abstract this)
        for (auto sample{get_sample(buffer)};
            not past_newest_sample(buffer, sample);
            advance_sample(buffer, sample)) {
          auto const &[first, last] = sample;
          if (ref_index >= 0 and
              static_cast<std::size_t>(ref_index) < buffer.n_channels()) {
            float const x_ref{*(first + ref_index)};
            std::transform(first, last, first, [=](float const x){
              return x - x_ref; });
            }
          if (use_highpass)
            dsp::tick(first, last, first, highpass_spec, highpass_state);
        }
      }

      compute(visualizer, buffer);
    }};

  std::thread t_compute;
  if (use_compute_thread)
    t_compute = std::thread{[&](){
    std::chrono::duration<double> update_interval{1 / update_rate};
    while (not graceful_exit_signal()) {
      compute_step(update_rate);
      std::this_thread::sleep_for(update_interval);
    }}};

  // Set up window
  sf::RenderWindow window{window_video_mode, window_title,
    static_cast<sf::Uint32>(window_style), window_context_settings};
  window.setFramerateLimit(frame_rate_limit);
  window.setVerticalSyncEnabled(vertical_sync);
  window.setView(default_view);
  sf::Event event{};
  sf::Vector2i mouse_pos_prev{};
  std::optional<sf::Vector2i> anchor_window{};

  while (window.isOpen()) {
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) graceful_exit_signal(true);
      else if (event.type == sf::Event::MouseWheelScrolled)
        zoom_anchored(window,
          event.mouseWheelScroll.delta, event.mouseWheelScroll.delta,
          {event.mouseWheelScroll.x, event.mouseWheelScroll.y},
          sf::Keyboard::isKeyPressed(key_fine));
      else if (event.type == sf::Event::KeyPressed and
          event.key.code == key_reset) window.setView(default_view);
    }
    if (graceful_exit_signal()) break;

    sf::Vector2i const mouse_pos{sf::Mouse::getPosition(window)};
    if (sf::Mouse::isButtonPressed(mouse_button_move)) {
      sf::Vector2i const mouse_delta{mouse_pos_prev - mouse_pos};
      pan(window, mouse_delta.x, mouse_delta.y, sf::Keyboard::isKeyPressed(key_fine));
    } else if (sf::Mouse::isButtonPressed(mouse_button_scale)) {
      if (not anchor_window.has_value()) anchor_window = mouse_pos_prev;
      sf::Vector2i const mouse_delta{mouse_pos_prev - mouse_pos};
      float const factor{.1f};
      zoom_anchored(window, factor * mouse_delta.x, factor * mouse_delta.y,
        anchor_window.value(), sf::Keyboard::isKeyPressed(key_fine));
    } else anchor_window = {};
    mouse_pos_prev = mouse_pos;

    window.clear();
    double const nominal_frame_rate{60.};
    if (not use_compute_thread) compute_step(nominal_frame_rate);
    draw(visualizer, buffer, window);
    window.display();
  }

  window.close();
  if (use_compute_thread && t_compute.joinable()) t_compute.join();
  return 0;
}

template<std::size_t mode_itr = 0, class... Ts>
int main_client_visualizer(std::size_t const mode_visualizer, Ts &... args) {
  if (mode_itr == mode_visualizer)
    return main_client_visualizer<mode_itr>(args...);
  if constexpr (mode_itr >= 0 and mode_itr < modes_visualizer.size())
    return main_client_visualizer<mode_itr + 1>(mode_visualizer, args...);
  return 1;
}

int main_client_visualizer(args::args_t::const_iterator &pos,
    args::args_t::const_iterator const &end, args::arg_t const &program_name) {
  auto flags_visualizer{args::init_flags(flag_descs_visualizer)};
  auto opts_visualizer{args::init_opts(opt_descs_visualizer)};
  args::parse_opts_and_flags(flags_visualizer, opts_visualizer, pos, end);

  std::size_t mode_visualizer{
    args::get_mode_index(modes_visualizer, pos, end, 0)};
  return main_client_visualizer(mode_visualizer, pos, end,
    program_name, flags_visualizer, opts_visualizer);
}