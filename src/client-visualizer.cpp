#include <type_traits>
#include <cstddef>
#include <cmath>
#include <complex>
#include <utility>
#include <tuple>
#include <memory>
#include <optional>
#include <array>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <mutex>
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
#include "number-printing.hpp"
#include "numbers.hpp"
#include "affine.hpp"
#include "dsp.hpp"
#include "plotting.hpp"
#include "lsl-customized.hpp"

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

struct VisualizerSingle2DFrameBase {
  Box<float> plot_box;
  PlotAnnotations pa;
  PlotAnnotationsStore pa_store;
};

Box<float> get_plot_to_frame_box(VisualizerSingle2DFrameBase const &self,
    Box<float> const &parent_box, sf::Vector2f const &view_scale) {
  return box_to_box(self.plot_box, get_frame_box(self.pa, parent_box,
    view_scale));
}

Box<float> get_frame_to_plot_box(VisualizerSingle2DFrameBase const &self,
    Box<float> const &parent_box, sf::Vector2f const &view_scale) {
  return box_to_box(get_frame_box(self.pa, parent_box, view_scale),
    self.plot_box);
}

void set_plot_ticks(VisualizerSingle2DFrameBase &self,
    Box<PlotAnnotations::majors_t> const &majors,
    Box<PlotAnnotations::minors_t> const &minors,
    Box<float> const &parent_box, sf::Vector2f const &view_scale) {
  std::array<float, 2> plot_to_frame;
  auto const set_edge{[&](auto &outs, auto const &ins, auto const &lerper){
      if (outs.size() != ins.size()) outs.resize(ins.size());
      std::transform(ins.begin(), ins.end(), outs.begin(), lerper);
    }};
  auto const lerp_minor{[&](PlotAnnotations::minor_t const &in){
      return lerp(in, plot_to_frame);
    }};
  auto const lerp_major{[&](PlotAnnotations::major_t const &in){
      auto const &[pos, label]{in};
      return std::make_pair(lerp_minor(pos), label);
    }};

  auto const plot_to_frame_box{
    get_plot_to_frame_box(self, parent_box, view_scale)};
  plot_to_frame = plot_to_frame_box.xs();
  set_edge(self.pa.majors.x0, majors.x0, lerp_major);
  set_edge(self.pa.majors.x1, majors.x1, lerp_major);
  set_edge(self.pa.minors.x0, minors.x0, lerp_minor);
  set_edge(self.pa.minors.x1, minors.x1, lerp_minor);
  plot_to_frame = plot_to_frame_box.ys();
  set_edge(self.pa.majors.y0, majors.y0, lerp_major);
  set_edge(self.pa.majors.y1, majors.y1, lerp_major);
  set_edge(self.pa.minors.y0, minors.y0, lerp_minor);
  set_edge(self.pa.minors.y1, minors.y1, lerp_minor);
}

void draw_underlay(VisualizerSingle2DFrameBase const &self,
    sf::RenderTarget &target, ViewTransform const &, sf::Vector2f const &) {
  draw_underlay(target, self.pa_store);
}

void draw_overlay(VisualizerSingle2DFrameBase const &self,
    sf::RenderTarget &target, ViewTransform const &vt,
    sf::Vector2f const &view_scale) {
  draw_overlay(target, self.pa_store);

  if constexpr (draw_bounding_boxes) {
    auto const target_box{get_render_target_box(target)},
      frame_box{get_frame_box(self.pa, target_box, view_scale)};
    draw_bounding_box(target, target_box, vt);
    draw_bounding_box(target, frame_box, vt);
  }
}

template<std::size_t mode>
struct Visualizer {};

template<std::size_t mode>
auto axes_pan(Visualizer<mode> &self, sf::Vector2f const &delta,
    bool const fine, Box<float> const &parent_box,
    sf::Vector2f const &view_scale){
  float const factor{fine ? .1f : 1.f};
  if constexpr (std::is_base_of_v<VisualizerSingle2DFrameBase,
      Visualizer<mode>>) {
    return self.plot_box -= scale(sf_vec_to_array(delta * factor),
      get_frame_to_plot_box(self, parent_box, view_scale));
  }
}

template<std::size_t mode>
auto axes_zoom(Visualizer<mode> &self, sf::Vector2f const &delta,
    array2d<float, 2, 1> const &anchor, bool const fine,
    Box<float> const &parent_box, sf::Vector2f const &view_scale) {
  if constexpr (std::is_base_of_v<VisualizerSingle2DFrameBase,
      Visualizer<mode>>) {
    float const factor{(fine ? .1f : 1.f) * -.01f};
    sf::Vector2f const alpha{std::exp2(delta.x * factor),
                             std::exp2(delta.y * factor)};
    auto const frame_to_plot_box{
      get_frame_to_plot_box(self, parent_box, view_scale)};
    auto const plot_anchor{lerp(anchor, frame_to_plot_box)};
    self.plot_box *= sf_vec_to_array(alpha);
    self.plot_box += sf_vec_to_array(sf::Vector2f{
      ((get<0, 0>(plot_anchor) - self.plot_box.x0) * (1.f - alpha.x)),
      ((get<1, 0>(plot_anchor) - self.plot_box.y0) * (1.f - alpha.y))});
    return self.plot_box;
  }
}

template<std::size_t mode>
void compute(Visualizer<mode> &, std::mutex &, Buffer &) {}

template<std::size_t mode>
void draw(Visualizer<mode> &, std::mutex &, Buffer const &, sf::RenderTarget &,
    ViewTransform const &, sf::Vector2f const &) {}

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
      sf::Color const &color, float const y_scale,
      bool const scrolling = true) :
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

void draw(TimeSeries &self, std::mutex &compute_draw_mutex,
    Buffer const &buffer, sf::RenderTarget &target, ViewTransform const &vt,
    sf::Vector2f const &view_scale) {
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
      args::get_mode_index(modes_visualizer, "short-time-spectrum")> : 
    public VisualizerSingle2DFrameBase {
  // Parameters
  std::size_t n_in;
  float curve_thickness;

  // State
  float * window, * in, * aggregate;
  fftwf_complex * out;
  fftwf_plan plan;
  //sf::VertexArray vertexes;
  Curve graph;

  std::size_t n_out() const { return this->n_in / 2u + 1u; }

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

  Visualizer(std::size_t const n_in, auto &&window_function,
      Box<float> const &plot_box, RGB8 const &color,
      float const curve_thickness, PlotAnnotations const &pa) :
      VisualizerSingle2DFrameBase{plot_box, pa, {}},
      n_in{n_in}, curve_thickness{curve_thickness},
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
      graph(this->n_out() - 1, color) {}

  ~Visualizer() {
    fftwf_destroy_plan(this->plan);
    fftwf_free(this->in);
    fftwf_free(this->out);
    fftwf_free(this->aggregate);
  }
};

using ShortTimeSpectrum =
  Visualizer<args::get_mode_index(modes_visualizer, "short-time-spectrum")>;

void compute(ShortTimeSpectrum &self, std::mutex &compute_draw_mutex,
    Buffer &buffer) {
  std::scoped_lock lock{compute_draw_mutex};
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

void draw(ShortTimeSpectrum &self, std::mutex &compute_draw_mutex,
    Buffer const &buffer, sf::RenderTarget &target, ViewTransform const &vt,
    sf::Vector2f const &view_scale) {
  auto render_target_box{get_render_target_box(target)};
  auto plot_to_frame_box{
    get_plot_to_frame_box(self, render_target_box, view_scale)};
  auto s{iterative_set_init(self.graph,
    {self.curve_thickness * .5f, self.curve_thickness * .5f})};
  {
    std::scoped_lock lock{compute_draw_mutex};
    for (std::size_t i{1}; i < self.n_out(); ++i) {
      float const i_l2Hz{dsp::to_l2Hz(lerp(static_cast<float>(i - 1),
        {0.f, static_cast<float>(self.n_out() - 2)},
        {self.f_min(buffer), self.f_max(buffer)}))};
      std::array<float, 2> const p{i_l2Hz, self.aggregate[i]};
      s = iterative_set_point(self.graph, s, vt,
        array_to_sf_vec(lerp(p, plot_to_frame_box)), {});
    }
  }
  draw_underlay(self, target, vt, view_scale);
  draw(target, self.graph);
  draw_overlay(self, target, vt, view_scale);
}

args::flag_descs_t const flag_descs_visualizer{
  {"no-compute-thread", "Do not separate the data processing from the graphical"
    " processing thread."},
  {"no-vertical-sync", nullptr},
  {"fullscreen", nullptr},
  {"no-highpass", "Disable preprocessing highpass filter."}};

args::opt_descs_t const opt_descs_visualizer{
  {"gui-scale", {"factor", "1", "Overall scaling factor for the graphical user "
    "interface.", {}}},
  {"window-width", {"n", "1024", "Window width in pixels.", {}}},
  {"window-height", {"n", "768", "Window height in pixels.", {}}},
  {"bits-per-pixel", {"n", "32", nullptr, {}}},
  {"antialiasing-level", {"n", "16", nullptr, {}}},
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
  std::string const lsl_property{opts_visualizer["lsl-property"].value()};
  std::string const lsl_value{opts_visualizer["lsl-value"].value()};
  double const timeout_resolve{args::parse_opt(
    args::parse_double, opts_visualizer["timeout-resolve"])};
  double const timeout_open{args::parse_opt(
    args::parse_double, opts_visualizer["timeout-open"])};

  double const buffer_duration{args::parse_opt(
    args::parse_double, opts_visualizer["buffer-duration"])};
  bool const use_compute_thread{not flags_visualizer["no-compute-thread"]};
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

  float const gui_scale{args::parse_opt(
    args::parse_float, opts_visualizer["gui-scale"])};

  sf::Vector2f const default_view_scale{1.f, 1.f},
    default_view_offset{0.f, 0.f};
  float const default_view_rotation{0.f};

  Box<float> default_plot_box{0.f, 0.f, 1.f, 1.f};
  Box<PlotAnnotations::majors_t> majors{{}, {}, {}, {}};
  Box<PlotAnnotations::minors_t> minors{{}, {}, {}, {}};

  sf::Color const background_color{0x00, 0x00, 0x00};

  auto const font{std::make_shared<sf::Font>()};
  font->loadFromFile("/System/Library/Fonts/Helvetica.ttc");
  //font->setSmooth(true);

  auto const mouse_button_mod_switch{sf::Mouse::Right};
  auto const key_left{sf::Keyboard::Left};
  auto const key_right{sf::Keyboard::Right};
  auto const key_down{sf::Keyboard::Down};
  auto const key_up{sf::Keyboard::Up};
  auto const key_reset{sf::Keyboard::Enter};
  auto const key_mod_view{sf::Keyboard::LShift};
  auto const key_mod_mode{sf::Keyboard::LControl};
  auto const key_mod_fine{sf::Keyboard::LAlt};
  auto const key_mod_switch{sf::Keyboard::RShift};
  auto const double_click_min_duration{std::chrono::milliseconds(0)};
  auto const double_click_max_duration{std::chrono::milliseconds(300)};

  // Set up LSL input and buffer
  auto lsl_attachment_var{lsl::attach(lsl_property, lsl_value,
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
    return TimeSeries(
      buffer.n_channels(), n_timepoints, color, y_scale, scrolling);
  } else if constexpr (mode_visualizer ==
      args::get_mode_index(modes_visualizer, "short-time-spectrum")) {
    double const fft_duration{args::parse_opt(
      args::parse_double, opts["fft-duration"])};
    std::size_t const n_in{static_cast<std::size_t>(
      std::round(buffer.sample_rate() * fft_duration))};
    auto const window_function{[](std::size_t const i, std::size_t const n){
      return dsp::window<float, dsp::WindowShape::blackman>(i, n, .16f); }};
    RGB8 const curve_color{0x55, 0x77, 0xff};
    float const x_min{.5f}, x_max{80.f}, y_min_dB{-24.f}, y_max_dB{24.f};
    default_plot_box = {dsp::to_l2Hz(x_min), y_min_dB,
                        dsp::to_l2Hz(x_max), y_max_dB};
    for (float const tick : {.2f, .3f, .4f, .5f, .6f, .7f, .8f, .9f, 2.f, 3.f,
        4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f,
        90.f})
      minors.x0.push_back(dsp::to_l2Hz(tick));
    for (float const tick : {.1f, 1.f, 10.f, 100.f})
      majors.x0.push_back(
        PlotAnnotations::major_t{dsp::to_l2Hz(tick), {fancy(tick)}});
    for (float const tick : {-21.f, -18.f, -15.f, -9.f, -6.f, -3.f, 3.f, 6.f,
        9.f, 15.f, 18.f, 21.f})
      minors.y0.push_back(tick);
    for (float const tick : {-24.f, -12.f, 0.f, 12.f, 24.f})
      majors.y0.push_back(
        PlotAnnotations::major_t{tick, {fancy(tick)}});
    float const curve_thickness{4.f};

    if (n_in > buffer.n_timepoints) throw(std::out_of_range{"Requested FFT "
      "size (" + std::to_string(n_in) + ") is bigger than available buffer "
      "size (" + std::to_string(buffer.n_timepoints) + ")"});

    PlotAnnotations const pa{
          {}, {},
          {{"Frequency [Hz]"}, {"Magnitude [dBFS]"}, {}, {}},
          {PlotAnnotations::on,
            PlotAnnotations::on | PlotAnnotations::vertical,
            PlotAnnotations::on | PlotAnnotations::alt | PlotAnnotations::nolabel,
            PlotAnnotations::on | PlotAnnotations::alt | PlotAnnotations::nolabel | PlotAnnotations::vertical},
          {PlotAnnotations::on,
            PlotAnnotations::on,
            PlotAnnotations::on | PlotAnnotations::alt,
            PlotAnnotations::on | PlotAnnotations::alt},
          {PlotAnnotations::on,
            PlotAnnotations::on,
            PlotAnnotations::off,
            PlotAnnotations::off},
          {PlotAnnotations::on,
            PlotAnnotations::on,
            PlotAnnotations::off,
            PlotAnnotations::off},
          {PlotAnnotations::on | PlotAnnotations::italic,
            PlotAnnotations::on | PlotAnnotations::italic | PlotAnnotations::vertical,
            PlotAnnotations::off | PlotAnnotations::italic,
            PlotAnnotations::off | PlotAnnotations::italic | PlotAnnotations::vertical},
          {0x7f, 0x7f, 0x7f}, 0xff,
          {0xff, 0xff, 0xff}, 0x37,
          {0xff, 0xff, 0xff}, 0x15,
          {0x7f, 0x7f, 0x7f}, 0xff,
          {0x7f, 0x7f, 0x7f}, 0xff,
          {0x00, 0x00, 0x00}, 0xcc,
          gui_scale * 7.f,
          gui_scale * 3.5f,
          gui_scale * 3.f,
          gui_scale * 2.f,
          gui_scale * 2.f,
          gui_scale * 30.f,
          gui_scale * 25.f,
          gui_scale * 10.f,
          gui_scale * 10.f,
          gui_scale * 10.f,
          gui_scale * 4096.f,
          font};

    return ShortTimeSpectrum(n_in, window_function, default_plot_box,
      curve_color, gui_scale * curve_thickness, pa);
  } else {
    return Visualizer<mode_visualizer>{};
  }}()};

  // Set up computation thread
  std::mutex compute_draw_mutex;
  
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

      compute(visualizer, compute_draw_mutex, buffer);
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

  sf::Cursor cursor_default;
  sf::Cursor cursor_cross;
  bool const cursor_default_supported{
    cursor_default.loadFromSystem(sf::Cursor::Arrow)};
  bool const cursor_cross_supported{
    cursor_cross.loadFromSystem(sf::Cursor::Cross)};

  sf::Vector2f view_offset, view_scale;
  float view_rotation;
  ViewTransform view_transform, view_transform_inv;
  auto const update_time_invariant_annotations{[&](){
      if constexpr (std::is_base_of_v<VisualizerSingle2DFrameBase,
          decltype(visualizer)>) {
        auto const render_target_box{get_render_target_box(window)};
        set_plot_ticks(visualizer, majors, minors, render_target_box,
          view_scale);
        lower(visualizer.pa_store, visualizer.pa, render_target_box,
          view_transform, view_scale, view_rotation);
      }
    }};
  auto const update_view_transform{[&](){
      std::tie(view_transform, view_transform_inv) =
        as_view_transform(view_scale, view_offset, view_rotation);
      update_time_invariant_annotations();
    }};
  auto const view_reset{[&](bool offset = true, bool scale = true,
        bool rotation = true){
      if (offset) view_offset = default_view_offset;
      if (scale) view_scale = default_view_scale;
      if (rotation) view_rotation = default_view_rotation;
      update_view_transform();
    }};
  view_reset();
  auto const view_zoom{[&](sf::Vector2f const &delta,
        array2d<float, 2, 1> const &anchor, bool const fine){
      float const factor{(fine ? .1f : 1.f) * .01f};
      auto const [r, r_inv]{as_rotation_matrix(view_rotation)};
      array2d<float, 2, 1> const r_inv_delta{
        muladd(r_inv, sf_vec_to_array2d(delta))};
      sf::Vector2f const alpha{std::exp2(get<0, 0>(r_inv_delta) * factor),
                               std::exp2(get<1, 0>(r_inv_delta) * factor)};
      array2d<float, 2, 1> const r_inv_anchor{muladd(r_inv, anchor)};
      view_offset += {
        ((get<0, 0>(r_inv_anchor) - view_offset.x) * (1.f - alpha.x)),
        ((get<1, 0>(r_inv_anchor) - view_offset.y) * (1.f - alpha.y))};
      view_scale = {view_scale.x * alpha.x, view_scale.y * alpha.y};
    }};
  auto const view_pan{[&](sf::Vector2f const &delta, bool const fine){
      float const factor{fine ? .1f : 1.f};
      auto const [r, r_inv]{as_rotation_matrix(view_rotation)};
      view_offset += array2d_to_sf_vec(
        muladd(r_inv, sf_vec_to_array2d(delta * factor)));
    }};
  auto const view_rotate{[&](std::optional<sf::Vector2f> const &delta,
        array2d<float, 2, 1> const &anchor, bool const fine){
      float const factor{(fine ? .1f : 1.f) *
        two_pi<float> / three_hundred_and_sixty<float>};
      auto const [r_pre, r_pre_inv]{as_rotation_matrix(view_rotation)};
      if (delta.has_value()) view_rotation += delta->x * factor;
      else view_rotation = 0.f;
      { // Wrap back to [0, 2π)
        view_rotation *= one_over_two_pi<float>;
        view_rotation -= std::floor(view_rotation);
        view_rotation *= two_pi<float>;
      }
      auto const [r_post, r_post_inv]{as_rotation_matrix(view_rotation)};
      view_offset += array2d_to_sf_vec(add(neg(muladd(r_pre_inv, anchor)),
        muladd(r_post_inv, anchor)));
    }};

  auto const axes_reset{[&](){
      if constexpr (std::is_base_of_v<VisualizerSingle2DFrameBase,
          decltype(visualizer)>) {
        visualizer.plot_box = default_plot_box;
        update_time_invariant_annotations();
      }
    }};
  axes_reset();

  auto const set_window_view{[&](){
      sf::Vector2f const sz{window.getSize()};
      window.setView(sf::View(sf::FloatRect(0.f, sz.y, sz.x, -sz.y)));
      update_time_invariant_annotations();
    }};
  set_window_view();

  sf::Event event{};
  sf::Vector2i mouse_pos_prev{};
  std::optional<sf::Vector2i> anchor_window{};
  sf::Vector2f anchor_window_center(window.getSize() / 2u);

  std::chrono::high_resolution_clock clock{};
  std::chrono::high_resolution_clock::time_point mouse_button_left_tic{
    std::chrono::high_resolution_clock::time_point::max()};
  std::chrono::high_resolution_clock::time_point mouse_button_left_toc{
    std::chrono::high_resolution_clock::time_point::max()};
  std::chrono::high_resolution_clock::time_point frame_tic{
    std::chrono::high_resolution_clock::time_point::max()};
  std::chrono::high_resolution_clock::time_point frame_toc{
    std::chrono::high_resolution_clock::time_point::max()};

  while (window.isOpen()) {
    frame_tic = frame_toc;
    frame_toc = std::chrono::high_resolution_clock::now();
    auto const frame_tictoc{frame_toc - frame_tic};

    std::optional<sf::Vector2f> delta{};
    auto const delta_plus{[&](float const x, float const y){
        // `-y` below to account for flipped `sf::View`
        return delta.value_or(sf::Vector2f{0.f, 0.f}) + sf::Vector2f{x, -y};
      }};
    bool reset{false};

    sf::Vector2i const mouse_pos{sf::Mouse::getPosition(window)};
    sf::Vector2i const mouse_delta{mouse_pos - mouse_pos_prev};
    bool const any_mouse_button_pressed{
      sf::Mouse::isButtonPressed(sf::Mouse::Left) or
      sf::Mouse::isButtonPressed(sf::Mouse::Right)};
    if (any_mouse_button_pressed) {
      delta = delta_plus(mouse_delta.x, mouse_delta.y);
      if (not anchor_window.has_value()) anchor_window = mouse_pos_prev;
      if (cursor_cross_supported) window.setMouseCursor(cursor_cross);
    } else {
      anchor_window = {};
      if (cursor_default_supported) window.setMouseCursor(cursor_default);
    }

    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed)
        graceful_exit_signal(true);
      else if (event.type == sf::Event::Resized) {
        set_window_view();
        anchor_window_center = static_cast<sf::Vector2f>(window.getSize() / 2u);
      } else if (event.type == sf::Event::MouseButtonPressed) {
        if (event.mouseButton.button == sf::Mouse::Button::Left) {
          mouse_button_left_tic = mouse_button_left_toc;
          mouse_button_left_toc = clock.now();
          auto const tictoc{mouse_button_left_toc - mouse_button_left_tic};
          if (tictoc >= double_click_min_duration and
              tictoc < double_click_max_duration)
            reset = true;
        }
      } else if (event.type == sf::Event::MouseWheelScrolled) {
        float const factor{-10.f};
        if (event.mouseWheelScroll.wheel == sf::Mouse::Wheel::HorizontalWheel)
          delta = delta_plus(factor * event.mouseWheelScroll.delta, 0.f);
        else
          delta = delta_plus(0.f, factor * event.mouseWheelScroll.delta);
        anchor_window = anchor_window.value_or(
          sf::Vector2i{event.mouseWheelScroll.x, event.mouseWheelScroll.y});
      } else if (event.type == sf::Event::KeyPressed) {
        if (event.key.code == key_reset)
          reset = true;
      }
    }
    if (graceful_exit_signal()) break;

    {
      float const factor{512.f};
      float const step_size{factor *
        std::chrono::duration_cast<std::chrono::duration<float>>(
          std::max(frame_tictoc, decltype(frame_tictoc)::zero())).count()};
      if (sf::Keyboard::isKeyPressed(key_left))
        delta = delta_plus(-step_size, 0.f);
      if (sf::Keyboard::isKeyPressed(key_right))
        delta = delta_plus(step_size, 0.f);
      if (sf::Keyboard::isKeyPressed(key_down))
        delta = delta_plus(0.f, step_size);
      if (sf::Keyboard::isKeyPressed(key_up))
        delta = delta_plus(0.f, -step_size);
    }

    auto const render_target_box{get_render_target_box(window)};
    if (reset) {
      if (sf::Keyboard::isKeyPressed(key_mod_view)) {
        if (not sf::Keyboard::isKeyPressed(key_mod_mode)) view_reset();
        else {
          if constexpr (std::is_base_of_v<VisualizerSingle2DFrameBase,
              decltype(visualizer)>)
            view_rotate({}, array_to_array2d(view_transform(get_frame_box(
              visualizer.pa, render_target_box, view_scale).center())),
                false);
          }
      } else
        axes_reset();
    } else if (delta.has_value()) {
      sf::Vector2f const anchor_window_or_center{anchor_window.has_value()
        ? static_cast<sf::Vector2f>(anchor_window.value())
        : anchor_window_center};
      array2d<float, 2, 1> const anchor{array_to_array2d(lerp(
        sf_vec_to_array(anchor_window_or_center),
        render_target_box, sf_view_to_box(window.getView())))};

      bool const mod_switch{sf::Mouse::isButtonPressed(mouse_button_mod_switch)
        != sf::Keyboard::isKeyPressed(key_mod_switch)};
      bool const mod_mode{sf::Keyboard::isKeyPressed(key_mod_mode)};
      bool const mod_fine{sf::Keyboard::isKeyPressed(key_mod_fine)};

      if (sf::Keyboard::isKeyPressed(key_mod_view)) {
        if (not mod_switch) {
          if (not mod_mode)
            view_pan(delta.value(), mod_fine);
          else
            view_rotate(delta.value(), anchor, mod_fine);
        } else {
          if (not mod_mode)
            view_zoom(delta.value(), anchor, mod_fine);
          else
            ; // Do nothing, for now
        }
      } else {
        // TODO: Make plot axes adjustments `view_transform`-aware
        //auto const [r, r_inv]{as_rotation_matrix(view_rotation)};
        //auto const view_anchor{view_transform_inv(anchor)};
        //auto const view_delta{array2d_to_sf_vec(muladd(view_transform_inv.a,
        //  sf_vec_to_array2d(delta.value())))};
        auto const view_anchor{anchor};
        auto const view_delta{delta.value()};
        if (not mod_switch) {
          if (not mod_mode)
            axes_pan(visualizer, view_delta, mod_fine, render_target_box,
              view_scale);
          else
            ; // Do nothing, for now
        } else {
          if (not mod_mode)
            axes_zoom(visualizer, view_delta, view_anchor, mod_fine,
              render_target_box, view_scale);
          else
            ; // Do nothing, for now
        }
      }
      update_view_transform();
    }


    mouse_pos_prev = mouse_pos;

    window.clear(background_color);
    double const nominal_frame_rate{60.};
    if (not use_compute_thread) compute_step(nominal_frame_rate);
    draw(visualizer, compute_draw_mutex, buffer, window, view_transform,
        view_scale);
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
