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
#include "fftw3.h"

#include "signal-handling.hpp"
#include "args.hpp"
#include "dsp.hpp"
#include "lsl-customized.hpp"
#include "numbers-and-math.hpp"

args::modes_t<4> constexpr visual_modes{"error", "help", "time-series",
  "short-time-spectrum"};

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

  std::size_t pull() {
      float * const buf{this->data.data()};
      std::size_t const size{this->data.size()};
      std::size_t const k{this->n_channels()};
      std::size_t i_in_k{this->i_in.load() * k};

      std::size_t n_pulled_k{this->inlet.pull_chunk_multiplexed(
        buf + i_in_k, nullptr, size - i_in_k, 0)};
      i_in_k += n_pulled_k;
      if (i_in_k >= size) {
        i_in_k = this->inlet.pull_chunk_multiplexed(buf, nullptr, size, 0);
        n_pulled_k += i_in_k;
      }

      std::size_t const i_in{i_in_k / k};
      this->n_pulled = n_pulled_k / k;
      this->i_in = i_in;
      this->i_out = i_in == 0 ? this->n_timepoints - 1 : i_in - 1;
      return this->n_pulled;
    }

    inline float & operator[](std::array<std::size_t, 2> const ji) {
    return this->data[get_scalar_index(*this, ji)];
  }

  inline float const & operator[](std::array<std::size_t, 2> const ji) const {
    return this->data[get_scalar_index(*this, ji)];
  }
};

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

struct VisualizerBase {
  virtual void compute(Buffer &) {}
  virtual void draw(Buffer const &, sf::RenderWindow &) { }

  virtual ~VisualizerBase() {}
};

template<std::size_t mode>
struct Visualizer : public VisualizerBase {};

template<>
struct Visualizer<args::get_mode_index(visual_modes, "time-series")> :
    public VisualizerBase {
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

  void draw(Buffer const &buffer, sf::RenderWindow &window) override {
    std::size_t const n_timepoints{
      std::min(this->n_timepoints(), buffer.n_timepoints)};
    std::size_t const n_channels{
      std::min(this->n_channels(), buffer.n_channels())};
    for (std::size_t i{0}; i < n_timepoints; ++i)
      for (std::size_t j{0}; j < n_channels; ++j) {
        float const s{this->scrolling
          ? buffer[{j, n_timepoints - 1 - i}]
          : buffer.data[i * n_channels + j]};
        float const x{i / static_cast<float>(n_timepoints - 1)};
        float const y{(s * (.5f * this->y_scale) + j + .5f) / n_channels};
        this->vertexess[j][i].position = sf::Vector2f{x, y};
      }
    for (auto const &vertexes : vertexess) window.draw(vertexes);
  }
};

template<>
struct Visualizer<args::get_mode_index(visual_modes, "short-time-spectrum")> :
    public VisualizerBase {
  // Parameters
  std::size_t const n_in;
  float x_min, x_max;
  float y_min_db, y_max_db;

  // State
  float * window, * in, * aggregate;
  fftwf_complex * out;
  fftwf_plan plan;
  sf::VertexArray vertexes;

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
      float const y_min_db, float const y_max_db) :
      n_in{n_in}, x_min{x_min}, x_max{x_max},
      y_min_db{y_min_db}, y_max_db{y_max_db},
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
      vertexes{[&](){
        sf::VertexArray vertexes(sf::LineStrip, this->n_out() - 1);
        for (std::size_t i{0}; i < vertexes.getVertexCount(); ++i)
          vertexes[i].color = color;
        return vertexes; }()} {}

  ~Visualizer() {
    fftwf_destroy_plan(this->plan);
    fftwf_free(this->in);
    fftwf_free(this->out);
    fftwf_free(this->aggregate);
  }

  virtual void compute(Buffer &buffer) override {
    std::fill(this->aggregate, this->aggregate + this->n_out(), 0.f);
    for (std::size_t j{0}; j < buffer.n_channels(); ++j) {
      for (std::size_t i{0}; i < this->n_in; ++i)
        this->in[i] = this->window[i] * buffer[{j, this->n_in - 1 - i}];
      fftwf_execute(this->plan);
      std::transform(this->out_std(), this->out_std() + this->n_out(),
        this->aggregate, this->aggregate,
        [](auto const &out, auto const &agg){
          auto const abs2{[](auto const &x){ return x * x; }};
          return abs2(std::real(out)) + abs2(std::imag(out)) + agg; });
    }
    std::transform(this->aggregate, this->aggregate + this->n_out(),
      this->aggregate, [&](auto const x){
        return 6.f * std::log2(std::sqrt(x * (1.f / buffer.n_channels()))); });
  }

  void draw(Buffer const &buffer, sf::RenderWindow &window) override {
    for (std::size_t i{1}; i < this->n_out(); ++i) {
      float const x{(std::log2(this->f_min(buffer) + this->f_range(buffer) *
          static_cast<float>(i - 1) / (this->n_out() - 2)) -
        this->x_min_log2()) / this->x_range_log2()};
      float const y{(this->aggregate[i] - this->y_min_db) / this->y_range_db()};
      this->vertexes[i - 1].position = sf::Vector2f{x, y};
    }
    window.draw(this->vertexes);
  }
};

int main_client_visual(args::args_t::const_iterator &pos,
    args::args_t::const_iterator const &end) {
  args::flags_t flags{{"use-compute-thread", false}};
  args::opts_t opts{{"frame-rate-limit", {}}};
  args::parse_opts_and_flags(flags, opts, pos, end);

  // Parameters
  std::string const property{"type"};
  std::string const value{"EEG"};
  double const timeout_resolve{5.};
  double const timeout_open{5.};
  double const buffer_duration{10.};

  bool const use_compute_thread{flags["use-compute-thread"]};
  double const update_rate{60.};

  std::optional<std::size_t> const ref_index_opt{};
  double const highpass_frequency{.25};
  double const highpass_quality_factor{.7071};

  std::string const window_title{"Visualizer"};
  auto const window_style{sf::Style::Resize | sf::Style::Close};
  sf::VideoMode const window_video_mode{
    1024, // width
    768,  // height
    32    // bits per pixel
  };
  sf::ContextSettings const window_context_settings{
    0, // depth bits
    0, // stencil bits
    8  // antialiasing level
  };
  unsigned int const frame_rate_limit(args::parse_opt(
    args::parse_unsigned_long_long, opts["frame-rate-limit"], 0ull));
  bool const vertical_sync{true};

  std::size_t visual_mode{
    args::get_mode_index(visual_modes, pos, end, 0)};

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
  std::unique_ptr<VisualizerBase> visualizer{
    (visual_mode == args::get_mode_index(visual_modes, "time-series"))
    ? new Visualizer<args::get_mode_index(visual_modes, "time-series")>(
      buffer.n_channels(), buffer.n_timepoints,
      sf::Color{0xff, 0xff, 0xaa, 0xaa}, 1.f, true)
    : (visual_mode == args::get_mode_index(visual_modes, "short-time-spectrum"))
    ? new Visualizer<args::get_mode_index(visual_modes, "short-time-spectrum")>(
      2048, [](std::size_t const i, std::size_t const n){
        return dsp::window<float, dsp::WindowShape::blackman>(i, n, .16f); },
      sf::Color{0xff, 0xff, 0xaa, 0xaa}, .5f, 80.f, -24.f, 48.f)
    : new VisualizerBase{}};

  // Set up computation thread
  auto const compute_step{[&](auto const /*update_rate*/){
      buffer.pull();

      { // Preprocessing (TODO: Abstract this)
        for (auto sample{get_sample(buffer)};
            not past_newest_sample(buffer, sample);
            advance_sample(buffer, sample)) {
          auto const &[first, last] = sample;
          if (ref_index_opt.has_value()) {
            float const x_ref{*(first + ref_index_opt.value())};
            std::transform(first, last, first, [=](float const x){
              return x - x_ref; });
            }
          dsp::tick(first, last, first, highpass_spec, highpass_state);
        }
      }

      visualizer->compute(buffer);
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
  sf::RenderWindow window{window_video_mode, window_title, window_style,
    window_context_settings};
  window.setFramerateLimit(frame_rate_limit);
  window.setVerticalSyncEnabled(vertical_sync);
  window.setView(sf::View(sf::Vector2f{.5f, .5f}, sf::Vector2f{1.f, -1.f}));
  sf::Event event{};

  while (window.isOpen()) {
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) graceful_exit_signal(true);
    }
    if (graceful_exit_signal()) break;

    window.clear();
    double const nominal_frame_rate{60.};
    if (not use_compute_thread) compute_step(nominal_frame_rate);
    visualizer->draw(buffer, window);
    window.display();
  }

  window.close();
  if (use_compute_thread && t_compute.joinable()) t_compute.join();
  return 0;
}
