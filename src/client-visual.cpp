#include <cstddef>
#include <cmath>
#include <vector>
//#include <iostream>

#include "lsl_cpp.h"
#include "SFML/Graphics.hpp"

#include "signal-handling.hpp"

// `constexpr` would be better, but `std::string` does not support it (yet).
std::vector<std::string> const visual_modes{"error", "help", "time-series"};

int main_client_visual(std::vector<std::string>::const_iterator pos,
    std::vector<std::string>::const_iterator const end) {
  // NOTE: I first tried to implement a threaded design, but along the way I
  // realized that it adds a lot of complexity, but only limited value.

  std::string const property{"type"};
  std::string const value{"EEG"};
  double timeout_open{5.};
  double timeout_resolve{5.};
  double const buffer_duration{10.};

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
  unsigned int const frame_rate_limit{60};
  bool vertical_sync{true};

  // Set up LSL input
  auto const streams{lsl::resolve_stream(property, value, 1, timeout_resolve)};
  if (graceful_exit_signal()) return 0;
  if (std::empty(streams)) return 1;

  auto const stream{streams.front()};
  auto const sample_rate{stream.nominal_srate()};
  if (sample_rate <= 0) return 1;

  lsl::stream_inlet inlet{stream};
  std::size_t const n_channels{static_cast<size_t>(inlet.get_channel_count())};

  std::size_t const n_timepoints{
    static_cast<std::size_t>(std::round(sample_rate * buffer_duration))};
  std::vector<float> buffer(n_channels * n_timepoints, 0.f);
  std::size_t i_marker{0}, n_pulled{0};
  auto const pull{[&](){
      std::size_t const n{inlet.pull_chunk_multiplexed(
        buffer.data() + i_marker, nullptr, buffer.size() - i_marker, 0)};
      i_marker += n;
      return n;
    }};

  try {
    inlet.open_stream(timeout_open);
  } catch (lsl::timeout_error const &e) {
    return 1;
  }

  // Set up view contents
  float const y_scale{1.f};
  sf::Color color{0xff, 0xff, 0xcc, 0xcc};
  std::vector<sf::VertexArray> vertexess(n_channels,
    sf::VertexArray{sf::LineStrip, n_timepoints});
  for (auto &vertexes : vertexess)
    for (std::size_t i{0}; i < n_timepoints; ++i)
      vertexes[i].color = color;

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

    n_pulled = pull();
    if (i_marker >= buffer.size()) { i_marker = 0; n_pulled += pull(); }

    for (std::size_t i{0}; i < n_timepoints; ++i)
      for (std::size_t j{0}; j < n_channels; ++j)
        vertexess[j][i].position = sf::Vector2f{
          i / static_cast<float>(n_timepoints - 1),
          (buffer[([&](auto const k){
              return k >= buffer.size() ? k - buffer.size() : k;
                }(i * n_channels + i_marker + j))] * (.5f * y_scale) +
            j + .5f) / n_channels};

    window.clear();
    for (auto const &vertexes : vertexess) window.draw(vertexes);
    window.display();
  }

  window.close();
  return 0;
}
