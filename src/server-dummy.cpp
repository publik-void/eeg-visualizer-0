#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <chrono>

#include "signal-handling.hpp"
#include "args.hpp"
#include "lsl-customized.hpp"

args::desc_t constexpr desc_dummy{
  "Send random data."};

args::flag_descs_t const flag_descs{
  {"print-data", "Print generated data to stdout."}};

args::opt_descs_t const opt_descs{
  {"sample-rate", {"f", "1000", "Generate data at a sampling frequency of <f> "
    "Hz.", {}}},
  {"amplitude", {"a", "1", "Scale the generated data by <a>.", {}}},
  {"dc-offset", {"b", "0", "Add <b> to the generated data.", {}}}};

int print_usage(args::arg_t const &program_name,
    args::UsageMessageMode const mode) {
  return args::print_usage(program_name, {"server", "dummy"}, desc_dummy,
    nullptr, nullptr, nullptr, nullptr, flag_descs, opt_descs, {}, mode);
}

int main_server_dummy(args::args_t::const_iterator &pos,
    args::args_t::const_iterator const &end, args::arg_t const &program_name) {
  auto flags{args::init_flags(flag_descs)};
  auto opts{args::init_opts(opt_descs)};
  args::parse_opts_and_flags(flags, opts, pos, end);
  if (flags["help"] or pos != end) return print_usage(
    program_name, flags["help"] ? args::UsageMessageMode::help :
      args::UsageMessageMode::error);

  std::string const name{"lsl-server-dummy"}, type{"EEG"},
    source_id{"lsl-server-dummy-id"};
  lsl::channel_format_t channel_format = lsl::cf_float32;
  double sample_rate{
    args::parse_last_optval(args::parse_double, opts, "sample-rate")};
  std::int32_t const n_channels = 8;

  lsl::stream_info info{name, type, n_channels, sample_rate, channel_format,
    source_id};

  lsl::xml_element info_channels = info.desc().append_child("channels");
  for (std::size_t i{0}; i < n_channels; ++i)
    info_channels.append_child("channel")
      .append_child_value("label", "channel-" + std::to_string(i))
      .append_child_value("unit", "microvolts")
      .append_child_value("type", type);

  lsl::stream_outlet outlet{info};

  bool const print_data{flags["print-data"]};
  if (print_data) std::cout.precision(2);

  std::vector sampling_buffer(n_channels, 0.f);
  std::random_device rd;
  std::mt19937 rng{rd()};
  std::uniform_real_distribution dist{-1.f, 1.f};
  float amplitude{
    args::parse_last_optval(args::parse_float, opts, "amplitude")};
  float dc_offset{
    args::parse_last_optval(args::parse_float, opts, "dc-offset")};

  auto const sampling_interval{
    std::chrono::duration_cast<std::chrono::high_resolution_clock::duration>(
      std::chrono::duration<decltype(sample_rate)>{1 / sample_rate})};
  auto t{std::chrono::high_resolution_clock::now()};

  while (true) {
    for (auto &x : sampling_buffer) x = dist(rng) * amplitude + dc_offset;

    if (wait_until(t += sampling_interval)) break;
    outlet.push_sample(sampling_buffer);

    if (print_data) for (auto const &x : sampling_buffer) std::cout
      << x << (&x == &sampling_buffer.back() ? "\n" : "\t") << std::flush;
  }

  return 0;
}

