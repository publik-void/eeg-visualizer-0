#include <utility>
#include <variant>
#include <string>
#include <vector>

#include "lsl_cpp.h"

#include "signal-handling.hpp"
#include "lsl-customized.hpp"

namespace lsl {

std::variant<attachment_t, int> attach(
    std::string const &property, std::string const &value,
    double const timeout_resolve, double const timeout_open) {
  auto const streams{resolve_stream(property, value, 1, timeout_resolve)};
  if (graceful_exit_signal()) return 0;
  if (std::empty(streams)) return 1;

  auto const stream{streams.front()};
  auto const sample_rate{stream.nominal_srate()};
  if (sample_rate <= 0) return 1;

  stream_inlet inlet{stream};
  try {
    inlet.open_stream(timeout_open);
  } catch (timeout_error const &e) {
    return 1;
  }
  if (graceful_exit_signal()) return 0;
  return std::make_pair(std::move(stream), std::move(inlet));
}

bool has_failed(std::variant<attachment_t, int> const &attachment_var) {
  return not std::holds_alternative<attachment_t>(attachment_var);
}

int exit_code(std::variant<attachment_t, int> const &attachment_var) {
  if (std::holds_alternative<int>(attachment_var))
    return std::get<int>(attachment_var);
  return 1;
}

attachment_t & get_attachment(attachment_t &attachment) {
  return attachment;
}

stream_info & get_stream(attachment_t &attachment) {
  return std::get<stream_info>(attachment);
}

stream_inlet & get_inlet(attachment_t &attachment) {
  return std::get<stream_inlet>(attachment);
}

attachment_t & get_attachment(std::variant<attachment_t, int> &attachment_var) {
  return std::get<attachment_t>(attachment_var);
}

stream_info & get_stream(std::variant<attachment_t, int> &attachment_var) {
  return std::get<stream_info>(get_attachment(attachment_var));
}

stream_inlet & get_inlet(std::variant<attachment_t, int> &attachment_var) {
  return std::get<stream_inlet>(get_attachment(attachment_var));
}

} // namespace lsl
