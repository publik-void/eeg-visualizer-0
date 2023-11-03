#pragma once

#include <utility>
#include <variant>
#include <string>

#include "lsl_cpp.h"

namespace lsl {

using attachment_t = std::pair<lsl::stream_info, lsl::stream_inlet>;

std::variant<attachment_t, int> attach(
    std::string const &, std::string const &, double const, double const);

bool has_failed(std::variant<attachment_t, int> const &);
int exit_code(std::variant<attachment_t, int> const &);

attachment_t & get_attachment(attachment_t &);
stream_info & get_stream(attachment_t &);
stream_inlet & get_inlet(attachment_t &);
attachment_t & get_attachment(std::variant<attachment_t, int> &);
stream_info & get_stream(std::variant<attachment_t, int> &);
stream_inlet & get_inlet(std::variant<attachment_t, int> &);

} // namespace lsl
