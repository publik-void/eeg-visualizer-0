#pragma once

#include <optional>
#include <array>
#include <vector>
#include <unordered_map>

namespace args {

using arg_t = char const *;
using args_t = std::vector<arg_t>;
using mode_t = char const *;
template <std::size_t n> using modes_t = std::array<mode_t, n>;
using key_t = arg_t;
using flag_t = bool;
using opt_t = std::optional<arg_t>;
using flags_t = std::unordered_map<key_t, flag_t>;
using opts_t = std::unordered_map<key_t, opt_t>;

args_t make_args(int const, char const * const []);
std::size_t parse_opts_and_flags(flags_t &, opts_t &,
    args_t::const_iterator &, args_t::const_iterator const &);

bool parse_bool(arg_t const &);
float parse_float(arg_t const &);
double parse_double(arg_t const &);
long long parse_long_long(arg_t const &);
unsigned long long parse_unsigned_long_long(arg_t const &);

auto parse_arg(auto &&parser, arg_t &arg, auto const &default_value) {
  try {
    return parser(arg);
  } catch (std::exception const &e) {
    return default_value;
  }
}

auto parse_opt(auto &&parser, opt_t &opt, auto const &default_value) {
  if (not opt.has_value()) return default_value;
  return parse_arg(parser, opt.value(), default_value);
}

// As a `constexpr`-qualified replacement for `std::strcmp`
bool constexpr is_equal(char const * const a, char const * const b) {
  //return std::string_view{a} == std::string_view{b};
  return *a == *b && (*a == '\0' || is_equal(a + 1, b + 1));
}

template <std::size_t n>
std::size_t constexpr get_mode_index(modes_t<n> const &modes,
    mode_t const &mode, std::size_t const error_index = -1) {
  // NOTE: Could be done with `std::find_if`, but would need heavy `<algorithm>`
  for (std::size_t i{0}; i < n; ++i)
    if (is_equal(modes[i], mode)) return i;
  return error_index;
}

template <std::size_t n>
std::size_t get_mode_index(modes_t<n> const &modes,
    args_t::const_iterator &pos, args_t::const_iterator const &end,
    std::size_t const error_index = -2) {
  if (pos >= end) return error_index;
  mode_t const mode{*(pos++)};
  return get_mode_index(modes, mode, error_index);
}

} // namespace args

