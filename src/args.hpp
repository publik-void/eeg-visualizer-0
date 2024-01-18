#pragma once

#include <cassert>
#include <tuple>
#include <optional>
#include <array>
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>

// Utilities for option parsing and usage messages.
// NOTE: Boost has a library that does this, but unfortunately, it looks to me
// like it is not designed with subcommands in mind. So once again, I'm
// reinventing the wheel. I should probably get into the habit of looking for
// good libraries instead of coding ahead immediately…
namespace args {

using arg_raw_t = char const *;
using arg_t = std::string;
using args_t = std::vector<arg_t>;
using mode_t = char const *;
template <std::size_t n> using modes_t = std::array<mode_t, n>;
using key_t = arg_t;
using flag_t = bool;
using optvals_t = std::vector<arg_t>;
using flags_t = std::unordered_map<key_t, flag_t>;
using opts_t = std::unordered_map<key_t, optvals_t>;
using desc_t = char const *;

// {key, description}
using flag_descs_t = std::unordered_map<key_t, desc_t>;

// {key, {variable, default, description, {{value, description}, ...}}}
using opt_descs_t = std::unordered_map<key_t, std::tuple<arg_raw_t, arg_raw_t,
  desc_t, std::vector<std::tuple<arg_raw_t, desc_t>>>>;

// {key, {alias_name, expanded_opts}}
// Alias keys must be a subset of the option keys. But it is okay if the option
// is not used for anything except aliases to other options.
using opt_aliases_t = std::unordered_map<key_t,
  std::unordered_map<arg_t, opts_t>>;

char constexpr const * const opt_prefix{"--"};
std::size_t constexpr opt_prefix_size{2u};
char constexpr opt_assignment{'='};
std::size_t constexpr opt_assignment_size{1u};
desc_t constexpr ind{"  "}, ind2{"    "}, ind3{"      "},
  normal{"\u001b[0m"}, bold{"\u001b[1m"};

args_t make_args(int const, char const * const []);
std::size_t parse_opts_and_flags(flags_t &, opts_t &,
    args_t::const_iterator &, args_t::const_iterator const &,
    opt_aliases_t const & = {});

bool parse_bool(arg_t const &);
float parse_float(arg_t const &);
double parse_double(arg_t const &);
int parse_int(arg_t const &);
long parse_long(arg_t const &);
long long parse_long_long(arg_t const &);
unsigned int parse_unsigned_int(arg_t const &);
unsigned long parse_unsigned_long(arg_t const &);
unsigned long long parse_unsigned_long_long(arg_t const &);
std::string parse_string(arg_t const &);

template <typename T>
std::vector<T> parse_optvals(T (*parser)(arg_t const &), opts_t &opts,
    key_t const &key, std::vector<T> const &defaults = {}) {
  optvals_t const &optvals{opts[key]};
  if (optvals.empty()) return defaults;
  std::vector<T> xs; xs.reserve(optvals.size());
  for (arg_t const &optval : optvals) {
    try {
      xs.push_back(parser(optval));
    } catch (std::exception const &e) {
      std::cerr << normal << "Failed to parse value " << bold << optval
        << normal << " for option " << bold << key << normal << " (" << e.what()
        << ")." << std::endl;
    }
  }
  return xs;
}

template <typename T>
T parse_last_optval(T (*parser)(arg_t const &), opts_t &opts, key_t const &key,
    std::optional<T> const &default_opt = {}) {
  auto xs{parse_optvals(parser, opts, key)};
  if (xs.empty()) {
    if (default_opt.has_value())
      return *default_opt;
    else
      throw(std::invalid_argument("No value for option \"" + key + "\"."));
  } else return xs.back();
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
  mode_t const mode{(*(pos++)).c_str()};
  auto const mode_index{get_mode_index(modes, mode, error_index)};
  assert(0 <= mode_index and mode_index < n);
  return mode_index;
}

enum struct UsageMessageMode { help, error };

int print_usage(arg_t const &, args_t const &, desc_t, mode_t const *,
    mode_t const * const, desc_t const *, desc_t const * const,
    flag_descs_t const &, opt_descs_t const &, opt_aliases_t const & = {},
    UsageMessageMode const = UsageMessageMode::error);

flags_t init_flags(flag_descs_t const &);
opts_t init_opts(opt_descs_t const &);

} // namespace args
