#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>

//#include <iostream>

#include "args.hpp"

namespace args {

args_t make_args(int const argc, char const * const argv[]) {
  args_t args{}; args.reserve(argc);
  for (int i{0}; i < argc; ++i) args.emplace_back(argv[i]);
  return args;
}

std::size_t parse_opts_and_flags(flags_t &flags, opts_t &opts,
    args_t::const_iterator &pos, args_t::const_iterator const &end) {
  char constexpr const * const opt_prefix{"--"};
  std::size_t const opt_prefix_size{std::strlen(opt_prefix)};
  auto const original_position{pos};
  while (pos != end) {
    auto const arg{std::string{*pos}};
    if (not arg.starts_with(opt_prefix)) break;
    if (arg == opt_prefix) { ++pos; break; }
    auto const equal_sign_pos{arg.find_first_of('=')};
    auto const arg_key{
      arg.substr(opt_prefix_size, equal_sign_pos - opt_prefix_size)};
    auto const is_keyed{[&](auto const &x){
      return not std::strcmp(x.first, arg_key.c_str()); }};
    if (equal_sign_pos == decltype(arg)::npos) { // arg is a flag
      auto flag_itr{std::find_if(std::begin(flags), std::end(flags), is_keyed)};
      if (flag_itr == std::end(flags)) break;
      flag_itr->second = true;
    } else { // arg is an opt
      auto opt_itr{std::find_if(std::begin(opts), std::end(opts), is_keyed)};
      if (opt_itr == std::end(opts)) break;
      opt_itr->second = (*pos) + equal_sign_pos + 1;
    }
    ++pos;
  }
  return original_position - pos;
}

bool parse_bool(arg_t const &str){
  for (auto const &token : {"true", "1", "on", "yes"})
    if (not std::strcmp(str, token)) return true;
  for (auto const &token : {"false", "0", "off", "no"})
    if (not std::strcmp(str, token)) return false;
  throw(std::invalid_argument(
    std::string{"\""} + str + "\" could not be parsed as bool"));
}

float parse_float(arg_t const &str){ return std::strtof(str, nullptr); }
double parse_double(arg_t const &str){ return std::strtod(str, nullptr); }
long long parse_long_long(arg_t const &str){
  return std::strtoll(str, nullptr, 10); }
unsigned long long parse_unsigned_long_long(arg_t const &str){
  return std::strtoull(str, nullptr, 10); }

} // namespace args

