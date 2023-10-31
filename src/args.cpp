#include <cstddef>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "args.hpp"

namespace args {

args_t make_args(int const argc, char const * const argv[]) {
  args_t args{}; args.reserve(argc);
  for (int i{0}; i < argc; ++i) args.emplace_back(argv[i]);
  return args;
}

// `constexpr` would be better, but `std::string` does not support it (yet).
std::size_t get_subcommand_index(
    std::vector<std::string> const &subcommands, std::string const &subcommand,
    std::size_t const error_index) {
  auto const index{
    std::find(std::cbegin(subcommands), std::cend(subcommands), subcommand)};
  return index == std::cend(subcommands)
    ? error_index
    : index - std::cbegin(subcommands);
}

std::size_t get_subcommand_index(std::vector<std::string> const &subcommands,
    args_t::const_iterator &pos, args_t::const_iterator const &end,
    std::size_t const error_index) {
  if (pos >= end) return error_index;
  std::string const subcommand{*(pos++)};
  return get_subcommand_index(subcommands, subcommand, error_index);
}

std::size_t parse_opts_and_flags(flags_t &flags, opts_t &opts,
    args_t::const_iterator &pos, args_t::const_iterator const &end) {
  char constexpr const * const opt_prefix{"--"};
  std::size_t constexpr opt_prefix_size{2};
  auto const original_position{pos};
  while (pos != end) {
    auto const arg{*pos};
    if (not arg.starts_with(opt_prefix)) break;
    if (arg == opt_prefix) { ++pos; break; }
    auto const equal_sign_pos{arg.find_first_of('=')};
    auto const arg_key{
      arg.substr(opt_prefix_size, equal_sign_pos - opt_prefix_size)};
    if (equal_sign_pos == arg_t::npos) { // arg is a flag
      auto flag_itr = flags.find(arg_key);
      if (flag_itr == std::end(flags)) break;
      flag_itr->second = true;
    } else { // arg is an opt
      auto opt_itr = opts.find(arg_key);
      if (opt_itr == std::end(opts)) break;
      opt_itr->second = arg.substr(equal_sign_pos + 1);
    }
    ++pos;
  }
  return original_position - pos;
}

} // namespace args

