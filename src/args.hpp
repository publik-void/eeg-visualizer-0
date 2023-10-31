#pragma once

#include <string>
#include <optional>
#include <vector>
#include <unordered_map>

namespace args {

using arg_t = std::string;
using args_t = std::vector<arg_t>;
using key_t = arg_t;
using flag_t = bool;
using opt_t = std::optional<arg_t>;
using flags_t = std::unordered_map<key_t, flag_t>;
using opts_t = std::unordered_map<key_t, opt_t>;

args_t make_args(int const, char const * const []);
std::size_t get_subcommand_index(
    std::vector<std::string> const &, std::string const &,
    std::size_t const = -1);
std::size_t get_subcommand_index(std::vector<std::string> const &,
    args_t::const_iterator &, args_t::const_iterator const &,
    std::size_t const = -2);
std::size_t parse_opts_and_flags(flags_t &, opts_t &,
    args_t::const_iterator &, args_t::const_iterator const &);

} // namespace args

