#include <cassert>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "args.hpp"

namespace args {

args_t make_args(int const argc, char const * const argv[]) {
  args_t args{}; args.reserve(argc);
  if (argv != nullptr)
    for (int i{0}; i < argc; ++i) {
      if (argv[i] == nullptr) break;
      args.emplace_back(argv[i]);
    }
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

template<typename To, typename From>
To safe_narrow(From const &x) {
  To min{std::numeric_limits<To>::min()};
  To max{std::numeric_limits<To>::max()};
  bool const underflow{x < min};
  bool const overflow{x > max};
  if (underflow || overflow) {
    auto what_str{std::to_string(x) + " is outside the permissible range [" +
      std::to_string(min) + ", " + std::to_string(max) + "]"};
    throw(std::underflow_error{what_str});
  }
  return static_cast<To>(x);
}

float parse_float(arg_t const &str){ return std::strtof(str, nullptr); }
double parse_double(arg_t const &str){ return std::strtod(str, nullptr); }
int parse_int(arg_t const &str){
  return safe_narrow<int>(std::strtol(str, nullptr, 10)); }
long parse_long(arg_t const &str){
  return std::strtol(str, nullptr, 10); }
long long parse_long_long(arg_t const &str){
  return std::strtoll(str, nullptr, 10); }
unsigned int parse_unsigned_int(arg_t const &str){
  return safe_narrow<unsigned int>(std::strtoul(str, nullptr, 10)); }
unsigned long parse_unsigned_long(arg_t const &str){
  return std::strtoul(str, nullptr, 10); }
unsigned long long parse_unsigned_long_long(arg_t const &str){
  return std::strtoull(str, nullptr, 10); }

bool constexpr is_empty(char const * const str) {
  return str == nullptr or *str == '\0';
}

int print_usage(arg_t const &program_name,
  args_t const &current_submodes, desc_t const desc,
  mode_t const * submodes_itr, mode_t const * const submodes_end,
  desc_t const * submodes_descs_itr, desc_t const * const submodes_descs_end,
  flag_descs_t const &flag_descs, opt_descs_t const &opt_descs,
  UsageMessageMode const mode) {
  std::ostream &out = mode == UsageMessageMode::help ? std::cout : std::cerr;
  desc_t constexpr ind{"  "}, ind2{"    "}, ind3{"      "},
    normal{"\u001b[0m"}, bold{"\u001b[1m"};

  out << normal << "Usage:\n";
  out << ind << bold << program_name;
  for (auto const &mode : current_submodes) out << " [--arg…] \\\n" << ind2
    << mode;
  for (auto const &[flag, desc] : flag_descs) out << " \\\n" << ind3 << "[--"
    << flag << "]";
  for (auto const &[opt, descs] : opt_descs) {
    auto const &[var, def, desc, val_descs]{descs};
    out << " \\\n" << ind3 << "[--" << opt << "=<";
    if (is_empty(var)) {
      if (val_descs.empty()) out << "…"; else {
        bool is_first_entry{true};
        for (auto const &[val, desc] : val_descs) {
          if (is_first_entry) is_first_entry = false; else out << "|";
          out << val;
        }
      }
    } else out << var;
    out << ">]";
  }
  out << " \\\n" << ind3 << "[--]";
  if (submodes_itr != nullptr and submodes_itr != submodes_end)
    out << " \\\n" << ind2 << "<mode> [--arg…]";
  out << normal << "\n";

  if (not is_empty(desc)) {
    out << "\n";
    out << ind << desc << "\n";
  }

  if (submodes_itr != nullptr and submodes_itr != submodes_end) {
    out << "\n" << "Modes:\n";
    bool is_first_entry{true};
    while (submodes_itr != submodes_end and
        submodes_descs_itr != submodes_descs_end) {
      assert(not is_empty(*submodes_itr));
      if (is_first_entry) is_first_entry = false; else out << "\n";
      out << ind << bold << *submodes_itr << normal << "\n";
      if (not std::strcmp(*submodes_itr, "help"))
        out << ind2 << "Print this usage message.\n";
      if (not is_empty(*submodes_descs_itr))
        out << ind2 << *submodes_descs_itr << "\n";
      ++submodes_itr; ++submodes_descs_itr;
    }
  }

  if (not flag_descs.empty()) {
    out << "\n" << "Flags:\n";
    bool is_first_entry{true};
    for (auto const &[flag, desc] : flag_descs) {
      if (is_first_entry) is_first_entry = false; else out << "\n";
      out << ind << bold << "--" << flag << normal << "\n";
      if (not std::strcmp(flag, "help"))
        out << ind2 << "Print this usage message.\n";
      if (not is_empty(desc)) out << ind2 << desc << "\n";
    }
  }

  if (not opt_descs.empty()) {
    out << "\n" << "Options:\n";
    bool is_first_entry{true};
    for (auto const &[opt, descs] : opt_descs) {
      auto const &[var, def, desc, val_descs]{descs};
      if (is_first_entry) is_first_entry = false; else out << "\n";
      out << ind << bold << "--" << opt << "=<";
      if (not is_empty(var)) out << var; else out << "…";
      out << ">" << normal;
      if (not is_empty(def)) out << " (default: " << def << ")";
      out << "\n";
      if (not is_empty(desc)) out << ind2 << desc << "\n";

      if (not val_descs.empty()) out << "\n";
      for (auto const &[val, desc] : val_descs) {
        out << ind2 << bold << "--" << opt << "=" << val << normal <<
          ": ";
        if (not is_empty(desc)) out << desc << "\n";
      }
    }
  }

  out << normal << std::flush;
  return mode == UsageMessageMode::help ? 0 : 1;
}

flags_t init_flags(flag_descs_t const &flag_descs) {
  flags_t flags{};
  for (auto const &[k, v] : flag_descs)
    flags.emplace(std::pair<key_t, flag_t>{k, false});
  return flags;
}

opts_t init_opts(opt_descs_t const &opt_descs) {
  opts_t opts{};
  for (auto const &[key, descs] : opt_descs) {
    auto const &[var, def, desc, val_descs]{descs};
    opts.emplace(
      std::pair<key_t, opt_t>{key, is_empty(def) ? opt_t{} : opt_t{def}});
  }
  return opts;
}

} // namespace args

