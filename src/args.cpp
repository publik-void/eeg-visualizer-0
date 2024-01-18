#include <cassert>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "args.hpp"

namespace args {

args_t make_args(int const argc, arg_raw_t const argv[]) {
  args_t args{}; args.reserve(argc);
  if (argv != nullptr)
    for (int i{0}; i < argc; ++i) {
      if (argv[i] == nullptr) break;
      args.emplace_back(argv[i]);
    }
  return args;
}

args_t make_args(flags_t const &flags, opts_t const &opts) {
  args_t args{};
  for (auto const &[k, v] : flags) if (v) args.emplace_back(opt_prefix + k);
  for (auto const &[k, optvals] : opts)
    for (arg_t const &optval : optvals)
      args.emplace_back(opt_prefix + k + opt_assignment + optval);
  return args;
}

bool expand(key_t const &key, arg_t const &val, opts_t &opts,
    opt_aliases_t const &opt_aliases,
    std::unordered_set<key_t> &opts_already_seen) {
  assert(opts.contains(key));

  if (not opts_already_seen.contains(key)) {
    opts[key] = {};
    opts_already_seen.insert(key);
  }

  bool is_alias{false};
  auto const opt_aliases_itr{opt_aliases.find(key)};
  if (opt_aliases_itr != std::cend(opt_aliases)) {
    auto const &aliases{opt_aliases_itr->second};
    auto const aliases_itr{aliases.find(val)};
    if (aliases_itr != std::cend(aliases)) {
      auto const &opts_ex{aliases_itr->second};
      is_alias = true;
      for (auto const &[key_ex, vals_ex] : opts_ex)
        for (auto const &val_ex : vals_ex)
          expand(key_ex, val_ex, opts, opt_aliases, opts_already_seen);
    }
  }

  if (not is_alias) opts[key].push_back(val);

  return is_alias;
}

std::size_t parse_opts_and_flags(flags_t &flags, opts_t &opts,
    args_t::const_iterator &pos, args_t::const_iterator const &end,
    opt_aliases_t const &opt_aliases) {
  auto const original_position{pos};
  std::unordered_set<key_t> opts_already_seen{};
  while (pos != end) {
    std::string const arg{*pos};
    if (not arg.starts_with(opt_prefix)) break;
    if (arg == opt_prefix) { ++pos; break; }
    auto const equal_sign_pos{arg.find_first_of(opt_assignment)};
    auto const arg_key{
      arg.substr(opt_prefix_size, equal_sign_pos - opt_prefix_size)};

    if (equal_sign_pos == decltype(arg)::npos) { // arg is a flag
      auto flag_itr{flags.find(arg_key)};
      if (flag_itr == std::end(flags)) break;
      flag_itr->second = true;
    } else { // arg is an opt
      arg_t const arg_val{arg.substr(equal_sign_pos + opt_assignment_size)};
      if (not opts.contains(arg_key)) break;
      expand(arg_key, arg_val, opts, opt_aliases, opts_already_seen);
    }
    ++pos;
  }

  // Go over unseen (default) opts, remove and expand
  // This assumes that there is not more than 1 (default) value per unseen opt
  auto const unseen{[&](opts_t::value_type const &opt){
      auto const &[key, val]{opt};
      return not opts_already_seen.contains(key);
    }};
  while (true) {
    auto opts_itr{std::find_if(std::begin(opts), std::end(opts), unseen)};
    if (opts_itr != std::end(opts)) {
      auto &[key, vals]{*opts_itr};
      if (vals.empty()) opts_already_seen.insert(key);
      else {
        arg_t const val{vals.back()};
        expand(key, val, opts, opt_aliases, opts_already_seen);
      }
    } else break;
  }

  return original_position - pos;
}

bool parse_bool(arg_t const &str){
  for (auto const &token : {"true", "1", "on", "yes"})
    if (str == token) return true;
  for (auto const &token : {"false", "0", "off", "no"})
    if (str == token) return false;
  throw(std::invalid_argument("\"" + str + "\" could not be parsed as bool"));
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

float parse_float(arg_t const &str){
  return std::strtof(str.c_str(), nullptr);
}
double parse_double(arg_t const &str){
  return std::strtod(str.c_str(), nullptr);
}
int parse_int(arg_t const &str){
  return safe_narrow<int>(std::strtol(str.c_str(), nullptr, 10)); }
long parse_long(arg_t const &str){
  return std::strtol(str.c_str(), nullptr, 10); }
long long parse_long_long(arg_t const &str){
  return std::strtoll(str.c_str(), nullptr, 10); }
unsigned int parse_unsigned_int(arg_t const &str){
  return safe_narrow<unsigned int>(std::strtoul(str.c_str(), nullptr, 10)); }
unsigned long parse_unsigned_long(arg_t const &str){
  return std::strtoul(str.c_str(), nullptr, 10); }
unsigned long long parse_unsigned_long_long(arg_t const &str){
  return std::strtoull(str.c_str(), nullptr, 10); }
std::string parse_string(arg_t const &str) { return str; }

bool constexpr is_empty(char const * const str) {
  return str == nullptr or *str == '\0';
}

int print_usage(arg_t const &program_name,
  args_t const &current_submodes, desc_t const desc,
  mode_t const * submodes_itr, mode_t const * const submodes_end,
  desc_t const * submodes_descs_itr, desc_t const * const submodes_descs_end,
  flag_descs_t const &flag_descs, opt_descs_t const &opt_descs,
  opt_aliases_t const &opt_aliases, UsageMessageMode const mode) {
  std::ostream &out = mode == UsageMessageMode::help ? std::cout : std::cerr;
  out << normal << "Usage:\n";
  out << ind << bold << program_name;
  for (auto const &mode : current_submodes) out << " [" << opt_prefix
    << "arg…] \\\n" << ind2 << mode;
  for (auto const &[flag, desc] : flag_descs) out << " \\\n" << ind3 << "["
    << opt_prefix << flag << "]";
  for (auto const &[opt, descs] : opt_descs) {
    auto const &[var, def, desc, val_descs]{descs};
    out << " \\\n" << ind3 << "[" << opt_prefix << opt << opt_assignment << "<";
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
  out << " \\\n" << ind3 << "[" << opt_prefix << "]";
  if (submodes_itr != nullptr and submodes_itr != submodes_end)
    out << " \\\n" << ind2 << "<mode> [" << opt_prefix << "arg…]";
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
      out << ind << bold << opt_prefix << flag << normal << "\n";
      if (flag == "help")
        out << ind2 << "Print this usage message.\n";
      if (not is_empty(desc)) out << ind2 << desc << "\n";
    }
  }

  if (not (opt_descs.empty() and opt_aliases.empty())) {
    out << "\n" << "Options:\n";
    bool is_first_entry{true};
    for (auto const &[opt, descs] : opt_descs) {
      auto const &[var, def, desc, val_descs]{descs};
      if (is_first_entry) is_first_entry = false; else out << "\n";
      out << ind << bold << opt_prefix << opt << "=<";
      if (not is_empty(var)) out << var; else out << "…";
      out << ">" << normal;
      if (not is_empty(def)) out << " (default: " << def << ")";
      out << "\n";
      if (not is_empty(desc)) out << ind2 << desc << "\n";

      if (not val_descs.empty()) out << "\n";
      for (auto const &[val, desc] : val_descs) {
        out << ind2 << bold << opt_prefix << opt << opt_assignment << val
          << normal;
        if (not is_empty(desc)) out << ": " << desc << "\n"; else out << "\n";
      }

      auto const opt_aliases_itr{opt_aliases.find(opt)};
      if (opt_aliases_itr != std::cend(opt_aliases) and
          not opt_aliases_itr->second.empty()) {
        out << "\n";
        for (auto const &[val, opts_ex] : opt_aliases_itr->second) {
          out << ind2 << bold << opt_prefix << opt << opt_assignment << val
            << normal << " expands into:\n" << bold;
          bool is_first_entry{true};
          for (auto const &[opt_ex, vals_ex] : opts_ex)
            for (arg_t const &val_ex : vals_ex) {
              if (is_first_entry) is_first_entry = false; else out << " \\\n";
              out << ind3 << opt_prefix << opt_ex << opt_assignment << val_ex;
            }
          out << normal << "\n";
        }
      }
    }
  }

  out << normal << std::flush;
  return mode == UsageMessageMode::help ? 0 : 1;
}

flags_t init_flags(flag_descs_t const &flag_descs) {
  flags_t flags{};
  for (auto const &[k, desc] : flag_descs)
    flags.emplace(std::pair<key_t, flag_t>{k, false});
  return flags;
}

opts_t init_opts(opt_descs_t const &opt_descs) {
  opts_t opts{};
  for (auto const &[key, descs] : opt_descs) {
    auto const &[var, def, desc, val_descs]{descs};
    opts.emplace(std::pair<key_t, optvals_t>{key,
      is_empty(def) ? optvals_t{} : optvals_t{def}});
  }
  return opts;
}

} // namespace args

