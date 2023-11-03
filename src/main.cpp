#include <cstddef>
#include <vector>
#include <iostream>

#include "signal-handling.hpp"
#include "args.hpp"
#include "server-dummy.hpp"
#include "client-visual.hpp"

args::modes_t<4> constexpr main_modes{"error", "help", "server-dummy",
  "client-visual"};

int print_usage(args::arg_t const &program_name, std::ostream &out = std::cout,
    int const exit_code = 0) {
  out
    << "Usage:\n"
    << "\t" << program_name << " …\n"
    << "\t\t(TODO)";
  return exit_code;
}

int main(int argc, char *argv[]) {
  set_signal_handlers();

  auto const cargs{args::make_args(argc, argv)};
  auto arg_itr{std::cbegin(cargs)};
  auto const arg_end{std::cend(cargs)};

  auto program_name{*(arg_itr++)};

  auto const main_mode{
    args::get_mode_index(main_modes, arg_itr, arg_end, 0)};

  if (main_mode == args::get_mode_index(main_modes, "help"))
    return print_usage(program_name);
  else if (main_mode == args::get_mode_index(main_modes, "server-dummy"))
    return main_server_dummy(arg_itr, arg_end);
  else if (main_mode == args::get_mode_index(main_modes, "client-visual"))
    return main_client_visual(arg_itr, arg_end);
  else // catch-all, includes `main_mode == error`
    return print_usage(program_name, std::cerr, 1);
}
