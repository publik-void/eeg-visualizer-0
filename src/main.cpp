#include <cstddef>
#include <vector>

#include "signal-handling.hpp"
#include "args.hpp"
#include "server-dummy.hpp"
#include "client-visualizer.hpp"

args::desc_t constexpr desc_main{nullptr};
args::desc_t constexpr desc_server{nullptr};
args::desc_t constexpr desc_client{nullptr};

args::modes_t<4> constexpr modes_main{"error", "help", "server", "client"};
args::modes_t<3> constexpr modes_server{"error", "help", "dummy"};
args::modes_t<3> constexpr modes_client{"error", "help", "visualizer"};

std::array<args::desc_t, 4> constexpr modes_descs_main{nullptr, nullptr,
  "Select a server (data sending) mode.",
  "Select a client (data receiving) mode."};

std::array<args::desc_t, 3> constexpr modes_descs_server{nullptr, nullptr,
  "Send random data."};

std::array<args::desc_t, 3> constexpr modes_descs_client{nullptr, nullptr,
  "Visualize incoming EEG data in real time."};

int print_usage_main(args::arg_t const &program_name,
    args::UsageMessageMode const mode) {
  return args::print_usage(program_name, {}, desc_main,
    std::next(std::cbegin(modes_main)), std::cend(modes_main),
    std::next(std::cbegin(modes_descs_main)),
    std::cend(modes_descs_main), {}, {}, mode);
}

int print_usage_server(args::arg_t const &program_name,
    args::UsageMessageMode const mode) {
  return args::print_usage(program_name, {"server"}, desc_server,
    std::next(std::cbegin(modes_server)), std::cend(modes_server),
    std::next(std::cbegin(modes_descs_server)), std::cend(modes_descs_server),
    {}, {}, mode);
}

int print_usage_client(args::arg_t const &program_name,
    args::UsageMessageMode const mode) {
  return args::print_usage(program_name, {"client"}, desc_client,
    std::next(std::cbegin(modes_client)), std::cend(modes_client),
    std::next(std::cbegin(modes_descs_client)), std::cend(modes_descs_client),
    {}, {}, mode);
}

int main(int argc, char *argv[]) {
  set_signal_handlers();

  auto const cargs{args::make_args(argc, argv)};
  auto arg_itr{std::cbegin(cargs)};
  auto const arg_end{std::cend(cargs)};

  auto program_name{arg_itr == arg_end ? "eeg-visualizer" : *(arg_itr++)};

  auto const mode_main{
    args::get_mode_index(modes_main, arg_itr, arg_end, 0)};

  if (mode_main == args::get_mode_index(modes_main, "help"))
    return print_usage_main(program_name, args::UsageMessageMode::help);

  else if (mode_main == args::get_mode_index(modes_main, "server")) {

    auto const mode_server{
      args::get_mode_index(modes_server, arg_itr, arg_end, 0)};

    if (mode_server == args::get_mode_index(modes_server, "help"))
      return print_usage_server(program_name, args::UsageMessageMode::help);

    else if (mode_server == args::get_mode_index(modes_server, "dummy"))
      return main_server_dummy(arg_itr, arg_end, program_name);

    else // catch-all, includes `mode_server == error`
      return print_usage_server(program_name, args::UsageMessageMode::error);

  } else if (mode_main == args::get_mode_index(modes_main, "client")) {

    auto const mode_client{
      args::get_mode_index(modes_client, arg_itr, arg_end, 0)};

    if (mode_client == args::get_mode_index(modes_client, "help"))
      return print_usage_client(program_name, args::UsageMessageMode::help);

    else if (mode_client == args::get_mode_index(modes_client, "visualizer"))
      return main_client_visualizer(arg_itr, arg_end, program_name);

    else // catch-all, includes `mode_client == error`
      return print_usage_client(program_name, args::UsageMessageMode::error);

  } else // catch-all, includes `mode_main == error`
    return print_usage_main(program_name, args::UsageMessageMode::error);
}
