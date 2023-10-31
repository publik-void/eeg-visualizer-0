#pragma once

#include <chrono>

void set_signal_handlers();
bool wait_until(std::chrono::high_resolution_clock::time_point);
bool graceful_exit_signal(bool const = false);
