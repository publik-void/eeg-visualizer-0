#include <csignal>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <chrono>

namespace {
  // Unnamed namespace for internal linkage
  // As far as I can tell, this basically makes definitions private to the
  // translation unit and is preferred over static variables.
  std::atomic_bool graceful_exit_signal_received{false};
  std::condition_variable cv;
  std::mutex cv_mutex;
  void graceful_exit(int const = 0) {
    graceful_exit_signal_received = true; cv.notify_all();
  }
}

void set_signal_handlers() {
  std::signal(SIGINT , graceful_exit);
  std::signal(SIGTERM, graceful_exit);
}

bool wait_until(std::chrono::high_resolution_clock::time_point t) {
  std::unique_lock<std::mutex> lk(cv_mutex);
  cv.wait_until(lk, t);
  return graceful_exit_signal_received;
}

bool graceful_exit_signal(bool const set) {
  if (set) graceful_exit_signal_received = true;
  return graceful_exit_signal_received;
}

