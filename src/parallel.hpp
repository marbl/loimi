// SPDX-License-Identifier: MIT-0
#pragma once
#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace hpc {
// Each action owns a distinct index; shared accumulators must be updated after
// this call. Join every worker before propagating the first unexpected failure.
// In particular, resource exhaustion must not silently discard read evidence.
template <class F> void parallel(size_t count, unsigned threads, F action) {
  if (!threads)
    throw std::invalid_argument("Worker count must be positive");
  std::atomic<size_t> next{0};
  std::atomic<bool> stop{false};
  std::exception_ptr failure;
  std::mutex failure_mutex;
  std::vector<std::thread> workers;
  auto join = [&] {
    for (auto &worker : workers)
      worker.join();
  };
  try {
    for (size_t t = 0; t < std::min<size_t>(count, threads); ++t)
      workers.emplace_back([&] {
        try {
          for (size_t i; !stop.load(std::memory_order_relaxed) &&
                         (i = next.fetch_add(1, std::memory_order_relaxed)) < count;)
            action(i);
        } catch (...) {
          std::lock_guard<std::mutex> lock(failure_mutex);
          if (!failure)
            failure = std::current_exception();
          stop.store(true, std::memory_order_relaxed);
        }
      });
  } catch (...) {
    // Thread creation can fail after earlier workers have started. Destructing
    // joinable threads would terminate the process instead of preserving logs.
    stop.store(true, std::memory_order_relaxed);
    join();
    throw;
  }
  join();
  if (failure)
    std::rethrow_exception(failure);
}
} // namespace hpc
