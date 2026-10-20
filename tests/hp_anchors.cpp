// SPDX-License-Identifier: MIT-0
#include "experiments/hp_anchors.hpp"
#include <stdexcept>

int main() {
  using hpc::experimental::has_run_flanks;
  auto check = [](bool v) {
    if (!v)
      throw std::runtime_error("HP flank contract failed");
  };
  check(!has_run_flanks(0, 5, 10, 100));
  check(!has_run_flanks(4, 5, 10, 100));
  check(has_run_flanks(1, 5, 10, 100));
  check(!has_run_flanks(0, 1, 10, 100));
  check(has_run_flanks(0, 5, 0, 100));
  check(has_run_flanks(4, 5, 99, 100));
  check(!has_run_flanks(5, 5, 10, 100));
  check(!has_run_flanks(0, 0, 0, 0));
  using hpc::experimental::has_span_flanks;
  check(has_span_flanks(1000, 1010, 2010, 1000, true, true, false, false));
  check(!has_span_flanks(999, 1010, 2010, 1000, true, true, false, false));
  check(!has_span_flanks(1000, 1010, 2009, 1000, true, true, false, false));
  check(!has_span_flanks(1000, 1010, 2010, 1000, false, true, false, false));
  check(!has_span_flanks(1000, 1010, 2010, 1000, true, false, false, false));
  check(has_span_flanks(0, 10, 1010, 1000, false, true, true, false));
  check(has_span_flanks(1000, 1010, 1010, 1000, true, false, false, true));
  check(!has_span_flanks(10, 9, 100, 1, true, true, false, false));
  check(!has_span_flanks(0, 101, 100, 1, true, true, true, true));
}
