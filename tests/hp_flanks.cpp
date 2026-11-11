// SPDX-License-Identifier: MIT-0
#include "hp_flanks.hpp"
#include <stdexcept>

int main() {
  auto check = [](bool v) {
    if (!v)
      throw std::runtime_error("HP flank contract failed");
  };
  using hpc::has_span_flanks;
  check(has_span_flanks(1000, 1010, 2010, 1000, true, true, false, false));
  check(!has_span_flanks(999, 1010, 2010, 1000, true, true, false, false));
  check(!has_span_flanks(1000, 1010, 2009, 1000, true, true, false, false));
  check(!has_span_flanks(1000, 1010, 2010, 1000, false, true, false, false));
  check(!has_span_flanks(1000, 1010, 2010, 1000, true, false, false, false));
  check(has_span_flanks(0, 10, 1010, 1000, false, true, true, false));
  check(has_span_flanks(1000, 1010, 1010, 1000, true, false, false, true));
  check(!has_span_flanks(10, 9, 100, 1, true, true, false, false));
  check(!has_span_flanks(0, 101, 100, 1, true, true, true, true));
  check(has_span_flanks(2, 4, 1004, 500, true, true, false, false, 498, 0));
  check(!has_span_flanks(2, 4, 1004, 500, true, true, false, false, 497, 0));
  check(has_span_flanks(1000, 1002, 1004, 500, true, true, false, false, 0, 498));
  check(!has_span_flanks(1000, 1002, 1004, 500, true, true, false, false, 0, 497));
  check(!has_span_flanks(2, 4, 1004, 500, false, true, false, false, 1000, 0));
  check(!has_span_flanks(2, 4, 1004, 500, true, true, false, false));
  check(!has_span_flanks(1000, 1002, 1004, 500, true, true, false, false));
  // Terminal neighborhoods waive outward distance only. Adjacent matches,
  // full inward context, and the original interior policy remain mandatory.
  check(has_span_flanks(1000, 1003, 1153, 500, true, true, false, false, 0, 0, false, true));
  check(has_span_flanks(150, 153, 1153, 500, true, true, false, false, 0, 0, true, false));
  check(!has_span_flanks(1000, 1003, 1153, 500, true, false, false, false, 0, 0, false, true));
  check(!has_span_flanks(150, 153, 1153, 500, false, true, false, false, 0, 0, true, false));
  check(!has_span_flanks(150, 153, 303, 500, true, true, false, false, 0, 0, true, true));
  check(!has_span_flanks(1000, 1003, 1153, 500, true, true, false, false));
}
