// SPDX-License-Identifier: MIT-0
#include "core.hpp"
#include <stdexcept>
using namespace hpc;
int main() {
  RunModel mode;
  auto check = [](bool b) {
    if (!b)
      throw std::runtime_error("Borrowing regression");
  };
  check(!mode.uncertain('A', {{12, 10}}));
  check(mode.borrow('A', {{12, 10}}, {{13, 100000}}, 12) == 12);
  check(mode.uncertain('A', {{12, 6}, {13, 4}}));
  check(mode.borrow('A', {{12, 6}, {13, 4}}, {{13, 100000}}, 12) == 13);
  check(mode.borrow('A', {{12, 6}, {13, 4}}, {{14, 100000}}, 12) == 12);
  check(mode.borrow('A', {{12, 6}, {13, 4}}, {{13, 4}}, 12) == 12);
  check(mode.borrow('A', {{12, 8}, {13, 2}}, {{13, 100000}}, 12) == 12);
}
