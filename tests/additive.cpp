// SPDX-License-Identifier: MIT-0
#include "core.hpp"
#include <stdexcept>

using namespace hpc;
static void check(bool ok) {
  if (!ok)
    throw std::runtime_error("Additive polishing regression");
}
int main() {
  const std::string ref = "ACGTAC", alt = "ACTTAC";
  std::vector<Alignment> reads;
  auto add = [&](int tier, int count, bool alternate) {
    for (int i = 0; i < count; ++i) {
      Alignment a;
      a.tier = tier;
      a.sequence = alternate ? alt : ref;
      a.end = ref.size();
      a.ops = align(a.sequence, ref, 0.5);
      a.reverse = i % 2;
      reads.push_back(a);
    }
  };
  add(0, 4, true);
  add(1, 5, true);
  add(1, 1, false);
  check(polish(ref, reads, 2, 5, .9, 500).empty());
  auto edits = polish(ref, reads, 2, 5, .9, 500, 2);
  check(edits.size() == 1 && edits[0].tier == 1 && edits[0].depth == 10 && edits[0].support == 9);
  reads.clear();
  add(0, 2, false);
  add(1, 6, true);
  check(polish(ref, reads, 2, 5, .9, 500).size() == 1);
  check(polish(ref, reads, 2, 5, .9, 500, 2).empty()); // Higher unique votes can veto.
  reads.clear();
  add(0, 5, false);
  add(1, 6, true);
  check(polish(ref, reads, 2, 5, .9, 500, 2).empty()); // Qualified first tier wins.
  reads.clear();
  add(0, 4, true);
  add(1, 4, true);
  check(polish(ref, reads, 2, 5, .9, 500, 2).empty()); // No cumulative qualification.
  reads.clear();
  add(0, 2, false);
  add(2, 6, true);
  check(polish(ref, reads, 3, 5, .9, 500, 2).size() == 1); // Ambiguous policy unchanged.
}
