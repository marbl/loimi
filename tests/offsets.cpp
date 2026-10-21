// SPDX-License-Identifier: MIT-0
// Deterministic, dataset-free contracts for the failure-only offset fallback.
#include "offsets.hpp"
#include <algorithm>
#include <random>
#include <stdexcept>

static void check(bool condition) {
  if (!condition)
    throw std::runtime_error("Internal offset regression");
}

int main() {
  using namespace hpc;
  std::mt19937 rng(913);
  std::string sequence;
  for (size_t i = 0; i < 12000; ++i) {
    char base;
    do {
      base = "ACGT"[rng() % 4];
    } while (!sequence.empty() && sequence.back() == base);
    sequence += base;
  }
  auto reference = sequence.substr(3000, 6000);
  check(discover_offsets("ACGT", "CG", 0, .15).empty());
  check(!validate_offset_overlap("ACGT", "CG", 0, 4, 1, .15));
  auto anchors = internal_offsets(sequence, reference, 100);
  check(!anchors.empty() &&
        std::all_of(anchors.begin(), anchors.end(), [](int64_t x) { return x == 3100; }));
  check(internal_offsets(sequence, reference + reference, 100).empty());
  check(discover_offsets(sequence, reference, 100, .15) == std::vector<int64_t>{3100});
  check(discover_offsets(reference, reference + reference, 100, .15).empty());
  check(discover_offsets(sequence.substr(4000, 1000), reference, 4100, .15) ==
        std::vector<int64_t>{3100});
  auto partial = discover_offsets(sequence.substr(0, 8000), reference, 100, .15);
  check(std::find(partial.begin(), partial.end(), 3100) != partial.end());
  check(validate_offset_overlap(sequence, reference, 100, 12100, 3100, .15));
  check(!validate_offset_overlap(sequence, reference, 100, 12100, 20000, .15));
  // The old free-end validator relocated this wrong origin to the right locus.
  check(!validate_offset_overlap(sequence, reference, 100, 12100, 1000, .15));
  const std::string uuid = "50cb72c7-194d-4194-9755-747e8f9c6e66";
  check(offset_molecule("hifi_" + uuid, 1) == offset_molecule("ont_" + uuid, 2));
  check(offset_molecule("hifi_m84031/1/ccs", 1) != offset_molecule("hifi_m84031/2/ccs", 2));
  auto joint = choose_joint_offset({3000, 36000}, {{"a", 0, true, {{20000, 3}, {}}},
                                                   {"b", 1, true, {{17000, 373}, {21000, 528}}}});
  check(joint.value == 3000 && joint.support == 2);
  // Sparse support can span tiers; multiple technologies cannot manufacture
  // additional independent molecules from a corrected/raw pair.
  joint = choose_joint_offset({3000}, {{"a", 0, true, {{10000, 5}}},
                                       {"a", 1, true, {{10000, 10}}},
                                       {"b", 1, false, {{10000, 20}}}});
  check(joint.support == 2 && joint.informative == 2);
  check(offset_candidates({{1, true, {59841}}, {2, true, {59891}}}, {0, 1}) ==
        std::vector<int64_t>{59841});
  bool joint_tie = false;
  try {
    (void)choose_joint_offset({3000, 36000}, {{"a", 0, true, {{10000, 5}, {10000, 5}}},
                                              {"b", 1, true, {{10000, 10}, {10000, 10}}}});
  } catch (const std::runtime_error &) {
    joint_tie = true;
  }
  check(joint_tie);
  // Five independent supporters and one outlier satisfy the 75% rule.
  check(internal_offset_cluster({{3100, 3101}, {3100}, {3099}, {3102}, {3100}, {8000}}, 3) == 3100);
  check(internal_offset_cluster({{-500}}, 1) == -500);
  for (auto evidence : std::vector<std::vector<std::vector<int64_t>>>{
           {},
           {{1, 1, 1, 1, 1, 1}},
           {{1}, {1}, {1}, {100}, {100}, {100}},
           {{1, 100}, {1, 100}, {1, 100}, {1, 100}, {1, 100}, {1, 100}}}) {
    bool rejected = false;
    try {
      (void)internal_offset_cluster(evidence, 3);
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    check(rejected);
  }
}
