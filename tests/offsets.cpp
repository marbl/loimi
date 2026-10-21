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
  std::vector<OffsetEvidence> rescue;
  for (uint32_t i = 0; i < 11; ++i)
    rescue.push_back({i, false, {i < 8 ? 201755 : i < 10 ? 234918 : 168588}});
  auto choice = choose_offset(rescue, 3);
  check(choice.value == 201755 && choice.support == 8 && choice.runner_up == 2);
  // Multiple placements of one ID are not independent confirmation.
  for (auto evidence : std::vector<std::vector<OffsetEvidence>>{
           {{1, false, {100}}, {1, false, {100}}, {1, false, {100}}},
           {{1, false, {100, 200}}, {2, false, {100, 200}}, {3, false, {100, 200}}},
           {{1, true, {200}},
            {2, true, {200}},
            {3, true, {200}},
            {4, false, {100}},
            {5, false, {100}},
            {6, false, {100}},
            {7, false, {100}},
            {8, false, {100}},
            {9, false, {100}},
            {10, false, {100}}}}) {
    bool rejected = false;
    try {
      (void)choose_offset(evidence, 3);
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    check(rejected);
  }
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
