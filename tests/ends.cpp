// SPDX-License-Identifier: MIT-0
// Terminal path repair must recover seed errors without changing the interior,
// violating tier precedence, or counting a clipped final run as full evidence.
#include "core.hpp"
#include <algorithm>
#include <random>
#include <stdexcept>
using namespace hpc;

static void check(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}

struct Input {
  std::string sequence;
  int tier = 0;
  bool unique = true;
};

static std::string extend(const std::vector<Input> &input, bool left = false, size_t limit = 1000,
                          bool pooled = true) {
  const std::string interior(100, 'C');
  std::vector<Alignment> alns(input.size());
  std::vector<EndRead> reads;
  for (size_t i = 0; i < input.size(); ++i) {
    auto &a = alns[i];
    auto tail = input[i].sequence;
    if (left)
      std::reverse(tail.begin(), tail.end());
    a.id = i + 1;
    a.tier = input[i].tier;
    a.sequence = left ? tail + interior : interior + tail;
    a.end = interior.size();
    a.ops = "100=";
    a.packed_ops = true;
    (left ? a.clip_left : a.clip_right) = tail.size();
    reads.push_back({&a, input[i].unique, a.clip_left, a.clip_right, true});
  }
  auto reference = interior;
  auto result = extend_ends(reference, reads, 2, {}, 5, limit, .15, pooled, left, !left);
  check(reference.substr(result.left, interior.size()) == interior, "Interior changed");
  check(result.left + result.right <= limit, "Extension exceeded cap");
  for (auto &a : alns) {
    size_t nm;
    md(a, reference, nm); // Validate query/reference consumption and match tags.
  }
  auto added = left ? reference.substr(0, result.left) : reference.substr(interior.size());
  if (left)
    std::reverse(added.begin(), added.end());
  return added;
}

int main() {
  // ACG has several equally good prefix endpoints in ACTG. These tails are
  // already anchored; rejecting them all would erase supported AC evidence.
  const std::string lead = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT";
  std::vector<Input> tied(6, {lead + "ACG"});
  tied[0].sequence = lead + "ACTG";
  check(extend(tied) == lead + "AC", "Tied terminal endpoints discarded anchored evidence");
  std::mt19937 rng(913);
  std::string truth;
  while (truth.size() < 250) {
    char c = "ACGT"[rng() % 4];
    if (truth.empty() || c != truth.back())
      truth += c;
  }
  const auto expected = truth.substr(0, truth.size() - 1);
  for (int kind = 0; kind < 3; ++kind) {
    auto damaged = truth;
    size_t p = 50;
    char alternate = 'A';
    for (char c : std::string("ACGT"))
      if (c != truth[p - 1] && c != truth[p] && c != truth[p + 1])
        alternate = c;
    if (kind == 0)
      damaged.insert(p, 1, alternate); // Extra seed run needs deletion.
    else if (kind == 1)
      damaged[p] = alternate; // Wrong seed base needs substitution.
    else {
      damaged.erase(p, 1);                // Missing seed run needs insertion.
      damaged.insert(10, 3, damaged[10]); // Make it longest without new HPC runs.
    }
    std::vector<Input> reads(7, {truth});
    reads[0].sequence = damaged;
    for (bool left : {false, true})
      check(extend(reads, left) == expected, "Supported seed repair failed");
    check(extend(reads, false, 80).size() <= 80, "Capped repair failed");
  }
  // More lower-tier reads cannot displace sufficient higher-tier evidence.
  auto other = truth;
  other.insert(50, 1, truth[52]);
  std::vector<Input> reads(6, {truth});
  reads.insert(reads.end(), 20, {other, 1});
  check(extend(reads) == expected, "Lower tier overrode unique evidence");
  // No unique tier qualifies: same-tier pooling remains available.
  reads.assign(7, {truth, 0, false});
  reads[0].sequence.insert(50, 1, truth[52]);
  check(extend(reads) == expected, "Pooled path repair failed");
  check(extend(reads, false, 1000, false).empty(), "Unique-only admitted ambiguous tails");
  reads.resize(4);
  check(extend(reads).empty(), "Sparse tails established extension");
}
