// SPDX-License-Identifier: MIT-0
#include "repeat_refinement.hpp"
#include <stdexcept>
#include <iostream>
#include <algorithm>
using namespace hpc;
static void check(bool value, int line = __builtin_LINE()) {
  if (!value)
    throw std::runtime_error("Repeat refinement regression at " + std::to_string(line));
}
static std::string random_bases(size_t n) {
  uint32_t state = 173;
  std::string s;
  for (size_t i = 0; i < n; ++i) {
    state = 1664525 * state + 1013904223;
    s += "ACGT"[(state >> 24) & 3];
  }
  return s;
}
int main() {
  for (auto q : {"A", "AC", "AAAAC", "ACAC", "ACCA", "CACA"})
    for (auto b : {"A", "AC", "ACA"})
      check(encode(project_repeat_backbone(q, b, 10000)).bases == b);
  bool limited = false;
  try {
    project_repeat_backbone("AAAA", "ACAC", 2);
  } catch (const std::length_error &) {
    limited = true;
  }
  check(limited);
  auto left = random_bases(1800), right = random_bases(1900);
  std::reverse(right.begin(), right.end());
  std::string correct, bad;
  for (size_t i = 0; i < 100; ++i) {
    correct += "AAAG";
    bad += "AG";
  }
  auto original = encode(left + bad + right), truth = encode(left + correct + right);
  check(original.bases == truth.bases);
  std::vector<Call> calls;
  for (auto n : original.lengths)
    calls.push_back({n, 6, 3, 0});
  auto regions = detect_repeats(original.bases, calls);
  check(!regions.empty());
  auto stable = calls;
  for (auto &c : stable)
    c.peak = c.depth;
  check(detect_repeats(original.bases, stable).empty());
  std::vector<Alignment> alignments(6);
  std::vector<RepeatRead> reads;
  for (size_t i = 0; i < alignments.size(); ++i) {
    auto &a = alignments[i];
    a.id = i + 1;
    a.tier = 0;
    a.sequence = left + correct + right;
    reads.push_back({&a, 0, original.bases.size(), true, 1});
  }
  RepeatConfig cfg;
  std::vector<RepeatDecision> reports;
  auto changed = refine_repeats(original.bases, calls, reads, 1, {}, cfg,
                                [&](auto d) { reports.push_back(d); });
  check(!changed.empty());
  size_t corrected = 0;
  for (size_t i = 0; i < calls.size(); ++i)
    if (original.lengths[i] != truth.lengths[i] && calls[i].length == truth.lengths[i])
      ++corrected;
  check(corrected >= 90);
  calls.clear();
  for (auto n : original.lengths)
    calls.push_back({n, 6, 3, 0});
  auto before = calls;
  cfg.max_reads = 2;
  reports.clear();
  refine_repeats(original.bases, calls, reads, 1, {}, cfg, [&](auto d) { reports.push_back(d); });
  for (size_t i = 0; i < calls.size(); ++i)
    check(calls[i].length == before[i].length);
  check(!reports.empty() && reports.front().status == "READ_BUDGET");
  // A sufficiently covered lower unique tier is used before ambiguous reads.
  calls.clear();
  for (auto n : original.lengths)
    calls.push_back({n, 6, 3, 0});
  for (auto &a : alignments)
    a.tier = 1;
  cfg.max_reads = 256;
  reports.clear();
  refine_repeats(original.bases, calls, reads, 2, {}, cfg, [&](auto d) { reports.push_back(d); });
  bool lower = false;
  for (auto &d : reports)
    lower |= d.status == "ACCEPTED" && d.tier == 1 && !d.pooled;
  check(lower);
  // Haplotype-conflict controls: a larger ambiguous or lower-priority pool
  // must not overwhelm adequately covered unique molecules in the first tier.
  // Both alleles share an HPC backbone, so only the raw run pattern differs.
  auto phase_case = [&](bool wrong_draft, bool ambiguous_competitors, bool all_ambiguous,
                        bool tied) {
    std::string other;
    for (size_t i = 0; i < 100; ++i)
      other += "AAAAAG";
    auto draft = encode(left + (wrong_draft ? other : correct) + right);
    std::vector<Call> phase_calls;
    for (auto n : draft.lengths)
      phase_calls.push_back({n, 30, 15, 0});
    std::vector<Alignment> pool(tied ? 6 : (all_ambiguous ? 6 : 30));
    std::vector<RepeatRead> phase_reads;
    for (size_t i = 0; i < pool.size(); ++i) {
      auto &a = pool[i];
      a.id = i + 1;
      bool competitor = !all_ambiguous && i >= (tied ? 3u : 6u);
      a.tier = competitor && !ambiguous_competitors && !tied ? 1 : 0;
      a.sequence = left + (competitor ? other : correct) + right;
      bool unique = !all_ambiguous && !(competitor && ambiguous_competitors);
      phase_reads.push_back({&a, 0, draft.bases.size(), unique, unique ? 1u : 2u});
    }
    std::vector<RepeatDecision> decisions;
    refine_repeats(draft.bases, phase_calls, phase_reads, 2, {}, cfg,
                   [&](auto d) { decisions.push_back(d); });
    for (size_t i = 0; i < phase_calls.size(); ++i)
      check(phase_calls[i].length == truth.lengths[i]);
    if (wrong_draft)
      check(std::any_of(decisions.begin(), decisions.end(), [&](const auto &d) {
        return d.status == "ACCEPTED" && d.tier == 0 && d.pooled == all_ambiguous;
      }));
  };
  phase_case(false, true, false, false);  // Keep correct unique allele vs 24 ambiguous.
  phase_case(true, true, false, false);   // Rescue unique allele despite ambiguous majority.
  phase_case(false, false, false, false); // Keep unique allele vs 24 lower-tier unique.
  phase_case(true, false, true, false);   // No unique evidence: pooled fallback works.
  phase_case(false, false, false, true);  // Balanced unique alleles: retain draft on ties.
  // Distinguish the draft-window bound from the anchored-read bound. Both
  // use the public limit, so a severely underexpanded draft can fit while its
  // correct supporting read intervals do not. These large synthetic repeats
  // verify that 15 kb actually admits work excluded by 10/12 kb.
  auto limit_case = [&](bool underexpanded) {
    std::string observed, initial;
    for (size_t i = 0; i < 3150; ++i) {
      observed += underexpanded ? "AAAG" : "AGGG";
      initial += underexpanded ? "AG" : "AAAG";
    }
    auto draft = encode(left + initial + right);
    auto expected = encode(left + observed + right);
    check(draft.bases == expected.bases);
    std::vector<Alignment> support(6);
    std::vector<RepeatRead> local_reads;
    for (size_t i = 0; i < support.size(); ++i) {
      support[i].id = i + 1;
      support[i].tier = 0;
      support[i].sequence = left + observed + right;
      local_reads.push_back({&support[i], 0, draft.bases.size(), true, 1});
    }
    for (size_t limit : {10000u, 12000u, 15000u}) {
      std::vector<Call> local_calls;
      for (auto n : draft.lengths)
        local_calls.push_back({n, 6, 3, 0});
      RepeatConfig config;
      config.max_window_bases = limit;
      std::vector<RepeatDecision> decisions;
      auto changes = refine_repeats(draft.bases, local_calls, local_reads, 1, {}, config,
                                    [&](auto d) { decisions.push_back(d); });
      if (limit < 15000) {
        check(changes.empty());
        for (size_t i = 0; i < local_calls.size(); ++i)
          check(local_calls[i].length == draft.lengths[i]);
        check(std::any_of(decisions.begin(), decisions.end(), [&](const auto &d) {
          return d.status == (underexpanded ? "LOW_ANCHORED_DEPTH" : "WINDOW_BUDGET");
        }));
      } else {
        check(!changes.empty());
        size_t corrected = 0;
        for (size_t i = 0; i < local_calls.size(); ++i)
          corrected += draft.lengths[i] != expected.lengths[i] &&
                       local_calls[i].length == expected.lengths[i];
        check(corrected >= (underexpanded ? 3100u : 6200u));
      }
    }
  };
  limit_case(false);
  limit_case(true);
  std::cout << "Repeat detection, projection, correction, tier and phase-conflict tests passed\n";
}
