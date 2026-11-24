// SPDX-License-Identifier: MIT-0
// Deterministic, dataset-free contracts for the failure-only offset fallback.
#include "offsets.hpp"
#include "core.hpp"
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
  // HG06807 part025/tig48: nearly periodic template with one distinguishing
  // base. Read-derived long chunks can all miss that marker, while short
  // template-derived anchors establish the correct copy without softening ties.
  auto unit = sequence.substr(0, 732);
  for (char base : std::string("ACGT"))
    if (base != unit.front() && base != unit[730]) {
      unit.back() = base;
      break;
    }
  std::string near_repeat;
  while (near_repeat.size() < 9958)
    near_repeat += unit;
  near_repeat.resize(9958);
  auto exact_repeat = near_repeat;
  for (char base : std::string("ACGT"))
    if (base != near_repeat[709] && base != near_repeat[710] && base != near_repeat[711]) {
      near_repeat[710] = base;
      break;
    }
  check(template_repeat_period(near_repeat) == 0);
  auto marker_query = sequence.substr(8000, 1024) + near_repeat + sequence.substr(9500, 1024);
  check(internal_offsets(marker_query, near_repeat, 3000).empty());
  auto marker_words = template_offset_anchors(near_repeat);
  check(!marker_words.empty() && marker_words.size() <= 128);
  check(template_offset_anchors(exact_repeat).empty());
  auto markers = template_anchor_offsets(marker_query, marker_words, 3000);
  check(markers == std::vector<int64_t>{4024});
  check(template_anchor_offsets(marker_query + marker_query, marker_words, 3000).empty());
  auto marker_choice = choose_template_anchor_offset({{"la", 0, markers}, {"ul", 2, {4013}}}, 2);
  check(marker_choice.value == 4024 && marker_choice.tier == 0 && marker_choice.support == 2);
  check(validate_offset_overlap(marker_query, near_repeat, 3000, 3000 + marker_query.size(),
                                marker_choice.value, .15));
  auto marker_rejected = [&](const std::vector<OffsetAnchorRow> &rows) {
    try {
      choose_template_anchor_offset(rows, 2);
      return false;
    } catch (const std::runtime_error &) {
      return true;
    }
  };
  check(marker_rejected({{"la", 0, {4024, 4024, 4024}}})); // Windows are not molecules.
  check(marker_rejected({{"alias", 1, {4024}}, {"alias", 2, {4013}}}));
  check(marker_rejected({{"la", 0, {4024}}, {"ul", 2, {4756}}})); // Wrong repeat copy.
  check(marker_rejected({{"la", 0, {4024, 4756}}, {"ul", 2, {4024}}}));
  check(marker_rejected({{"la", 0, {4024}}, {"ul1", 2, {4756}}, {"ul2", 2, {4756}}}));
  check(choose_template_anchor_offset({{"la", 0, {4024}}}, 1).value == 4024);
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
  // A partial overlap with isolated mismatches has no long exact chain, but
  // smaller well-separated anchors still propose the same validated origin.
  auto noisy_partial = sequence.substr(0, 8000);
  for (size_t i = 3500; i < noisy_partial.size(); i += 1200)
    noisy_partial[i] = noisy_partial[i] == 'A' ? 'C' : 'A';
  auto short_anchors = discover_offsets(noisy_partial, reference, 100, .15);
  check(std::find(short_anchors.begin(), short_anchors.end(), 3100) != short_anchors.end());
  check(validate_offset_overlap(noisy_partial, reference, 100, 8100, 3100, .15));
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
  // Three-molecule sparse layout: LAO and one UL agree; another UL prefers a
  // distant copy. Two independent supporters suffice without changing tiers.
  joint = choose_joint_offset({1030, 66363}, {{"lao", 0, false, {{}, {30724, 52}}},
                                              {"ul1", 1, true, {{41045, 371}, {41045, 508}}},
                                              {"ul2", 1, true, {{}, {41045, 1696}}}});
  check(joint.value == 66363 && joint.support == 2 && joint.informative == 3);
  // Part008/tig562: competing offsets look good to individual reads,
  // but bounded overlaps across all molecules resolve the later candidate.
  const auto rescue_policy = OffsetSupportPolicy::SparseRescue;
  const std::vector<int64_t> rescue_loci = {11170, 13360, 14823};
  const std::vector<OffsetRow> rescue_rows = {{"a", 0, true, {{3621, 1}, {3621, 2}, {3067, 4}}},
                                              {"b", 0, true, {{3621, 1}, {3621, 1}, {3621, 4}}},
                                              {"c", 0, true, {{3621, 13}, {3621, 4}, {3621, 5}}}};
  joint = choose_joint_offset(rescue_loci, rescue_rows, {}, rescue_policy);
  check(joint.value == 13360 && joint.support == 3 && joint.runner_up == 2);
  auto sparse_rescue = [&](const std::vector<OffsetRow> &rows, size_t support) {
    try {
      return choose_joint_offset({3000, 36000}, rows, {}, rescue_policy).support == support;
    } catch (const std::runtime_error &) {
      return support == 0;
    }
  };
  check(sparse_rescue({{"a", 0, true, {{10000, 5}, {}}}, {"b", 1, true, {{10000, 5}, {}}}}, 2));
  check(sparse_rescue({{"a", 0, true, {{10000, 5}, {}}}}, 1));
  check(sparse_rescue({{"a", 0, false, {{10000, 5}, {}}}}, 0));
  check(sparse_rescue({{"a", 0, true, {{10000, 5}, {10000, 5}}}}, 0));
  check(sparse_rescue({{"a", 0, true, {{10000, 5}, {}}},
                       {"a", 1, true, {{10000, 5}, {}}},
                       {"b", 1, true, {{}, {10000, 5}}}},
                      0));
  check(sparse_rescue(
      {{"a", 0, true, {{10000, 5}, {10000, 5}}}, {"b", 0, true, {{10000, 5}, {10000, 5}}}}, 0));
  // A nonperiodic template may occur several times in an enclosing read.
  // Ordinary discovery rejects that read. Shared discovery keeps both loci.
  auto nonperiodic_repeat = reference + sequence.substr(0, 1000) + reference;
  check(template_repeat_period(reference) == 0);
  check(discover_offsets(nonperiodic_repeat, reference, 0, .15).empty());
  check(discover_rescue_offsets(nonperiodic_repeat, reference, {0}, nonperiodic_repeat.size(), 0,
                                .15) ==
        std::vector<int64_t>({0, int64_t(reference.size() + 1000)}));
  // The short-anchor retry can propose a partial overlap when the whole
  // shorter sequence has no admissible match.
  auto short_overlap_read = sequence.substr(0, 5400);
  check(
      discover_rescue_offsets(short_overlap_read, reference, {0}, short_overlap_read.size(), 0, .15)
          .empty());
  check(discover_rescue_offsets(short_overlap_read, reference, {0}, short_overlap_read.size(), 0,
                                .15, true) == std::vector<int64_t>{3000});
  // Do not truncate the competitor set to make an over-budget query resolve.
  std::vector<int64_t> many_starts;
  for (int64_t begin = 0; begin <= 4096; ++begin)
    many_starts.push_back(begin);
  bool budget_rejected = false;
  try {
    (void)discover_rescue_offsets(nonperiodic_repeat, reference, many_starts,
                                  nonperiodic_repeat.size() + 4096, 0, .15);
  } catch (const std::runtime_error &) {
    budget_rejected = true;
  }
  check(budget_rejected);
  // Uninformative molecules cannot veto a unanimous informative tier.
  for (size_t n = 1; n <= 8; ++n)
    for (size_t supporters = 1; supporters <= n; ++supporters) {
      std::vector<OffsetRow> rows;
      for (size_t i = 0; i < n; ++i)
        rows.push_back(
            {std::to_string(i), i < supporters ? 0 : 1, true,
             i < supporters ? std::vector<OffsetScore>{{10000, 5}} : std::vector<OffsetScore>{{}}});
      bool accepted = false;
      OffsetChoice result;
      try {
        result = choose_joint_offset({3000}, rows);
        accepted = true;
      } catch (const std::runtime_error &) {
      }
      check(accepted);
      if (accepted)
        check(result.value == 3000 && result.support == supporters);
    }
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
  // Partition 026/tig74: a split higher tier is not decisive even when its
  // unique-support tie break favors the wrong candidate. UL can resolve it.
  std::vector<OffsetRow> tier_rows = {{"lao1", 0, true, {{10000, 5}, {}}},
                                      {"lao2", 0, false, {{}, {10000, 5}}}};
  for (size_t i = 0; i < 8; ++i)
    tier_rows.push_back({"ul" + std::to_string(i), 1, true, {{}, {10000, 5}}});
  joint = choose_joint_offset({-810, 59199}, tier_rows);
  check(joint.value == 59199 && joint.support == 9 && joint.tier == 1);
  // A decisive higher tier still wins even against many lower-tier reads.
  tier_rows.push_back({"lao3", 0, true, {{10000, 5}, {}}});
  tier_rows.push_back({"lao4", 0, true, {{10000, 5}, {}}});
  joint = choose_joint_offset({-810, 59199}, tier_rows);
  check(joint.value == -810 && joint.support == 3 && joint.tier == 0);
  // Sparse unique-overlap fallback: two placements of one non-unique LA
  // molecule conflict with a clearly better unique origin. The mixed vote
  // fails, but the same selector can use the sparse unique subset. Keeping
  // that selector also means a tied unique read must remain unresolved.
  std::vector<OffsetRow> sparse_conflict = {
      {"unique", 0, true, {{6278, 10}, {6278, 48}, {1760, 36}}},
      {"ambiguous", 0, false, {{5137, 124}, {6278, 0}, {2823, 13}}},
      {"ambiguous", 0, false, {{}, {6278, 30}, {6278, 0}}}};
  bool mixed_rejected = false;
  try {
    choose_joint_offset({208, 10842, 18826}, sparse_conflict);
  } catch (const std::runtime_error &) {
    mixed_rejected = true;
  }
  check(mixed_rejected);
  check(choose_joint_offset({208, 10842, 18826}, {sparse_conflict.front()}).value == 208);
  auto tied_unique = sparse_conflict.front();
  tied_unique.scores = {{6278, 10}, {6278, 10}, {}};
  bool unique_rejected = false;
  try {
    choose_joint_offset({208, 10842, 18826}, {tied_unique});
  } catch (const std::runtime_error &) {
    unique_rejected = true;
  }
  check(unique_rejected);
  // A unanimous informative higher tier wins, including a single molecule.
  joint = choose_joint_offset({3000, 36000}, {{"a", 0, true, {{10000, 5}, {}}},
                                              {"b", 1, true, {{}, {10000, 5}}},
                                              {"c", 1, true, {{}, {10000, 5}}},
                                              {"d", 1, true, {{}, {10000, 5}}}});
  check(joint.value == 3000 && joint.tier == 0);
  // Reported 1/1 UL case: three sampled molecules are uninformative and
  // none is uniquely placed. The one resolved informative vote is sufficient.
  joint = choose_joint_offset({16888}, {{"a", 0, false, {{}}},
                                        {"b", 1, false, {{}}},
                                        {"c", 1, false, {{}}},
                                        {"d", 2, false, {{10000, 10}}}});
  check(joint.value == 16888 && joint.support == 1 && joint.tier == 2);
  // A higher tier can also be inconclusive because every read supports both
  // loci. Lower-tier resolution is allowed; unresolved ties must still fail.
  tier_rows = {{"a", 0, true, {{10000, 5}, {10000, 5}}},
               {"b", 0, true, {{10000, 5}, {10000, 5}}},
               {"c", 1, true, {{}, {10000, 5}}}};
  check(choose_joint_offset({3000, 36000}, tier_rows).tier == 1);
  tier_rows.push_back({"d", 1, true, {{10000, 5}, {}}});
  bool tiers_unresolved = false;
  try {
    (void)choose_joint_offset({3000, 36000}, tier_rows);
  } catch (const std::runtime_error &) {
    tiers_unresolved = true;
  }
  check(tiers_unresolved);
  // Nov-2026 part026/tig399: a whole-template tandem repeat makes two
  // near-period-separated origins indistinguishable. The lower origin is a
  // final tie-break, never a replacement for a stronger evidence decision.
  auto motif = sequence.substr(0, 971);
  std::string periodic;
  while (periodic.size() < 5536)
    periodic += motif;
  periodic.resize(5536);
  check(template_repeat_period(periodic) == 971);
  // Fully periodic templates previously lost every tied proposal before the
  // repeat-aware selector ran. Retain only actual optimal origins inside the
  // known layout, not invented phases or unsupported extensions.
  std::string repeat_query;
  for (int i = 0; i < 10; ++i)
    repeat_query += motif;
  check(discover_offsets(repeat_query, periodic, 0, .15).empty());
  auto periodic_origins =
      discover_rescue_offsets(repeat_query, periodic, {0}, repeat_query.size(), 971, .15);
  check(periodic_origins == std::vector<int64_t>({0, 3884}));
  auto short_query = motif + motif;
  check(discover_rescue_offsets(short_query, periodic, {0}, periodic.size(), 971, .15) ==
        std::vector<int64_t>{0});
  check(discover_rescue_offsets(short_query, periodic, {5}, periodic.size(), 971, .15).empty());
  check(discover_rescue_offsets(short_query, periodic, {0}, periodic.size() - 1, 971, .15).empty());
  auto bad_periodic = periodic;
  bad_periodic[2700] = bad_periodic[2700] == 'A' ? 'C' : 'A';
  bool bad_period_rejected = false;
  try {
    discover_rescue_offsets(short_query, bad_periodic, {0}, periodic.size(), 971, .15);
  } catch (const std::invalid_argument &) {
    bad_period_rejected = true;
  }
  check(bad_period_rejected);
  // Cacheable contained placements have identical target strings after an
  // exact-period shift. Clipped edge placements must retain separate scores.
  auto repeat_read = periodic.substr(100, 2000);
  for (char base : std::string("ACGT"))
    if (base != repeat_read[499] && base != repeat_read[500] && base != repeat_read[501]) {
      repeat_read[500] = base;
      break;
    }
  auto first_score = score_offset_overlap_local(repeat_read, periodic, 100, 2100, 0, .15);
  auto shifted_score = score_offset_overlap_local(repeat_read, periodic, 1071, 3071, 0, .15);
  check(first_score.overlap && first_score.overlap == shifted_score.overlap &&
        first_score.edits == shifted_score.edits);
  auto edge_score = score_offset_overlap_local(repeat_read, periodic, 3984, 5984, 0, .15);
  check(edge_score.overlap != first_score.overlap);

  std::vector<OffsetRow> repeated = {
      {"a", 0, true, {{5536, 11}, {5273, 10}}},  {"b", 0, false, {{5536, 10}, {5536, 8}}},
      {"c", 0, true, {{5536, 14}, {5536, 12}}},  {"d", 0, false, {{5536, 8}, {5536, 10}}},
      {"e", 0, true, {{5536, 20}, {5536, 0}}},   {"f", 0, false, {{5536, 8}, {5536, 10}}},
      {"g", 0, true, {{5536, 12}, {5536, 12}}},  {"h", 0, false, {{5536, 8}, {5536, 10}}},
      {"i", 0, false, {{5536, 10}, {5536, 10}}}, {"j", 0, true, {{5536, 4}, {5536, 16}}}};
  auto repeat_choice = choose_joint_offset({10254, 12187}, repeated, periodic);
  check(repeat_choice.value == 10254 && repeat_choice.support == 9 &&
        repeat_choice.repeat_period == 971);
  for (auto &r : repeated)
    std::reverse(r.scores.begin(), r.scores.end());
  check(choose_joint_offset({12187, 10254}, repeated, periodic).value == 10254);
  // Greater later-offset support wins regardless of candidate enumeration.
  std::vector<OffsetRow> better_later;
  for (int i = 0; i < 5; ++i)
    better_later.push_back({std::to_string(i), 0, true, {{5536, 30}, {5536, 0}}});
  auto better = choose_joint_offset({10254, 12187}, better_later, periodic);
  check(better.value == 12187 && better.repeat_period == 0);
  // A decisive lower tier also precedes the otherwise tied higher-tier repeat
  // fallback, preserving the established tier-resolution policy.
  better_later = {{"a", 0, true, {{5536, 0}, {5536, 0}}}, {"b", 1, true, {{5536, 30}, {5536, 0}}}};
  check(choose_joint_offset({10254, 12187}, better_later, periodic).value == 12187);
  auto rejects = [&](const std::string &ref, int64_t later) {
    try {
      choose_joint_offset({10254, later}, {{"a", 0, true, {{5536, 0}, {5536, 0}}}}, ref);
      return false;
    } catch (const std::runtime_error &) {
      return true;
    }
  };
  auto broken = periodic;
  broken[2700] = broken[2700] == 'A' ? 'C' : 'A';
  check(template_repeat_period(broken) == 0);
  check(rejects(broken, 12187));               // Local repeats alone cannot prove equivalence.
  check(rejects(periodic, 10254 + 971 + 400)); // Wrong repeat phase.
  check(rejects("", 12187));                   // No template: preserve the original ambiguity veto.

  for (size_t period : {127u, 971u, 2003u}) {
    auto unit = sequence.substr(0, period);
    std::string query;
    for (int i = 0; i < 10; ++i)
      query += unit;
    auto ref = query.substr(0, 3 * period);
    for (auto sides : {std::pair{true, true}, std::pair{true, false}, std::pair{false, true}}) {
      size_t predicted = sides.first ? (sides.second ? 2 : 7) * period : 0;
      auto a = align_repeat_boundary(query, ref, sides.first, sides.second, predicted, period, .15);
      check(a.query_begin == predicted && a.query_end == predicted + ref.size());
      check(a.ops == std::string(ref.size(), '='));
      auto anchored = align_marker_boundary(query, ref, sides.first, sides.second, predicted, .15);
      check(anchored.query_begin == predicted && anchored.query_end == predicted + ref.size());
      check(anchored.ops == std::string(ref.size(), '='));
    }
    // Sequencing indels may adjust the predicted endpoint without moving to
    // an adjacent repeat copy. Validate consumption as well as start position.
    auto noisy = query;
    noisy.insert(2 * period + 50, 1, 'A');
    auto a = align_repeat_boundary(noisy, ref, true, true, 2 * period, period, .15);
    check(a.query_begin == 2 * period);
    check(size_t(std::count_if(a.ops.begin(), a.ops.end(), [](char c) { return c != 'I'; })) ==
          ref.size());
    check(size_t(std::count_if(a.ops.begin(), a.ops.end(), [](char c) { return c != 'D'; })) ==
          a.query_end - a.query_begin);
    // An exact match in a distant copy must NOT rescue an incompatible local
    // placement. A zero error budget makes this negative selection explicit.
    ref[30] = ref[30] == 'A' ? 'C' : 'A';
    query[6 * period + 30] = ref[30];
    bool refused = false;
    try {
      align_repeat_boundary(query, ref, true, true, 2 * period, period, 0);
    } catch (const std::runtime_error &) {
      refused = true;
    }
    check(refused);
  }
  // Part012/tig357: fixed endpoints split four LAO reads across three repeat
  // copies. Bounded endpoint freedom removes layout drift without lowering the
  // existing tier-majority or independent-molecule requirements.
  std::vector<int64_t> drift_candidates{10609, 42360, 74512};
  std::vector<OffsetRow> fixed_rows{{"a", 1, true, {{42974, 366}, {42305, 438}, {10153, 168}}},
                                    {"b", 1, true, {{17261, 689}, {42974, 122}, {42974, 267}}},
                                    {"c", 1, true, {{14761, 497}, {42974, 280}, {42974, 126}}},
                                    {"d", 1, false, {{2873, 347}, {18502, 246}, {}}},
                                    {"e", 2, false, {{}, {18541, 348}, {}}},
                                    {"f", 2, true, {{}, {29807, 602}, {29170, 398}}}};
  bool drift_rejected = false;
  try {
    choose_joint_offset(drift_candidates, fixed_rows);
  } catch (const std::runtime_error &) {
    drift_rejected = true;
  }
  check(drift_rejected);
  auto local_rows = fixed_rows;
  local_rows[0].scores = {{42462, 238}, {41793, 79}, {9641, 2}};
  local_rows[1].scores = {{16749, 141}, {42462, 80}, {42462, 111}};
  local_rows[2].scores = {{14249, 41}, {42462, 86}, {42462, 76}};
  local_rows[3].scores = {{2361, 0}, {17990, 27}, {}};
  local_rows[4].scores = {{895, 8}, {18029, 179}, {}};
  local_rows[5].scores = {{}, {29295, 355}, {28658, 324}};
  check(choose_joint_offset(drift_candidates, local_rows).value == 42360);
  // Better evidence at a later origin still wins; never default to the first.
  for (auto &row : local_rows)
    row.scores = {{4000, 100}, {4000, 50}, {4000, 0}};
  check(choose_joint_offset(drift_candidates, local_rows).value == 74512);
  for (auto &row : local_rows)
    row.scores = {{4000, 0}, {4000, 0}, {4000, 0}};
  bool local_tied = false;
  try {
    choose_joint_offset(drift_candidates, local_rows);
  } catch (const std::runtime_error &) {
    local_tied = true;
  }
  check(local_tied);
  // Sequence-level contract: recover endpoint drift in either direction but
  // do not reach a distant exact copy, or accept too little aligned sequence.
  for (int delta : {-200, 0, 200}) {
    auto local = score_offset_overlap_local(sequence, reference, 100, 12100, 3100 + delta, 0);
    check(local.overlap == reference.size() - 512 && local.edits == 0);
    auto boundary = align_layout_boundary(sequence, reference, true, true, 3000 + delta, 256, 0);
    check(boundary.query_begin == 3000 && boundary.query_end == 9000);
  }
  check(!score_offset_overlap_local(sequence, reference, 100, 12100, 3700, 0).overlap);
  check(!score_offset_overlap_local(sequence, reference, 100, 12100, 12100, .15).overlap);
  check(!score_offset_overlap_local(sequence, reference.substr(0, 999), 100, 12100, 3100, .15)
             .overlap);
  bool distant = false;
  try {
    align_layout_boundary(sequence, reference, true, true, 3600, 256, 0);
  } catch (const std::runtime_error &) {
    distant = true;
  }
  check(distant);
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
