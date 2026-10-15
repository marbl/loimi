// SPDX-License-Identifier: MIT-0
// Small deterministic contracts for coordinate transforms, evidence selection,
// and packed traces. Throwing checks remain active in Release builds (unlike assert).
#include "core.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <random>
#include <sstream>
using namespace hpc;
static void check(bool v, const char *msg) {
  if (!v)
    throw std::runtime_error(msg);
}
int main() {
  {
    std::string interior;
    for (int i = 0; i < 50; ++i)
      interior += "ACGT";
    std::vector<Alignment> alns(5);
    std::vector<EndRead> reads;
    for (size_t i = 0; i < alns.size(); ++i) {
      auto &a = alns[i];
      a.id = i + 1;
      a.sequence = "CCAATTGGCC" + interior + "GGTTAACCGG";
      a.clip_left = a.clip_right = 10;
      a.begin = 0;
      a.end = interior.size();
      a.ops = cigar(std::string(interior.size(), '='));
      a.packed_ops = true;
      reads.push_back({&a, true, 10, 10});
    }
    auto reference = interior;
    auto added = extend_ends(reference, reads, 1, {}, 5, 100, .15, true, true, true);
    check(added.left == 8 && added.right == 8, "Both ends stop before truncated terminal runs");
    check(reference == "AATTGGCC" + interior + "GGTTAACC", "End extension preserves interior");
    for (auto &a : alns) {
      size_t nm;
      md(a, reference, nm);
      check(nm == 0 && a.clip_left == 2 && a.clip_right == 2, "End traces and clips updated");
    }
    // Four molecules cannot establish a five-molecule extension, even with
    // repeated placement records of one molecule.
    for (auto &a : alns) {
      a.clip_left = a.clip_right = 10;
      a.begin = 0;
      a.end = interior.size();
      a.ops = cigar(std::string(interior.size(), '='));
    }
    alns.back().id = alns.front().id;
    reference = interior;
    added = extend_ends(reference, reads, 1, {}, 5, 100, .15, true, true, true);
    check(added.left == 0 && added.right == 0 && reference == interior,
          "Distinct molecule depth required");
    alns.back().id = 5;
    reads.back().eligible = false;
    added = extend_ends(reference, reads, 1, {}, 5, 100, .15, true, true, true);
    check(added.left == 0 && added.right == 0, "Unknown placements cannot establish ends");
  }
  check(local_hpc_votes(std::string(2000, '='), 500, .15).empty(), "Exact trace fast path");
  check(local_hpc_votes("", 500, .15).empty(), "Empty trace");
  std::string damaged(2000, '=');
  for (size_t i = 1400; i < damaged.size(); i += 3)
    damaged[i] = "XID"[(i / 3) % 3];
  auto local = local_hpc_votes(damaged, 100, .15);
  check(local[500] && !local[1900], "Preserve clean prefix; reject damaged tail");
  auto reversed = damaged;
  std::reverse(reversed.begin(), reversed.end());
  auto backward = local_hpc_votes(reversed, 100, .15);
  std::reverse(backward.begin(), backward.end());
  check(local == backward, "Local quality is orientation symmetric");
  // Compare the rolling implementation with a direct count at every position,
  // including short traces, insertions, deletions, and truncated terminal windows.
  for (size_t radius : {1u, 7u, 5000u}) {
    auto observed = local_hpc_votes(damaged, radius, .15);
    for (size_t i = 0; i < damaged.size(); ++i) {
      size_t b = i > radius ? i - radius : 0;
      size_t e = i + std::min(radius, damaged.size() - i - 1) + 1;
      auto errors =
          std::count_if(damaged.begin() + b, damaged.begin() + e, [](char c) { return c != '='; });
      check(bool(observed[i]) == (damaged[i] == '=' && errors <= .15 * (e - b)),
            "Rolling local quality equals reference implementation");
    }
  }
  std::istringstream model_input("HPC_RUN_MODEL_V1 3 5\n"
                                 "A 1 100 .01 .96 .02 .01\n"
                                 "A 2 100 .01 .25 .4 .34\n"
                                 "A 3 100 .01 .02 .8 .17\n");
  auto model = RunModel::read(model_input);
  check(model.estimate('A', {{2, 5}}, 2) == 3, "Calibrated single-allele likelihood");
  check(model.estimate_single('A', 2, 5) == 3, "Allocation-free single-allele likelihood");
  check(model.estimate_single('A', 2, 4) == 2, "Single-allele low-depth fallback");
  check(model.estimate('A', {{2, 4}, {3, 1}}, 2) == 3, "Calibrated mixed likelihood");
  check(model.estimate('A', {{2, 4}}, 2) == 2, "Low-depth model fallback");
  check(model.estimate('C', {{2, 5}}, 2) == 2, "Unmodeled base fallback");
  check(model.estimate('A', {{4, 5}}, 4) == 4, "Unmodeled true-length fallback");
  for (const std::string &bad :
       {"HPC_RUN_MODEL_V1 3 5\nA 2 1 .1 .2 .3 .4\n", "HPC_RUN_MODEL_V1 3 5\nA 1 1 0 .2 .3 .5\n",
        "HPC_RUN_MODEL_V1 3 5\nA 1 1 .1 .2 .3 .5\n"}) {
    bool invalid = false;
    try {
      std::istringstream input(bad);
      RunModel::read(input);
    } catch (const std::runtime_error &) {
      invalid = true;
    }
    check(invalid, "Malformed run model rejected");
  }
  // Exercise each free-overhang direction independently: swapping edlib query
  // and target must not invert insertion/deletion semantics in our read trace.
  auto left = align_boundary("TGCACGT", "ACGT", true, false, 0.15);
  check(left.query_begin == 3 && left.query_end == 7 && left.ops == "====",
        "Left overhang clipping");
  auto right = align_boundary("ACGTCAC", "ACGT", false, true, 0.15);
  check(right.query_begin == 0 && right.query_end == 4 && right.ops == "====",
        "Right overhang clipping");
  auto both = align_boundary("TGCACGTCAC", "ACGT", true, true, 0.15);
  check(both.query_begin == 3 && both.query_end == 7 && both.ops == "====", "Both overhangs");
  check(locate("ACGT", "TTACGTAA", 0.15) == 2, "Infix location");
  check(locate("ACGT", "ACGTACGT", 0.15) == -1, "Repeated anchor rejected");
  check(offset_cluster({5608, 5609, 5608, 5608}) == 5608, "Offset cluster");
  bool rejected = false;
  try {
    offset_cluster({0, 100, 200});
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  check(rejected, "Conflicting offsets rejected");
  auto r = encode("AAAAACCGTTT");
  check(r.bases == "ACGT" && r.lengths == std::vector<uint32_t>({5, 2, 1, 3}), "RLE");
  check(encode(std::string(1000, 'A')).lengths[0] == 1000, "Long run");
  check(reverse_complement(reverse_complement("AACGTN")) == "AACGTN", "Reverse complement");
  check(tier_of("abc/ccs", {{"hifi", "*/ccs"}, {"other", "*"}}) == 0, "Glob slash matching");
  check(tier_of("abc", {{"first", "a*"}, {"second", "*"}}) == 0, "First rule");
  check(tier_of("xyz", {{"first", "a*"}}) == -1, "Unmatched");
  // The empirical mode is neither an arithmetic average nor pooled tier voting.
  check(call({{{5, 3}, {6, 2}, {20, 1}}}, 1).length == 5, "Mode");
  check(call({{{8, 1}}, {{5, 100}}}, 1).length == 8, "Tier precedence");
  check(call({{}, {{5, 100}}}, 1).tier == 1, "Fallback");
  check(call({{{8, 1}}, {{5, 100}}}, 2).tier == 1, "Support threshold");
  check(call({{{4, 2}, {6, 2}}}, 1).length == 4, "Deterministic tie");
  check(call({{{25, 12}, {26, 9}, {27, 4}, {24, 3}, {23, 1}, {29, 1}}}, 1, true).length == 26,
        "Experimental rounded mean uses the full distribution");
  check(call({{{8, 1}}, {{5, 100}}}, 1, true).length == 8,
        "Experimental mean preserves tier precedence");
  auto insertion = align("ACGAT", "ACGT", 1);
  check(insertion == "===I=", "Edlib insertion orientation");
  auto deletion = align("ACGT", "ACGAT", 1);
  check(deletion == "===D=", "Edlib deletion orientation");
  Alignment a;
  a.sequence = "ACGAT";
  a.ops = insertion;
  a.end = 4;
  size_t nm = 0;
  check(md(a, "ACGT", nm) == "4" && nm == 1, "MD insertion");
  // Crop inside a read, keeping an insertion at the left edge, but not at the
  // right edge. Reference coordinates and query clips must remain independent.
  Alignment cropped;
  cropped.sequence = "ACGAT";
  cropped.ops = "===I=";
  cropped.end = 4;
  check(crop_alignment(cropped, 3, 4) && cropped.begin == 0 && cropped.end == 1 &&
            cropped.clip_left == 3 && cropped.clip_right == 0,
        "Crop insertion at left boundary");
  check(md(cropped, "T", nm) == "1" && nm == 1, "Cropped insertion NM/MD");
  cropped = a;
  check(crop_alignment(cropped, 0, 3) && cropped.clip_right == 2 && md(cropped, "ACG", nm) == "3" &&
            nm == 0,
        "Crop insertion at right boundary");
  cropped.sequence = "ACT";
  cropped.shared_sequence.reset();
  cropped.ops = "==D=";
  cropped.packed_ops = false;
  cropped.begin = 0;
  cropped.end = 4;
  cropped.clip_left = cropped.clip_right = 0;
  check(!crop_alignment(cropped, 2, 3), "Do not export deletion-only crop");
  // Six spanning reads give a high-confidence allele; alternate strands ensure
  // packing the trace does not accidentally change strand-support decisions.
  std::vector<Alignment> reads(6);
  for (size_t i = 0; i < reads.size(); ++i) {
    reads[i].sequence = "ACAT";
    reads[i].ops = "==X=";
    reads[i].end = 4;
    reads[i].reverse = i % 2;
  }
  auto edits = polish("ACGT", reads, 1, 5, 0.9, 20);
  check(edits.size() == 1 && hpc::apply("ACGT", edits) == "ACAT", "Substitution polishing");
  std::vector<const Alignment *> packed;
  for (auto &x : reads) {
    x.ops = cigar(x.ops);
    x.packed_ops = true;
    packed.push_back(&x);
  }
  check(polish("ACGT", packed, 1, 5, 0.9, 20).size() == 1, "Packed polishing without copies");
  for (auto &x : reads) {
    check(md(x, "ACGT", nm) == "2G1" && nm == 1, "Packed MD");
    x.packed_ops = false;
  }
  check(expand_cigar("10=2I3D1X") == std::string(10, '=') + "IIDDDX", "Packed trace roundtrip");
  for (auto &x : reads) {
    x.sequence = "ACCGT";
    x.ops = "==I==";
  }
  check(polish("ACGT", reads, 1, 5, 0.9, 20).empty(), "Do not repolish homopolymers");
  for (auto &x : reads) {
    x.sequence = "ACAGT";
    x.ops = "==I==";
  }
  edits = polish("ACGT", reads, 1, 5, 0.9, 20);
  check(edits.size() == 1 && hpc::apply("ACGT", edits) == "ACAGT" && lift(4, edits) == 5,
        "Non-HPC insertion");
  std::cout << "core tests passed\n";
  // Check both accepted large edits and the first rejected size at the limit.
  for (size_t length : {size_t(200), size_t(500), size_t(501)}) {
    std::string insert;
    for (size_t i = 0; i < length; ++i)
      insert += (i % 2 ? 'T' : 'A');
    for (auto &x : reads) {
      x.sequence = "AC" + insert + "GT";
      x.ops = "==" + std::string(length, 'I') + "==";
      x.end = 4;
    }
    auto large = polish("ACGT", reads, 1, 5, 0.9, 500);
    check(large.size() == (length <= 500 ? 1 : 0), "Large insertion size limit");
    if (!large.empty())
      check(hpc::apply("ACGT", large) == reads[0].sequence, "Large insertion application");
    std::string ref = "AC" + insert + "GT";
    for (auto &x : reads) {
      x.sequence = "ACGT";
      x.ops = "==" + std::string(length, 'D') + "==";
      x.end = ref.size();
    }
    large = polish(ref, reads, 1, 5, 0.9, 500);
    check(large.size() == (length <= 500 ? 1 : 0), "Large deletion size limit");
    if (!large.empty())
      check(hpc::apply(ref, large) == "ACGT", "Large deletion application");
  }
  // Differential depth check: indexed spanning counts must match a brute-force scan.
  std::mt19937 rng(123);
  std::string reference;
  for (size_t i = 0; i < 200; ++i)
    reference += "ACGT"[i % 4];
  for (size_t trial = 0; trial < 100; ++trial) {
    std::vector<Alignment> sample;
    unsigned depth[2][2]{}, votes[2][2]{};
    for (size_t i = 0; i < 40; ++i) {
      Alignment a;
      a.begin = rng() % 105;
      a.end = a.begin + 1 + rng() % (200 - a.begin);
      a.tier = rng() % 2;
      a.reverse = rng() % 2;
      a.sequence = reference.substr(a.begin, a.end - a.begin);
      a.ops = std::string(a.sequence.size(), '=');
      bool supports = rng() % 5 != 0;
      if (a.begin <= 100 && a.end > 100 && supports) {
        a.sequence[100 - a.begin] = 'G';
        a.ops[100 - a.begin] = 'X';
      }
      if (a.begin < 100 && a.end > 101) {
        ++depth[a.tier][a.reverse];
        if (supports)
          ++votes[a.tier][a.reverse];
      }
      a.ops = cigar(a.ops);
      a.packed_ops = true;
      sample.push_back(std::move(a));
    }
    bool expect = false;
    int tier = -1;
    for (int t = 0; t < 2; ++t)
      if (depth[t][0] + depth[t][1] >= 3) {
        auto n = votes[t][0] + votes[t][1];
        expect = n >= 3 && double(n) / (depth[t][0] + depth[t][1]) >= 0.6 &&
                 (!depth[t][0] || !depth[t][1] || (votes[t][0] && votes[t][1]));
        tier = t;
        break;
      }
    auto result = polish(reference, sample, 2, 3, 0.6, 500);
    check(result.size() == size_t(expect), "Indexed depth agrees with brute force");
    if (expect)
      check(result[0].tier == tier && result[0].depth == depth[tier][0] + depth[tier][1],
            "Indexed tier depth");
  }
}
