// SPDX-License-Identifier: MIT-0
#include "hp_refinement.hpp"
#include <sstream>
#include <stdexcept>
using namespace hpc;
static void check(bool ok) {
  if (!ok)
    throw std::runtime_error("Refinement regression");
}
int main() {
  std::istringstream calibration("HPC_RUN_MODEL_V1 4 5\nG 1 1 .01 .96 .01 .01 .01\n"
                                 "G 2 1 .01 .01 .55 .42 .01\nG 3 1 .01 .01 .42 .55 .01\n");
  auto model = RunModel::read(calibration);
  check(model.uncertain('G', {{2, 3}, {3, 2}}));
  check(!model.uncertain('G', {{2, 50}}));
  check(RunModel{}.uncertain('A', {{2, 3}, {3, 2}}));
  check(!RunModel{}.uncertain('A', {{2, 5}}));
  std::string left, right;
  for (size_t i = 0; i < 25; ++i) {
    left += "ACTG";
    right += "TCAG";
  }
  left.back() = 'C'; // immediate anchor must not join the G run
  std::string ref = left + "GG" + right;
  std::vector<Alignment> reads(10);
  std::vector<const Alignment *> ptrs;
  for (size_t i = 0; i < reads.size(); ++i) {
    auto &a = reads[i];
    a.id = i + 1;
    a.begin = 0;
    a.end = ref.size();
    a.sequence = left + "GGG" + right;
    a.ops = cigar(align(a.sequence, ref, .15));
    a.packed_ops = true;
    a.reverse = i % 2;
    ptrs.push_back(&a);
  }
  std::vector<HpSite> sites{{100, 102, 'G'}};
  auto call_hp = [&] { return realign_homopolymers(ref, ptrs, sites, 2, 5, {}, true, 500); };
  auto edits = call_hp();
  check(edits.size() == 1 && edits[0].edit.alt == "GGG" && edits[0].tier == 0);
  check(hpc::apply(ref, edits) == left + "GGG" + right);
  // The flank rule applies to rebuilt evidence, not only initial voting.
  check(realign_homopolymers(ref, ptrs, sites, 2, 5, {}, true, 500, 500).empty());
  // A qualifying lower UNIQUE tier outranks non-unique high-tier observations.
  for (size_t i = 0; i < 10; ++i) {
    reads[i].tier = i < 5 ? 2 : 1;
    if (i >= 5) {
      reads[i].sequence = ref;
      reads[i].ops = std::to_string(ref.size()) + "=";
    }
  }
  check(call_hp().empty());
  for (auto &a : reads)
    a.tier = 2;
  check(call_hp().empty()); // 5/5 mode tie chooses existing length 2
  for (auto &a : reads) {
    a.sequence = left + "GGG" + right;
    a.ops = cigar(align(a.sequence, ref, .15));
  }
  check(call_hp().size() == 1); // pooled fallback is permitted
  check(realign_homopolymers(ref, ptrs, sites, 2, 5, {}, false, 500).empty());
  ptrs.resize(4);
  check(call_hp().empty()); // no sparse post-polishing edits
  ptrs.clear();
  check(call_hp().empty());
  // Read clipping deprives the target of its required local context.
  reads[0].begin = 50;
  reads[0].end = ref.size();
  reads[0].clip_left = 50;
  reads[0].ops = cigar(align(reads[0].sequence.substr(50), ref.substr(50), .15));
  ptrs.assign(5, &reads[0]);
  check(call_hp().empty());
}
