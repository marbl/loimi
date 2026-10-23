// SPDX-License-Identifier: MIT-0
#include "core.hpp"
#include <random>
#include <stdexcept>

static void check(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
int main() {
  using namespace hpc;
  const std::string ref = "GATGCC", truth = "GACTATTGCC";
  std::vector<Alignment> reads(6);
  for (size_t i = 0; i < reads.size(); ++i) {
    auto &a = reads[i];
    a.sequence = truth;
    a.end = ref.size();
    a.reverse = i % 2;
    // Same molecule sequence, represented as one insertion or two insertions
    // separated by a matching T. Both pieces normalize onto the same site.
    a.ops = i < 3 ? "2=1I1=3I3=" : "2=4I4=";
    a.packed_ops = true;
    size_t nm;
    md(a, ref, nm);
    auto es = events(a, ref);
    check(es.size() == 1 && es[0].begin == 2 && es[0].end == 2 && es[0].alt == "CTAT",
          "Equivalent linked and contiguous insertions have one canonical allele");
  }
  auto edits = polish(ref, reads, 1, 5, .9, 500);
  check(polish(ref, reads, 1, 5, .9, 500, 0, true).empty(),
        "Frozen master comparison must not silently gain linked-event support");
  check(edits.size() == 1 && edits[0].support == 6 && hpc::apply(ref, edits) == truth,
        "Linked support survives exact allele voting and conflict filtering");
  check(polish(ref, reads, 1, 5, .9, 3).empty(), "Combined allele obeys size cap");
  for (size_t i = 0; i < reads.size(); ++i) {
    reads[i].sequence = i < 3 ? "GACTGCC" : "GATATTGCC";
    reads[i].ops = i < 3 ? "2=1I4=" : "3=3I3=";
  }
  check(polish(ref, reads, 1, 5, .9, 500).empty(),
        "Do not assemble haplotypes from unrelated reads carrying different pieces");
  Alignment independent;
  independent.sequence = "TCGTTCGT";
  independent.ops = "X===X===";
  independent.end = 8;
  check(events(independent, "ACGTACGT").size() == 2,
        "Nearby changes without a normalization collision remain independent");

  // Property test: normalized, linked replacements must reproduce the EXACT
  // sequence implied by a trace, including noisy repeat runs and replacements.
  std::mt19937 rng(8191);
  for (size_t trial = 0; trial < 2000; ++trial) {
    std::string reference;
    for (size_t p = 0; p < 100; ++p)
      reference += "ACGT"[rng() % 4];
    Alignment a;
    a.end = reference.size();
    for (char base : reference) {
      unsigned choice = rng() % 8;
      if (choice == 0) {
        a.sequence += "ACGT"[rng() % 4];
        a.ops += 'I';
      }
      if (choice == 1)
        a.ops += 'D';
      else if (choice == 2) {
        char alt = base == 'A' ? 'C' : 'A';
        a.sequence += alt;
        a.ops += 'X';
      } else {
        a.sequence += base;
        a.ops += '=';
      }
    }
    if (trial % 2) {
      a.ops = cigar(a.ops);
      a.packed_ops = true;
    }
    std::vector<Accepted> replacements;
    for (auto e : events(a, reference))
      replacements.push_back({e, 0, 1, 1});
    check(hpc::apply(reference, replacements) == a.sequence,
          "Linked normalization preserves the complete read haplotype");
  }
}
