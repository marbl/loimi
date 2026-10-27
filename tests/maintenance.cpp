// SPDX-License-Identifier: MIT-0
// Differential tests keep optimized code tied to simple reference algorithms.
#include "core.hpp"
#include "options.hpp"
#include "io.hpp"
#include "parallel.hpp"
#include <random>
#include <sstream>
#include <stdexcept>

static void check(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}
template <class F> static void rejects(F action) {
  bool failed = false;
  try {
    action();
  } catch (const std::exception &) {
    failed = true;
  }
  check(failed, "Expected invalid input to be rejected");
}
int main() {
  using namespace hpc;
  auto options = [](std::vector<std::string> extra) {
    std::vector<std::string> arguments{"hpc-consensus", "--import", "unused.cnspack", "--output",
                                       "unused",        "--tier",   "all=*"};
    arguments.insert(arguments.end(), extra.begin(), extra.end());
    std::vector<char *> argv;
    for (auto &arg : arguments)
      argv.push_back(arg.data());
    return parse_options(int(argv.size()), argv.data());
  };
  auto fast = options({});
  check(fast.hp_flank_bases == 500 && !fast.hp_realign && !fast.polish && !fast.extend,
        "Fast mode filters flanks without enabling raw realignment or extension");
  auto defaults = options({"--polish-non-hpc"});
  check(defaults.hp_realign && defaults.hp_flank_bases == 500 && defaults.linked_events &&
            !defaults.no_extend,
        "Polishing defaults include HP refinement, linked events and extension");
  auto disabled = options({"--polish-non-hpc", "--no-hp-realign", "--no-hp-flank-filter",
                           "--no-linked-events", "--no-extend-ends"});
  check(!disabled.hp_realign && !disabled.hp_flank_bases && !disabled.linked_events &&
            disabled.no_extend,
        "Each promoted feature can be independently disabled");
  auto explicit_flags =
      options({"--polish-non-hpc", "--no-hp-realign", "--hp-realign", "--no-hp-flank-filter",
               "--hp-flank-bases", "300", "--no-linked-events", "--linked-events"});
  check(explicit_flags.hp_realign && explicit_flags.hp_flank_bases == 300 &&
            explicit_flags.linked_events,
        "Explicit enable switches and custom flank values override earlier disables");
  for (const auto &removed : {"--experiment-run-mean", "--experiment-consistent-ambiguous",
                              "--experiment-retain-higher-unique",
                              "--experiment-terminal-hp-filter", "--experiment-hp-realign",
                              "--experiment-hp-flank-bases", "--experiment-ambiguous-fallback"})
    rejects([&] { options({removed}); });
  std::vector<unsigned> visited(1000);
  parallel(visited.size(), 4, [&](size_t i) { ++visited[i]; });
  for (auto count : visited)
    check(count == 1, "Workers process each slot exactly once");
  bool allocation_failed = false;
  try {
    parallel(100, 4, [](size_t i) {
      if (i == 3)
        throw std::bad_alloc();
    });
  } catch (const std::bad_alloc &) {
    allocation_failed = true;
  }
  check(allocation_failed, "Worker allocation failures reach caller after joining");
  parallel(0, 4, [](size_t) { throw std::runtime_error("No empty work"); });
  std::mt19937 random(271828);
  for (size_t trial = 0; trial < 1000; ++trial) {
    // Periodic references exercise many rotations, including a shift that is
    // not a multiple of the inserted motif length and boundary position zero.
    std::string motif, reference;
    for (size_t j = 0, n = 1 + random() % 15; j < n; ++j)
      motif += "ACGT"[random() % 4];
    for (size_t j = 0; j < 100; ++j)
      reference += motif;
    size_t position = random() % (reference.size() + 1);
    std::string allele = reference.substr(position, motif.size());
    if (allele.empty())
      allele = motif;
    Alignment alignment;
    alignment.sequence = reference.substr(0, position) + allele + reference.substr(position);
    alignment.ops = std::string(position, '=') + std::string(allele.size(), 'I') +
                    std::string(reference.size() - position, '=');
    alignment.end = reference.size();
    size_t expected = position;
    while (expected && reference[expected - 1] == allele.back()) {
      allele = allele.back() + allele.substr(0, allele.size() - 1);
      --expected;
    }
    if (trial % 2) {
      alignment.ops = cigar(alignment.ops);
      alignment.packed_ops = true;
    }
    auto observed = events(alignment, reference);
    check(observed.size() == 1 && observed[0].begin == expected && observed[0].end == expected &&
              observed[0].alt == allele,
          "Allocation-free normalization agrees with legacy rotation");
    std::vector<Accepted> edits;
    size_t cursor = 0;
    std::string legacy;
    while (cursor + 20 < reference.size()) {
      size_t begin = cursor + random() % 10;
      size_t end = begin + random() % 10;
      std::string alt(random() % 10, 'G');
      edits.push_back({{begin, end, alt}, 0, 5, 5});
      legacy += reference.substr(cursor, begin - cursor) + alt;
      cursor = end;
      // Ensure progress even for an insertion at cursor.
      if (cursor == begin && begin == 0)
        break;
      if (edits.size() > 200)
        break;
    }
    legacy += reference.substr(cursor);
    check(hpc::apply(reference, edits) == legacy,
          "Direct append agrees with legacy edit application");
    check(hpc::apply(reference, {}) == reference, "Empty edits preserve reference");
  }
  rejects([] { hpc::apply("AC", {{{2, 1, ""}, 0, 1, 1}}); });
  rejects([] { hpc::apply("AC", {{{0, 3, ""}, 0, 1, 1}}); });
  check(parse_coordinate("0") == 0 && parse_positive_integer("5") == 5, "Numeric parsers");
  for (auto s : {"-1", "1x", "", "18446744073709551616"})
    rejects([&] { parse_coordinate(s); });
  rejects([] { parse_positive_integer("0"); });
  rejects([] { run({}); });
  std::ostringstream fasta_output;
  write_fasta(fasta_output, "test", std::string(81, 'A'));
  check(fasta_output.str() == ">test\n" + std::string(80, 'A') + "\nA\n", "FASTA wrapping");
}
