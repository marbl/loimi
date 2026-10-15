// SPDX-License-Identifier: MIT-0
#pragma once
#include "core.hpp"

namespace hpc {
// Shared strict numeric readers for CLI arguments and diagnostic input tables.
// Positive integers fit uint32_t; coordinates fit both size_t and int64_t.
unsigned parse_positive_integer(const std::string &text);
size_t parse_coordinate(const std::string &text);
struct Options {
  // These defaults are part of the CLI contract; update help, README, and tests
  // together whenever changing an algorithmic threshold or resource policy.
  std::string input, templates, tig_map, output, cram;
  std::string evidence_sites, evidence_exclude_reads;
  std::map<std::string, std::string> model_paths;
  std::vector<Tier> tiers;
  unsigned threads = 1, min_support = 5, polish_support = 5;
  size_t batch_bases = 8000000;
  double max_error = 0.15, polish_fraction = 0.9;
  double local_hpc_error = 0.15;
  size_t local_hpc_radius = 500;
  bool filter_hpc = true;
  size_t max_indel = 500;
  bool polish = false, drop_unmatched = false, missing_error = false;
  bool extend = false, no_extend = false, extend_left = true, extend_right = true;
  unsigned end_support = 5;
  size_t max_end_extension = 10000;
  bool experiment_mean = false, ambiguous_fallback = true;
  bool region = false;
  uint32_t region_tig = 0;
  size_t region_begin = 0, region_end = 0, region_flank = 100000;
};
// Parse and validate CLI configuration. --help prints usage and exits successfully.
// All other failures throw before any package is read or output is created.
Options parse_options(int argc, char **argv);
} // namespace hpc
