// SPDX-License-Identifier: MIT-0
#include "core.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <iostream>
#include <cstdlib>
#include "options.hpp"

namespace hpc {
static void require(bool b, const std::string &s) {
  if (!b)
    throw std::runtime_error(s);
}
unsigned parse_positive_integer(const std::string &s) {
  size_t used = 0;
  auto n = std::stoul(s, &used);
  require(used == s.size() && n > 0 && n <= UINT32_MAX, "Expected positive integer: " + s);
  return n;
}
size_t parse_coordinate(const std::string &s) {
  require(!s.empty() && s.find_first_not_of("0123456789") == std::string::npos,
          "Expected nonnegative coordinate: " + s);
  auto n = std::stoull(s);
  require(n <= SIZE_MAX && n <= INT64_MAX, "Coordinate out of range: " + s);
  return size_t(n);
}
static void parse_region(Options &o, const std::string &s) {
  require(!o.region, "Specify --region only once");
  auto colon = s.find(':');
  auto dash = s.find('-', colon == std::string::npos ? 0 : colon + 1);
  require(s.rfind("tig", 0) == 0 && colon != std::string::npos && dash != std::string::npos,
          "Expected --region tig00000109:0-10000000 (draft bases, zero-based half-open)");
  auto id = parse_coordinate(s.substr(3, colon - 3));
  require(id <= UINT32_MAX, "Region tig ID out of range");
  o.region_tig = uint32_t(id);
  o.region_begin = parse_coordinate(s.substr(colon + 1, dash - colon - 1));
  o.region_end = parse_coordinate(s.substr(dash + 1));
  require(o.region_begin < o.region_end, "Region must have begin < end");
  o.region = true;
}
static double probability(const std::string &s) {
  size_t used = 0;
  double x = std::stod(s, &used);
  require(used == s.size() && std::isfinite(x) && x > 0 && x <= 1,
          "Expected probability in (0,1]: " + s);
  return x;
}
Options parse_options(int argc, char **argv) {
  Options o;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto value = [&]() {
      require(i + 1 < argc, "Missing value for " + a);
      return std::string(argv[++i]);
    };
    if (a == "call")
      continue;
    if (a == "--help" || a == "-h") {
      std::cout
          << "hpc-consensus 0.1.0\nUsage: hpc-consensus --import part.cnspack --template-hpc "
             "paths.fa\n"
             "  --tier 'hifi=hifi_*' --tier 'other=*' --output PREFIX [options]\n"
             "--template-hpc FILE   Required when package has no embedded HPC bases\n"
             "Template offsets are inferred per tig; singleton reads pass through directly.\n"
             "--tig-map FILE        Template association: internal_ID FASTA_name per line\n"
             "Output names are tig + eight-digit packaged tig ID (e.g. tig00000109).\n"
             "--region TIG:BEGIN-END  One tig, zero-based half-open uncompressed draft bases\n"
             "--region-flank N      Context in bases/runs (default 100000; region mode only)\n"
             "--experiment-run-mean  Diagnostic rounded-mean run caller (default: mode)\n"
             "--evidence-sites FILE  Export all-tier HPC histograms at TIG RUN rows (zero-based)\n"
             "--evidence-exclude-reads FILE  Read IDs excluded only from diagnostic histograms\n"
             "--unique-only        Disable default fallback to multiply placed reads\n"
             "--ambiguous-fallback Pool within tiers if no unique tier qualifies (default)\n"
             "--tier [LABEL=]GLOB   Repeated in priority order; first matching glob wins\n"
             "--tier-model LABEL=FILE  Optional calibrated run-length model for one tier\n"
             "--unmatched error|drop (default error)\n--min-run-support N (default 5)\n"
             "Try unique tiers, then pooled tiers; sparse evidence is retained as last resort.\n"
             "--missing-run one|error (default one)\n--max-error F (default 0.15)\n"
             "--hpc-filter-radius N  Local vote window radius in HPC operations (default 500)\n"
             "--hpc-filter-error F   Maximum local X/I/D fraction (default 0.15)\n"
             "--no-hpc-filter        Disable local vote filtering for comparisons\n"
             "--threads N (default 1; parallel alignments within each tig)\n"
             "--batch-bases N (default 8000000; raw bases per HPC alignment batch)\n"
             "--polish-non-hpc      Optional conservative correction\n"
             "--extend-ends         Coverage-supported ends (also enabled by polishing)\n"
             "--no-extend-ends      Preserve template boundaries, including with polishing\n"
             "--min-end-support N  Minimum molecules at every added run (default 5)\n"
             "--max-end-extension N  Maximum raw bases per end (default 10000)\n"
             "--end-extension-side both|left|right (default both)\n"
             "--min-polish-support N (default 5)\n--min-polish-fraction F (default 0.9)\n"
             "--max-polish-indel N (default 500)\n"
             "--cram FILE          Optional embedded-reference CRAM; requires samtools on PATH\n"
             "Input packages must match the linked Verkko store ABI. Outputs must not exist.\n";
      std::exit(0);
    } else if (a == "--import")
      o.input = value();
    else if (a == "--template-hpc")
      o.templates = value();
    else if (a == "--tig-map")
      o.tig_map = value();
    else if (a == "--region")
      parse_region(o, value());
    else if (a == "--region-flank")
      o.region_flank = parse_coordinate(value());
    else if (a == "--output")
      o.output = value();
    else if (a == "--cram")
      o.cram = value();
    else if (a == "--extend-ends")
      o.extend = true;
    else if (a == "--no-extend-ends")
      o.no_extend = true;
    else if (a == "--min-end-support")
      o.end_support = parse_positive_integer(value());
    else if (a == "--max-end-extension")
      o.max_end_extension = parse_positive_integer(value());
    else if (a == "--end-extension-side") {
      auto side = value();
      require(side == "both" || side == "left" || side == "right", "Invalid extension side");
      o.extend_left = side != "right";
      o.extend_right = side != "left";
    } else if (a == "--polish-non-hpc")
      o.polish = true;
    else if (a == "--experiment-run-mean")
      o.experiment_mean = true;
    else if (a == "--evidence-sites")
      o.evidence_sites = value();
    else if (a == "--evidence-exclude-reads")
      o.evidence_exclude_reads = value();
    else if (a == "--tier-model") {
      auto v = value();
      auto equal = v.find('=');
      require(equal != std::string::npos && equal > 0 && equal + 1 < v.size(),
              "Expected --tier-model LABEL=FILE");
      require(o.model_paths.emplace(v.substr(0, equal), v.substr(equal + 1)).second,
              "Duplicate tier model");
    } else if (a == "--ambiguous-fallback" || a == "--experiment-ambiguous-fallback")
      o.ambiguous_fallback = true;
    else if (a == "--unique-only")
      o.ambiguous_fallback = false;
    else if (a == "--threads")
      o.threads = parse_positive_integer(value());
    else if (a == "--batch-bases")
      o.batch_bases = parse_positive_integer(value());
    else if (a == "--min-run-support")
      o.min_support = parse_positive_integer(value());
    else if (a == "--min-polish-support")
      o.polish_support = parse_positive_integer(value());
    else if (a == "--max-polish-indel")
      o.max_indel = parse_positive_integer(value());
    else if (a == "--hpc-filter-radius")
      o.local_hpc_radius = parse_positive_integer(value());
    else if (a == "--hpc-filter-error")
      o.local_hpc_error = probability(value());
    else if (a == "--no-hpc-filter")
      o.filter_hpc = false;
    else if (a == "--max-error")
      o.max_error = probability(value());
    else if (a == "--min-polish-fraction")
      o.polish_fraction = probability(value());
    else if (a == "--unmatched") {
      auto v = value();
      require(v == "error" || v == "drop", "Invalid unmatched policy");
      o.drop_unmatched = v == "drop";
    } else if (a == "--missing-run") {
      auto v = value();
      require(v == "one" || v == "error", "Invalid missing-run policy");
      o.missing_error = v == "error";
    } else if (a == "--tier") {
      auto v = value();
      auto p = v.find('=');
      Tier t = p == std::string::npos ? Tier{"tier" + std::to_string(o.tiers.size()), v}
                                      : Tier{v.substr(0, p), v.substr(p + 1)};
      require(!t.label.empty() && !t.pattern.empty(), "Empty tier label/pattern");
      require(t.label.find_first_not_of(
                  "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-") ==
                  std::string::npos,
              "Invalid tier label");
      for (const auto &old : o.tiers)
        require(old.label != t.label, "Duplicate tier label");
      o.tiers.push_back(t);
    } else
      throw std::runtime_error("Unknown option: " + a);
  }
  require(!o.input.empty() && !o.output.empty() && !o.tiers.empty(),
          "Require --import, --output and at least one --tier; see --help");
  require(o.polish_fraction > 0.5, "Polishing fraction must exceed 0.5");
  require(o.model_paths.empty() || !o.experiment_mean,
          "--tier-model cannot be combined with --experiment-run-mean");
  for (const auto &[label, path] : o.model_paths)
    require(std::any_of(o.tiers.begin(), o.tiers.end(),
                        [&](const auto &tier) { return tier.label == label; }),
            "Model refers to unknown tier: " + label);
  require(o.evidence_exclude_reads.empty() || !o.evidence_sites.empty(),
          "--evidence-exclude-reads requires --evidence-sites");
  require(
      o.evidence_sites.empty() || !o.region,
      "Evidence export uses complete template coordinates; extract a fixture instead of --region");
  require(!o.region || o.region_end <= size_t(INT64_MAX) - o.region_flank,
          "Region plus flank overflows coordinate range");
  return o;
}
} // namespace hpc
