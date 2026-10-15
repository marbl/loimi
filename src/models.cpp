// SPDX-License-Identifier: MIT-0
#include "core.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace hpc {
RunModel RunModel::read(std::istream &input) {
  RunModel model;
  std::string line;
  bool header = false;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream row(line);
    if (!header) {
      std::string magic, extra;
      if (!(row >> magic >> model.observed_max >> model.minimum_depth) ||
          magic != "HPC_RUN_MODEL_V1" || model.observed_max < 1 || model.observed_max > 1000 ||
          model.minimum_depth < 1 || model.minimum_depth > 1000000 || (row >> extra))
        throw std::runtime_error("Invalid run model header");
      header = true;
      continue;
    }
    char base;
    size_t truth, loci;
    if (!(row >> base >> truth >> loci) || std::string("ACGT").find(base) == std::string::npos ||
        truth < 1 || truth > 1000 || loci < 1)
      throw std::runtime_error("Invalid run model row");
    auto &rows = model.log_probabilities[base];
    if (truth != rows.size() + 1)
      throw std::runtime_error("Run model lengths must be contiguous from one per base");
    std::vector<double> probabilities(model.observed_max + 1);
    double sum = 0;
    for (double &p : probabilities) {
      if (!(row >> p) || !std::isfinite(p) || p <= 0 || p > 1)
        throw std::runtime_error("Run model probabilities must be finite and positive");
      sum += p;
      p = std::log(p);
    }
    std::string extra;
    if ((row >> extra) || std::abs(sum - 1) > 1e-6)
      throw std::runtime_error("Invalid run model probability row size or normalization");
    rows.push_back(std::move(probabilities));
  }
  if (!input.eof() || !header || model.log_probabilities.empty())
    throw std::runtime_error("Empty or unreadable run model");
  for (const auto &[base, rows] : model.log_probabilities) {
    auto &calls = model.single_allele_calls[base];
    calls.resize(model.observed_max + 1);
    for (size_t observed = 0; observed <= model.observed_max; ++observed) {
      double best = -std::numeric_limits<double>::infinity();
      for (size_t length = 0; length < rows.size(); ++length)
        if (rows[length][observed] > best) {
          best = rows[length][observed];
          calls[observed] = length + 1;
        }
    }
  }
  return model;
}

uint32_t RunModel::estimate_single(char base, uint32_t observed, uint32_t depth) const {
  // Avoid constructing an allocating map for the common single-allele run.
  auto found = log_probabilities.find(base);
  if (depth < minimum_depth || found == log_probabilities.end() || observed < 1 ||
      observed > found->second.size())
    return observed;
  return single_allele_calls.at(base)[std::min<size_t>(observed, observed_max)];
}

uint32_t RunModel::estimate(char base, const Histogram &observations, uint32_t fallback) const {
  auto found = log_probabilities.find(base);
  if (found == log_probabilities.end() || fallback < 1 || fallback > found->second.size())
    return fallback;
  uint64_t depth = 0;
  for (auto [length, count] : observations)
    depth += count;
  if (depth < minimum_depth)
    return fallback;
  if (observations.size() == 1)
    return single_allele_calls.at(
        base)[std::min<size_t>(observations.begin()->first, observed_max)];
  double best = -std::numeric_limits<double>::infinity();
  uint32_t result = fallback;
  for (size_t length = 0; length < found->second.size(); ++length) {
    double score = 0;
    for (auto [observed, count] : observations)
      score += count * found->second[length][std::min<size_t>(observed, observed_max)];
    // Rows are ascending: exact score ties deterministically favor shorter runs.
    if (score > best) {
      best = score;
      result = length + 1;
    }
  }
  return result;
}

} // namespace hpc
