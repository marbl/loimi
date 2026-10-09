// SPDX-License-Identifier: MIT-0
// Replay histograms through the production predictor without repeating alignments.
// Input rows: BASE LENGTH:COUNT [LENGTH:COUNT ...]. Output: MODE CALIBRATED.
#include "core.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

int main(int argc, char **argv) {
  try {
    if (argc != 2)
      throw std::runtime_error("Usage: apply-run-model MODEL < histograms.txt");
    std::ifstream input(argv[1]);
    if (!input)
      throw std::runtime_error("Cannot open model");
    auto model = hpc::RunModel::read(input);
    std::string line;
    while (std::getline(std::cin, line)) {
      std::istringstream row(line);
      char base;
      if (!(row >> base))
        throw std::runtime_error("Missing base");
      hpc::Histogram histogram;
      std::string allele;
      while (row >> allele) {
        auto colon = allele.find(':');
        if (colon == std::string::npos)
          throw std::runtime_error("Expected LENGTH:COUNT");
        size_t used = 0;
        auto length = std::stoull(allele.substr(0, colon), &used);
        if (used != colon || length < 1 || length > UINT32_MAX)
          throw std::runtime_error("Invalid length");
        auto count = std::stoull(allele.substr(colon + 1), &used);
        if (used != allele.size() - colon - 1 || count < 1 || count > UINT32_MAX ||
            !histogram.emplace(length, count).second)
          throw std::runtime_error("Invalid count or duplicate length");
      }
      if (histogram.empty())
        throw std::runtime_error("Empty histogram");
      auto modal = hpc::call({histogram}, 1).length;
      std::cout << modal << '\t' << model.estimate(base, histogram, modal) << '\n';
    }
    if (!std::cin.eof())
      throw std::runtime_error("Cannot read histograms");
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
