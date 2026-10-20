// SPDX-License-Identifier: MIT-0
// Diagnostic parity interface: BASE CURRENT UNIQUE_HIST ; AMBIGUOUS_HIST.
#include "core.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
int main(int argc, char **argv) {
  try {
    if (argc != 2)
      throw std::runtime_error("Usage: apply-borrow MODEL_OR_DASH");
    hpc::RunModel model;
    if (std::string(argv[1]) != "-") {
      std::ifstream in(argv[1]);
      model = hpc::RunModel::read(in);
    }
    std::string line;
    while (std::getline(std::cin, line)) {
      std::istringstream in(line);
      char base;
      uint32_t current;
      in >> base >> current;
      hpc::Histogram u, a;
      bool ambiguous = false;
      std::string token;
      while (in >> token) {
        if (token == ";") {
          ambiguous = true;
          continue;
        }
        auto colon = token.find(':');
        if (colon == std::string::npos)
          throw std::runtime_error("Expected length:count");
        (ambiguous ? a : u)[std::stoul(token.substr(0, colon))] +=
            std::stoul(token.substr(colon + 1));
      }
      std::cout << model.uncertain(base, u) << ' ' << model.borrow(base, u, a, current) << '\n';
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
