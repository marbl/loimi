// SPDX-License-Identifier: MIT-0
#pragma once
#include "core.hpp"
namespace hpc::experimental {
// Frozen master event extraction for controlled comparisons. Unlike the branch
// extractor, independently normalized neighboring events are not joined.
std::vector<Edit> master_events(const Alignment &alignment, const std::string &reference);
} // namespace hpc::experimental
