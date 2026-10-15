// SPDX-License-Identifier: MIT-0
#pragma once
#include <cstdint>
#include <map>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

namespace hpc {
// Plain FASTA, uppercase bases, first header token as ID. Reject duplicate IDs.
// An empty optional path returns an empty map.
std::map<std::string, std::string> fasta(const std::string &path);
// Authoritative internal ID -> template name, accepting Verkko's table header.
std::map<uint32_t, std::string> template_names(const std::string &path);
// Execute literal argv (never a shell); throw if the process fails.
void run(std::vector<std::string> args);
// Wrapped at 80 bases; sequence is borrowed, not copied.
void write_fasta(std::ostream &output, const std::string &name, std::string_view sequence);
} // namespace hpc
