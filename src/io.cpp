// SPDX-License-Identifier: MIT-0
// Text formats and subprocess boundaries; no consensus policy belongs here.
#include "io.hpp"
#include <cctype>
#include <cerrno>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>

namespace hpc {
namespace {
void require(bool value, const std::string &message) {
  if (!value)
    throw std::runtime_error(message);
}
} // namespace
std::map<std::string, std::string> fasta(const std::string &path) {
  // Read plain FASTA only. Header descriptions are ignored, while duplicate IDs
  // are errors: silently replacing a template could select the wrong chromosome.
  std::map<std::string, std::string> out;
  if (path.empty())
    return out;
  std::ifstream in(path);
  require(bool(in), "Cannot open template " + path);
  std::string line, key;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    if (line[0] == '>') {
      std::istringstream fields(line.substr(1));
      fields >> key;
      require(!key.empty() && !out.count(key), "Duplicate/empty template name");
      out[key] = {};
    } else {
      require(!key.empty(), "Sequence before FASTA header");
      for (unsigned char c : line)
        if (!std::isspace(c))
          out[key] += std::toupper(c);
    }
  }
  return out;
}
// Keep template names independent of output names. A supplied map is authoritative:
// never infer an ID from a piece suffix when the package has its own numbering.
std::map<uint32_t, std::string> template_names(const std::string &path) {
  std::map<uint32_t, std::string> names;
  if (path.empty())
    return names;
  std::ifstream in(path);
  require(bool(in), "Cannot open tig map " + path);
  std::set<std::string> seen;
  std::string line;
  size_t line_number = 0;
  while (std::getline(in, line)) {
    ++line_number;
    std::istringstream fields(line);
    std::string id, name, extra;
    if (!(fields >> id) || id[0] == '#')
      continue;
    const auto context = " in tig map at line " + std::to_string(line_number);
    // Verkko emits a human-readable table heading. Recognize its exact tokens
    // only before data, so malformed mapping rows are still reported as errors.
    if (names.empty() && id == "tigID") {
      std::string trailing;
      require(bool(fields >> name >> extra) && name == "tig" && extra == "name" &&
                  !(fields >> trailing),
              "Invalid tig map heading" + context);
      continue;
    }
    if (names.empty() && id.find_first_not_of('-') == std::string::npos) {
      require(bool(fields >> name) && name.find_first_not_of('-') == std::string::npos &&
                  !(fields >> extra),
              "Invalid tig map separator" + context);
      continue;
    }
    require(bool(fields >> name) && !(fields >> extra), "Expected ID NAME" + context);
    require(id.find_first_not_of("0123456789") == std::string::npos,
            "Invalid internal tig ID" + context);
    auto value = std::stoull(id);
    require(value <= UINT32_MAX, "Internal tig ID out of range" + context);
    require(names.emplace(uint32_t(value), name).second, "Duplicate internal tig ID" + context);
    require(seen.insert(name).second, "Duplicate template name" + context);
  }
  return names;
}
void run(std::vector<std::string> args) {
  require(!args.empty(), "Cannot execute an empty command");
  // No shell interpolation: filenames and tier-derived strings remain arguments.
  std::vector<char *> ptr;
  for (auto &a : args)
    ptr.push_back(a.data());
  ptr.push_back(nullptr);
  pid_t pid = fork();
  require(pid >= 0, "fork failed");
  if (pid == 0) {
    execvp(ptr[0], ptr.data());
    perror(ptr[0]);
    _exit(127);
  }
  int status = 0;
  pid_t waited;
  do {
    waited = waitpid(pid, &status, 0);
  } while (waited < 0 && errno == EINTR);
  require(waited == pid, "waitpid failed");
  require(WIFEXITED(status) && WEXITSTATUS(status) == 0, "Command failed: " + args[0]);
}
void write_fasta(std::ostream &f, const std::string &name, std::string_view s) {
  f << '>' << name << '\n';
  for (size_t p = 0; p < s.size(); p += 80)
    f.write(s.data() + p, std::min<size_t>(80, s.size() - p)) << '\n';
}
} // namespace hpc
