// SPDX-License-Identifier: MIT-0
// Text formats and subprocess boundaries; no consensus policy belongs here.
#include "io.hpp"
#include <cctype>
#include <cerrno>
#include <fstream>
#include <spawn.h>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <sys/wait.h>

extern char **environ;

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
void run(std::vector<std::string> args) {
  require(!args.empty(), "Cannot execute an empty command");
  // No shell interpolation: filenames and tier-derived strings remain arguments.
  std::vector<char *> ptr;
  for (auto &a : args)
    ptr.push_back(a.data());
  ptr.push_back(nullptr);
  // The caller can retain a chromosome-sized address space even after its
  // working buffers are freed. Avoid fork's page-table/commit overhead merely
  // to exec samtools. Spawn preserves PATH lookup, environment and literal argv.
  pid_t pid;
  const int error = posix_spawnp(&pid, ptr[0], nullptr, nullptr, ptr.data(), environ);
  if (error)
    throw std::system_error(error, std::generic_category(), "Cannot launch command: " + args[0]);
  int status = 0;
  pid_t waited;
  do {
    waited = waitpid(pid, &status, 0);
  } while (waited < 0 && errno == EINTR);
  if (waited < 0)
    throw std::system_error(errno, std::generic_category(), "Cannot wait for command: " + args[0]);
  if (WIFSIGNALED(status))
    throw std::runtime_error("Command killed by signal " + std::to_string(WTERMSIG(status)) + ": " +
                             args[0]);
  require(WIFEXITED(status) && WEXITSTATUS(status) == 0,
          "Command failed: " + args[0] +
              (WIFEXITED(status) ? " (exit status " + std::to_string(WEXITSTATUS(status)) + ")"
                                 : " (unexpected process status)"));
}
void write_fasta(std::ostream &f, const std::string &name, std::string_view s) {
  f << '>' << name << '\n';
  for (size_t p = 0; p < s.size(); p += 80)
    f.write(s.data() + p, std::min<size_t>(80, s.size() - p)) << '\n';
}
} // namespace hpc
