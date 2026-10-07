// SPDX-License-Identifier: MIT-0
// Pipeline: package prepass -> per-tig offset -> batched HPC votes -> expansion
// -> optional raw alignment/polish -> optional final-reference CRAM export.
#include "core.hpp"
#include "sqStore.H"
#include "tgTig.H"
#include "files.H"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <sys/wait.h>
#include <sys/resource.h>
#include <unistd.h>

namespace fs = std::filesystem;
using namespace hpc;
struct Options {
  // These defaults are part of the CLI contract; update help, README, and tests
  // together whenever changing an algorithmic threshold or resource policy.
  std::string input, templates, tig_map, output, cram;
  std::vector<Tier> tiers;
  unsigned threads = 1, min_support = 1, polish_support = 5;
  size_t batch_bases = 8000000;
  double max_error = 0.15, polish_fraction = 0.9;
  size_t max_indel = 500;
  bool polish = false, drop_unmatched = false, missing_error = false;
};
static void require(bool b, const std::string &s) {
  if (!b)
    throw std::runtime_error(s);
}
static unsigned number(const std::string &s) {
  size_t used = 0;
  auto n = std::stoul(s, &used);
  require(used == s.size() && n > 0 && n <= UINT32_MAX, "Expected positive integer: " + s);
  return n;
}
static double probability(const std::string &s) {
  size_t used = 0;
  double x = std::stod(s, &used);
  require(used == s.size() && std::isfinite(x) && x > 0 && x <= 1,
          "Expected probability in (0,1]: " + s);
  return x;
}
static Options parse(int argc, char **argv) {
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
             "--tier [LABEL=]GLOB   Repeated in priority order; first matching glob wins\n"
             "--unmatched error|drop (default error)\n--min-run-support N (default 1)\n"
             "--missing-run one|error (default one)\n--max-error F (default 0.15)\n"
             "--threads N (default 1; parallel alignments within each tig)\n"
             "--batch-bases N (default 8000000; raw bases per HPC alignment batch)\n"
             "--polish-non-hpc      Optional conservative correction\n"
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
    else if (a == "--output")
      o.output = value();
    else if (a == "--cram")
      o.cram = value();
    else if (a == "--polish-non-hpc")
      o.polish = true;
    else if (a == "--threads")
      o.threads = number(value());
    else if (a == "--batch-bases")
      o.batch_bases = number(value());
    else if (a == "--min-run-support")
      o.min_support = number(value());
    else if (a == "--min-polish-support")
      o.polish_support = number(value());
    else if (a == "--max-polish-indel")
      o.max_indel = number(value());
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
  return o;
}
static std::map<std::string, std::string> fasta(const std::string &path) {
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
static std::map<uint32_t, std::string> template_names(const std::string &path) {
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
static void run(std::vector<std::string> args) {
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
  require(waitpid(pid, &status, 0) == pid, "waitpid failed");
  require(WIFEXITED(status) && WEXITSTATUS(status) == 0,
          "Command failed: " + args[0] + " " + args[1]);
}
struct ReadOwner {
  // Placements share immutable sequence storage; construct reverse orientation
  // once per read only when needed. Work objects keep it alive after this map clears.
  struct Read {
    std::string name;
    std::shared_ptr<const std::string> sequence, reverse;
  };
  std::map<uint32, Read> reads;
};
// Use native binary readers but retain only decoded sequences, not packed blobs too.
static bool load_package(readBuffer &input, tgTig &tig, ReadOwner *owner) {
  // Native packages contain children+1 read records, including a redundant first
  // record. Consume every record, but retain each ID only once. A null owner is
  // the metadata prepass; the reusable sqRead then holds just one decoded record.
  if (!tig.loadFromBuffer(&input))
    return false;
  sqRead read;
  for (uint32_t i = 0; i <= tig.numberOfChildren(); ++i) {
    sqStore::sqStore_loadReadFromBuffer(&input, &read);
    if (owner && !owner->reads.count(read.sqRead_readID()))
      owner->reads.emplace(read.sqRead_readID(),
                           ReadOwner::Read{read.sqRead_name(),
                                           std::make_shared<const std::string>(
                                               read.sqRead_sequence(sqRead_raw | sqRead_normal)),
                                           {}});
  }
  return true;
}
static void memory_stage(const std::string &stage) {
  // ru_maxrss is a process high-water mark, not current usage. macOS reports
  // bytes; Linux reports KiB. This does not include simultaneous child RSS.
  rusage u{};
  getrusage(RUSAGE_SELF, &u);
  double bytes = u.ru_maxrss;
#ifndef __APPLE__
  bytes *= 1024;
#endif
  std::cerr << "[memory] " << stage << " peak_rss_GiB=" << bytes / (1024.0 * 1024 * 1024) << '\n';
}
struct Work {
  // hbegin/hend start in packaged HPC layout space, then become template-relative.
  // rle is batch-local scratch; trace changes from expanded ops to packed CIGAR
  // after voting. Alignment::ops has its own explicit packed_ops discriminator.
  Alignment a;
  size_t hbegin = 0, hend = 0;
  Rle rle;
  std::string trace, reason;
  bool boundary_left = false, boundary_right = false;
  bool unique = false;
  bool raw_alignment_attempted = false;
  unsigned placements = 0;
};
static auto placement_key(tgPosition *p) {
  // Same ID alone is not a duplicate: distinct loci/strands must survive for CRAM.
  return std::make_tuple(p->ident(), p->min(), p->max(), p->isReverse(), p->_askip, p->_bskip);
}
static unsigned placement_mapq(unsigned n) {
  // Equal alternatives imply P(wrong)=1-1/n. The unique cap of 60 is a user
  // convention, not a confidence estimate calibrated from alignment scores.
  if (n == 0)
    return 0;
  if (n == 1)
    return 60;
  return unsigned(std::lround(-10 * std::log10(1.0 - 1.0 / n)));
}
template <class F> static void parallel(size_t n, unsigned nt, F action) {
  // Each worker owns one Work slot; histogram accumulation is serial afterward.
  // Actions must catch per-read failures themselves: exceptions cannot cross threads.
  std::atomic<size_t> next{0};
  std::vector<std::thread> workers;
  for (size_t t = 0; t < std::min<size_t>(n, nt); ++t)
    workers.emplace_back([&]() {
      for (size_t i; (i = next.fetch_add(1)) < n;)
        action(i);
    });
  for (auto &t : workers)
    t.join();
}
static void write_fasta(std::ostream &f, const std::string &name, const std::string &s) {
  f << '>' << name << '\n';
  for (size_t p = 0; p < s.size(); p += 80)
    f << s.substr(p, 80) << '\n';
}
int main(int argc, char **argv) {
  fs::path staging;
  try {
    Options o = parse(argc, argv);
    auto templates = fasta(o.templates);
    auto mapped_names = template_names(o.tig_map);
    std::vector<std::string> suffixes = {".fasta", ".run-calls.tsv", ".read-decisions.tsv",
                                         ".summary.json", ".template-offsets.tsv"};
    if (o.polish) {
      suffixes.push_back(".non-hpc-edits.vcf");
      suffixes.push_back(".draft.fasta");
    }
    for (const auto &s : suffixes)
      require(!fs::exists(o.output + s), "Output exists: " + o.output + s);
    if (!o.cram.empty()) {
      require(!fs::exists(o.cram) && !fs::exists(o.cram + ".crai"), "CRAM output exists");
      require(!fs::exists(o.output + ".fasta.fai"), "FASTA index exists");
    }
    // Stage beside the prefix for rename-based publication. Failures retain the
    // directory for diagnosis; publication across multiple files is not atomic.
    std::string temp = o.output + ".work.XXXXXX";
    std::vector<char> tmp(temp.begin(), temp.end());
    tmp.push_back(0);
    require(mkdtemp(tmp.data()) != nullptr, "Cannot create staging directory beside output");
    staging = tmp.data();
    std::map<std::string, std::ofstream> outputs;
    for (auto &s : suffixes) {
      outputs[s].open(staging / ("assembly" + s));
      require(bool(outputs[s]), "Cannot create output");
      outputs[s].exceptions(std::ios::badbit | std::ios::failbit);
    }
    auto &fa = outputs[".fasta"];
    auto &calls = outputs[".run-calls.tsv"];
    auto &decisions = outputs[".read-decisions.tsv"];
    calls << "tig\trun\tdraft_begin\tdraft_end\tbase\tlength\ttier\tdepth\tpeak\n";
    decisions << "tig\tread_id\tname\ttier\tstatus\n";
    outputs[".template-offsets.tsv"]
        << "tig\toffset\tattempted\tplaced\tvalidation_support\ttolerance\tmethod\n";
    if (o.polish)
      outputs[".non-hpc-edits.vcf"]
          << "##fileformat=VCFv4.2\n##reference=uncompressed_draft\n##INFO=<ID=TIER,Number=1,Type="
             "String,Description=\"Winning "
             "tier\">\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Eligible spanning "
             "depth\">\n##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Alternate "
             "support\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    std::ofstream body;
    std::vector<std::pair<std::string, size_t>> refs;
    if (!o.cram.empty()) {
      body.open(staging / "body.sam");
      body.exceptions(std::ios::badbit | std::ios::failbit);
    }
    uint64_t tigs = 0, read_count = 0, used = 0, dropped = 0, no_evidence = 0, n_edits = 0,
             n_cram = 0, exact = 0, base_alignments = 0;
    std::set<uint32_t> tig_ids, emitted_ids;
    std::map<uint32_t, unsigned> occurrences, declared;
    sqRead_setDefaultVersion(sqRead_raw | sqRead_normal);
    memory_stage("starting placement prepass");
    // Count distinct placements across the whole package before assigning MAPQ.
    // Counts in the package may also describe alternatives in other packages.
    {
      readBuffer scan(o.input.c_str());
      while (!scan.eof()) {
        tgTig t;
        if (!load_package(scan, t, nullptr))
          break;
        std::set<decltype(placement_key(t.getChild(0)))> seen;
        for (uint32_t i = 0; i < t.numberOfChildren(); ++i) {
          auto *p = t.getChild(i);
          if (seen.insert(placement_key(p)).second)
            ++occurrences[p->ident()];
          declared[p->ident()] = std::max(declared[p->ident()], unsigned(p->_numPlacement));
        }
      }
    }
    memory_stage("placement prepass complete");
    readBuffer input(o.input.c_str());
    while (!input.eof()) {
      tgTig tig;
      ReadOwner owner;
      memory_stage("loading next tig");
      if (!load_package(input, tig, &owner))
        break;
      memory_stage("package loaded");
      ++tigs;
      require(tig_ids.insert(tig.tigID()).second, "Repeated tig ID");
      std::ostringstream label;
      label << "tig" << std::setfill('0') << std::setw(8) << tig.tigID();
      std::string name = label.str();
      require(name.find_first_of("\t\r\n ") == std::string::npos, "Invalid contig name");
      std::string ref(tig.bases() ? tig.bases() : "", tig._basesLen);
      auto it = templates.end();
      if (!o.tig_map.empty()) {
        auto mapped = mapped_names.find(tig.tigID());
        if (mapped != mapped_names.end())
          it = templates.find(mapped->second);
      } else {
        it = templates.find(name);
        if (it == templates.end())
          it = templates.find(std::to_string(tig.tigID()));
        if (it == templates.end()) {
          std::ostringstream legacy;
          legacy << "piece" << std::setfill('0') << std::setw(6) << tig.tigID();
          it = templates.find(legacy.str());
        }
      }
      if (it != templates.end()) {
        require(ref.empty() || ref == it->second, "Embedded/external template mismatch");
        ref = it->second;
      }
      std::vector<Work> work;
      std::set<decltype(placement_key(tig.getChild(0)))> local;
      for (uint32_t i = 0; i < tig.numberOfChildren(); ++i) {
        auto *p = tig.getChild(i);
        if (!local.insert(placement_key(p)).second)
          continue;
        auto r = owner.reads.find(p->ident());
        require(r != owner.reads.end(), "Missing packaged read");
        auto &rd = r->second;
        Work w;
        w.a.id = p->ident();
        w.a.name = rd.name;
        ++read_count;
        require(!w.a.name.empty() && w.a.name.find_first_of("\t\n\r ") == std::string::npos,
                "Read name cannot be represented as SAM QNAME");
        w.a.tier = tier_of(w.a.name, o.tiers);
        w.placements = std::max(occurrences[w.a.id], declared[w.a.id]);
        if (declared[w.a.id] == 0 && occurrences[w.a.id] == 1)
          w.placements = 0;
        w.unique = p->_numPlacement == 1 && w.placements == 1;
        // Ambiguous reads may be exported, but never seed the offset or vote in
        // this caller. Counts include alternatives even if later alignment fails.
        if (!w.unique && o.cram.empty())
          w.reason = "AMBIGUOUS_PLACEMENT";
        // Legacy skipConsensus encodes template-building priority, not exclusion.
        // Explicit tiers replace that policy; placement count determines uniqueness.
        else if (w.a.tier < 0) {
          require(o.drop_unmatched, "Unmatched read name: " + w.a.name);
          w.reason = "UNMATCHED_TIER";
        } else {
          require(p->min() >= 0 && p->max() > p->min() && size_t(p->max()) <= tig.length(),
                  "Invalid layout interval");
          w.hbegin = p->min();
          w.hend = p->max();
          w.a.reverse = p->isReverse();
          if (w.a.reverse && !rd.reverse)
            rd.reverse = std::make_shared<const std::string>(reverse_complement(*rd.sequence));
          w.a.shared_sequence = w.a.reverse ? rd.reverse : rd.sequence;
          require(p->_askip >= 0 && p->_bskip >= 0 &&
                      size_t(p->_askip) + size_t(p->_bskip) < w.a.bases().size(),
                  "Invalid raw clipping");
          w.a.clip_left = p->_askip;
          w.a.clip_right = p->_bskip;
        }
        work.push_back(std::move(w));
      }
      owner.reads.clear();
      if (work.size() == 1) {
        // Passthrough has no graph-coordinate inference or polishing. Preserve
        // oriented full read storage for SAM, emitting only its unclipped interval.
        auto &w = work.front();
        require(w.reason.empty() && w.unique,
                "Singleton tig " + name + " has no eligible unique read: " + w.reason);
        auto seq =
            w.a.bases().substr(w.a.clip_left, w.a.bases().size() - w.a.clip_left - w.a.clip_right);
        auto rle = encode(seq);
        ++used;
        write_fasta(fa, name, seq);
        refs.emplace_back(name, seq.size());
        if (o.polish)
          write_fasta(outputs[".draft.fasta"], name, seq);
        for (size_t i = 0; i < rle.bases.size(); ++i)
          calls << name << '\t' << i << '\t' << rle.offsets[i] << '\t' << rle.offsets[i + 1] << '\t'
                << rle.bases[i] << '\t' << rle.lengths[i] << '\t' << o.tiers[w.a.tier].label
                << "\t1\t1\n";
        outputs[".template-offsets.tsv"] << name << "\t.\t0\t0\t0\t0\tSINGLE_READ\n";
        decisions << name << '\t' << w.a.id << '\t' << w.a.name << '\t' << o.tiers[w.a.tier].label
                  << "\tACCEPTED\n";
        if (!o.cram.empty()) {
          std::string c;
          if (w.a.clip_left)
            c = std::to_string(w.a.clip_left) + "S";
          c += std::to_string(seq.size()) + "=";
          if (w.a.clip_right)
            c += std::to_string(w.a.clip_right) + "S";
          bool primary = emitted_ids.insert(w.a.id).second;
          body << w.a.name << '\t' << ((w.a.reverse ? 16 : 0) | (primary ? 0 : 256)) << '\t' << name
               << "\t1\t" << placement_mapq(w.placements) << '\t' << c << "\t*\t0\t0\t"
               << w.a.bases() << "\t*\tRG:Z:" << o.tiers[w.a.tier].label
               << "\tNM:i:0\tMD:Z:" << seq.size() << "\tAS:i:0\n";
          ++n_cram;
        }
        std::cerr << "tig " << name << ": single-read passthrough, " << seq.size() << " bases\n";
        continue;
      }
      require(!ref.empty(),
              "Tig " + name +
                  " has no HPC template; check --template-hpc and --tig-map associations");
      for (size_t i = 0; i < ref.size(); ++i)
        require(std::string("ACGTN").find(ref[i]) != std::string::npos &&
                    (!i || ref[i] != ref[i - 1]),
                "Template is not homopolymer compressed: " + name);
      // Prefer high-tier reads away from possible overhangs, then try boundary reads.
      std::vector<size_t> seeds;
      for (size_t i = 0; i < work.size(); ++i)
        if (work[i].reason.empty() && work[i].unique)
          seeds.push_back(i);
      size_t margin = tig.length() > ref.size() ? tig.length() - ref.size() : 0;
      auto rank = [&](size_t i) {
        auto &w = work[i];
        return std::make_tuple(w.a.tier, w.hbegin < margin || w.hend > ref.size(), w.hbegin,
                               w.a.id);
      };
      std::sort(seeds.begin(), seeds.end(), [&](size_t a, size_t b) { return rank(a) < rank(b); });
      std::vector<int64_t> estimates;
      int64_t offset = 0;
      size_t attempted = 0, validated = 0;
      bool inferred = false;
      // Disjoint discovery/validation sets. Never choose an offset from lengths alone.
      for (size_t i : seeds) {
        if (attempted == 32)
          break;
        ++attempted;
        auto &w = work[i];
        auto q = encode(w.a.bases().substr(w.a.clip_left,
                                           w.a.bases().size() - w.a.clip_left - w.a.clip_right))
                     .bases;
        auto pos = locate(q, ref, o.max_error);
        if (pos < 0)
          continue;
        estimates.push_back(int64_t(w.hbegin) - pos);
        if (seeds.size() < 6 || estimates.size() < 6)
          continue;
        size_t split = estimates.size() / 2;
        try {
          auto first = offset_cluster({estimates.begin(), estimates.begin() + split});
          auto second = offset_cluster({estimates.begin() + split, estimates.end()});
          if (std::abs(first - second) > 20)
            continue;
          offset = offset_cluster(estimates);
          validated = 0;
          for (size_t j = split; j < estimates.size(); ++j)
            if (std::abs(estimates[j] - offset) <= 20)
              ++validated;
          if (validated >= 3) {
            inferred = true;
            break;
          }
        } catch (const std::runtime_error &) {
        }
      }
      if (seeds.size() < 6 && !estimates.empty()) {
        // Align every available read. Sparse evidence has no held-out validation.
        try {
          offset = offset_cluster(estimates, 1);
          inferred = true;
        } catch (const std::runtime_error &) {
        }
      }
      require(inferred, "Cannot infer reliable template offset for tig " + name + " from " +
                            std::to_string(attempted) + " reads (" +
                            std::to_string(estimates.size()) + " unique placements)");
      outputs[".template-offsets.tsv"] << name << '\t' << offset << '\t' << attempted << '\t'
                                       << estimates.size() << '\t' << validated << "\t20\t"
                                       << (seeds.size() < 6 ? "SPARSE" : "VALIDATED") << '\n';
      std::cerr << "tig " << name << ": inferred template offset " << offset << " from "
                << estimates.size() << " reads\n";
      for (auto &w : work)
        if (w.reason.empty()) {
          // Recover overlap in template space; coordinate clipping here does not
          // yet clip the read. Semi-global HPC alignment determines its true bounds.
          int64_t b = int64_t(w.hbegin) - offset, e = int64_t(w.hend) - offset;
          if (e <= 0 || b >= int64_t(ref.size()))
            w.reason = "OUTSIDE_TEMPLATE";
          else {
            w.boundary_left = b < 0;
            w.boundary_right = e > int64_t(ref.size());
            w.hbegin = std::max<int64_t>(0, b);
            w.hend = std::min<int64_t>(ref.size(), e);
          }
        }
      // Histograms persist; per-read RLE and expanded traces exist only inside a bounded batch.
      struct RunBin {
        uint32_t length = 0, count = 0;
        int tier = -1;
      };
      // One inline allele handles the common case. Alternative lengths use
      // sparse overflow. Higher support thresholds must retain every tier until
      // its total depth is known, otherwise fallback would depend on read order.
      size_t histogram_tiers = o.min_support == 1 ? 1 : o.tiers.size();
      require(ref.size() <= SIZE_MAX / histogram_tiers, "Histogram allocation overflow");
      std::vector<RunBin> bins(ref.size() * histogram_tiers);
      std::map<size_t, Histogram> overflow;
      memory_stage("starting batched HPC alignment");
      for (size_t batch = 0; batch < work.size();) {
        // This bounds scratch, not total RSS. Always admit one read even if it
        // exceeds the base budget; also cap placement count for tiny/filtered reads.
        size_t end = batch, bases = 0;
        do {
          bases += work[end].a.bases().size();
          ++end;
        } while (end < work.size() && end - batch < 4096 && bases < o.batch_bases);
        parallel(end - batch, o.threads, [&](size_t j) {
          size_t i = batch + j;
          auto &w = work[i];
          if (!w.reason.empty())
            return;
          try {
            auto seq = w.a.bases().substr(w.a.clip_left,
                                          w.a.bases().size() - w.a.clip_left - w.a.clip_right);
            w.rle = encode(seq);
            auto target = ref.substr(w.hbegin, w.hend - w.hbegin);
            if (w.boundary_left || w.boundary_right) {
              auto a = align_boundary(w.rle.bases, target, w.boundary_left, w.boundary_right,
                                      o.max_error);
              size_t b = w.rle.offsets.at(a.query_begin), e = w.rle.offsets.at(a.query_end);
              // Convert aligned run bounds to raw bases before dropping the RLE
              // scratch. Existing raw clips are additive and retain SAM orientation.
              require(e > b, "Empty boundary overlap");
              w.a.clip_left += b;
              w.a.clip_right += seq.size() - e;
              w.rle = encode(seq.substr(b, e - b));
              w.trace = std::move(a.ops);
            } else
              w.trace = align(w.rle.bases, target, o.max_error);
          } catch (const std::exception &e) {
            w.reason = "HPC_ALIGNMENT_FAILED";
          }
        });
        for (size_t i = batch; i < end; ++i) {
          auto &w = work[i];
          if (w.reason.empty() && w.unique) {
            ++used;
            if (w.rle.bases == ref.substr(w.hbegin, w.hend - w.hbegin))
              ++exact;
            size_t q = 0, r = w.hbegin;
            for (char op : w.trace) {
              // Original raw clipping may cut a run. Do not vote with a truncated length.
              bool partial = (q == 0 && w.a.clip_left > 0 &&
                              w.a.bases()[w.a.clip_left - 1] == w.a.bases()[w.a.clip_left]) ||
                             (q + 1 == w.rle.bases.size() && w.a.clip_right > 0 &&
                              w.a.bases()[w.a.bases().size() - w.a.clip_right - 1] ==
                                  w.a.bases()[w.a.bases().size() - w.a.clip_right]);
              if (op == '=' && ref[r] != 'N' && !partial) {
                auto len = w.rle.lengths[q];
                {
                  size_t key = r * histogram_tiers + (o.min_support == 1 ? 0 : w.a.tier);
                  auto &bin = bins[key];
                  if (bin.tier < 0 || w.a.tier < bin.tier) {
                    bin = {len, 1, w.a.tier};
                    overflow.erase(key);
                  } else if (w.a.tier == bin.tier) {
                    if (len == bin.length)
                      ++bin.count;
                    else
                      ++overflow[key][len];
                  }
                }
              }
              if (op != 'D')
                ++q;
              if (op != 'I')
                ++r;
            }
          }
          w.rle = {};
          // clear() alone retains string capacity; swap releases no-longer-needed
          // storage. Never keep per-base traces for the entire chromosome.
          if ((o.polish || !o.cram.empty()) && w.reason.empty())
            w.trace = cigar(w.trace);
          else {
            std::string().swap(w.trace);
            std::string().swap(w.a.sequence);
            w.a.shared_sequence.reset();
          }
        }
        batch = end;
      }
      memory_stage("HPC evidence complete");
      std::vector<size_t> prefix(ref.size() + 1);
      std::string draft;
      // prefix is the exact HPC-run -> draft-base coordinate transform; the
      // sentinel supports reads ending at the final template run boundary.
      for (size_t i = 0; i < ref.size(); ++i) {
        Call c;
        for (size_t t = 0; t < histogram_tiers; ++t) {
          auto &bin = bins[i * histogram_tiers + t];
          if (bin.tier < 0)
            continue;
          auto it = overflow.find(i * histogram_tiers + t);
          Call candidate;
          if (it == overflow.end()) {
            if (bin.count >= o.min_support)
              candidate = {bin.length, bin.count, bin.count, bin.tier};
          } else {
            it->second[bin.length] = bin.count;
            candidate = call({it->second}, o.min_support);
            if (candidate.tier >= 0)
              candidate.tier = bin.tier;
            overflow.erase(it);
          }
          if (c.tier < 0 && candidate.tier >= 0)
            c = candidate;
        }
        if (c.tier < 0) {
          ++no_evidence;
          require(!o.missing_error, "No evidence at " + name + ":" + std::to_string(i));
        }
        prefix[i] = draft.size();
        draft.append(c.length, ref[i]);
        prefix[i + 1] = draft.size();
        calls << name << '\t' << i << '\t' << prefix[i] << '\t' << prefix[i + 1] << '\t' << ref[i]
              << '\t' << c.length << '\t' << (c.tier < 0 ? "NO_EVIDENCE" : o.tiers[c.tier].label)
              << '\t' << c.depth << '\t' << c.peak << '\n';
      }
      std::vector<RunBin>().swap(bins);
      overflow.clear();
      memory_stage("draft called");
      if (o.polish || !o.cram.empty()) {
        parallel(work.size(), o.threads, [&](size_t i) {
          auto &w = work[i];
          if (!w.reason.empty())
            return;
          try {
            w.a.begin = prefix[w.hbegin];
            w.a.end = prefix[w.hend];
            w.rle = encode(w.a.bases().substr(w.a.clip_left,
                                              w.a.bases().size() - w.a.clip_left - w.a.clip_right));
            // Re-encoding is linear and worker-local; it avoids retaining every
            // read's run arrays between stages without repeating HPC alignment.
            auto trace = expand_cigar(w.trace);
            // Build a valid raw trace directly for runs, including unequal-length mismatches.
            size_t q = 0, r = w.hbegin;
            std::string ops;
            for (char op : trace) {
              size_t nq = op == 'D' ? 0 : w.rle.lengths[q];
              size_t nr = op == 'I' ? 0 : prefix[r + 1] - prefix[r];
              size_t both = std::min(nq, nr);
              ops.append(both, op == '=' ? '=' : 'X');
              ops.append(nq - both, 'I');
              ops.append(nr - both, 'D');
              if (op != 'D')
                ++q;
              if (op != 'I')
                ++r;
            }
            w.a.ops = std::move(ops);
            // Prototype: realign the layout interval when its HPC trace has errors.
            if (w.trace.find_first_of("XID") != std::string::npos) {
              w.raw_alignment_attempted = true;
              w.a.ops = align(w.a.bases().substr(w.a.clip_left, w.a.bases().size() - w.a.clip_left -
                                                                    w.a.clip_right),
                              draft.substr(w.a.begin, w.a.end - w.a.begin), o.max_error);
            }
            size_t nm;
            md(w.a, draft, nm);
            w.a.ops = cigar(w.a.ops);
            w.a.packed_ops = true;
          } catch (const std::exception &e) {
            w.reason = "FINAL_ALIGNMENT_FAILED";
            std::string().swap(w.a.ops);
            w.a.shared_sequence.reset();
          }
          w.rle = {};
        });
        for (auto &w : work)
          if (w.raw_alignment_attempted)
            ++base_alignments;
      }
      for (auto &w : work)
        std::string().swap(w.trace);
      std::vector<size_t>().swap(prefix);
      memory_stage("raw alignments complete");
      std::vector<Accepted> edits;
      if (o.polish) {
        write_fasta(outputs[".draft.fasta"], name, draft);
        std::vector<const Alignment *> alns;
        for (auto &w : work)
          if (w.reason.empty() && w.unique)
            alns.push_back(&w.a);
        edits =
            polish(draft, alns, o.tiers.size(), o.polish_support, o.polish_fraction, o.max_indel);
        n_edits += edits.size();
        memory_stage("polishing complete");
        for (const auto &e : edits) {
          // VCF anchor indels with the preceding reference base, in draft coordinates.
          size_t p = e.edit.begin;
          std::string rb = draft.substr(p, e.edit.end - p), ab = e.edit.alt;
          if (rb.empty() || ab.empty()) {
            --p;
            rb = draft[p] + rb;
            ab = draft[p] + ab;
          }
          outputs[".non-hpc-edits.vcf"] << name << '\t' << p + 1 << "\t.\t" << rb << '\t' << ab
                                        << "\t.\tPASS\tTIER=" << o.tiers[e.tier].label
                                        << ";DP=" << e.depth << ";AC=" << e.support << '\n';
        }
      }
      std::string final = hpc::apply(draft, edits);
      write_fasta(fa, name, final);
      refs.emplace_back(name, final.size());
      if (!o.cram.empty()) {
        for (auto &w : work)
          if (w.reason.empty()) {
            bool touches = false;
            for (auto &e : edits)
              if (e.edit.begin <= w.a.end && e.edit.end >= w.a.begin)
                touches = true;
            w.a.begin = lift(w.a.begin, edits);
            w.a.end = lift(w.a.end, edits);
            // Coordinate shifts alone suffice only outside edits. A touched read
            // needs a new trace and fresh NM/MD against the final reference.
            try {
              if (touches) {
                ++base_alignments;
                w.a.ops =
                    align(w.a.bases().substr(w.a.clip_left,
                                             w.a.bases().size() - w.a.clip_left - w.a.clip_right),
                          final.substr(w.a.begin, w.a.end - w.a.begin), o.max_error);
                w.a.packed_ops = false;
              }
              size_t nm;
              auto mdtag = md(w.a, final, nm);
              std::string c;
              if (w.a.clip_left)
                c = std::to_string(w.a.clip_left) + "S";
              c += w.a.packed_ops ? w.a.ops : cigar(w.a.ops);
              if (w.a.clip_right)
                c += std::to_string(w.a.clip_right) + "S";
              bool primary = emitted_ids.insert(w.a.id).second;
              // Primary selection follows successful emission order, not the
              // coordinate order introduced later by samtools sort.
              body << w.a.name << '\t' << ((w.a.reverse ? 16 : 0) | (primary ? 0 : 256)) << '\t'
                   << name << '\t' << w.a.begin + 1 << '\t' << placement_mapq(w.placements) << '\t'
                   << c << "\t*\t0\t0\t" << w.a.bases() << "\t*\tRG:Z:" << o.tiers[w.a.tier].label
                   << "\tNM:i:" << nm << "\tMD:Z:" << mdtag << "\tAS:i:" << -int64_t(nm) << '\n';
              ++n_cram;
            } catch (const std::exception &e) {
              w.reason = "FINAL_ALIGNMENT_FAILED";
            }
          }
      }
      for (auto &w : work) {
        if (!w.reason.empty())
          ++dropped;
        decisions << name << '\t' << w.a.id << '\t' << w.a.name << '\t'
                  << (w.a.tier < 0 ? "." : o.tiers[w.a.tier].label) << '\t'
                  << (w.reason.empty() ? "ACCEPTED" : w.reason) << '\n';
      }
      std::cerr << "tig " << name << ": " << ref.size() << " runs -> " << final.size() << " bases, "
                << edits.size() << " edits\n";
    }
    require(tigs > 0, "Empty package");
    outputs[".summary.json"] << "{\n  \"tigs\": " << tigs << ",\n  \"reads\": " << read_count
                             << ",\n  \"consensus_reads\": " << used
                             << ",\n  \"dropped\": " << dropped
                             << ",\n  \"exact_hpc_reads\": " << exact
                             << ",\n  \"no_evidence_runs\": " << no_evidence
                             << ",\n  \"non_hpc_edits\": " << n_edits
                             << ",\n  \"base_alignments\": " << base_alignments
                             << ",\n  \"cram_records\": " << n_cram << "\n}\n";
    for (auto &[s, f] : outputs)
      f.close();
    if (!o.cram.empty()) {
      // Contig lengths are now final: assemble the SAM header, then sort/index
      // CRAM with embedded references so the staging FASTA is not a dependency.
      body.close();
      std::ofstream sam(staging / "alignments.sam");
      sam.exceptions(std::ios::badbit | std::ios::failbit);
      sam << "@HD\tVN:1.6\tSO:unsorted\n";
      for (auto &[n, l] : refs)
        sam << "@SQ\tSN:" << n << "\tLN:" << l << '\n';
      for (auto &t : o.tiers)
        sam << "@RG\tID:" << t.label << "\tSM:assembly\n";
      sam << "@PG\tID:hpc-consensus\tPN:hpc-consensus\tVN:0.1.0\n";
      std::ifstream in(staging / "body.sam");
      if (in.peek() != EOF)
        sam << in.rdbuf();
      sam.close();
      run({"samtools", "faidx", (staging / "assembly.fasta").string()});
      run({"samtools", "sort", "-@", std::to_string(o.threads), "-O", "CRAM", "--reference",
           (staging / "assembly.fasta").string(), "--output-fmt-option", "embed_ref=1", "-o",
           (staging / "assembly.cram").string(), (staging / "alignments.sam").string()});
      run({"samtools", "index", (staging / "assembly.cram").string()});
      run({"samtools", "quickcheck", (staging / "assembly.cram").string()});
    }
    for (auto &s : suffixes)
      fs::rename(staging / ("assembly" + s), o.output + s);
    if (!o.cram.empty()) {
      fs::rename(staging / "assembly.fasta.fai", o.output + ".fasta.fai");
      fs::rename(staging / "assembly.cram", o.cram);
      fs::rename(staging / "assembly.cram.crai", o.cram + ".crai");
    }
    fs::remove_all(staging);
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "hpc-consensus: " << e.what() << '\n';
    if (!staging.empty())
      std::cerr << "Incomplete outputs retained in " << staging << '\n';
    return 1;
  }
}
