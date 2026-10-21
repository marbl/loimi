// SPDX-License-Identifier: MIT-0
// Pipeline: package prepass -> per-tig offset -> batched HPC votes -> expansion
// -> optional raw alignment/polish -> optional final-reference CRAM export.
#include "core.hpp"
#include "offsets.hpp"
#include "options.hpp"
#include "io.hpp"
#include "parallel.hpp"
#include "sqStore.H"
#include "tgTig.H"
#include "files.H"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
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
static void require(bool b, const std::string &s) {
  if (!b)
    throw std::runtime_error(s);
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
using ReadFilter = std::function<bool(tgTig &, sqRead &, ReadOwner &)>;
static bool load_package(readBuffer &input, tgTig &tig, ReadOwner *owner,
                         const ReadFilter &retain = {}) {
  // Native packages contain children+1 read records, including a redundant first
  // record. Consume every record, but retain each ID only once. A null owner is
  // the metadata prepass; the reusable sqRead then holds just one decoded record.
  if (!tig.loadFromBuffer(&input))
    return false;
  sqRead read;
  for (uint32_t i = 0; i <= tig.numberOfChildren(); ++i) {
    sqStore::sqStore_loadReadFromBuffer(&input, &read);
    if (owner && !owner->reads.count(read.sqRead_readID()) &&
        (!retain || retain(tig, read, *owner)))
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
  size_t original_clip_left = 0, original_clip_right = 0;
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
// Shared lookup for region planning and the full caller; the name map remains
// authoritative even when the chosen region is only a prefix of a template.
static std::string template_for(tgTig &tig, const Options &o,
                                const std::map<std::string, std::string> &templates,
                                const std::map<uint32_t, std::string> &names) {
  std::string ref(tig.bases() ? tig.bases() : "", tig._basesLen);
  auto it = templates.end();
  if (!o.tig_map.empty()) {
    auto n = names.find(tig.tigID());
    if (n != names.end())
      it = templates.find(n->second);
  } else {
    for (const auto &prefix : {std::string("tig"), std::string(""), std::string("piece")}) {
      std::ostringstream key;
      key << prefix << std::setfill('0')
          << std::setw(prefix == "tig"     ? 8
                       : prefix == "piece" ? 6
                                           : 0)
          << tig.tigID();
      it = templates.find(key.str());
      if (it != templates.end())
        break;
    }
  }
  if (it != templates.end()) {
    require(ref.empty() || ref == it->second, "Embedded/external template mismatch");
    ref = it->second;
  }
  return ref;
}
struct Offset {
  int64_t value = 0;
  size_t attempted = 0, placed = 0, validated = 0;
  bool sparse = false;
  bool internal = false;
  bool layout_rescue = false;
};
static Offset infer_offset(std::vector<Work> &work, size_t layout_length, const std::string &ref,
                           double max_error, size_t eligible_count = SIZE_MAX,
                           bool rescue_only = false) {
  std::vector<size_t> seeds;
  for (size_t i = 0; i < work.size(); ++i)
    if (work[i].reason.empty() && work[i].unique)
      seeds.push_back(i);
  size_t margin = layout_length > ref.size() ? layout_length - ref.size() : 0;
  auto rank = [&](size_t i) {
    auto &w = work[i];
    return std::make_tuple(w.a.tier, w.hbegin < margin || w.hend > ref.size(), w.hbegin, w.a.id);
  };
  std::sort(seeds.begin(), seeds.end(), [&](size_t a, size_t b) { return rank(a) < rank(b); });
  Offset out;
  out.sparse = (eligible_count == SIZE_MAX ? seeds.size() : eligible_count) < 6;
  std::vector<int64_t> estimates;
  if (!rescue_only) {
    for (size_t i : seeds) {
      if (out.attempted == 32)
        break;
      ++out.attempted;
      auto &w = work[i];
      auto q = encode(w.a.bases().substr(w.a.clip_left,
                                         w.a.bases().size() - w.a.clip_left - w.a.clip_right))
                   .bases;
      auto pos = locate(q, ref, max_error);
      if (pos < 0)
        continue;
      estimates.push_back(int64_t(w.hbegin) - pos);
      out.placed = estimates.size();
      if (out.sparse || estimates.size() < 6)
        continue;
      size_t split = estimates.size() / 2;
      try {
        auto first = offset_cluster({estimates.begin(), estimates.begin() + split});
        auto second = offset_cluster({estimates.begin() + split, estimates.end()});
        if (std::abs(first - second) > 20)
          continue;
        out.value = offset_cluster(estimates);
        out.validated = 0;
        for (size_t j = split; j < estimates.size(); ++j)
          out.validated += std::abs(estimates[j] - out.value) <= 20;
        if (out.validated >= 3)
          return out;
      } catch (const std::runtime_error &) {
      }
    }
    if (out.sparse && !estimates.empty()) {
      try {
        out.value = offset_cluster(estimates, 1);
        return out;
      } catch (const std::runtime_error &) {
        // Sparse layouts can also contain reads overhanging both template ends.
      }
    }
    // Preserve the established fast path exactly. Only unresolved offsets use
    // internal anchors; retain bounded scratch and count each read once.
    std::vector<std::vector<int64_t>> anchors;
    size_t limit = std::min<size_t>(32, seeds.size());
    for (size_t j = 0; j < limit; ++j) {
      auto &w = work[seeds[j]];
      auto q = encode(w.a.bases().substr(w.a.clip_left,
                                         w.a.bases().size() - w.a.clip_left - w.a.clip_right))
                   .bases;
      anchors.push_back(internal_offsets(q, ref, int64_t(w.hbegin)));
    }
    try {
      // Small informative samples must tolerate the outliers allowed by the
      // 75% agreement rule. Count independent reads, never individual chunks.
      size_t informative = std::count_if(anchors.begin(), anchors.end(),
                                         [](const auto &read) { return !read.empty(); });
      size_t required =
          out.sparse ? 1 : std::max<size_t>(3, std::min<size_t>(6, (3 * informative + 3) / 4));
      auto candidate = internal_offset_cluster(anchors, required);
      size_t validated = 0;
      for (size_t j = 0; j < limit; ++j) {
        if (std::none_of(anchors[j].begin(), anchors[j].end(),
                         [&](int64_t x) { return std::abs(x - candidate) <= 20; }))
          continue;
        auto &w = work[seeds[j]];
        auto q = encode(w.a.bases().substr(w.a.clip_left,
                                           w.a.bases().size() - w.a.clip_left - w.a.clip_right))
                     .bases;
        if (validate_offset_overlap(q, ref, w.hbegin, w.hend, candidate, max_error))
          ++validated;
      }
      if (validated >= (out.sparse ? 1u : 3u)) {
        out.value = candidate;
        out.internal = true;
        out.validated = validated;
        out.placed = std::count_if(anchors.begin(), anchors.end(), [&](const auto &read) {
          return std::any_of(read.begin(), read.end(),
                             [&](int64_t x) { return std::abs(x - candidate) <= 20; });
        });
        std::cerr << "offset rescue: internal anchors, offset " << candidate
                  << ", validated overlaps " << validated << '\n';
        return out;
      }
    } catch (const std::runtime_error &) {
      // Do not guess an offset when internal evidence is weak or contradictory.
    }
  }
  // Final rescue: known placements from each technology, including non-unique
  // ones, can anchor the coordinate transform without gaining consensus votes.
  std::map<int, std::vector<size_t>> tiers;
  for (size_t i = 0; i < work.size(); ++i) {
    const auto &w = work[i];
    if ((w.reason.empty() || w.reason == "AMBIGUOUS_PLACEMENT") && w.a.tier >= 0 &&
        w.placements > 0 && !w.a.bases().empty()) {
      tiers[w.a.tier].push_back(i);
    }
  }
  std::vector<OffsetEvidence> evidence;
  std::vector<int> evidence_tiers;
  std::vector<size_t> selected;
  for (auto &[tier, indices] : tiers) {
    std::set<std::pair<uint64_t, uint32_t>> sample;
    for (auto i : indices) {
      auto id = work[i].a.id;
      sample.emplace(offset_sample_key(id), id);
      if (sample.size() > 64)
        sample.erase(std::prev(sample.end()));
    }
    for (auto i : indices) {
      auto &w = work[i];
      if (!sample.count({offset_sample_key(w.a.id), w.a.id}))
        continue;
      auto q = encode(w.a.bases().substr(w.a.clip_left,
                                         w.a.bases().size() - w.a.clip_left - w.a.clip_right))
                   .bases;
      evidence.push_back({w.a.id, w.unique, discover_offsets(q, ref, w.hbegin, max_error)});
      evidence_tiers.push_back(tier);
      selected.push_back(i);
    }
  }
  // Score the SAME candidates against every sampled placement, not just reads
  // that proposed that candidate. This allows sparse cross-tier corroboration.
  // Only one read's compressed sequence exists at a time; the persistent matrix
  // is bounded by sampled placements x 128 candidate loci, not chromosome size.
  auto candidates = offset_candidates(evidence, evidence_tiers);
  std::vector<OffsetRow> rows;
  for (auto i : selected) {
    auto &w = work[i];
    auto q = encode(w.a.bases().substr(w.a.clip_left,
                                       w.a.bases().size() - w.a.clip_left - w.a.clip_right))
                 .bases;
    OffsetRow row{offset_molecule(w.a.name, w.a.id), w.a.tier, w.unique, {}};
    for (auto candidate : candidates)
      row.scores.push_back(score_offset_overlap(q, ref, w.hbegin, w.hend, candidate, max_error));
    rows.push_back(std::move(row));
  }
  try {
    auto choice = choose_joint_offset(candidates, rows);
    out.value = choice.value;
    out.layout_rescue = true;
    out.attempted = selected.size();
    out.placed = out.validated = choice.support;
    std::cerr << "offset joint rescue: offset=" << choice.value << " support=" << choice.support
              << " informative=" << choice.informative << " runner_up=" << choice.runner_up
              << " candidates=" << candidates.size() << " tier=" << choice.tier << '\n';
    return out;
  } catch (const std::runtime_error &e) {
    std::cerr << "offset joint rescue: " << e.what() << '\n';
  }
  throw std::runtime_error("Cannot infer reliable template offset from " +
                           std::to_string(out.attempted) + " reads (" +
                           std::to_string(estimates.size()) + " unique placements)");
}
// Region mode normally retains only its established unique fast-path anchors.
// On failure, reread the package with a bounded per-tier sample of known reads.
// This avoids loading a whole chromosome merely to test a small region, and
// uses the same ID sampling and placement multiplicities as whole-tig rescue.
static std::vector<Work> region_rescue_reads(const Options &o,
                                             const std::map<uint32_t, unsigned> &occurrences,
                                             const std::map<uint32_t, unsigned> &declared) {
  readBuffer input(o.input.c_str());
  while (!input.eof()) {
    tgTig t;
    ReadOwner owner;
    std::map<int, std::set<std::pair<uint64_t, uint32_t>>> samples;
    std::set<uint32_t> seen;
    auto retain = [&](tgTig &tig, sqRead &read, ReadOwner &storage) {
      if (tig.tigID() != o.region_tig)
        return false;
      auto id = read.sqRead_readID();
      if (!seen.insert(id).second || declared.at(id) == 0)
        return false;
      int tier = tier_of(read.sqRead_name(), o.tiers);
      if (tier < 0)
        return false;
      auto &sample = samples[tier];
      sample.emplace(offset_sample_key(id), id);
      if (sample.size() > 64) {
        auto last = std::prev(sample.end());
        uint32_t removed = last->second;
        sample.erase(last);
        storage.reads.erase(removed);
        if (removed == id)
          return false;
      }
      return true;
    };
    if (!load_package(input, t, &owner, retain))
      break;
    if (t.tigID() != o.region_tig)
      continue;
    std::vector<Work> result;
    for (uint32_t i = 0; i < t.numberOfChildren(); ++i) {
      auto p = t.getChild(i);
      auto it = owner.reads.find(p->ident());
      if (it == owner.reads.end() || p->_numPlacement == 0)
        continue;
      auto &rd = it->second;
      Work w;
      w.a.id = p->ident();
      w.a.name = rd.name;
      w.a.tier = tier_of(rd.name, o.tiers);
      w.placements = std::max(occurrences.at(w.a.id), declared.at(w.a.id));
      w.unique = p->_numPlacement == 1 && w.placements == 1;
      w.a.reverse = p->isReverse();
      if (w.a.reverse && !rd.reverse)
        rd.reverse = std::make_shared<const std::string>(reverse_complement(*rd.sequence));
      w.a.shared_sequence = w.a.reverse ? rd.reverse : rd.sequence;
      require(p->min() >= 0 && p->max() > p->min() && size_t(p->max()) <= t.length(),
              "Invalid layout interval in rescue anchor");
      require(p->_askip >= 0 && p->_bskip >= 0 &&
                  size_t(p->_askip) + p->_bskip < w.a.bases().size(),
              "Invalid clipping in rescue anchor");
      w.hbegin = p->min();
      w.hend = p->max();
      w.a.clip_left = p->_askip;
      w.a.clip_right = p->_bskip;
      result.push_back(std::move(w));
    }
    return result;
  }
  throw std::runtime_error("Missing region tig during offset rescue");
}
int main(int argc, char **argv) {
  fs::path staging;
  try {
    Options o = parse_options(argc, argv);
    auto templates = fasta(o.templates);
    auto mapped_names = template_names(o.tig_map);
    std::map<int, RunModel> models;
    for (size_t i = 0; i < o.tiers.size(); ++i) {
      auto path = o.model_paths.find(o.tiers[i].label);
      if (path != o.model_paths.end()) {
        std::ifstream input(path->second);
        require(bool(input), "Cannot open tier model: " + path->second);
        models.emplace(i, RunModel::read(input));
      }
    }
    // Calibration observes the exact accepted HPC votes, before priority drops
    // lower tiers. Exclusions affect this diagnostic only, never consensus.
    std::map<std::string, std::set<size_t>> evidence_sites;
    std::set<uint32_t> evidence_excluded;
    if (!o.evidence_sites.empty()) {
      std::ifstream in(o.evidence_sites);
      require(bool(in), "Cannot open evidence sites");
      std::string line;
      while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#')
          continue;
        std::istringstream row(line);
        std::string name, pos, extra;
        require(bool(row >> name >> pos) && !(row >> extra), "Expected TIG RUN in evidence sites");
        require(evidence_sites[name].insert(parse_coordinate(pos)).second,
                "Duplicate evidence site");
      }
      require(in.eof(), "Cannot read evidence sites");
    }
    if (!o.evidence_exclude_reads.empty()) {
      std::ifstream in(o.evidence_exclude_reads);
      require(bool(in), "Cannot open evidence read exclusions");
      std::string id;
      while (in >> id)
        evidence_excluded.insert(parse_positive_integer(id));
      require(in.eof(), "Cannot read evidence exclusions");
    }
    std::vector<std::string> suffixes = {".fasta", ".run-calls.tsv", ".read-decisions.tsv",
                                         ".summary.json", ".template-offsets.tsv"};
    const bool extend = (o.polish || o.extend) && !o.no_extend;
    if (extend)
      suffixes.push_back(".end-extensions.tsv");
    if (o.region)
      suffixes.push_back(".region.tsv");
    if (o.consistent_ambiguous)
      suffixes.push_back(".borrowing.tsv");
    if (!o.evidence_sites.empty())
      suffixes.push_back(".run-evidence.tsv");
    if (o.evidence_read_votes)
      suffixes.push_back(".run-read-votes.tsv");
    if (!models.empty())
      suffixes.push_back(".models.tsv");
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
    if (extend)
      outputs[".end-extensions.tsv"] << "tig\toriginal_draft_length\tleft_added\tright_added\tleft_"
                                        "candidates\tright_candidates\n";
    if (!models.empty()) {
      outputs[".models.tsv"] << "tier\tmodel_path\n";
      for (const auto &[label, path] : o.model_paths)
        outputs[".models.tsv"] << label << '\t' << path << '\n';
    }
    calls << "tig\trun\tdraft_begin\tdraft_end\tbase\tlength\ttier\tdepth\tpeak\tplacement_class\n";
    if (!o.evidence_sites.empty())
      outputs[".run-evidence.tsv"] << "tig\trun\tbase\ttier\tplacement_class\tlength\tcount\n";
    decisions << "tig\tread_id\tname\ttier\tstatus\n";
    if (o.evidence_read_votes)
      outputs[".run-read-votes.tsv"]
          << "tig\trun\tread_id\tname\ttier\tplacements\treverse\thpc_begin\thpc_end"
             "\top\tlength\tstatus\tquery_raw_begin\tquery_raw_end\n";
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
             n_cram = 0, exact = 0, base_alignments = 0, ambiguous_run_calls = 0,
             model_changed_runs = 0, filtered_hpc_votes = 0, locally_filtered_reads = 0,
             left_added = 0, right_added = 0;
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
    Offset region_offset;
    size_t region_hpc_end = 0;
    std::set<uint32_t> region_reads;
    if (o.region) {
      // Region mode needs one additional sequential pass: keep only the first
      // 32 ranked offset anchors, infer against the complete HPC template, then
      // retain just prefix-overlapping reads in the normal processing pass.
      // Multiplicity was counted across the COMPLETE package above.
      readBuffer scan(o.input.c_str());
      bool found = false;
      while (!scan.eof()) {
        tgTig t;
        ReadOwner anchors;
        std::map<uint32_t, tgPosition *> positions;
        using Rank = std::tuple<int, bool, size_t, uint32_t>;
        std::set<Rank> ranks;
        std::set<uint32_t> considered;
        size_t eligible = 0;
        std::string template_seq;
        bool prepared = false;
        auto retain = [&](tgTig &tig, sqRead &read, ReadOwner &owner) {
          if (tig.tigID() != o.region_tig)
            return false;
          if (!prepared) {
            prepared = true;
            template_seq = template_for(tig, o, templates, mapped_names);
            for (uint32_t i = 0; i < tig.numberOfChildren(); ++i) {
              auto p = tig.getChild(i);
              if (p->_numPlacement == 1 && occurrences[p->ident()] == 1 &&
                  declared[p->ident()] == 1)
                positions.emplace(p->ident(), p);
            }
          }
          uint32_t id = read.sqRead_readID();
          auto it = positions.find(id);
          if (it == positions.end() || !considered.insert(id).second)
            return false;
          int tier = tier_of(read.sqRead_name(), o.tiers);
          if (tier < 0)
            return false;
          ++eligible;
          auto p = it->second;
          require(p->min() >= 0 && p->max() > p->min() && size_t(p->max()) <= tig.length(),
                  "Invalid layout interval in region anchor");
          size_t margin =
              tig.length() > template_seq.size() ? tig.length() - template_seq.size() : 0;
          Rank rank{tier, size_t(p->min()) < margin || size_t(p->max()) > template_seq.size(),
                    size_t(p->min()), id};
          ranks.insert(rank);
          if (ranks.size() > 32) {
            auto worst = std::prev(ranks.end());
            uint32_t removed = std::get<3>(*worst);
            ranks.erase(worst);
            owner.reads.erase(removed);
            if (removed == id)
              return false;
          }
          return true;
        };
        if (!load_package(scan, t, &anchors, retain))
          break;
        if (t.tigID() != o.region_tig)
          continue;
        require(!found, "Repeated region tig ID");
        found = true;
        // The existing singleton path crops the oriented read directly.
        if (t.numberOfChildren() == 1) {
          region_reads.insert(t.getChild(0)->ident());
          continue;
        }
        require(!template_seq.empty(), "Region tig has no HPC template");
        std::vector<Work> seeds;
        for (auto &[id, rd] : anchors.reads) {
          auto p = positions.at(id);
          Work w;
          w.unique = true;
          w.a.id = id;
          w.a.tier = tier_of(rd.name, o.tiers);
          w.a.sequence = p->isReverse() ? reverse_complement(*rd.sequence) : *rd.sequence;
          require(p->_askip >= 0 && p->_bskip >= 0 &&
                      size_t(p->_askip) + size_t(p->_bskip) < w.a.sequence.size(),
                  "Invalid raw clipping in region anchor");
          w.a.clip_left = p->_askip;
          w.a.clip_right = p->_bskip;
          w.hbegin = p->min();
          w.hend = p->max();
          seeds.push_back(std::move(w));
        }
        std::cerr << "tig " << t.tigID() << ": planning region offset\n";
        try {
          region_offset = infer_offset(seeds, t.length(), template_seq, o.max_error, eligible);
        } catch (const std::runtime_error &) {
          seeds.clear();
          auto rescue = region_rescue_reads(o, occurrences, declared);
          region_offset = infer_offset(rescue, t.length(), template_seq, o.max_error, eligible, true);
        }
        // At least one raw base per run guarantees this prefix reaches END.
        // Keep full overlapping placements as context so artificial boundaries
        // do not change their run correspondence. Do not recursively extend.
        size_t limit = std::min(template_seq.size(), o.region_end + o.region_flank);
        region_hpc_end = limit;
        for (uint32_t i = 0; i < t.numberOfChildren(); ++i) {
          auto p = t.getChild(i);
          int64_t b = int64_t(p->min()) - region_offset.value;
          int64_t e = int64_t(p->max()) - region_offset.value;
          if (b < int64_t(limit) && e > 0) {
            region_reads.insert(p->ident());
            region_hpc_end =
                std::max(region_hpc_end, size_t(std::min<int64_t>(e, template_seq.size())));
          }
        }
        std::cerr << "Region: retaining " << region_reads.size() << " read IDs for "
                  << region_hpc_end << " HPC prefix runs\n";
      }
      require(found, "Requested region tig is absent from package");
      outputs[".region.tsv"]
          << "tig\tdraft_begin\tdraft_end\tfinal_begin\tfinal_end\tprocessed_hpc_end\n";
    }
    readBuffer input(o.input.c_str());
    while (!input.eof()) {
      tgTig tig;
      ReadOwner owner;
      memory_stage("loading next tig");
      auto retain_region = [&](tgTig &t, sqRead &r, ReadOwner &) {
        return !o.region || (t.tigID() == o.region_tig && region_reads.count(r.sqRead_readID()));
      };
      if (!load_package(input, tig, &owner, retain_region))
        break;
      if (o.region && tig.tigID() != o.region_tig)
        continue;
      memory_stage("package loaded");
      ++tigs;
      require(tig_ids.insert(tig.tigID()).second, "Repeated tig ID");
      std::ostringstream label;
      label << "tig" << std::setfill('0') << std::setw(8) << tig.tigID();
      std::string name = label.str();
      require(name.find_first_of("\t\r\n ") == std::string::npos, "Invalid contig name");
      std::string ref = template_for(tig, o, templates, mapped_names);
      const size_t complete_template_size = ref.size();
      if (o.region && tig.numberOfChildren() != 1)
        ref.resize(region_hpc_end);
      std::vector<Work> work;
      std::set<decltype(placement_key(tig.getChild(0)))> local;
      for (uint32_t i = 0; i < tig.numberOfChildren(); ++i) {
        auto *p = tig.getChild(i);
        if (!local.insert(placement_key(p)).second)
          continue;
        if (o.region && tig.numberOfChildren() != 1) {
          int64_t b = int64_t(p->min()) - region_offset.value;
          int64_t e = int64_t(p->max()) - region_offset.value;
          if (b >= int64_t(std::min(ref.size(), o.region_end + o.region_flank)) || e <= 0)
            continue;
        }
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
        // Known non-unique reads can supply fallback evidence. Unknown placement
        // counts remain export-only. Known placements can still anchor offsets
        // when --unique-only excludes them from consensus voting.
        if (!w.unique && o.cram.empty() && !(o.ambiguous_fallback && w.placements > 1))
          w.reason = "AMBIGUOUS_PLACEMENT";
        // Legacy skipConsensus encodes template-building priority, not exclusion.
        // Explicit tiers replace that policy; placement count determines uniqueness.
        else if (w.a.tier < 0) {
          require(o.drop_unmatched, "Unmatched read name: " + w.a.name);
          w.reason = "UNMATCHED_TIER";
        }
        if (w.a.tier >= 0 && (w.reason.empty() || w.placements > 0)) {
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
          w.original_clip_left = w.a.clip_left;
          w.original_clip_right = w.a.clip_right;
        }
        work.push_back(std::move(w));
      }
      owner.reads.clear();
      if (work.size() == 1 && (!o.region || tig.numberOfChildren() == 1)) {
        // Passthrough has no graph-coordinate inference or polishing. Preserve
        // oriented full read storage for SAM, emitting only its unclipped interval.
        auto &w = work.front();
        require(w.reason.empty() && (w.unique || (o.ambiguous_fallback && w.placements > 1)),
                "Singleton tig " + name + " has no eligible read: " + w.reason);
        auto seq =
            w.a.bases().substr(w.a.clip_left, w.a.bases().size() - w.a.clip_left - w.a.clip_right);
        if (o.region) {
          require(o.region_begin < seq.size(), "Region starts beyond singleton length");
          size_t end = std::min(o.region_end, seq.size());
          outputs[".region.tsv"] << name << '\t' << o.region_begin << '\t' << end << '\t'
                                 << o.region_begin << '\t' << end << "\t.\n";
          w.a.clip_left += o.region_begin;
          w.a.clip_right += seq.size() - end;
          seq = seq.substr(o.region_begin, end - o.region_begin);
        }
        auto rle = encode(seq);
        if (!w.unique)
          ambiguous_run_calls += rle.bases.size();
        ++used;
        write_fasta(fa, name, seq);
        refs.emplace_back(name, seq.size());
        if (o.polish)
          write_fasta(outputs[".draft.fasta"], name, seq);
        for (size_t i = 0; i < rle.bases.size(); ++i)
          calls << name << '\t' << i << '\t' << rle.offsets[i] << '\t' << rle.offsets[i + 1] << '\t'
                << rle.bases[i] << '\t' << rle.lengths[i] << '\t' << o.tiers[w.a.tier].label
                << "\t1\t1\t" << (w.unique ? "UNIQUE" : "AMBIGUOUS") << '\n';
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
        if (std::string_view("ACGTN").find(ref[i]) == std::string_view::npos ||
            (i && ref[i] == ref[i - 1]))
          throw std::runtime_error("Template is not homopolymer compressed: " + name);
      std::cerr << "tig " << name << ": resolving template offset\n";
      auto inferred = o.region ? region_offset : infer_offset(work, tig.length(), ref, o.max_error);
      // Offset-only reads must not remain resident through consensus calling.
      for (auto &w : work)
        if (!w.reason.empty() && o.cram.empty())
          w.a.shared_sequence.reset();
      int64_t offset = inferred.value;
      outputs[".template-offsets.tsv"]
          << name << '\t' << offset << '\t' << inferred.attempted << '\t' << inferred.placed << '\t'
          << inferred.validated << '\t' << (inferred.layout_rescue ? offset_local_span : 20) << '\t'
          << (inferred.layout_rescue ? "LAYOUT_OVERLAP"
              : inferred.internal    ? "INTERNAL_OVERLAP"
              : inferred.sparse      ? "SPARSE"
                                     : "VALIDATED")
          << '\n';
      std::cerr << "tig " << name << ": inferred template offset " << offset << " from "
                << inferred.placed << " reads\n";
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
      };
      // Process unique tiers before pooled tiers without multiplying dense
      // chromosome storage by tier count. Below-threshold unique observations
      // survive sparsely until their same-tier ambiguous observations arrive.
      // Only the processing index changes: work remains in original order for
      // placement bookkeeping, primary/secondary assignment and CRAM output.
      const size_t nt = o.tiers.size();
      std::vector<std::vector<size_t>> groups(2 * nt);
      for (size_t i = 0; i < work.size(); ++i)
        groups[(work[i].unique ? 0 : nt) + std::max(0, work[i].a.tier)].push_back(i);
      std::vector<RunBin> bins(ref.size());
      std::vector<uint8_t> vote_classes(ref.size()); // 1 unique, 2 ambiguous
      std::map<size_t, Histogram> overflow;
      struct SparseVote {
        size_t run;
        uint32_t length, count;
      };
      std::vector<std::vector<SparseVote>> pending(nt);
      // Retain only uncertain qualified histograms, not dense all-tier evidence.
      // Release each tier's records after its ambiguous phase has completed.
      std::vector<std::map<size_t, Histogram>> uncertain(o.consistent_ambiguous ? nt : 0);
      const RunModel unmodeled;
      size_t active_tier = 0;
      std::vector<Call> selected(ref.size());
      // Bits: 1 threshold-qualified/final, 2 model changed length, 4 pooled,
      // 8 additive higher-tier unique evidence (experimental mode calling),
      // 16 same-tier ambiguous borrowing changed the final length.
      std::vector<uint8_t> selected_flags(ref.size());
      auto add_vote = [&](size_t run, uint32_t length, uint32_t count, uint8_t kind) {
        if ((selected_flags[run] & 1) &&
            !(o.consistent_ambiguous && uncertain[active_tier].count(run)))
          return;
        auto &bin = bins[run];
        vote_classes[run] |= kind;
        if (!bin.count)
          bin = {length, count};
        else if (bin.length == length)
          bin.count += count;
        else
          overflow[run][length] += count;
      };
      std::vector<uint8_t> evidence_mask;
      std::map<std::pair<size_t, int>, Histogram> evidence;
      auto selected_sites = evidence_sites.find(name);
      if (selected_sites != evidence_sites.end()) {
        evidence_mask.resize(ref.size());
        for (size_t pos : selected_sites->second) {
          require(pos < ref.size(), "Evidence site outside template: " + name);
          evidence_mask[pos] = 1;
        }
      }
      memory_stage("starting batched HPC alignment");
      for (size_t phase = 0; phase < groups.size(); ++phase) {
        bool pooled_phase = phase >= nt;
        size_t tier = phase % nt;
        active_tier = tier;
        if (pooled_phase && o.ambiguous_fallback) {
          for (const auto &vote : pending[tier])
            add_vote(vote.run, vote.length, vote.count, 1);
          std::vector<SparseVote>().swap(pending[tier]);
        }
        const auto &order = groups[phase];
        for (size_t batch = 0; batch < order.size();) {
          // This bounds scratch, not total RSS. Always admit one read even if it
          // exceeds the base budget; also cap placement count for tiny/filtered reads.
          size_t end = batch, bases = 0;
          do {
            bases += work[order[end]].a.bases().size();
            ++end;
          } while (end < order.size() && end - batch < 4096 && bases < o.batch_bases);
          parallel(end - batch, o.threads, [&](size_t j) {
            size_t i = order[batch + j];
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
            } catch (const std::bad_alloc &) {
              throw;
            } catch (const std::exception &e) {
              w.reason = "HPC_ALIGNMENT_FAILED";
            }
          });
          for (size_t i = batch; i < end; ++i) {
            auto &w = work[order[i]];
            if (w.reason.empty() && (w.unique || (o.ambiguous_fallback && w.placements > 1))) {
              ++used;
              if (ref.compare(w.hbegin, w.hend - w.hbegin, w.rle.bases) == 0)
                ++exact;
              auto local_votes =
                  o.filter_hpc ? local_hpc_votes(w.trace, o.local_hpc_radius, o.local_hpc_error)
                               : std::vector<uint8_t>{};
              size_t q = 0, r = w.hbegin, operation = 0;
              bool filtered_read = false;
              for (char op : w.trace) {
                // Original raw clipping may cut a run. Do not vote with a truncated length.
                bool partial = (q == 0 && w.a.clip_left > 0 &&
                                w.a.bases()[w.a.clip_left - 1] == w.a.bases()[w.a.clip_left]) ||
                               (q + 1 == w.rle.bases.size() && w.a.clip_right > 0 &&
                                w.a.bases()[w.a.bases().size() - w.a.clip_right - 1] ==
                                    w.a.bases()[w.a.bases().size() - w.a.clip_right]);
                bool unflanked = o.terminal_hp_filter;
                // Stream only requested sites, retaining no per-read trace table.
                // A row reports evidence eligibility, not whether this tier won.
                if (o.evidence_read_votes && op != 'I' && r < evidence_mask.size() &&
                    evidence_mask[r]) {
                  const bool eligible = local_votes.empty() || local_votes[operation];
                  const char *status = op != '='                         ? "NONMATCH"
                                       : ref[r] == 'N'                   ? "N_TEMPLATE"
                                       : partial                         ? "PARTIAL_RUN"
                                       : unflanked                       ? "UNFLANKED_RUN"
                                       : !eligible                       ? "LOCAL_FILTER"
                                       : evidence_excluded.count(w.a.id) ? "DIAGNOSTIC_EXCLUDED"
                                                                         : "ELIGIBLE";
                  size_t qb = w.a.clip_left + w.rle.offsets.at(q);
                  size_t qe = op == 'D' ? qb : w.a.clip_left + w.rle.offsets.at(q + 1);
                  outputs[".run-read-votes.tsv"]
                      << name << '\t' << r << '\t' << w.a.id << '\t' << w.a.name << '\t'
                      << o.tiers[w.a.tier].label << '\t' << w.placements << '\t' << w.a.reverse
                      << '\t' << w.hbegin << '\t' << w.hend << '\t' << op << '\t'
                      << (op == 'D' ? 0 : w.rle.lengths[q]) << '\t' << status << '\t' << qb << '\t'
                      << qe << '\n';
                }
                if (op == '=' && ref[r] != 'N' && !partial) {
                  bool eligible = !unflanked && (local_votes.empty() || local_votes[operation]);
                  if (!eligible) {
                    ++filtered_hpc_votes;
                    filtered_read = true;
                  }
                  auto len = w.rle.lengths[q];
                  if (eligible) {
                    int vote_tier = w.a.tier + (w.unique ? 0 : int(o.tiers.size()));
                    if (!evidence_mask.empty() && evidence_mask[r] &&
                        !evidence_excluded.count(w.a.id))
                      ++evidence[{r, vote_tier}][len];
                    add_vote(r, len, 1, w.unique ? 1 : 2);
                  }
                }
                if (op != 'D')
                  ++q;
                if (op != 'I')
                  ++r;
                ++operation;
              }
              locally_filtered_reads += filtered_read;
            }
            w.rle = {};
            // clear() alone retains string capacity; swap releases no-longer-needed
            // storage. Never keep per-base traces for the entire chromosome.
            if ((o.polish || extend || !o.cram.empty()) && w.reason.empty())
              w.trace = cigar(w.trace);
            else {
              std::string().swap(w.trace);
              std::string().swap(w.a.sequence);
              w.a.shared_sequence.reset();
            }
          }
          batch = end;
        }
        // Sorted sparse histories let additive fallback retain higher-tier
        // votes without another dense chromosome-sized histogram allocation.
        std::vector<size_t> higher_cursors(tier);
        // A qualifying unique tier seals the call before any ambiguous tier
        // is considered. Otherwise retain its histogram for same-tier pooling.
        for (size_t run = 0; run < ref.size(); ++run) {
          auto &bin = bins[run];
          if (!bin.count)
            continue;
          auto it = overflow.find(run);
          Call candidate{bin.length, bin.count, bin.count, int(tier)};
          if (it != overflow.end()) {
            it->second[bin.length] = bin.count;
            candidate = call({it->second}, 1, o.experiment_mean);
          }
          if (o.consistent_ambiguous && pooled_phase && (selected_flags[run] & 1)) {
            auto found = uncertain[tier].find(run);
            require(found != uncertain[tier].end(), "Missing uncertain unique histogram");
            Histogram a = it == overflow.end() ? Histogram{{bin.length, bin.count}} : it->second;
            auto m = models.find(int(tier));
            const auto &model = m == models.end() ? unmodeled : m->second;
            auto before = selected[run].length;
			/*
            selected[run].length = model.borrow(ref[run], found->second, a, before);
            if (selected[run].length != before)
              selected_flags[run] |= 16;
            outputs[".borrowing.tsv"] << name << '\t' << run << '\t' << ref[run] << '\t'
                                      << o.tiers[tier].label << '\t' << selected[run].depth << '\t'
                                      << candidate.depth << '\t' << before << '\t'
                                      << selected[run].length << '\n';
									  */
            bin = {};
            vote_classes[run] = 0;
            continue;
          }
          bool qualified = candidate.depth >= o.min_support;
          if (!pooled_phase && !qualified && (o.ambiguous_fallback || o.retain_higher_unique)) {
            pending[tier].push_back({run, bin.length, bin.count});
            if (it != overflow.end())
              for (auto [length, count] : it->second)
                if (length != bin.length)
                  pending[tier].push_back({run, length, count});
          }
          if ((pooled_phase || qualified || !o.ambiguous_fallback) &&
              (qualified || selected[run].tier < 0)) {
            bool additive = false;
            Histogram mixed;
            if (o.retain_higher_unique && !pooled_phase && qualified) {
              // Qualification uses this tier alone, never accumulated depth.
              // All earlier unique tiers failed that threshold at this run.
              for (size_t higher = 0; higher < tier; ++higher) {
                auto &cursor = higher_cursors[higher];
                const auto &history = pending[higher];
                while (cursor < history.size() && history[cursor].run < run)
                  ++cursor;
                while (cursor < history.size() && history[cursor].run == run) {
                  const auto &v = history[cursor++];
                  mixed[v.length] += v.count;
                  additive = true;
                }
              }
              if (additive) {
                if (it == overflow.end())
                  mixed[bin.length] += bin.count;
                else
                  for (auto [length, count] : it->second)
                    mixed[length] += count;
                candidate = call({mixed}, 1);
              }
            }
            // Sparse calls are tentative: a later qualifying tier overrides
            // them, but another sparse tier cannot displace the first one.
            candidate.tier = int(tier + ((vote_classes[run] & 2) ? nt : 0));
            uint32_t uncalibrated = candidate.length;
            auto model = models.find(int(tier));
            // A single technology's likelihood is invalid for mixed evidence.
            // This experiment deliberately uses raw mode at mixed sites only.
            if (model != models.end() && !additive)
              candidate.length =
                  it == overflow.end()
                      ? model->second.estimate_single(ref[run], bin.length, bin.count)
                      : model->second.estimate(ref[run], it->second, candidate.length);
            selected[run] = candidate;
            selected_flags[run] = (qualified ? 1 : 0) | (candidate.length != uncalibrated ? 2 : 0) |
                                  (vote_classes[run] == 3 ? 4 : 0) | (additive ? 8 : 0);
          }
          bin = {};
          vote_classes[run] = 0;
        }
        overflow.clear();
        if (o.consistent_ambiguous && pooled_phase)
          uncertain[tier].clear();
      }
      // Release accumulation storage before allocating raw-coordinate prefixes.
      std::vector<RunBin>().swap(bins);
      std::vector<uint8_t>().swap(vote_classes);
      std::vector<std::vector<SparseVote>>().swap(pending);
      std::vector<std::vector<size_t>>().swap(groups);
      memory_stage("HPC evidence complete");
      for (const auto &[key, histogram] : evidence) {
        auto [run, tier] = key;
        for (auto [length, count] : histogram)
          outputs[".run-evidence.tsv"] << name << '\t' << run << '\t' << ref[run] << '\t'
                                       << o.tiers[tier % o.tiers.size()].label << '\t'
                                       << (tier < int(o.tiers.size()) ? "UNIQUE" : "AMBIGUOUS")
                                       << '\t' << length << '\t' << count << '\n';
      }
      evidence.clear();
      std::vector<uint8_t>().swap(evidence_mask);
      std::vector<size_t> prefix(ref.size() + 1);
      std::string draft;
      // prefix is the exact HPC-run -> draft-base coordinate transform; the
      // sentinel supports reads ending at the final template run boundary.
      for (size_t i = 0; i < ref.size(); ++i) {
        Call c = selected[i];
        model_changed_runs += bool(selected_flags[i] & 2);
        if (c.tier < 0) {
          ++no_evidence;
          require(!o.missing_error, "No evidence at " + name + ":" + std::to_string(i));
        }
        bool ambiguous_call = c.tier >= int(o.tiers.size());
        ambiguous_run_calls += ambiguous_call;
        prefix[i] = draft.size();
        draft.append(c.length, ref[i]);
        prefix[i + 1] = draft.size();
        calls << name << '\t' << i << '\t' << prefix[i] << '\t' << prefix[i + 1] << '\t' << ref[i]
              << '\t' << c.length << '\t'
              << (c.tier < 0 ? "NO_EVIDENCE" : o.tiers[c.tier % o.tiers.size()].label) << '\t'
              << c.depth << '\t' << c.peak << '\t'
              << (c.tier < 0                 ? "NONE"
                  : (selected_flags[i] & 16) ? "UNIQUE_BORROWED"
                  : (selected_flags[i] & 8)  ? "UNIQUE_ADDITIVE"
                  : (selected_flags[i] & 4)  ? "POOLED"
                  : ambiguous_call           ? "AMBIGUOUS"
                                             : "UNIQUE")
              << '\n';
      }
      std::vector<Call>().swap(selected);
      std::vector<uint8_t>().swap(selected_flags);
      memory_stage("draft called");
      size_t crop_begin = 0, crop_end = draft.size();
      if (o.region) {
        require(o.region_begin < draft.size(), "Region starts beyond draft length");
        crop_begin = o.region_begin;
        crop_end = std::min(o.region_end, draft.size());
        // Prefix calling establishes raw coordinates. Expensive raw alignment
        // and polishing only need reads overlapping the requested raw context.
        size_t context_begin = crop_begin > o.region_flank ? crop_begin - o.region_flank : 0;
        size_t context_end = std::min(draft.size(), o.region_end + o.region_flank);
        for (auto &w : work)
          if (w.reason.empty() &&
              (prefix[w.hend] <= context_begin || prefix[w.hbegin] >= context_end)) {
            w.reason = "OUTSIDE_REGION_CONTEXT";
            std::string().swap(w.trace);
            w.a.shared_sequence.reset();
          }
      }
      if (o.polish || extend || !o.cram.empty()) {
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
          } catch (const std::bad_alloc &) {
            throw;
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
      if (extend) {
        std::vector<EndRead> terminal;
        for (auto &w : work)
          if (w.reason.empty())
            terminal.push_back({&w.a, w.unique, w.a.clip_left - w.original_clip_left,
                                w.a.clip_right - w.original_clip_right,
                                w.unique || w.placements > 1});
        size_t before = draft.size();
        // Regional crops are not biological ends. Right extension is allowed
        // only if processing reached the complete original template boundary.
        auto added = extend_ends(
            draft, terminal, o.tiers.size(), models, o.end_support, o.max_end_extension,
            o.max_error, o.ambiguous_fallback, o.extend_left && crop_begin == 0,
            o.extend_right && crop_end == before && ref.size() == complete_template_size);
        left_added += added.left;
        right_added += added.right;
        outputs[".end-extensions.tsv"] << name << '\t' << before << '\t' << added.left << '\t'
                                       << added.right << '\t' << added.left_candidates << '\t'
                                       << added.right_candidates << '\n';
        crop_end += added.left + added.right;
        // Interior traces were already validated; a uniform coordinate shift
        // does not require another chromosome-wide decode/MD pass.
        for (const auto &r : terminal)
          if (r.left || r.right) {
            size_t nm;
            md(*r.alignment, draft, nm);
          }
        std::cerr << "tig " << name << ": extended left " << added.left << ", right " << added.right
                  << " bases (anchored molecules " << added.left_candidates << ", "
                  << added.right_candidates << ")\n";
      }
      std::vector<Accepted> edits;
      if (o.polish) {
        write_fasta(outputs[".draft.fasta"], name,
                    std::string_view(draft).substr(crop_begin, crop_end - crop_begin));
        std::vector<const Alignment *> alns;
        for (auto &w : work)
          if (w.reason.empty() && (w.unique || (o.ambiguous_fallback && w.placements > 1))) {
            // Borrow alignments without copying their traces. Temporary tier
            // ranks preserve the unique-first policy used by run voting.
            if (!w.unique)
              w.a.tier += o.tiers.size();
            alns.push_back(&w.a);
          }
        edits =
            polish(draft, alns, o.tiers.size() * (o.ambiguous_fallback ? 2 : 1), o.polish_support,
                   o.polish_fraction, o.max_indel);
        for (auto &w : work)
          if (w.a.tier >= int(o.tiers.size()))
            w.a.tier -= o.tiers.size();
        // A crop boundary must not split a replacement allele or lose its VCF
        // anchor. Leave such boundary events unchanged, rather than inventing a
        // partial allele. Context is still used to evaluate interior candidates.
        if (o.region)
          edits.erase(std::remove_if(edits.begin(), edits.end(),
                                     [&](const Accepted &e) {
                                       return e.edit.begin <= crop_begin || e.edit.end >= crop_end;
                                     }),
                      edits.end());
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
          outputs[".non-hpc-edits.vcf"]
              << name << '\t' << p - crop_begin + 1 << "\t.\t" << rb << '\t' << ab
              << "\t.\tPASS\tTIER=" << o.tiers[e.tier % o.tiers.size()].label << ";DP=" << e.depth
              << ";AC=" << e.support << '\n';
        }
      }
      // No later phase needs the old draft. Fast/CRAM-only runs and unchanged
      // polished tigs can transfer its storage instead of copying a chromosome.
      std::string final = edits.empty() ? std::move(draft) : hpc::apply(draft, edits);
      size_t final_begin = lift(crop_begin, edits), final_end = lift(crop_end, edits);
      std::string region_reference;
      if (o.region) {
        outputs[".region.tsv"] << name << '\t' << crop_begin << '\t' << crop_end << '\t'
                               << final_begin << '\t' << final_end << '\t' << ref.size() << '\n';
        region_reference = std::move(final);
        final = region_reference.substr(final_begin, final_end - final_begin);
      }
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
                // Realign against the polished context before rebasing to the crop.
                auto target =
                    (o.region ? region_reference : final).substr(w.a.begin, w.a.end - w.a.begin);
                w.a.ops =
                    align(w.a.bases().substr(w.a.clip_left,
                                             w.a.bases().size() - w.a.clip_left - w.a.clip_right),
                          target, o.max_error);
                w.a.packed_ops = false;
              }
              if (o.region && !crop_alignment(w.a, final_begin, final_end))
                continue;
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
            } catch (const std::bad_alloc &) {
              throw;
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
                             << ",\n  \"ambiguous_run_calls\": " << ambiguous_run_calls
                             << ",\n  \"model_changed_runs\": " << model_changed_runs
                             << ",\n  \"filtered_hpc_votes\": " << filtered_hpc_votes
                             << ",\n  \"locally_filtered_reads\": " << locally_filtered_reads
                             << ",\n  \"non_hpc_edits\": " << n_edits
                             << ",\n  \"left_extension_bases\": " << left_added
                             << ",\n  \"right_extension_bases\": " << right_added
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
