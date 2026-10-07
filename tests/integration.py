# SPDX-License-Identifier: MIT-0
"""Generate packages with real Verkko tools, then verify the standalone CLI.

All generated data stays in a temporary directory. Native package construction
exercises ABI compatibility as well as CLI behavior; it is not a performance or
biological-accuracy benchmark. CRAM checks require samtools on PATH.
"""
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile
import random

exe, verkko = map(pathlib.Path, sys.argv[1:])
tools = verkko / "build/bin"

def run(args, ok=True):
    """Assert success (or an expected failure), preserving tool logs on failure."""
    result = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    if (result.returncode == 0) != ok:
        raise AssertionError(f"{args}\n{result.stdout}\n{result.stderr}")
    return result

def sequence(path):
    """Join wrapped FASTA lines for the single-contig fixtures below."""
    return "".join(x.strip() for x in path.read_text().splitlines() if not x.startswith(">"))

with tempfile.TemporaryDirectory(prefix="hpc-consensus-test-") as td:
    root = pathlib.Path(td)
    # High tier A lengths 5,5,5,6,6,20: mode 5. Lower tier must not outvote.
    records = [(f"high_{i}", "A" * length + "CCGTTT") for i, length in enumerate([5,5,5,6,6,20])]
    records += [(f"low_{i}", "AACCGTTT") for i in range(10)]
    records += [("ambiguous", "A" * 40 + "CCGTTT"), ("reverse", "AAACGGTTTTT")]
    (root / "reads.fa").write_text("".join(f">{n}\n{s}\n" for n,s in records))
    rows = [f"reads_{n}\t0\t4\t{2 if n == 'ambiguous' else 1}\t{int(n.startswith('low_'))}" for n,s in records if n != "reverse"]
    rows += ["reads_reverse\t4\t0\t1"]
    # Two distinct placements plus an identical duplicate: export alternatives
    # once each, but never let either contribute unique consensus evidence.
    rows += ["reads_ambiguous\t4\t0\t2", "reads_ambiguous\t0\t4\t2\t0"]
    (root / "layout").write_text(f"tig\tpiece\nlen\t4\nrds\t{len(rows)}\n" + "\n".join(rows) + "\nend\n")
    (root / "template.fa").write_text(">1\nACGT\n")
    run([tools / "sqStoreCreate", "-o", root / "reads.seqStore", "-minlength", "1", "-raw", "-pacbio", "reads", root / "reads.fa"])
    run([tools / "layoutToPackage", "-layout", root / "layout", "-store", root / "reads.seqStore", "-output", root / "part#.cnspack", "-idmap", root / "maps"])
    packs = list(root.glob("*.cnspack"))
    assert len(packs) == 1
    base = [exe, "--import", packs[0], "--template-hpc", root / "template.fa", "--tier", "high=reads_high_*", "--tier", "low=*", "--threads", "2"]
    run(base + ["--output", root / "fast"])
    assert sequence(root / "fast.fasta") == "AAAAACCGTTT"
    assert (root / "fast.fasta").read_text().startswith(">tig00000001\n")
    # Template numbering need not match native package numbering. A conflicting
    # numeric header must not override the explicitly mapped template name.
    (root / "mapped.fa").write_text(">piece000111\nACGT\n>1\nTGCA\n")
    (root / "tig.map").write_text("# internal ID, template name\n\n1 piece000111\n")
    mapped = base + ["--template-hpc", root / "mapped.fa", "--tig-map", root / "tig.map"]
    run(mapped + ["--output", root / "mapped"])
    assert (root / "mapped.fasta").read_bytes() == (root / "fast.fasta").read_bytes()
    # Accept the native Verkko table heading as well as headerless maps.
    (root / "tig.map").write_text("     tigID  tig name\n----------  --------------------\n         1  piece000111\n         2  piece000002\n         3  piece000003\n")
    run(mapped + ["--output", root / "mapped-header"])
    assert (root / "mapped-header.fasta").read_bytes() == (root / "fast.fasta").read_bytes()
    for label, contents in [("duplicate-id", "1 piece000111\n1 other\n"),
                            ("duplicate-name", "1 piece000111\n2 piece000111\n"),
                            ("invalid-id", "-1 piece000111\n"),
                            ("extra-column", "1 piece000111 extra\n"),
                            ("bad-header", "tigID wrong name\n1 piece000111\n"),
                            ("bad-row-after-header", "tigID tig name\n---------- --------------------\nnot-an-id piece000111\n"),
                            ("missing-entry", "2 piece000111\n")]:
        (root / "bad.map").write_text(contents)
        run(mapped + ["--tig-map", root / "bad.map", "--output", root / label], ok=False)
    metrics = json.loads((root / "fast.summary.json").read_text())
    assert metrics["base_alignments"] == 0 and metrics["cram_records"] == 0
    assert metrics["dropped"] == 2
    assert "AMBIGUOUS_PLACEMENT" in (root / "fast.read-decisions.tsv").read_text()
    run(base + ["--output", root / "serial", "--threads", "1"])
    assert (root / "serial.fasta").read_bytes() == (root / "fast.fasta").read_bytes()
    assert (root / "serial.run-calls.tsv").read_bytes() == (root / "fast.run-calls.tsv").read_bytes()
    # Force one-read batches to catch accidental dependence on batch boundaries.
    run(base + ["--output", root / "tiny-batch", "--batch-bases", "1"])
    assert (root / "tiny-batch.run-calls.tsv").read_bytes() == (root / "fast.run-calls.tsv").read_bytes()
    for label,batch in [("support-small", "1"), ("support-large", "8000000")]:
        run(base + ["--output", root / label, "--batch-bases", batch, "--min-run-support", "7"])
        assert sequence(root / f"{label}.fasta") == "AACCGTTT"
    assert (root / "support-small.run-calls.tsv").read_bytes() == (root / "support-large.run-calls.tsv").read_bytes()
    run(base + ["--output", root / "fast"], ok=False)
    run([exe, "--import", packs[0], "--template-hpc", root / "template.fa", "--tier", "absent_*", "--output", root / "unmatched"], ok=False)
    run([exe, "--import", packs[0], "--template-hpc", root / "template.fa", "--tier", "low=reads_low_*", "--tier", "other=*", "--output", root / "reordered"])
    assert sequence(root / "reordered.fasta") == "AACCGTTT"
    # Third tier is reached when first two tiers have no eligible evidence.
    run([exe, "--import", packs[0], "--template-hpc", root / "template.fa", "--tier", "ambiguous=reads_ambiguous", "--tier", "absent=absent_*", "--tier", "ont=reads_low_*", "--unmatched", "drop", "--output", root / "third"])
    assert sequence(root / "third.fasta") == "AACCGTTT"
    (root / "interior.fa").write_text(">1\nCG\n")
    run([exe, "--import", packs[0], "--template-hpc", root / "interior.fa", "--tier", "*", "--output", root / "boundary"], ok=False)
    run(base + ["--output", root / "polished", "--polish-non-hpc"])
    assert sequence(root / "polished.fasta") == "AAAAACCGTTT"
    # Compare decoded records, not compressed file bytes, whose headers/encoding
    # may vary. Primary/secondary flags and MAPQ are behavioral contracts.
    if shutil.which("samtools"):
        for polish in (False, True):
            label = "both" if polish else "cram"
            args = base + ["--output", root / label, "--cram", root / f"{label}.cram"]
            if polish:
                args += ["--polish-non-hpc"]
            run(args)
            sam = run(["samtools", "view", root / f"{label}.cram"]).stdout.splitlines()
            assert len(sam) == len(records) + 1
            names = [line.split("\t")[0] for line in sam]
            ambiguous = [line.split("\t") for line in sam if line.split("\t")[0] == "reads_ambiguous"]
            assert len(ambiguous)==2
            assert {int(row[1]) for row in ambiguous} == {0,272}
            assert all(int(row[4])==3 for row in ambiguous)
            assert all(int(line.split("\t")[4])==60 for line in sam if not line.startswith("reads_ambiguous\t"))
            assert sum(n.startswith("reads_low_") for n in names) == 10
            assert all(int(line.split("\t")[1]) & 2048 == 0 for line in sam)
            assert (root / f"{label}.cram.crai").exists()
    else:
        print("samtools unavailable: CRAM integration skipped")
    # An incorrect non-HPC template base must change only in polishing mode.
    (root / "fix.fa").write_text("".join(f">r{i}\nACATACGT\n" for i in range(6)))
    (root / "fix.layout").write_text("tig\tfix\nlen\t8\nrds\t6\n" + "".join(f"fix_r{i}\t0\t8\t1\n" for i in range(6)) + "end\n")
    (root / "fix.template.fa").write_text(">1\nACGTACGT\n")
    run([tools / "sqStoreCreate", "-o", root / "fix.seqStore", "-minlength", "1", "-raw", "-pacbio", "fix", root / "fix.fa"])
    run([tools / "layoutToPackage", "-layout", root / "fix.layout", "-store", root / "fix.seqStore", "-output", root / "fix#.pack"])
    fixbase = [exe, "--import", next(root.glob("*.pack")), "--template-hpc", root / "fix.template.fa", "--tier", "*"]
    run(fixbase + ["--output", root / "fix-fast"])
    assert sequence(root / "fix-fast.fasta") == "ACGTACGT"
    args = fixbase + ["--output", root / "fix-polish", "--polish-non-hpc"]
    if shutil.which("samtools"):
        args += ["--cram", root / "fix.cram"]
    run(args)
    assert sequence(root / "fix-polish.fasta") == "ACATACGT"
    assert sequence(root / "fix-polish.draft.fasta") == "ACGTACGT"
    assert json.loads((root / "fix-polish.summary.json").read_text())["non_hpc_edits"] == 1
    if shutil.which("samtools"):
        rows = run(["samtools", "view", root / "fix.cram"]).stdout.splitlines()
        assert len(rows) == 6 and all("NM:i:0" in r and "MD:Z:8" in r for r in rows)
        checked = run(["samtools", "calmd", root / "fix.cram", root / "fix-polish.fasta"])
        assert "different" not in checked.stderr
    # One package mixes independent origins, sparse evidence, boundary overlaps,
    # and singletons. Deterministic nonrepetitive templates avoid accidental ties.
    rng = random.Random(42)
    records, layouts, templates = [], [], []
    for tig_id, offset in [(1, 37), (2, 113), (3, 29), (4, 41), (5, 53), (6, 67)]:
        template = "A"
        for _ in range(399):
            template += rng.choice([b for b in "ACGT" if b != template[-1]])
        templates.append(f">{tig_id}\n{template}\n")
        rows = []
        # This endpoint-spanning read preserves the original layout extent.
        # Legacy skipConsensus is ignored; count=2 excludes sparse-fixture anchors
        # from inference while still permitting their CRAM export.
        records.append((f"edge{tig_id}", template))
        rows.append(f"multi_edge{tig_id}\t0\t{offset+400}\t{1 if tig_id<3 else 2}\t1")
        for i in range(8 if tig_id<3 else tig_id-1):
            start = 130 + 10*i
            seq = template[start:start+120]
            b, e = start+offset, start+offset+120
            if i%2:
                seq = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                b, e = e, b
            name = f"t{tig_id}r{i}"
            records.append((name, seq))
            rows.append(f"multi_{name}\t{b}\t{e}\t1\t0")
        prefix = "".join("ACGT"[j%4] for j in range(offset))
        if prefix[-1] == template[0]:
            prefix = prefix[:-1] + "T"
        records.append((f"boundary{tig_id}", prefix + template[:120]))
        rows.append(f"multi_boundary{tig_id}\t0\t{offset+120}\t{1 if tig_id<3 else 2}\t0")
        layouts.append(f"tig\tt{tig_id}\nlen\t{offset+400}\nrds\t{len(rows)}\n" + "\n".join(rows) + "\nend\n")
    records += [("single7", "AACCCGT"), ("single8", "ACGGGTT")]
    layouts += ["tig\tsingle7\nlen\t4\nrds\t1\nmulti_single7\t0\t4\t1\t1\nend\n", "tig\tsingle8\nlen\t4\nrds\t1\nmulti_single8\t4\t0\t1\t1\nend\n"]
    (root / "multi.fa").write_text("".join(f">{n}\n{s}\n" for n,s in records))
    (root / "multi.layout").write_text("".join(layouts))
    (root / "multi.template.fa").write_text("".join(templates))
    run([tools / "sqStoreCreate", "-o", root / "multi.seqStore", "-minlength", "1", "-raw", "-pacbio", "multi", root / "multi.fa"])
    run([tools / "layoutToPackage", "-layout", root / "multi.layout", "-store", root / "multi.seqStore", "-partition", "100", "1", "2", "-output", root / "multi#.pkg"])
    packs = list(root.glob("*.pkg"))
    assert len(packs) == 1
    multiargs = [exe, "--import", packs[0], "--template-hpc", root / "multi.template.fa", "--tier", "*", "--polish-non-hpc", "--output", root / "multi"]
    if shutil.which("samtools"):
        multiargs += ["--cram", root / "multi.cram"]
    run(multiargs)
    offsets = [line.split("\t") for line in (root / "multi.template-offsets.tsv").read_text().splitlines()[1:]]
    assert {int(row[0].removeprefix("tig")): int(row[1]) for row in offsets if row[1] != "."} == {1:37, 2:113, 3:29, 4:41, 5:53, 6:67}
    for row in offsets:
        if 3 <= int(row[0].removeprefix("tig")) <= 6:
            assert int(row[2]) == int(row[0].removeprefix("tig"))-1 and row[6] == "SPARSE"
    fa = (root / "multi.fasta").read_text()
    assert ">tig00000007\nAACCCGT\n" in fa and ">tig00000008\nAACCCGT\n" in fa
    if shutil.which("samtools"):
        sam = run(["samtools", "view", root / "multi.cram"]).stdout.splitlines()
        assert all(row.split("\t")[2].startswith("tig000000") for row in sam)
        singles = [line.split("\t") for line in sam if "multi_single" in line]
        assert len(singles) == 2 and {int(row[1]) for row in singles} == {0,16}
        assert all("NM:i:0" in row and "MD:Z:7" in row for row in singles)
    decisions = (root / "multi.read-decisions.tsv").read_text()
    assert "TEMPLATE_BOUNDARY" not in decisions
    for tid in (1,2):
        assert any(f"multi_boundary{tid}\t" in line and line.endswith("ACCEPTED") for line in decisions.splitlines())
    if shutil.which("samtools"):
        boundaries = [line.split("\t") for line in sam if "multi_boundary" in line]
        assert len(boundaries)==6 and all("S" in row[5] and row[3]=="1" for row in boundaries)
        checked=run(["samtools", "calmd", root / "multi.cram", root / "multi.fasta"])
        assert "different" not in checked.stderr
print("native package integration passed")
