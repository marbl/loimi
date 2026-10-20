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
    run(base + ["--no-hpc-filter", "--output", root / "unfiltered-exact"])
    assert sequence(root / "unfiltered-exact.fasta") == sequence(root / "fast.fasta")
    assert json.loads((root / "fast.summary.json").read_text())["filtered_hpc_votes"] == 0
    for flag, value in [("--hpc-filter-radius", "0"), ("--hpc-filter-error", "1.1")]:
        run(base + [flag, value, "--output", root / "invalid-filter"], ok=False)
    # A long accurate prefix lets a non-collinear tail pass the whole-read
    # cutoff. Its isolated '=' operations must not contaminate run lengths.
    rng_filter = random.Random(139)
    def random_hpc(length, previous=""):
        bases = []
        for _ in range(length):
            previous = rng_filter.choice([b for b in "ACGT" if b != previous])
            bases.append(previous)
        return "".join(bases)
    target = random_hpc(10000)
    # Native boundary clipping must remain recoverable at BOTH biological ends.
    core = "T" + random_hpc(498, "T")
    core += next(b for b in "ACT" if b != core[-1])
    left_tail, right_tail = "CCAATTGGCC", "GGTTAACCGG"
    terminal = [(f"anchor{i}",core,10,510) for i in range(6)]
    terminal += [(f"left{i}",left_tail+core,5,510) for i in range(5)]
    terminal += [(f"right{i}",core+right_tail,10,515) for i in range(5)]
    # Pack both orientations; extension must operate on reference-oriented
    # sequence, preserving original read identity and SAM strand flags.
    terminal = [(n, s.translate(str.maketrans('ACGT','TGCA'))[::-1],e,b)
                if n in ('left1','left3','right1','right3') else (n,s,b,e)
                for n,s,b,e in terminal]
    (root/'ends.fa').write_text(''.join(f'>{n}\n{s}\n' for n,s,_,_ in terminal))
    (root/'ends.template.fa').write_text('>1\n'+core+'\n')
    (root/'ends.layout').write_text('tig\tends\nlen\t515\nrds\t16\n'+''.join(f'ends_{n}\t{b}\t{e}\t1\n' for n,s,b,e in terminal)+'end\n')
    run([tools/'sqStoreCreate','-o',root/'ends.store','-minlength','1','-raw','-pacbio','ends',root/'ends.fa'])
    run([tools/'layoutToPackage','-layout',root/'ends.layout','-store',root/'ends.store','-output',root/'ends#.endpack'])
    endbase=[exe,'--import',next(root.glob('*.endpack')),'--template-hpc',root/'ends.template.fa','--tier','anchors=ends_anchor*','--tier','tails=*']
    expected='AATTGGCC'+core+'GGTTAACC'
    for label,flags,expected_seq in [
        ('ends-fast',[],core),('ends-enabled',['--extend-ends'],expected),
        ('ends-polish',['--polish-non-hpc'],expected),
        ('ends-disabled',['--polish-non-hpc','--no-extend-ends'],core),
        ('ends-low-depth',['--extend-ends','--min-end-support','6'],core),
        ('ends-right',['--extend-ends','--end-extension-side','right'],core+'GGTTAACC')]:
        run(endbase+flags+['--output',root/label])
        actual_seq = sequence(root/(label+'.fasta'))
        assert actual_seq==expected_seq,(label,len(actual_seq),len(expected_seq))
    if shutil.which('samtools'):
        run(endbase+['--polish-non-hpc','--cram',root/'ends.cram','--output',root/'ends-cram'])
        checked=run(['samtools','calmd',root/'ends.cram',root/'ends-cram.fasta'])
        assert 'different' not in checked.stderr
        assert len(checked.stdout.splitlines()) >= 16
    run(endbase+['--extend-ends','--region','tig00000001:100-400','--output',root/'ends-region'])
    assert sequence(root/'ends-region.fasta')==core[100:400]
    corrupt = target[:9000] + random_hpc(1000, target[8999])
    double = lambda s: "".join(b * 2 for b in s)
    local_reads = [(f"anchor{i}", double(target[:1500])) for i in range(6)]
    local_reads += [("good", double(target)), ("bad", double(corrupt[:9000]) + corrupt[9000:])]
    (root / "local.fa").write_text("".join(f">{n}\n{s}\n" for n,s in local_reads))
    local_rows = [f"local_{n}\t0\t{1500 if n.startswith('anchor') else 10000}\t1" for n,s in local_reads]
    (root / "local.layout").write_text("tig\tlocal\nlen\t10000\nrds\t8\n" + "\n".join(local_rows) + "\nend\n")
    (root / "local.template.fa").write_text(">1\n" + target + "\n")
    run([tools / "sqStoreCreate", "-o", root / "local.seqStore", "-minlength", "1", "-raw", "-pacbio", "local", root / "local.fa"])
    run([tools / "layoutToPackage", "-layout", root / "local.layout", "-store", root / "local.seqStore", "-output", root / "local#.cnspack", "-idmap", root / "local.maps"])
    local_base = [exe, "--import", next(root.glob("local*.cnspack")), "--template-hpc", root / "local.template.fa", "--tier", "reads=*", "--threads", "2"]
    (root / "local.sites").write_text("tig00000001\t9500\n")
    run(local_base + ["--evidence-sites", root / "local.sites", "--output", root / "local-filtered"])
    run(local_base + ["--no-hpc-filter", "--output", root / "local-unfiltered"])
    def bad_lengths(prefix):
        rows = (root / (prefix + ".run-calls.tsv")).read_text().splitlines()[1:]
        return sum(int(row.split("\t")[5]) != 2 for row in rows)
    assert bad_lengths("local-unfiltered") > 100
    assert bad_lengths("local-filtered") < bad_lengths("local-unfiltered") / 4
    metrics = json.loads((root / "local-filtered.summary.json").read_text())
    assert metrics["filtered_hpc_votes"] > 100 and metrics["locally_filtered_reads"] == 1
    local_evidence = (root / "local-filtered.run-evidence.tsv").read_text().splitlines()[1:]
    assert len(local_evidence) == 1 and local_evidence[0].split("\t")[5:] == ["2", "1"]
    assert sequence(root / "local-filtered.fasta").startswith(double(target[:8500]))
    if shutil.which("samtools"):
        run(local_base + ["--output", root / "local-cram", "--cram", root / "local.cram"])
        assert sequence(root / "local-cram.fasta") == sequence(root / "local-filtered.fasta")
        run(["samtools", "quickcheck", root / "local.cram"])
        records_out = run(["samtools", "view", root / "local.cram"]).stdout.splitlines()
        assert any(row.split("\t")[0] == "local_bad" for row in records_out)
    # Calibration exports every tier before priority selection, but only at
    # requested HPC indexes. Read exclusions must not alter the called FASTA.
    (root / "sites.tsv").write_text("tig00000001\t0\n")
    first_high = next(row.split("\t")[1] for row in
                      (root / "fast.read-decisions.tsv").read_text().splitlines()[1:]
                      if row.split("\t")[2] == "reads_high_0")
    (root / "exclude.ids").write_text(first_high + "\n")
    run(base + ["--evidence-sites", root / "sites.tsv", "--evidence-exclude-reads",
                root / "exclude.ids", "--output", root / "evidence"])
    assert (root / "evidence.fasta").read_bytes() == (root / "fast.fasta").read_bytes()
    evidence = [r.split("\t") for r in (root / "evidence.run-evidence.tsv").read_text().splitlines()[1:]]
    assert all(r[1] == "0" for r in evidence)
    observed = {(r[3], r[4], int(r[5])): int(r[6]) for r in evidence}
    assert observed[("high", "UNIQUE", 5)] == 2
    assert observed[("high", "UNIQUE", 6)] == 2
    assert observed[("high", "UNIQUE", 20)] == 1
    assert observed[("low", "UNIQUE", 2)] == 10
    assert any(r[4] == "AMBIGUOUS" for r in evidence)
    # An optional per-tier model changes length, not tier eligibility. Missing
    # base rows retain mode, and malformed/unknown model assignments fail early.
    model_rows = ["HPC_RUN_MODEL_V1 20 5"]
    for length in range(1, 7):
        weights = [1.0] * 21
        for observed in ([5, 6, 20] if length == 6 else [length]):
            weights[observed] = 100.0
        total = sum(weights)
        model_rows.append(f"A {length} 100 " + " ".join(str(x/total) for x in weights))
    (root / "high.model").write_text("\n".join(model_rows) + "\n")
    run(base + ["--tier-model", f"high={root / 'high.model'}", "--output", root / "modeled"])
    assert sequence(root / "modeled.fasta") == "AAAAAACCGTTT"
    assert json.loads((root / "modeled.summary.json").read_text())["model_changed_runs"] == 1
    if shutil.which("samtools"):
        run(base + ["--tier-model", f"high={root / 'high.model'}", "--output", root / "modeled-cram",
                    "--cram", root / "modeled.cram"])
        assert sequence(root / "modeled-cram.fasta") == "AAAAAACCGTTT"
        checked = run(["samtools", "calmd", root / "modeled.cram", root / "modeled-cram.fasta"])
        assert "different" not in checked.stderr
    run(base + ["--tier-model", f"unknown={root / 'high.model'}", "--output", root / "bad-model-tier"], ok=False)
    run(base + ["--tier-model", f"high={root / 'high.model'}", "--experiment-run-mean",
                "--output", root / "bad-model-mean"], ok=False)
    # Raw crop boundaries can fall within a homopolymer. Polishing and CRAM are
    # independent switches; region names stay canonical and SAM is crop-relative.
    for mode in ("fast", "polish", "cram", "both"):
        if mode in ("cram", "both") and not shutil.which("samtools"):
            continue
        label = f"region-{mode}"
        args = base + ["--region", "tig00000001:2-9", "--output", root / label]
        if mode in ("polish", "both"):
            args += ["--polish-non-hpc"]
        if mode in ("cram", "both"):
            args += ["--cram", root / f"{label}.cram"]
        run(args)
        assert sequence(root / f"{label}.fasta") == "AAACCGT"
        if mode in ("cram", "both"):
            rows = run(["samtools", "view", root / f"{label}.cram"]).stdout.splitlines()
            ambiguous = [r.split("\t") for r in rows if r.startswith("reads_ambiguous\t")]
            assert len(ambiguous) == 2 and all(int(r[4]) == 3 for r in ambiguous)
            assert all(r.split("\t")[3] == "1" and "S" in r.split("\t")[5] for r in rows)
            checked = run(["samtools", "calmd", root / f"{label}.cram", root / f"{label}.fasta"])
            assert "different" not in checked.stderr
    for i, region in enumerate(("tig00000001:9-2", "tig00000001:-1-3", "tig00000002:0-10",
                                "tig00000001:100-200")):
        run(base + ["--region", region, "--output", root / f"bad-region-{i}"], ok=False)
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
    assert metrics["dropped"] == 0
    # Unique evidence wins even over a higher-quality ambiguous tier; ambiguous
    # evidence is used only in the uncovered tail and does not require CRAM.
    fallback_records = [(n, "A" * 40 + "CCGTTTGGAACCTT" if n == "ambiguous" else s)
                        for n, s in records]
    (root / "fallback.fa").write_text("".join(f">{n}\n{s}\n" for n, s in fallback_records))
    fallback_rows = [f"fallback_{n}\t0\t{8 if n == 'ambiguous' else 4}\t{2 if n == 'ambiguous' else 1}"
                     for n, s in fallback_records if n != "reverse"]
    (root / "fallback.layout").write_text(f"tig\tfallback\nlen\t8\nrds\t{len(fallback_rows)}\n" + "\n".join(fallback_rows) + "\nend\n")
    (root / "fallback.template.fa").write_text(">1\nACGTGACT\n")
    run([tools / "sqStoreCreate", "-o", root / "fallback.seqStore", "-minlength", "1", "-raw", "-pacbio", "fallback", root / "fallback.fa"])
    run([tools / "layoutToPackage", "-layout", root / "fallback.layout", "-store", root / "fallback.seqStore", "-output", root / "fallback#.fpack"])
    fallback_base = [exe, "--import", next(root.glob("*.fpack")), "--template-hpc", root / "fallback.template.fa",
                     "--tier", "amb=fallback_ambiguous", "--tier", "high=fallback_high_*", "--tier", "low=*"]
    run(fallback_base + ["--unique-only", "--output", root / "fallback-default"])
    run(fallback_base + ["--output", root / "fallback-enabled"])
    assert sequence(root / "fallback-default.fasta") == "AAAAACCGTTTGACT"
    assert sequence(root / "fallback-enabled.fasta") == "AAAAACCGTTTGGAACCTT"
    fallback_metrics = json.loads((root / "fallback-enabled.summary.json").read_text())
    assert fallback_metrics["no_evidence_runs"] == 0 and fallback_metrics["ambiguous_run_calls"] == 4
    # Six prefix-only anchors establish offset independently of the tail's
    # evidence. Tail tuples are (tier, unique count, ambiguous count, G length).
    def check_tier_policy(label, tiers, expected, selected_tier, depth, placement, extra=()):
        reads = [(f"anchor{i}", "ACGT", 4, 1) for i in range(6)]
        for tier, unique, ambiguous, length in tiers:
            for kind, count, copies in [("u", unique, 1), ("a", ambiguous, 2)]:
                reads += [(f"{tier}{kind}{i}", "ACGT" + "G" * length + "ACT", 8, copies)
                          for i in range(count)]
        stem = root / label
        stem.with_suffix(".fa").write_text("".join(f">{n}\n{s}\n" for n,s,_,_ in reads))
        stem.with_suffix(".layout").write_text(f"tig\t{label}\nlen\t8\nrds\t{len(reads)}\n" +
            "".join(f"{label}_{n}\t0\t{end}\t{copies}\n" for n,s,end,copies in reads) + "end\n")
        run([tools / "sqStoreCreate", "-o", str(stem) + ".store", "-minlength", "1", "-raw", "-pacbio", label, stem.with_suffix(".fa")])
        run([tools / "layoutToPackage", "-layout", stem.with_suffix(".layout"), "-store", str(stem) + ".store", "-output", str(stem) + "#.tierpack"])
        args = [exe, "--import", next(root.glob(label + "*.tierpack")), "--template-hpc", root / "fallback.template.fa"]
        for tier in ("LA", "LAO", "UL"):
            # Explicit separator avoids LA matching LAO names.
            args += ["--tier", f"{tier}={label}_{tier}[ua]*"]
        args += ["--tier", "anchors=*", *extra]
        for suffix, batch, threads in [("small", "1", "1"), ("large", "8000000", "2")]:
            output = str(stem) + "-" + suffix
            run(args + ["--output", output, "--batch-bases", batch, "--threads", threads])
            assert sequence(pathlib.Path(output + ".fasta")) == "ACGT" + "G" * expected + "ACT"
            row = pathlib.Path(output + ".run-calls.tsv").read_text().splitlines()[5].split("\t")
            assert row[6:8] == [selected_tier, str(depth)], row
            assert row[-1] == placement, row
        assert pathlib.Path(str(stem) + "-small.run-calls.tsv").read_bytes() == pathlib.Path(str(stem) + "-large.run-calls.tsv").read_bytes()
    check_tier_policy("unique-next", [("LA",4,6,2),("LAO",5,0,3)], 3,"LAO",5,"UNIQUE")
    check_tier_policy("unique-third", [("LA",2,8,2),("LAO",4,0,3),("UL",5,0,4)], 4,"UL",5,"UNIQUE")
    check_tier_policy("pooled-first", [("LA",2,3,2),("LAO",4,4,3)], 2,"LA",5,"POOLED")
    check_tier_policy("pooled-next", [("LA",2,2,2),("LAO",3,2,3)], 3,"LAO",5,"POOLED")
    check_tier_policy("sparse-first", [("LA",1,1,2),("LAO",2,2,3)], 2,"LA",2,"POOLED")
    check_tier_policy("unique-sparse", [("LA",1,5,2),("LAO",4,4,3)], 2,"LA",1,"UNIQUE", ["--unique-only"])
    check_tier_policy("unique-only-next", [("LA",4,5,2),("LAO",5,4,3)], 3,"LAO",5,"UNIQUE", ["--unique-only"])
    check_tier_policy("borrow-confident", [("LA",5,8,2),("LAO",5,0,3)], 2,"LA",5,"UNIQUE", ["--experiment-consistent-ambiguous"])
    check_tier_policy("borrow-pooled", [("LA",2,3,2),("LAO",4,4,3)], 2,"LA",5,"POOLED", ["--experiment-consistent-ambiguous"])
    additive = ["--experiment-retain-higher-unique"]
    check_tier_policy("additive-next", [("LA",4,6,2),("LAO",5,0,3)], 3,"LAO",9,"UNIQUE_ADDITIVE", additive)
    check_tier_policy("additive-third", [("LA",2,8,2),("LAO",4,0,2),("UL",5,0,4)], 2,"UL",11,"UNIQUE_ADDITIVE", additive)
    check_tier_policy("additive-unique-only", [("LA",4,5,2),("LAO",5,4,3)], 3,"LAO",9,"UNIQUE_ADDITIVE", additive+["--unique-only"])
    check_tier_policy("additive-first", [("LA",5,0,2),("LAO",8,0,3)], 2,"LA",5,"UNIQUE", additive)
    check_tier_policy("additive-pooled-unchanged", [("LA",2,3,2),("LAO",4,4,3)], 2,"LA",5,"POOLED", additive)
    check_tier_policy("additive-sparse-unchanged", [("LA",1,0,2),("LAO",4,0,3)], 2,"LA",1,"UNIQUE", additive+["--unique-only"])
    # Calibration sees combined depth five, not two separate subthreshold
    # histograms, and uses the winning tier's model only.
    pooled_rows = ["HPC_RUN_MODEL_V1 20 5"]
    for length in range(1, 7):
        weights = [1.0] * 21
        weights[5 if length == 6 else (4 if length == 5 else length)] = 100.0
        pooled_rows.append(f"G {length} 100 " + " ".join(str(x/sum(weights)) for x in weights))
    (root / "pooled.model").write_text("\n".join(pooled_rows) + "\n")
    check_tier_policy("additive-no-mixed-model", [("LA",1,0,2),("LAO",5,0,5)], 5,"LAO",6,"UNIQUE_ADDITIVE",
                      additive+["--tier-model", f"LAO={root / 'pooled.model'}"])
    check_tier_policy("additive-pure-model", [("LA",5,0,5),("LAO",5,0,2)], 6,"LA",5,"UNIQUE",
                      additive+["--tier-model", f"LA={root / 'pooled.model'}"])
    check_tier_policy("pooled-model", [("LA",2,3,5),("LAO",4,0,3)], 6,"LA",5,"POOLED",
                      ["--tier-model", f"LA={root / 'pooled.model'}"])
    check_tier_policy("model-lower-unique", [("LA",2,3,5),("LAO",5,0,3)], 3,"LAO",5,"UNIQUE",
                      ["--tier-model", f"LA={root / 'pooled.model'}"])
    # Polishing uses the same unique-first groups, without leaking temporary
    # group ranks into exported read-group labels.
    pol_reads = [(f"u{i}", "ACGT", 4, 1) for i in range(6)] + [(f"a{i}", "ACGTGAGT", 8, 2) for i in range(6)]
    (root / "polfallback.fa").write_text("".join(f">{n}\n{s}\n" for n, s, _, _ in pol_reads))
    (root / "polfallback.layout").write_text("tig\tpolfallback\nlen\t8\nrds\t12\n" + "".join(f"polfallback_{n}\t0\t{e}\t{copies}\n" for n, s, e, copies in pol_reads) + "end\n")
    run([tools / "sqStoreCreate", "-o", root / "polfallback.seqStore", "-minlength", "1", "-raw", "-pacbio", "polfallback", root / "polfallback.fa"])
    run([tools / "layoutToPackage", "-layout", root / "polfallback.layout", "-store", root / "polfallback.seqStore", "-output", root / "polfallback#.pfpack"])
    args = [exe, "--import", next(root.glob("*.pfpack")), "--template-hpc", root / "fallback.template.fa", "--tier", "high=polfallback_u*", "--tier", "low=*", "--polish-non-hpc", "--output", root / "polfallback"]
    if shutil.which("samtools"):
        args += ["--cram", root / "polfallback.cram"]
    run(args)
    assert sequence(root / "polfallback.fasta") == "ACGTGAGT"
    assert "TIER=low" in (root / "polfallback.non-hpc-edits.vcf").read_text()
    if shutil.which("samtools"):
        checked = run(["samtools", "calmd", root / "polfallback.cram", root / "polfallback.fasta"])
        assert "different" not in checked.stderr
    assert "AMBIGUOUS_PLACEMENT" in (root / "fallback-default.read-decisions.tsv").read_text()
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
    args = fixbase + ["--region", "tig00000001:1-7", "--polish-non-hpc", "--output", root / "fix-region"]
    if shutil.which("samtools"):
        args += ["--cram", root / "fix-region.cram"]
    run(args)
    assert sequence(root / "fix-region.fasta") == "CATACG"
    assert sequence(root / "fix-region.draft.fasta") == "CGTACG"
    vcf_rows = [r for r in (root / "fix-region.non-hpc-edits.vcf").read_text().splitlines() if not r.startswith("#")]
    assert len(vcf_rows) == 1 and vcf_rows[0].split("\t")[1] == "2"
    if shutil.which("samtools"):
        checked = run(["samtools", "calmd", root / "fix-region.cram", root / "fix-region.fasta"])
        assert "different" not in checked.stderr
    # A regional non-HPC insertion changes final length, not the draft-defined
    # biological interval. Its CRAM must align against the lifted crop boundaries.
    (root / "ins.fa").write_text("".join(f">r{i}\nACGATACGT\n" for i in range(6)))
    (root / "ins.layout").write_text("tig\tins\nlen\t8\nrds\t6\n" + "".join(f"ins_r{i}\t0\t8\t1\n" for i in range(6)) + "end\n")
    run([tools / "sqStoreCreate", "-o", root / "ins.seqStore", "-minlength", "1", "-raw", "-pacbio", "ins", root / "ins.fa"])
    run([tools / "layoutToPackage", "-layout", root / "ins.layout", "-store", root / "ins.seqStore", "-output", root / "ins#.ipack"])
    args = [exe, "--import", next(root.glob("*.ipack")), "--template-hpc", root / "fix.template.fa", "--tier", "*", "--region", "tig00000001:1-7", "--polish-non-hpc", "--output", root / "ins-region"]
    if shutil.which("samtools"):
        args += ["--cram", root / "ins-region.cram"]
    run(args)
    assert sequence(root / "ins-region.draft.fasta") == "CGTACG"
    assert sequence(root / "ins-region.fasta") == "CGATACG"
    if shutil.which("samtools"):
        checked = run(["samtools", "calmd", root / "ins-region.cram", root / "ins-region.fasta"])
        assert "different" not in checked.stderr
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
    # A native fixture preserves serialized read names, orientation and declared
    # ambiguity. Full-span extraction must reproduce the same tig sequence.
    (root / "extract.tsv").write_text(f"{root / 'extracted'}\t0\t400\n")
    run([exe.parent / "extract-regions", packs[0], "1", "37", root / "extract.tsv"])
    extracted_args = [exe, "--import", root / "extracted.cnspack", "--template-hpc", root / "multi.template.fa", "--tier", "*", "--output", root / "extracted-result"]
    if shutil.which("samtools"):
        extracted_args += ["--cram", root / "extracted-result.cram"]
    run(extracted_args)
    assert sequence(root / "extracted-result.fasta") == fa.split(">tig00000001\n")[1].split(">")[0].replace("\n", "")
    multi_region_base = [exe, "--import", packs[0], "--template-hpc", root / "multi.template.fa", "--tier", "*"]
    for tid in (1, 3, 7, 8):
        label = f"multi-region-{tid}"
        end = 60 if tid < 7 else 5
        args = multi_region_base + ["--region", f"tig{tid:08d}:2-{end}", "--region-flank", "10", "--output", root / label]
        if shutil.which("samtools"):
            args += ["--cram", root / f"{label}.cram"]
        run(args)
        full = fa.split(f">tig{tid:08d}\n")[1].split(">")[0].replace("\n", "")
        assert sequence(root / f"{label}.fasta") == full[2:end]
        assert (root / f"{label}.fasta").read_text().count(">") == 1
        if shutil.which("samtools"):
            checked = run(["samtools", "calmd", root / f"{label}.cram", root / f"{label}.fasta"])
            assert "different" not in checked.stderr
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
