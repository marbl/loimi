# SPDX-License-Identifier: MIT-0
"""Run two builds on identical inputs; compare artifacts and decoded CRAM.

Use a fresh --work directory. Inputs and common CLI flags follow '--'; omit
--output/--cram (this script assigns them). Timing is sequential, not concurrent.
OS time reports include subprocesses; logs retain the caller's own RSS telemetry.
"""
import argparse
import hashlib
import json
import pathlib
import platform
import re
import subprocess
import time


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=pathlib.Path, required=True)
    parser.add_argument('--after', type=pathlib.Path, required=True)
    parser.add_argument('--work', type=pathlib.Path, required=True)
    parser.add_argument('--cram', action='store_true')
    parser.add_argument('--after-first', action='store_true', help='Reverse execution order to check cache effects')
    parser.add_argument('args', nargs=argparse.REMAINDER)
    options = parser.parse_args()
    args = options.args[1:] if options.args[:1] == ['--'] else options.args
    if '--output' in args or '--cram' in args:
        parser.error('Output paths are assigned by this script')
    options.work.mkdir(parents=True, exist_ok=False)
    report = {'common_args': args, 'platform': platform.platform(), 'runs': {}, 'matching_artifacts': []}
    builds = [('before', options.before), ('after', options.after)]
    for label, executable in (reversed(builds) if options.after_first else builds):
        prefix = options.work / label
        command = [str(executable.resolve()), *args, '--output', str(prefix)]
        if options.cram:
            command += ['--cram', str(prefix) + '.cram']
        timing = ['/usr/bin/time', '-l' if platform.system() == 'Darwin' else '-v']
        started = time.monotonic()
        with (options.work / (label + '.log')).open('w') as log:
            subprocess.run(timing + command, stdout=log, stderr=log, check=True)
        report['runs'][label] = {'elapsed_seconds': time.monotonic() - started,
                                 'command': command, 'executable_sha256': digest(executable)}
        log_text = (options.work / (label + '.log')).read_text()
        rss = re.search(r'(\d+)\s+maximum resident set size', log_text)
        if not rss:
            rss = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', log_text)
        if rss:
            report['runs'][label]['os_peak_rss_bytes'] = int(rss[1]) * (1 if platform.system() == 'Darwin' else 1024)
        stages = re.findall(r'peak_rss_GiB=([0-9.]+)', log_text)
        if stages:
            report['runs'][label]['last_caller_peak_rss_gib'] = float(stages[-1])
        if options.cram:
            subprocess.run(['samtools', 'quickcheck', str(prefix) + '.cram'], check=True)
            with (options.work / (label + '.decoded.sam')).open('w') as sam:
                subprocess.run(['samtools', 'view', str(prefix) + '.cram'], stdout=sam, check=True)
    # CRAM headers contain paths and compression can vary; compare SAM records
    # in their original order instead, retaining flags, names, bases and tags.
    ignored = {'.log', '.cram', '.cram.crai'}
    suffixes = lambda label: {p.name[len(label):] for p in options.work.glob(label + '.*')} - ignored
    if suffixes('before') != suffixes('after'):
        raise AssertionError('The builds produced different artifact sets')
    for before in sorted(options.work.glob('before.*')):
        suffix = before.name[len('before'):]
        if suffix in ignored:
            continue
        after = options.work / ('after' + suffix)
        if not after.exists() or digest(before) != digest(after):
            raise AssertionError(f'Artifacts differ: {before} / {after}')
        report['matching_artifacts'].append(suffix)
    (options.work / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
