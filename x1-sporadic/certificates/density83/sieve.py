"""Verify or replay the factor-degree sieve for the signed fibres on X_1(83).

Run `python3 sieve.py WORKDIR --replay` after creating the three case files,
series.dat and the norm2, norm3, norm4 executables there. The recorded manifest
uses a fixed chunk size. Factor and output ordering do not affect its hashes.
"""
import argparse
import gzip
import hashlib
import json
import subprocess
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from itertools import islice
from pathlib import Path

GROUPS = (
    (5, 7, 11, 13, 17, 19, 23, 29),
    (31, 37, 41, 43, 47, 53, 59, 61),
    (67, 71, 73, 79, 89, 97, 101),
    (103, 107, 109, 113, 127, 131, 137),
    (139, 149, 151, 157, 163, 167, 173),
    (179, 181, 191, 193, 197, 199, 211),
    (223, 227, 229, 233, 239, 241, 251),
    (257, 263, 269, 271, 277, 281, 283),
)
CHUNK = 32000
LIMIT = 56
ALL = (1 << (LIMIT + 1)) - 1


def audit(inputs, log, primes):
    """Match each residual divisor to its norm polynomial factor-degree rows."""
    rows = {row[:-2]: row[-2:] for row in inputs}
    assert len(rows) == len(inputs)
    recorded = {}
    pending = {}
    with (gzip.open(log, 'rt') if log.suffix == '.gz' else log.open()) as stream:
        for line in stream:
            head, tail = line.split(' left')
            data = head.split()
            n = len(inputs[0])-2
            key = tuple(map(int, data[:n]))
            degree = int(data[n])
            assert key in rows and key not in recorded and degree == rows[key][0]
            mask = (1 << (min(LIMIT, degree//2)+1)) - 2
            terms = []
            previous = -1
            for token in data[n+1:]:
                prime, factors = token.split(':')
                prime = int(prime)
                i = primes.index(prime)
                assert i > previous
                previous = i
                degrees = Counter()
                for term in factors.rstrip(',').split(','):
                    f, e = map(int, term.split('^'))
                    assert f > 0 and e > 0
                    degrees[f] += e
                assert sum(f*e for f,e in degrees.items()) == degree
                possible = 1
                for f, e in degrees.items():
                    for _ in range(e):
                        possible |= (possible << f) & ALL
                mask &= possible
                terms.append(str(prime)+'/'+','.join(f'{f}^{degrees[f]}' for f in sorted(degrees)))
            assert [j for j in range(1,LIMIT+1) if mask & (1 << j)] == list(map(int, tail.split()))
            if mask:
                pending[key] = mask
            recorded[key] = ' '.join(map(str, (*key,*rows[key]))) + ' ' + ' '.join(terms)
    assert recorded.keys() == rows.keys()
    h = hashlib.sha256()
    for key in sorted(recorded):
        h.update((recorded[key]+'\n').encode('ascii'))
    return {'rows': len(rows), 'sha256': h.hexdigest()}, pending


def run_engine(root, k, inputs, path, log, group, threads):
    with path.open('w') as out:
        for row in inputs:
            print(*row, file=out)
    primes = ','.join(map(str, GROUPS[group]))
    args = [str(root/f'norm{k}'), str(path), str(root/'series.dat'), '0',
            str(len(inputs)), str(threads), primes, str(LIMIT)]
    tmp = log.with_suffix('.tmp')
    with gzip.open(tmp, 'wt', compresslevel=1) as out, subprocess.Popen(
            args, stdout=subprocess.PIPE, text=True) as process:
        for line in process.stdout:
            out.write(line)
        assert process.wait() == 0
    tmp.replace(log)


def check_block(root, k, name, inputs, initial, threads, replay):
    work = root/'sieve'
    remaining = {row[:-2]: (1 << (min(LIMIT,row[-2]//2)+1))-2 for row in inputs}
    rows = {row[:-2]: row for row in inputs}
    report = []
    for g, primes in enumerate(GROUPS):
        if not remaining:
            break
        todo = [rows[key] for key in sorted(remaining)]
        log = work/f'{name}-{g}.gz'
        if g == 0 and initial is not None:
            log = initial
        elif replay or not log.exists():
            run_engine(root, k, todo, work/f'{name}-{g}.in', log, g, threads)
        result, possible = audit(todo, log, primes)
        remaining = {key: mask & remaining[key] for key, mask in possible.items()
                     if mask & remaining[key]}
        result['pending'] = len(remaining)
        report.append(result)
    residuals = [[*rows[key], [i for i in range(1,LIMIT+1) if mask & (1<<i)]]
                 for key, mask in sorted(remaining.items())]
    return report, residuals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdir', type=Path)
    parser.add_argument('--replay', action='store_true', help='compute every factorization anew')
    parser.add_argument('--write-manifest', type=Path)
    parser.add_argument('--await-chunks', action='store_true',
                        help='consume first-pass output from a running four-factor search')
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    root = args.workdir.resolve()
    (root/'sieve').mkdir(exist_ok=True)
    report = {}
    unresolved = []

    def run(k, name, block, initial):
        result, left = check_block(root, k, name, block, initial, args.threads, args.replay)
        return name, result, left

    for k in (2,3,4):
        with (root/f'cases{k}.dat').open() as f, ThreadPoolExecutor(args.jobs) as pool:
            # Only a bounded number of case blocks are held in memory.
            jobs = []
            def collect(job):
                name, result, left = job.result()
                report[name] = result
                unresolved.extend(left)
            index = 0
            while True:
                count = CHUNK if k == 4 else 50000
                block = [tuple(map(int,s.split())) for s in islice(f,count)]
                if not block:
                    break
                block = [row for row in block if row[-2]]
                name = f'{k}-{index:03}'
                first = None
                if not args.replay:
                    first = (root/'chunks'/f'chunk-{index:03}.gz' if k==4 else root/f'pass{k}-0.log')
                    if k==4 and args.await_chunks:
                        while not first.exists():
                            time.sleep(2)
                    if not first.exists():
                        first = None
                jobs.append(pool.submit(run,k,name,block,first))
                if len(jobs) >= args.jobs:
                    collect(jobs.pop(0))
                index += 1
            for job in jobs:
                collect(job)
    manifest = dict(blocks=dict(sorted(report.items())), remaining=sorted(unresolved))
    if args.write_manifest:
        args.write_manifest.write_text(json.dumps(manifest, indent=1)+'\n')
    else:
        expected = json.loads(Path(__file__).with_name('manifest.json').read_text())
        assert manifest == expected
    for k in (2,3,4):
        blocks = [stages for name,stages in report.items() if name.startswith(str(k)+'-')]
        print(k, 'checked', sum(s[0]['rows'] for s in blocks), 'pending by stage',
              [sum(s[g]['pending'] for s in blocks if len(s)>g)
               for g in range(max(map(len,blocks)))], flush=True)
    print('Remaining fibres:', len(unresolved))


if __name__ == '__main__':
    main()
