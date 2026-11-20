"""Replay or verify the three-over-three signed fibre search at level 131."""
import argparse
import gzip
import os
from functools import lru_cache
from pathlib import Path
import subprocess

from sympy import isprime
from vanhoeij131_triples import enumerate_fibres

SEARCH_BOUND = 142
STAGES = (
    (129, (5,7,11,13,17,19,23,29)),
    (142, (31,37,41,43,47,53,59,61)),
)


def initial_mask(degree, limit):
    return (1 << (min(degree//2,limit)+1))-2


def case_data():
    rows = {(*u,s):(d,pole) for u,s,d,pole in enumerate_fibres()}
    assert len(rows) == 235516
    positive = {k:v for k,v in rows.items() if v[0]}
    assert len(positive) == 234213
    assert min(d for d,_ in positive.values()) == 143
    return rows,positive


def save_cases(path, cases):
    with Path(path).open('w') as out:
        for k,(d,pole) in cases.items():
            print(*k,d,pole,file=out)


@lru_cache(maxsize=None)
def good_prime(p):
    return p>3 and p!=131 and isprime(p)


def audit(line, stage_limit, data, masks):
    lhs, tails = line.rstrip().split(' left')
    entries = lhs.split()
    key = tuple(map(int,entries[:7]))
    d = int(entries[7])
    assert key in data and data[key][0] == d
    possible = initial_mask(d, SEARCH_BOUND)
    primes = set()
    for entry in entries[8:]:
        prime, vals = entry.split(':')
        prime = int(prime)
        assert good_prime(prime) and prime not in primes
        primes.add(prime)
        subsets, total = 1, 0
        for item in vals.rstrip(',').split(','):
            n,m = map(int,item.split('^'))
            assert n>0 and m>0
            total += n*m
            for _ in range(m):
                subsets |= subsets << n
        assert total == d
        possible &= subsets
    short = possible & initial_mask(d,stage_limit)
    assert [i for i in range(1,stage_limit+1) if short & (1<<i)] == list(map(int,tails.split()))
    masks[key] &= possible
    return key


def audit_log(stream, limit, data, masks, expected):
    seen = set()
    for line in stream:
        k = audit(line,limit,data,masks)
        assert k not in seen
        seen.add(k)
    assert seen == set(expected), (len(seen),len(expected))


def verify_archive(path, rows):
    masks = {k:initial_mask(d,SEARCH_BOUND) for k,(d,_) in rows.items()}
    with gzip.open(path,'rt') as archive:
        seen, expected = set(), set()
        for line in archive:
            if line.startswith('# norm '):
                assert seen == expected
                limit = int(line.split()[2])
                expected = {k for k,m in masks.items() if m}
                seen = set()
                continue
            key = audit(line,limit,rows,masks)
            assert key not in seen
            seen.add(key)
        assert seen == expected
    return masks


def report(masks):
    open_keys = sum(bool(v) for v in masks.values())
    small = sum(bool(v & ((1<<130)-2)) for v in masks.values())
    print('positive signed fibres',len(masks),'still compatible with a proper factor through 129',
          small,'through 142',open_keys,flush=True)
    return open_keys


def run(args, rows, positive):
    from triple_series131 import write_series
    dest = args.directory
    dest.mkdir(parents=True,exist_ok=True)
    series = dest/'triple_series131.txt'
    write_series(series)
    masks = {k:initial_mask(d,SEARCH_BOUND) for k,(d,_) in positive.items()}
    logs = []
    for stage,(limit,primes) in enumerate(STAGES,1):
        pending = {k:v for k,v in positive.items() if masks[k]}
        if not pending:
            break
        cases = dest/('triple_fibres.dat' if stage==1 else f'triple_open{stage}.dat')
        save_cases(cases, rows if stage==1 else pending)
        log = dest/f'triple_norm_pass{stage}.log'
        if not log.exists():
            tmp = log.with_suffix('.tmp')
            with tmp.open('w') as out:
                cmd = [args.engine,str(cases),str(series),'0',str(len(rows)),
                       str(args.threads),','.join(map(str,primes)),str(limit)]
                subprocess.run(cmd,stdout=out,check=True)
            tmp.replace(log)
        with log.open() as inp:
            audit_log(inp,limit,positive,masks,pending)
        logs.append((limit,log))
        print('norm pass',stage,flush=True)
        report(masks)
    if args.archive:
        with gzip.open(args.archive,'wt',compresslevel=6) as out:
            for limit,log in logs:
                out.write(f'# norm {limit}\n')
                with log.open() as inp:
                    for line in inp:
                        out.write(line)
    remaining = {k:positive[k] for k,m in masks.items() if m}
    if remaining:
        save_cases(dest/'triple_unresolved.dat', remaining)
    return masks


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--verify',type=Path)
    mode.add_argument('--directory',type=Path)
    parser.add_argument('--archive',type=Path)
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--engine',default=os.environ.get('X1_TRIPLE_NORM',
                        str(Path(__file__).resolve().with_name('triple_norm131'))))
    args = parser.parse_args()
    rows, positive = case_data()
    masks = verify_archive(args.verify,positive) if args.verify else run(args,rows,positive)
    if report(masks):
        for k,m in masks.items():
            if m:
                print('remaining',k,[i for i in range(1,143) if m & (1<<i)])
        raise SystemExit(1)
