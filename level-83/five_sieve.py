"""Run the finite-field norm sieve on all retained five-factor fibres."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from itertools import islice
import json
from pathlib import Path
import subprocess


GROUPS = ["5,7,11,13,17,19,23,29", "31,37,41,43,47,53,59,61",
          "67,71,73,79,89,97,101", "103,107,109,113,127,131,137"]


def sieve_block(block, engine, series, threads, primes, path, factors, terms):
    nkeys = terms*(factors+1) if factors else 11
    pos = -1 if factors else -2
    keys = {tuple(row[:nkeys]): (row, mask) for row, mask in block}
    assert len(keys) == len(block)
    path.write_text("".join(" ".join(map(str, row))+"\n" for row, _ in block))
    command = [str(engine), str(path), str(series), "0", str(len(keys)),
               str(threads), primes, "56"]
    if factors:
        command.extend([str(factors),str(terms)])
    p = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    patterns, remaining = [], []
    for line in p.stdout:
        fields, left = line.split(" left")
        fields = fields.split()
        key = tuple(map(int, fields[:nkeys]))
        row, old = keys.pop(key)
        d = int(fields[nkeys])
        assert d == row[pos] and d>56
        mask = (1 << (min(56, d//2)+1))-2
        for term in fields[nkeys+1:]:
            prime, pattern = term.split(":")
            assert prime in primes.split(",")
            total = 0
            sums = 1
            for item in pattern.rstrip(",").split(","):
                degree, count = map(int, item.split("^"))
                total += degree*count
                for _ in range(count):
                    sums |= sums << degree
            assert total == d
            mask &= sums
        assert list(map(int, left.split())) == [i for i in range(1, 57) if (mask>>i)&1]
        if mask & old:
            remaining.append((row, mask & old))
        patterns.append(line)
    assert p.wait() == 0 and not keys
    return sha256("".join(sorted(patterns)).encode()).hexdigest(), remaining


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("cases", type=Path)
    cli.add_argument("series", type=Path)
    cli.add_argument("work", type=Path)
    cli.add_argument("--engine", type=Path, default=Path("./norm_five"))
    cli.add_argument("--jobs", type=int, default=3)
    cli.add_argument("--threads", type=int, default=4)
    cli.add_argument("--factors", type=int, default=0)
    cli.add_argument("--terms", type=int, default=3)
    a = cli.parse_args()
    a.work.mkdir(exist_ok=True, parents=True)
    engine, series = a.engine.resolve(), a.series.resolve()
    def rows():
        with a.cases.open() as f:
            for line in f:
                row = list(map(int,line.split()))
                d = row[-1] if a.factors else row[-2]
                assert d>56
                yield row, (1 << (min(56,d//2)+1))-2
    source = rows()
    report = []
    for stage, primes in enumerate(GROUPS):
        pending, hashes, batches = [], [], []
        total = 0
        while batch := list(islice(source, 12000)):
            batches.append(batch)
            total += len(batch)
        if not total:
            break
        def run(item):
            i, batch = item
            return sieve_block(batch, engine, series, a.threads, primes,
                               a.work/f"stage{stage}-block{i:03}.dat", a.factors, a.terms)
        with ThreadPoolExecutor(a.jobs) as pool:
            for hash_, survivors in pool.map(run, enumerate(batches)):
                hashes.append(hash_)
                pending.extend(survivors)
        pending.sort()
        report.append(dict(primes=primes, tested=total, survivors=len(pending),
                           checksums=hashes))
        (a.work/"record.json").write_text(json.dumps(report, indent=2)+"\n")
        (a.work/"remaining.json").write_text(json.dumps(pending)+"\n")
        print("stage",stage,"tested",total,"survived",len(pending),flush=True)
        source = iter(pending)


if __name__ == "__main__":
    main()
