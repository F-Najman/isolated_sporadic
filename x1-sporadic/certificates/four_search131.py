"""Verify the four-factor cusp census and the irreducibility witnesses through total degree 200."""
import argparse
import gzip
import subprocess
import tempfile
from pathlib import Path

BASE = Path(__file__).resolve().parent
PRIMES = {5,7,11,13,17,19,23,29}


def verify(rows, witnesses):
    seen = set()
    for line in witnesses:
        entries = line.split()
        key = tuple(map(int,entries[:9]))
        d = int(entries[9])
        assert key in rows and rows[key] == d and key not in seen
        seen.add(key)
        assert entries[-1] == "left"
        possible = (1 << (d//2+1))-2
        used = set()
        for token in entries[10:-1]:
            p, factors = token.split(":")
            p = int(p)
            assert p in PRIMES and p not in used
            used.add(p)
            degrees, total = 1, 0
            for factor in factors.rstrip(",").split(","):
                k, multiplicity = map(int,factor.split("^"))
                assert k > 0 and multiplicity > 0
                total += k*multiplicity
                for _ in range(multiplicity):
                    degrees |= degrees << k
            assert total == d
            possible &= degrees
        assert possible == 0
    assert seen == set(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path)
    parser.add_argument("--witnesses",type=Path,default=BASE/"four131-witnesses.txt.gz")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        executable = args.census.resolve() if args.census else tmp/"census"
        if not args.census:
            subprocess.run(["c++","-O3","-std=c++17",str(BASE/"four_cusps131.cpp"),"-o",str(executable)],check=True)
        selected = tmp/"four.dat"
        output = subprocess.check_output([str(executable),"200",str(selected)],text=True)
        head,*lines = output.splitlines()
        assert head.startswith("units 30734768 signed 61469536 ")
        counts = dict(tuple(map(int,s.split())) for s in lines)
        assert sum(counts.values()) == 61469536 and counts[0] == 10447
        assert (min(d for d in counts if d),max(counts)) == (143,674)
        rows = {}
        for line in selected.read_text().splitlines():
            *key,d,pole = map(int,line.split())
            assert tuple(key) not in rows and 0<d<=200 and d<=pole
            rows[tuple(key)] = d
        assert len(rows) == sum(n for d,n in counts.items() if 0<d<=200) == 15980
        with gzip.open(args.witnesses,"rt") as inp:
            verify(rows,inp)
        print("61469536 signed fibres; minimum positive residual degree 143; all 15980 remainders through degree 200 irreducible.")


if __name__ == "__main__":
    main()
