"""Replay the Magma tests and check coverage against the degree-14 to 18 residue census."""
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

DIR = Path(__file__).resolve().parent
STAGES = ("reductions", "places", "pairs", "products", "cusps")


def replay(stage):
    result = subprocess.run(["magma", "-b", f"higher131-{stage}.m"], cwd=DIR,
                            check=True, stdout=subprocess.PIPE, text=True)
    return result.stdout.splitlines()


def audit(logs):
    totals = {}
    cusp_free = set()
    mixed_cusp = set()
    singles, pairs, repeated = set(), set(), set()
    categories = defaultdict(Counter)
    for line in logs["reductions"]:
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "total":
            b, n = map(int, fields[1:])
            totals[b] = n
            continue
        if fields[0] != "row":
            continue
        b = int(fields[1])
        c, ci = map(int, fields[5:7])
        data = list(map(int, fields[8:]))
        assert len(data) % 4 == 0
        components = tuple((f, m, j) for f, m, j, dual in zip(*[iter(data)]*4))
        assert all(j.bit_length()-1 == f and 1 <= m <= 2 for f, m, j in components)
        assert c + ci + sum(f*m for f,m,j in components) == b
        assert ci in (0, 13)
        signature = tuple(sorted(components))
        if not components:
            categories[b]["cusp"] += 1
        elif c or ci:
            mixed_cusp.add(signature)
            categories[b]["mixed"] += 1
        else:
            cusp_free.add(signature)
            if len(components) == 1:
                f,m,j = components[0]
                assert f == b and m == 1
                singles.add((f,j))
                categories[b]["single"] += 1
            else:
                assert len(components) == 2
                f,m,j = components[0]
                g,n,l = components[1]
                assert n == 1 and j != l
                categories[b]["multiple"] += 1
                if m == 1:
                    pairs.add((f,j,g,l))
                else:
                    assert (f,m,j,g) == (3,2,13,12)
                    repeated.add(l)
    assert list(totals.values()) == [241,465,945,1905,3825]
    assert [[categories[b][t] for t in ("cusp","mixed","single","multiple")]
            for b in range(14,19)] == [[2,2,5,0],[2,7,5,0],[2,12,12,0],[2,24,15,5],[2,44,20,9]]
    degree13 = {((13,1,j),) for j in (12037,12901)}
    assert mixed_cusp <= cusp_free | degree13

    tested, exceptional = set(), set()
    for line in logs["places"]:
        fields = line.split()
        if not fields or fields[0] != "case":
            continue
        b, j = map(int,fields[1:3])
        assert fields[8] == "false" and int(fields[-1]) == 25+7*b
        if int(fields[-3]) < int(fields[-1]):
            exceptional.add(j)
        tested.add((b,j))
    assert tested == singles and exceptional == {378895,524145}

    tested_pairs, hard_pairs = set(), set()
    for line in logs["pairs"]:
        if not line.startswith("pair "):
            continue
        main, rest = line.split(" ranks ")
        key = tuple(map(int,main.split()[1:]))
        ranks, want = rest.split(" expected ")
        values = list(map(int,re.findall(r"\d+",ranks)))
        assert len(values) == 13 and int(want) == 7*(key[0]+key[2])
        assert max(values) <= int(want)
        if min(values) < int(want):
            hard_pairs.add(key)
        tested_pairs.add(key)
    assert tested_pairs == pairs
    assert hard_pairs == {(3,13,15,35825),(9,617,9,865),(8,419,10,1591)}
    assert repeated == {4411,5451}

    info = next(s for s in logs["products"] if s.startswith("hard singles "))
    a, tail = info.split(" pairs ")
    b, c = tail.split(" repeated ")
    integers = lambda s: list(map(int,re.findall(r"\d+",s)))
    assert set(integers(a)) == exceptional
    assert set(zip(*[iter(integers(b))]*4)) == hard_pairs
    assert set(integers(c)) == repeated
    product_line = "canonical product degree 126 dimensions [ 5, 5 ] product rank 25 symmetric bound 21"
    assert logs["products"].count(product_line) == 2
    assert "confirmed 223 of 223" in logs["products"]
    assert {int(s.split()[3]) for s in logs["cusps"] if s.startswith("cuspidal cotangent")} == set(range(14,19))
    print("All reductions excluded in base degrees 14 through 18.")


if __name__ == "__main__":
    audit({stage: replay(stage) for stage in STAGES})
