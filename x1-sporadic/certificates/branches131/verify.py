"""Audit the recorded Magma calculations, or regenerate them with --replay."""
import argparse
import gzip
import json
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRANSCRIPT = HERE / "transcript.json.gz"
COUNTS = {3: {5: 1, 10: 2, 11: 2, 12: 5, 13: 7, 14: 14, 15: 34,
              16: 105, 17: 357, 18: 1001, 19: 2853},
          2: dict(zip(range(19, 26), (62, 121, 217, 407, 793, 1531, 3001)))}
SIDE = {3: ("ternary", "three"), 2: ("binary25", "two")}
SMALL = {"cusp25_two": "inert orbit 1 cap 6 rectangles 1333 failures 0",
         "cusp25_three": "inert orbit 1 cap 6 rectangles 1333 failures 0",
         "extra25_three": "nonvanishing anti-Fricke differential at both cusps",
         "cusp25_extra": "empty shifted classes excluded",
         "shift25_two": "residual degree 12, rank 84"}

def clean(output):
    assert "error" not in output.lower() and "failed" not in output.lower(), output[-1000:]
    return output.splitlines()

def row(values):
    return " ".join(map(str, values))

def census_jobs(p):
    for shift in range(3 if p == 3 else 2):
        lo, hi = (1, 12) if shift else ((1, 19) if p == 3 else (19, 25))
        for cls in range(25 if p == 3 else 15):
            yield f"{p}/base/{shift}-{cls}", (shift, cls, lo, hi)

def census(p, logs):
    cases = set()
    invariant = []
    for name, (shift, cls, lo, hi) in census_jobs(p):
        nmod = 5 if p == 3 else 3
        h, n = 13*(cls//nmod-2), cls % nmod
        found, reported = Counter(), Counter()
        lines = clean(logs[name])
        assert lines[-1] == f"complete {lo} {hi} {cls} {shift}"
        for line in lines:
            fields = line.split()
            if not fields or fields[0] not in ("case", "summary"):
                continue
            b, hh, nn = map(int, fields[1:4])
            assert (hh, nn) == (h, n) and lo <= b <= hi
            if fields[0] == "summary":
                assert b not in reported
                reported[b] = int(fields[-1])
                continue
            found[b] += 1
            data = list(map(int, fields[4:]))
            assert data and len(data) % 4 == 0
            factors = sorted(tuple(data[i:i+4]) for i in range(0, len(data), 4))
            assert sum(f*m for f,m,j,dual in factors) == b
            assert all(m >= 1 and p**f <= j < 2*p**f for f,m,j,dual in factors)
            if p == 2 and shift:
                invariant.append((b,h,n,factors))
            else:
                flag = [] if p == 3 else [int(n != 0)]
                cases.add(row(flag + [v for term in factors for v in term]))
        assert all(reported.get(b, 0) == found[b] for b in range(lo, hi+1))
    if p == 2:
        assert invariant == [(12,13,1,[(12,1,7963,6179)])]
    return cases

def primary(p, logs, cases):
    seen = {}
    degrees = Counter()
    for log in logs:
        local = 0
        for line in clean(log):
            v = line.split()
            if not v:
                continue
            if v[0] == "hard":
                local += 1
            if v[0] == "checked":
                b, choices, hard = map(int,v[1:4])
                case = " ".join(v[4:])
                assert case not in seen and local == hard and choices >= hard
                degrees[b] += 1
                seen[case] = choices, hard
                local = 0
        assert local == 0
    assert set(seen) == cases and degrees == COUNTS[p], (p, degrees)
    assert sum(n for n,h in seen.values()) == (20486 if p == 3 else 116663)
    prefix = "0 " if p == 3 else ""
    hard = {prefix+s: n for s,(n,f) in seen.items() if f}
    assert len(hard) == (37 if p == 3 else 353)
    return hard

def exact(logs, expected, moving=False):
    seen, survivors = {}, {}
    moved = 0
    for log in logs:
        if moving:
            assert "ordinary fibre 12707 rank 25" in clean(log)
        remaining = Counter()
        for line in clean(log):
            if line.startswith("remaining "):
                remaining[line.split(" row ")[1]] += 1
            if line.startswith("moved "):
                assert moving and line.endswith("sign 1 degree 100 rank 125"), line
                moved += 1
            if line.startswith("exact "):
                header, case = line.split(" row ")
                _, n, failed = header.split()
                assert case not in seen and remaining[case] == int(failed)
                seen[case] = int(n)
                if int(failed):
                    survivors[case] = int(failed)
        assert all(case in seen for case in remaining)
    assert seen == expected
    if moving:
        assert moved == 41
    return survivors

def stage_logs(archive, p, key):
    prefix = f"{p}/{key}/"
    return [v for k,v in archive.items() if k.startswith(prefix)]

def audit(archive):
    for name, ending in SMALL.items():
        ls = clean(archive[name])
        assert ls[-1].endswith(ending)
        if name in ("cusp25_two","cusp25_three"):
            assert len([v for v in ls if "failures 0" in v]) == 2
    assert clean(archive["prym_at_two"])[-1].endswith("773727027031129171413080804164258787077")
    for p in (3,2):
        cases = census(p,archive)
        hard = primary(p,stage_logs(archive,p,"primary"),cases)
        left = exact(stage_logs(archive,p,"exact"),hard)
        if p == 3:
            assert not left
        else:
            assert len(left) == 20 and sum(left.values()) == 41
            assert not exact(stage_logs(archive,p,"moving"),{s:hard[s] for s in left},True)
        print(f"Characteristic {p}: {len(cases)} base configurations, all markings excluded.")

def replay(dest, workers):
    dest = dest.resolve()
    dest.mkdir(parents=True,exist_ok=True)
    logs = {}
    def launch(job):
        name, argv, rows = job
        key = name.replace("/","_")
        if rows is not None:
            data = dest/(key+".dat")
            data.write_text("\n".join(rows)+"\n")
            argv = [*argv[:-1], "file:="+str(data), argv[-1]]
        out = subprocess.run(argv,cwd=HERE,stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,text=True,check=True).stdout
        (dest/(key+".log")).write_text(out)
        clean(out)
        if "/base/" in name:
            assert out.splitlines()[-1].startswith("complete ")
        return name,out
    def run(jobs):
        with ThreadPoolExecutor(workers) as pool:
            logs.update(pool.map(launch,jobs))
    def marked_jobs(p, key, script, cases, extra=()):
        cases = sorted(cases)
        count = min(workers, len(cases))
        return [(f"{p}/{key}/{i}",["magma","-b",*extra,script+".m"],cases[i::count])
                for i in range(count)]
    run([(name,["magma","-b",name+".m"],None) for name in SMALL] +
        [("prym_at_two",["sage","-python","prym_at_two.py"],None)])
    base = []
    for p in (3,2):
        for name,(shift,cls,lo,hi) in census_jobs(p):
            base.append((name,["magma","-b",f"lo:={lo}",f"hi:={hi}",f"cls:={cls}",
                               f"shift:={shift}",SIDE[p][0]+"_census.m"],None))
    run(base)
    for p in (3,2):
        cases = census(p,logs)
        run(marked_jobs(p,"primary",SIDE[p][0]+"_hecke",cases))
        hard = primary(p,stage_logs(logs,p,"primary"),cases)
        run(marked_jobs(p,"exact","exact25_"+SIDE[p][1],hard))
        left = exact(stage_logs(logs,p,"exact"),hard)
        if p == 2:
            run(marked_jobs(p,"moving","exact25_two",left,("move:=1",)))
    with gzip.open(dest/"transcript.json.gz","wt") as out:
        json.dump(logs,out)
    return logs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay",type=Path,metavar="OUTPUT_DIRECTORY",
                        help="recompute all stages with Magma and Sage")
    parser.add_argument("--jobs",type=int,default=4)
    args = parser.parse_args()
    if args.replay:
        assert args.jobs >= 1
        data = replay(args.replay,args.jobs)
    else:
        with gzip.open(TRANSCRIPT,"rt") as inp:
            data = json.load(inp)
    audit(data)
