"""Verify the full-cover census and marked Hecke tests for X_1(131)."""
import argparse
import gzip
import json
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from math import prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
RANGES = {"small": (1,25,2), "lower": (26,31,4), "upper": (32,33,5), "ss95": (1,30,4)}

def valid(text):
    assert "error" not in text.lower() and "assertion failed" not in text.lower(), text[-600:]
    return text.splitlines()

def jobs(segment):
    low,high,bits = RANGES[segment]
    for shift in (0,1,2):
        cut = 0 if shift else bits
        for n in (0,1,2):
            for part in range(2**cut):
                yield f"{segment}/base/{n}_{shift}_{part}", (low,high,n,shift,part,cut)

def parse_census(segment, logs):
    cases = set()
    counts = Counter()
    for name,tag in jobs(segment):
        low,high,n,shift,shard,cut = tag
        lines = valid(logs[name]); assert lines[-1] == "complete "+" ".join(map(str,tag))
        found,reported = Counter(),{}
        for s in lines:
            v=s.split()
            if not v or v[0] not in ("case","summary"):
                continue
            e,cl,a=map(int,v[1:4]);assert low <= e <= high and (cl,a)==(n,shift)
            if v[0]=="summary":
                assert int(v[4])==shard and e not in reported
                reported[e]=int(v[-1])
                continue
            data=list(map(int,v[4:]));assert data and len(data)%5==0
            terms=sorted(tuple(data[i:i+5]) for i in range(0,len(data),5))
            assert all(j.bit_length()-1==f and m>0 for f,m,r,j,jw in terms)
            assert all(r==1 or (r==5 and f==3 and j==13) for f,m,r,j,jw in terms)
            assert sum(m for f,m,r,j,jw in terms if r==5)==a
            assert sum(f*m*r for f,m,r,j,jw in terms)==e
            case=" ".join(str(x) for term in terms for x in term)
            assert case not in cases
            cases.add(case); found[e]+=1; counts[e]+=1
        for e in range(low,high+1):
            if e != 15*shift:
                assert reported.get(e,0)==found[e]
    return cases, counts

def stage(logs,segment,key):
    return [t for n,t in logs.items() if n.startswith(f"{segment}/{key}/")]

def parse_primary(logs,cases):
    results={}; hard={}; tally=Counter()
    for t in logs:
        ls=valid(t)
        assert "forms 779 cusp rank 129 regular 650" in ls
        for s in ls:
            if not s.startswith("checked "):
                continue
            hd,row=s.split(" row "); _,e,n,failed=hd.split();e,n,failed=map(int,(e,n,failed))
            assert row not in results and 0 <= failed <= n
            v=list(map(int,row.split())); terms=[v[i:i+5] for i in range(0,len(v),5)]
            assert e==sum(f*m*r for f,m,r,j,jw in terms)
            maxdeg=max(f*r for f,m,r,j,jw in terms)
            combos=prod((65//r)**m for f,m,r,j,jw in terms)
            assert any(n==combos//(65//r) for f,m,r,j,jw in terms if f*r==maxdeg)
            results[row]=n;tally[e]+=n
            if failed:
                hard[row]=n
    assert set(results)==cases
    return hard,tally

def parse_exact(logs,hard):
    checked={}
    for t in logs:
        ls=valid(t)
        assert "forms 779 cusp rank 129 regular 650" in ls and not any(l.startswith("survivor ") for l in ls)
        for l in ls:
            if l.startswith("exact "):
                head,row=l.split(" row ");_,count,miss=head.split()
                assert row not in checked and int(miss)==0
                checked[row]=int(count)
    assert checked==hard

def audit(logs):
    assert valid(logs["arithmetic"])[-1]=="full Mordell-Weil certificate complete"
    cusp_log=valid(logs["cusps"])
    for row in (
        "three split cusps: rectangles 4032 failures 0",
        "inert orbit and three split cusps: rectangles 4032 failures 0",
        "six jets on all split cusps: rank 384 with inert differences 448",
        "six jets and three higher multiplicities: rectangles 4032 failures 0",
        "six jets, inert orbit and three higher multiplicities: rectangles 4032 failures 0",
        "all split jets through 9 and inert differences: rank 640",
    ):
        assert row in cusp_log
    assert "trace differential split vanishing 11 nonsplit vanishing 0" in valid(logs["trace_differential"])
    for segment in RANGES:
        cases,counts=parse_census(segment,logs)
        hard,tally=parse_primary(stage(logs,segment,"primary"),cases)
        parse_exact(stage(logs,segment,"exact"),hard)
        print(segment, "base configurations", dict(sorted(counts.items())))
        print(segment, "marked choices", sum(tally.values()), "exact rechecks", len(hard), "survivors 0")

def schedule(cases, workers):
    groups=[[] for _ in range(min(len(cases),workers))]
    load=[0]*len(groups)
    weights=[]
    for row in cases:
        data=list(map(int,row.split()));terms=[data[i:i+5] for i in range(0,len(data),5)]
        terms.sort(key=lambda v:-v[0]*v[2])
        sizes=[65//r for f,m,r,j,jw in terms for _ in range(m)]
        weights.append((prod(sizes[1:])+32,row))
    for w,row in sorted(weights,reverse=True):
        i=min(range(len(groups)),key=load.__getitem__);load[i]+=w;groups[i].append(row)
    return groups

def replay(dest,workers):
    dest=dest.resolve();dest.mkdir(parents=True,exist_ok=True); logs={}
    def run(job):
        name,cmd,data=job;prefix=name.replace("/","_")
        if data is not None:
            inp=dest/(prefix+".dat");inp.write_text("\n".join(data)+"\n")
            cmd=cmd[:-1]+["file:="+str(inp),cmd[-1]]
        text=subprocess.check_output(cmd,cwd=HERE,stderr=subprocess.STDOUT,text=True)
        (dest/(prefix+".log")).write_text(text);valid(text);return name,text
    def batch(js):
        with ThreadPoolExecutor(workers) as ex:logs.update(ex.map(run,js))
    batch([("arithmetic",["sage","-python","mordell_weil.py"],None),
           ("cusps",["magma","-b","cusps.m"],None),
           ("trace_differential",["magma","-b","trace_differential.m"],None)])
    for seg in RANGES:
        batch([(key,["magma","-b","ss:="+str(int(seg=="ss95")),*(f"{a}:={v}" for a,v in zip(("lo","hi","cl","shift","part","cut"),args)),"census.m"],None)
               for key,args in jobs(seg)])
        cases,_=parse_census(seg,logs)
        batch([(f"{seg}/primary/{i}",["magma","-b","hecke.m"],rows) for i,rows in enumerate(schedule(cases,workers))])
        hard,_=parse_primary(stage(logs,seg,"primary"),cases)
        batch([(f"{seg}/exact/{i}",["magma","-b","exact.m"],rows) for i,rows in enumerate(schedule(hard,workers))])
    with gzip.open(dest/"full-transcript.json.gz","wt") as f:json.dump(logs,f)
    return logs

if __name__=="__main__":
    args=argparse.ArgumentParser(description=__doc__)
    args.add_argument("--replay",type=Path,metavar="DIRECTORY",help="rerun all Magma and Sage computations")
    args.add_argument("--jobs",type=int,default=4)
    opt=args.parse_args()
    if opt.replay:
        assert opt.jobs>=1
        data=replay(opt.replay,opt.jobs)
    else:
        with gzip.open(HERE/"full-transcript.json.gz","rt") as f:data=json.load(f)
    audit(data)
