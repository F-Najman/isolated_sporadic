"""Verify the degree-95 sieve on X_1(131)/H_5 over F_3."""
import argparse
import gzip
import json
import os
import subprocess
from collections import Counter, defaultdict
from math import comb
from pathlib import Path

from exceptional import pencil_layout, cusp_and_ordinary, certify_boxes

HERE = Path(__file__).resolve().parent
LAST_SUPPORT = 8


def lines(text,prefix):
    return [line for line in text.splitlines() if line.startswith(prefix)]


def integers(line):
    return tuple(map(int,line.split()))


def eigentests(norm_text,hecke_text):
    data = defaultdict(list)
    tested = {}
    for line in lines(norm_text,"data "):
        _,i,f,j,jw = line.split()
        data[int(i)].append(int(f))
    for line in lines(hecke_text,"full Hecke fibre "):
        v = line.split()
        i,n,rank = int(v[3]),int(v[6]),int(v[9])
        assert v[8]=="[" and v[10]=="]" and len(v)==11
        assert i not in tested
        tested[i] = (n,rank)
    assert set(tested) == set(data)
    for i,degs in data.items():
        assert tested[i] == (65**(len(degs)-1),5*sum(degs))

    accepted = {}
    for line in lines(norm_text,"norm "):
        left,right = line.split(" places")
        v=left.split()
        d=int(v[2]); idx=int(v[1])
        assert idx not in accepted and len(v)==7
        a=integers(right)
        terms=list(zip(a[::3],a[1::3],a[2::3]))
        assert sum(f*m for f,m,r in terms)==d
        possible=all(m%r==0 for f,m,r in terms)
        assert possible==(v[4]=="true")
        if possible:
            assert v[6]=="true"
            ordinary=sorted(f for f,m,r in terms if m==r==1)
            inert=[t for t in terms if t[1:]!=(1,1)]
            assert all(t==(1,65,65) for t in inert)
            assert ordinary==sorted(data.get(idx,[]))
            assert ordinary or inert
        accepted[idx]=possible
    return accepted


def audit(data):
    assert "trace-zero space 130, full space 131, jets 13 by 260" in data["jets"]
    assert "eigenvalue 51, ranks 65 30 95" in data["nonrepresentable"]
    assert len(lines(data["invariant"],"invariant quotient "))==5
    invariant=[s for s in lines(data["invariant"],"invariant quotient ") if "lifts true" in s]
    assert len(invariant)==2 and {s.split(" places ")[1] for s in invariant}=={"1 5 65 true","5 1 13 false"}
    unit_status=eigentests(data["unit_factors"],data["unit_hecke"])
    assert set(unit_status)==set(range(36)) and {i for i,v in unit_status.items() if v}=={0,6,8,25}
    assert "norm 6 95 pass true unamb true places 95 1 1" in data["unit_factors"]
    assert "data 6 95 4175000167681667948618744480333547672455993093 2554018553896954323719512172264650245764679020" in data["unit_factors"]

    ex=data["exceptions"]
    assert len(set(tuple(v[1:]) for v in ex))==len(ex)
    expected=Counter(sum(bool(x) for x in v[1:]) for v in ex)
    for v in ex:
        assert len(v)==14 and sum(v[1:])==95
    cov={}
    for line in lines(data["cover"],"support "):
        v=line.split()
        k,n,tests,passed,pencils,other=map(int,[v[1],v[4],v[7],v[10],v[12],v[14]])
        assert k not in cov and comb(13,k)==13*n and tests>0
        assert passed+pencils+other==comb(13,k)*comb(94,k-1)
        assert other==13*expected[k]
        cov[k]=passed+pencils+other
    assert set(cov)==set(range(1,LAST_SUPPORT+1))
    primitive,directions=pencil_layout(ex,data["ranks"])
    pure,ordinary=cusp_and_ordinary(ex,primitive,directions,data["pencils"])
    boxes=data["boxes"]; refs=data["box_refs"]
    assert len(refs)==len(pure) and set(refs)==set(range(len(boxes)))
    for v,j in zip(pure,refs):
        assert all(a<=b for a,b in zip(v,boxes[j]))
    assert len(lines(data["anti"],"anti "))==len(boxes)
    for i,line in enumerate(lines(data["anti"],"anti ")):
        d=sum(boxes[i]);assert line==f"anti {i} {d} {4*d}"
    floors=lines(data["floors"],"fibre ")
    assert len(floors)==len(boxes)
    for i,line in enumerate(floors):
        v=line.split()
        assert int(v[1])==i and int(v[3])==sum(t//5 for t in boxes[i]) and v[6]==v[9]
    assert "absolute trace 6, projective marking orders 5 and 65" in data["frobenius17"]
    results=eigentests(data["other_factors"],data["other_hecke"])
    assoc={int(k):v for k,v in data["norm_ids"].items()}
    assert set(assoc)==set(ordinary) and set(assoc.values())==set(results)
    for s in lines(data["other_factors"],"norm "):
        v=s.split()
        assert all(ordinary[i]==int(v[2]) for i,n in assoc.items() if n==int(v[1]))
    total=sum(cov.values())
    print(f"All {total:,} effective split-cusp divisors of degree 95 with support at most {LAST_SUPPORT} are excluded as class representatives.")
    print("Unit fibres",len(unit_status),"split ordinary configurations",sum(unit_status.values()))
    print("Exceptional cusp-divisor orbits",len(ex),"distinct exceptional residual norms",len(results))


def write_transcript(path,data):
    with gzip.open(path,"wt") as f:
        json.dump(data,f,separators=(",",":"))


def build(dest,jobs,saved_cover=None):
    from units import Norms,prepare,write_rows,valuation_matrix
    from numpy import array
    dest=dest.resolve();dest.mkdir(parents=True,exist_ok=True)
    engine=Norms()
    def run(cmd,log,cwd=HERE):
        output=subprocess.check_output(cmd,cwd=cwd,stderr=subprocess.STDOUT,text=True,
                                       env=dict(os.environ,OMP_NUM_THREADS=str(jobs)))
        (dest/(log+".log")).write_text(output)
        assert "error" not in output.lower() and "assertion failed" not in output.lower(),output[-500:]
        return output
    def magma(script,name,**files):
        return run(["magma","-b",*(f"{k}:={dest/v}" for k,v in files.items()),script+".m"],name)
    data={}
    # Regenerate the integral Fricke expansions independently.
    run(["sage","-python",str(HERE.parent/"triple_series131.py"),str(dest/"fricke-series.txt")],"fricke")
    assert (dest/"fricke-series.txt").read_text()==(HERE/"fricke-series.txt").read_text()
    prepare(dest,engine)
    e=(-1,0,0,0,-1,0,0,1,0,0,0,1,0)
    unit_order=valuation_matrix @ e
    assert tuple(max(0,-v) for v in unit_order)==(4,32,1,0,23,29,0,0,2,0,2,0,2)
    ternary_fibre,_=engine.signed_fibre(e,1)
    saved=integers((dest/"unit-norms.dat").read_text().splitlines()[6])
    assert saved==(6,*ternary_fibre)
    f,_=Norms(17).signed_fibre(e,1)
    write_rows(dest/"mod17.dat",[(17,*f)])
    data["frobenius17"]=magma("frobenius17","frobenius17",file="mod17.dat")
    data["jets"]=magma("quotient_jets","jets",file="jets260.dat")
    raw=(dest/"jets260.dat").read_text().splitlines()
    assert len(raw)==13*260 and all(len(s)==131 for s in raw)
    (dest/"jets100.dat").write_text("".join(s+"\n" for i,s in enumerate(raw) if i%260<100))
    data["invariant"]=magma("invariant","invariant")
    data["nonrepresentable"]=magma("nonrepresentable","nonrepresentable")
    data["unit_factors"]=magma("factor_norms","unit_factors",file="unit-norms.dat")
    data["unit_hecke"]=magma("full_hecke","unit_hecke",file="unit_factors.log")
    if saved_cover is None:
        run(["g++","-O3","-std=c++17","-fopenmp","cusp_classes.cpp","-o",str(dest/"cover")],"compile")
        data["cover"]=run([str(dest/"cover"),str(dest/"jets100.dat"),str(dest/"unit-divisors.dat"),
                           "1",str(LAST_SUPPORT),str(dest/"exceptional.dat")],"cover")
        ex=sorted((integers(s) for s in (dest/"exceptional.dat").read_text().splitlines()),key=lambda v:v[1:])
    else:
        data["cover"],ex=saved_cover
    data["exceptions"]=ex
    write_rows(dest/"rank-input.dat",[(i,*v[1:]) for i,v in enumerate(ex)])
    data["ranks"]=magma("principal_parts","ranks",file="jets100.dat",patterns="rank-input.dat")
    primitive,directions=pencil_layout(ex,data["ranks"])
    write_rows(dest/"pencil-input.dat",[(i,*v) for i,v in enumerate(primitive)])
    data["pencils"]=magma("pencil","pencils",file="jets260.dat",patterns="pencil-input.dat",result="norm-series.dat")
    pure,_=cusp_and_ordinary(ex,primitive,directions,data["pencils"])
    data["boxes"],data["box_refs"],data["anti"],data["floors"]=certify_boxes(pure,dest,magma,write_rows)
    # Equal normalized series have the same rational function on B.
    series={}
    for s in (dest/"norm-series.dat").read_text().splitlines():
        idx,d,*v=integers(s); a=array(v,dtype=int)
        assert len(v)==d+1 and a[0] in (1,2)
        key=tuple(a*a[0]%3)
        series.setdefault(key,[]).append(idx)
    polynomials=sorted((engine.recover(key),indices) for key,indices in series.items())
    write_rows(dest/"other-norms.dat",[(j,*v) for j,(v,_) in enumerate(polynomials)])
    data["norm_ids"]={str(idx):j for j,(_,indices) in enumerate(polynomials) for idx in indices}
    data["other_factors"]=magma("factor_norms","other_factors",file="other-norms.dat")
    data["other_hecke"]=magma("full_hecke","other_hecke",file="other_factors.log")
    audit(data)
    write_transcript(dest/"transcript.json.gz",data)
    return data


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay",type=Path,metavar="DIRECTORY",help="repeat the Magma, Sage, Python and C++ calculations")
    parser.add_argument("--jobs",type=int,default=6,help="threads for the finite rank enumeration")
    opt=parser.parse_args()
    if opt.replay:
        assert opt.jobs>0
        build(opt.replay,opt.jobs)
    else:
        with gzip.open(HERE/"transcript.json.gz","rt") as f:audit(json.load(f))

