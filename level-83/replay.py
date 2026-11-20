"""Rebuild the arithmetic certificates for the search below degree 57 on X_1(83)."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work",type=Path)
    parser.add_argument("--full",action="store_true",help="repeat all larger enumerations and sieves")
    parser.add_argument("--hecke",action="store_true",help="also verify Hecke correspondences through 31")
    parser.add_argument("--rebuild-jets",action="store_true")
    parser.add_argument("--threads",type=int,default=8)
    a=parser.parse_args()
    src=Path(__file__).resolve().parent
    data=src/"inputs"
    ref=src/"records"
    w=a.work.resolve()
    w.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,OMP_NUM_THREADS=str(a.threads))
    def run(args,log=None):
        with (w/log).open("w") if log else open(os.devnull,"w") as f:
            subprocess.run(list(map(str,args)),check=True,cwd=w,env=env,stdout=f if log else None)
    def matched(name):
        if name.endswith(".json"):
            assert json.loads((w/name).read_text())==json.loads((ref/name).read_text()),name
        else:
            assert (w/name).read_bytes()==(ref/name).read_bytes(),name
    def py(name,*args):
        run([sys.executable,src/name,*args])
    def sage(name,*args):
        run(["sage","-python",src/name,*args])
    for bound in (89,113):
        name=f"linear{bound}.json"
        py("linear.py",bound,name)
        matched(name)
    ser=data/"series.dat"
    sage("exact.py","linear89.json",ser,"linear89-norms.json")
    matched("linear89-norms.json")
    sage("classes.py","linear89-norms.json","class-witnesses.json")
    matched("class-witnesses.json")
    jets=data/"eis83.npz"
    if a.rebuild_jets:
        sage("eisenjets.py",83,83)
        jets=w/"eis83.npz"
    sage("second.py",jets,ref/"minima-all.json","second.json")
    matched("second.json")
    from exact import section_norm,make_basis
    basis,anti=make_basis(ser,57)
    norms=json.loads((ref/"minima-all.json").read_text())["x_norms"]
    for i,(s,weights) in enumerate([(((3,34),(35,40)),[1,-1]),(((2,9),(18,33)),[1,1])]):
        _,_,f=section_norm(s,weights,57,basis,anti)
        assert list(f)==norms[i]
    if not (a.full or a.hecke):
        print("Exact norms, unrestricted cusp search and isogeny certificates verified.")
        return
    try:
        flags=shlex.split(subprocess.check_output(["pkg-config","--cflags","--libs","flint"],text=True,
                                                   stderr=subprocess.DEVNULL))
    except (OSError,subprocess.CalledProcessError):
        flags=shlex.split(subprocess.check_output(["sage","-sh","-c","pkg-config --cflags --libs flint"],
                                                  text=True))
    def compile_(name,flint=False,extra=()):
        run(["g++","-O3","-std=c++17","-fopenmp","-I",src,*extra,src/f"{name}.cpp","-o",w/name,
             *(flags+["-lgmpxx","-lgmp"] if flint else [])])
    if a.full:
        for tool in ("five","three","fourterms","nonsplit"):
            compile_(tool)
        compile_("trinorm",True)
        compile_("norm_five",True,["-DUNIT_FACTORS=5"])
        def enum(name,tool,args):
            run([w/tool,*args,f"{name}.dat"],f"{name}-hist.txt")
            matched(f"{name}-hist.txt")
        enum("five","five",[115])
        enum("three2","three",[2,4,280])
        enum("three3","three",[3,2,115])
        enum("fourterms","fourterms",[115])
        enum("nonsplit4","nonsplit",[4,4,300])
        enum("nonsplit5","nonsplit",[5,2,300])
        py("pencils.py","linear113.json",100,"height100")
        for key,fac,terms in [("five",0,2),("three2",2,3),("three3",3,3),("fourterms",3,4),
                              ("nonsplit5",5,3),("height100-2",2,2),("height100-3",2,3),
                              ("height100-4",2,4)]:
            py("five_sieve.py",f"{key}.dat",ser,key+"-run","--engine",w/("trinorm" if fac else "norm_five"),
               "--factors",fac,"--terms",terms,"--threads",a.threads,"--jobs",1)
            matched(key+"-run/record.json")
            assert json.loads((w/(key+"-run/remaining.json")).read_text())==[]
        sage("minima.py","linear89.json","five.dat","three3.dat",ser,"minima-all.json",
             "--fourterms","fourterms.dat")
        matched("minima-all.json")
    if a.hecke:
        compile_("hecke",True)
        sage("series_hecke.py",1850,"hecke-series.dat")
        original=ser.read_text().splitlines()
        longer=(w/"hecke-series.dat").read_text().splitlines()
        assert all(u.split()==v.split()[:int(original[0])] for u,v in zip(original[1:],longer[1:]))
        py("isogenies.py",w/"hecke","hecke-series.dat",ref/"minima-all.json","isogenies.json",
           "--jobs",min(a.threads,4))
        matched("isogenies.json")
    print("Checked every requested record.")


if __name__=="__main__":
    main()
