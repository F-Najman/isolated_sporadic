"""Audit the archived level-83 search. Add --factor to recompute all norms."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdir',type=Path)
    parser.add_argument('--factor',action='store_true')
    parser.add_argument('--spaces',action='store_true')
    parser.add_argument('--threads',type=int,default=8)
    args=parser.parse_args()
    src=Path(__file__).resolve().parent
    root=args.workdir.resolve()
    root.mkdir(exist_ok=True,parents=True)
    with gzip.open(src/'catalogue.txt.gz','rb') as inp,(root/'units.txt').open('wb') as out:
        shutil.copyfileobj(inp,out)
    subprocess.run([sys.executable,str(src/'weighted.py'),str(root),'--from-units',
                    str(root/'units.txt')],check=True)
    shutil.copyfile(src/'record.json',root/'record.json')
    series=src/'series.dat'
    if args.factor:
        subprocess.run(['sage','-python',str(src.parent/'density83'/'series.py'),
                        '--units',str(root)],check=True)
        assert series.read_bytes()==(root/'series.dat').read_bytes()
        flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
        subprocess.run(['g++','-O3','-std=c++17','-fopenmp',str(src/'search.cpp'),
                        '-o',str(root/'search'),*flags,'-lgmpxx','-lgmp'],check=True)
    else:
        for name in ('stage0.gz','stage1.gz'):
            shutil.copyfile(src/name,root/name)
    cmd=[sys.executable,str(src/'verify.py'),str(root),str(series),'--threads',str(args.threads)]
    subprocess.run(cmd+(['--replay'] if args.factor else []),check=True)
    if args.spaces:
        subprocess.run(['sage','-python',str(src.parent/'eisenjets.py'),'83','83'],cwd=root,check=True)
        report=root/'cusp-spaces.json'
        subprocess.run(['sage','-python',str(src/'cusp_spaces.py'),str(root/'eis83.npz'),
                        str(report)],check=True)
        assert json.loads(report.read_text())==json.loads((src/'spaces.json').read_text())


if __name__=='__main__':
    main()
