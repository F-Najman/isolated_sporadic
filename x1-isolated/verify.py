"""Replay the isolation certificates. Use --rebuild-jets to regenerate the 131 cusp expansions."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--rebuild-jets', action='store_true')
parser.add_argument('--exact-points', action='store_true', help='also evaluate the nine small Riemann–Roch dimensions over Q')
args = parser.parse_args()


def run(command, expected=None, **kwargs):
    result = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kwargs)
    print(result.stdout, end='', flush=True)
    result.check_returncode()
    if expected and expected not in result.stdout:
        raise RuntimeError(f'missing success message from {command}')


for n, order in [(13, 19), (16, 20), (18, 21), (20, 60)]:
    run(['magma', '-b', f'n:={n}', 'negative.m'], f'total {order}')
run(['sage', '-python', 'bounds.py'], '20 60')
for n, degree in [(17, 5), (19, 7), (22, 5), (23, 10), (24, 5), (25, 6), (26, 6), (27, 6), (32, 9)]:
    message = f'level {n} degree {degree} l(P) 1 marked order {n}'
    run(['magma', '-b', f'n:={n}', 'prime:=101', 'point.m'], message)
    if args.exact_points:
        run(['magma', '-b', f'n:={n}', 'point.m'], message)

with tempfile.TemporaryDirectory() as tmp:
    engine = str(Path(tmp) / 'unit_traces')
    flags = subprocess.check_output(['pkg-config', '--cflags', '--libs', 'flint'], text=True)
    run(['c++', '-std=c++17', '-O3', 'fibre/unit_traces.cpp', '-o', engine, *shlex.split(flags)])
    run(['sage', '-python', 'irreducible131.py'], 'one closed point of degree 143', env={**os.environ, 'X1_UNIT_TRACES': engine})
run(['sage', '-python', 'winding131.py'], 'elliptic quotient at level 131 has rank one')
if args.rebuild_jets:
    run(['sage', '-python', 'eisenjets.py', '131', '130'], 'constant rank 129 / 129')
run(['sage', '-python', 'level131.py', 'eis131.npz'], 'effective degree 357 rank 356')
print('All isolation certificates passed.')
