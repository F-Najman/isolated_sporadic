"""Rebuild the CGPS, census and van Hoeij extract from local checkouts."""
from pathlib import Path
import json, re, sys
if len(sys.argv) != 3:
    raise SystemExit("Usage: python import_inputs.py /path/to/least-cm-degree /path/to/x1mn-sporadic.github.io")
base=Path(__file__).resolve().parent/'data'
cgps=Path(sys.argv[1])/'Least Degrees'/'X1'
web=Path(sys.argv[2])/'data'
def levels(f):return sorted(map(int,re.findall(r'\d+', (cgps/f).read_text().split(':=')[1])))
cm={int(n):int(d) for n,d in re.findall(r'\[\*\s*(\d+),\s*(\d+),', (cgps/'dcm_list_all_min_orders_X1_10k.m').read_text()) if int(n)<=720}
assert len(cm)==720
p={'cm_minimum':[cm[n] for n in range(1,721)],'cm_no':levels('no_sporadic_CM_X1.m'),'cm_undecided':levels('unknown_X1.m')}
curves=json.loads((web/'curves.json').read_text())['curves'];p['census_gonality_low']={x['n']:x['gonality']['lb'] for x in curves if x['m']==1 and x['gonality']['lb'] is not None}
pts=json.loads((web/'points.json').read_text())['points'];w={}; proven={}
for x in pts:
 n,d=x['n'],x['degree'];item=(d,'census '+x['id']);
 if n not in w or d<w[n][0]:w[n]=item
 if x['m']==1 and x['sporadic']=='yes' and (n not in proven or d>proven[n][0]):proven[n]=item
for path in (web/'knowledge/sources').glob('vanhoeij_LowDegreePlaces*.txt'):
 for ns,ds in re.findall(r'^N\s*=\s*(\d+),\s*degv\s*=\s*(\d+)',path.read_text(),re.M):
  n,d=int(ns),int(ds)
  if n not in w or d<w[n][0]:w[n]=(d,'van Hoeij level %d degree %d'%(n,d))
p['noncuspidal']=w;p['proved_sporadic']=proven
(base/'level_inputs.json').write_text(json.dumps(p,indent=1)+'\n')
