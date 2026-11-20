import numpy as np,sys,itertools,collections,random
from sage.all import matrix,GF
D=np.load('eis180.npz');p=int(D['prime']);R=matrix(GF(p),D['rows']);O=[]
for w,m,start in D['profile']:
 w,m,start=map(int,(w,m,start));O.append((w,[list(range(start+j*m,start+(j+1)*m)) for j in range(w)]))
def ix(o,u):
 w,points=o;return [row for seq in points for row in seq[:u//w+(w==1)]]
if len(sys.argv)==1:
 for s in range(4,10):
  rows=sum([ix(o,s) for o in O],[]);print(s,len(rows),R.matrix_from_rows(rows).rank(),flush=True)
else:
 G=GF(p);k=81;lookup={v:(i,w,ell+1-(w==1)) for i,(w,pts) in enumerate(O) for seq in pts for ell,v in enumerate(seq)}
 def loww(z,ids):
  choices={}
  for j,y in enumerate(z):
   if y:
    i,degree,extra=lookup[ids[j]];weight=degree*extra
    if i not in choices or weight>choices[i][0]:choices[i]=(weight,j)
  return sorted(choices.values(),reverse=True)
 def exclude(K,ids):
  h=K.nrows()
  if not h:return True
  chosen=[];total=0
  for x,j in loww(K.row(0),ids):
   chosen.append(j);total+=x
   if total>k:break
  if total<=k:return False
  if h==1:return True
  gone=set()
  for j in chosen:
   col=K.column(j);position=next(i for i,y in enumerate(col) if y); key=tuple(col/col[position]);
   if key in gone:continue
   gone.add(key);sub=matrix(G,[col]).right_kernel_matrix()*K
   if not exclude(sub,ids):return False
  return True
 s,cap,t=map(int,sys.argv[1:4]);base=sum([ix(o,s) for o in O],[]);S=set(base);B=R.matrix_from_rows(base);print('base',s,len(base),B.rank(),flush=True);T=R*B.right_kernel_matrix().transpose();ext=[tuple(j for j in ix(o,cap) if j not in S) for o in O]; bad=0;tested=0; hist=collections.Counter()
 from sage.all import identity_matrix
 extmat=B.stack(identity_matrix(G,R.ncols()));fr=extmat.matrix_from_rows(extmat.transpose().pivots());Z=(R*fr.inverse()).matrix_from_columns(range(len(base)))
 for which in itertools.combinations(range(len(O)),t):
  U=sum([ext[i] for i in which],());tested+=1
  M=T.matrix_from_rows(U);r=M.rank();hist[len(U)-r]+=1
  if r<len(U):
   K=M.left_kernel_matrix();K=K.augment(-K*Z.matrix_from_rows(U));
   if not exclude(K,list(U)+base):bad+=1;print("possible",which,flush=True)
  if tested%500==0:print(tested,hist,flush=True)
 print('TAIL',s,cap,t,tested,bad,flush=True)
 assert (s,cap,t,tested,bad)==(7,81,2,4095,0)
 all8=sum([ix(o,8) for o in O],[]);K=R.matrix_from_rows(all8).left_kernel_matrix();assert K.nrows()==8 and exclude(K,all8);print('PASS uniform cap 8',flush=True)
