import math,sys,itertools,collections
import numpy as np
import sympy as sp
from sage.all import matrix, GF
from pathlib import Path
N,d=map(int,sys.argv[1:3]); base=Path('.'); q=next(n for n in range(1+N,100000,N) if sp.isprime(n) and n>2*N); z=pow(int(sp.primitive_root(q)),(q-1)//N,q); roots=[pow(z,i,q) for i in range(N)];
from extra_eis import extra_eis
extraform=extra_eis(N,q,roots)
cusps=[]; orbs=[]
for c in range(N//2+1):
 h=math.gcd(c,N); ds=[i for i in range(h) if math.gcd(i,h)==1];
 if 2*c%N==0:ds=[i for i in ds if i<=(h-i)%h];
 inds=list(range(len(cusps),len(cusps)+len(ds)));orbs.append(inds);cusps.extend((c,b,h) for b in ds)
r=sum(len(x)==1 for x in orbs); extra=d-r; print('N, p, r, extra',N,q,r,extra,flush=True); M=(N-1)//2; aidx,bidx=np.triu_indices(M); C=[]; J=[]; group=[]; 
for O in orbs:
 mult=(1 if len(O)==1 else 0)+extra//len(O); start=len(J); group.append((len(O),mult,start));
 for ip in O:
  c,b,h=cusps[ip]; terms=np.zeros((M,mult+1),np.int64)
  for i in range(1,M+1):
   a=c*i%N; beta=b*i%N;
   if a: terms[i-1,0]=(pow(2,-1,q)-a*pow(N,-1,q))%q
   elif beta:terms[i-1,0]=(1+roots[beta])*pow(2*(1-roots[beta]),-1,q)%q
   for v in range(1,mult+1):
    s=0
    for e0 in sp.divisors(v*h):
     k=v*h//e0
     if e0%N==a: s+=roots[beta*k%N]
     if (-e0)%N==a: s-=roots[-beta*k%N]
    terms[i-1,v]=s%q
  for k in range(mult+1):
   row=np.zeros(len(aidx),np.int64)
   for l in range(k+1):row+=terms[aidx,l]*terms[bidx,k-l];
   if k:Ck=J
   else:Ck=C
   Ck.append((row%q).tolist()+extraform(c,b,h,k))
# Sparse linear combinations retain valid weight two forms and reduce the dimensions of the kernel computation.
mu=N*N
for P in sp.factorint(N):mu=mu*(P*P-1)//(P*P)
g=1+mu//24-len(C)//2; nb=len(C)+g+40; nc=len(C[0]); rng=np.random.default_rng(20261118);
if nc>nb:
 slots=rng.integers(0,nb,size=nc); coeff=rng.integers(1,q,size=nc,dtype=np.int64)
 def compress(rows):return [(np.bincount(slots,weights=(np.asarray(r)*coeff)%q,minlength=nb).astype(np.int64)%q).tolist() for r in rows]
 C=compress(C);J=compress(J)
A=matrix(GF(q),C); rank=A.rank();print('constant rank',rank,'/',len(C)-1,flush=True);assert rank==len(C)-1
K=A.right_kernel_matrix();J=matrix(GF(q),J)*K.transpose();pivs=J.pivots();print('jet rank',len(pivs),flush=True);V=J.matrix_from_columns(pivs);np.savez(base/f'eis{N}.npz',rows=np.asarray(V,dtype=np.int32),profile=np.array(group),prime=q)
