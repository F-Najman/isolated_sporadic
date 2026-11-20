from pathlib import Path
import ast,math
from sage.all import matrix,ZZ,QQ
n=83;d=64;B=matrix(ZZ,ast.literal_eval(Path(__file__).with_name('valuations83.txt').read_text()));B=B.matrix_from_columns([0]).left_kernel_matrix()*B;assert not any(B.column(0));B=B.matrix_from_columns(range(1,42));print(B.dimensions(),flush=True)
L=B.BKZ(block_size=22);assert L.row_module()==B.row_module();G,M=L.change_ring(QQ).gram_schmidt();assert M*G==L
N,m=L.nrows(),L.ncols();S=1<<15;V=S;U=1<<12;rho=[g.dot_product(g) for g in G];P=[];done=set();units=range(1,83);fold=lambda j:min(j%n,(-j)%n)-1
for a in units:
 perm=[fold(a*(j+1)) for j in range(m)];assert B.row_module()==B.matrix_from_columns(perm).row_module()
for i in range(m):
 for j in range(i+1,m):
  if (i,j) in done:continue
  P.append((i,j));done|={tuple(sorted((fold((i+1)*a),fold((j+1)*a)))) for a in units}
assert len(done)==m*(m-1)//2;print('pairs',len(P),'min LLL L1',min(sum(map(abs,row)) for row in L),flush=True)
with open('sphere83.in','w') as f:
 print(N,m,d,S,V,U,len(P),file=f)
 for row in L:print(*row,file=f)
 for row in M:print(*[int((v*S).round()) for v in row],file=f)
 print(*[int((r*U).floor()) for r in rho],file=f)
 for i,j in P:print(i,j,*[int((QQ(d)*(g[i]-g[j])*S/(2*r)).round()) for g,r in zip(G,rho)],file=f)
