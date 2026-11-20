SetQuitOnError(true);SetColumns(0);SetSeed(13195);
k:=GF(3);S<t>:=LaurentSeriesRing(k,320);
fd:=Open(file,"r");J:=[];
for i in [0..12] do
    V:=[];
    for n in [1..260] do Append(~V,[StringToInteger(a):a in Eltseq(Gets(fd))]);end for;
    Append(~J,Matrix(k,V));
end for;
function combi(M,a) return Eltseq(a*Transpose(M));end function;
function productMap(a)
    A:=ZeroMatrix(k,0,131);
    for M in J do
        U:=Submatrix(M,1,1,42,131);g:=combi(U,a);
        T:=Matrix(k,42,42,[[(i-j+1 ge 1) select g[i-j+1] else 0:j in [1..42]]:i in [1..42]]);
        A:=VerticalJoin(A,T*U);
    end for;
    return A;
end function;
function series(M,a)
    v:=combi(M,a);
    return &+[S|v[i]*t^i:i in [1..#v]]+O(t^(#v+1));
end function;
function functions(D)
    A:=ZeroMatrix(k,0,131);
    for i in [1..13] do if D[i] gt 0 then A:=VerticalJoin(A,Submatrix(J[i],1,1,D[i],131));end if;end for;
    W:=Nullspace(Transpose(A));ell:=1+&+D-Rank(A);
    if ell eq 1 then return [],ell;end if;
    // A pair with no additional common base point gives the known dimension.
    repeat
        a:=Random(W);b:=Random(W);
        M:=HorizontalJoin(productMap(b),-productMap(a));
        N:=Nullspace(Transpose(M));
    until Dimension(N) eq ell;
    H:=[series(J[i],a):i in [1..13]];
    fs:=[];
    for h in Basis(N) do
        v:=Vector(k,Eltseq(h)[1..131]);
        Append(~fs,[series(J[i],v)/H[i]:i in [1..13]]);
    end for;
    assert &and[AbsolutePrecision(f[i]) ge 115:f in fs,i in [1..13]];
    return fs,ell;
end function;
src:=Split(Read(patterns),"\n"); out:=Open(result,"w"); id:=0;
for line in src do
   v:=[StringToInteger(a):a in Split(line," ")];ix:=v[1];D:=v[2..14];
   fs,ell:=functions(D);printf "cuspidal representative %o degree %o ell %o\n",ix,&+D,ell;
   if ell eq 1 then continue;end if;
   for h in Points(ProjectiveSpace(k,ell-1)) do
       a:=Eltseq(h);
       gs:=[&+[a[j]*fs[j][i]:j in [1..ell]]:i in [1..13]];
       // A nonzero section cannot vanish identically as a rational function.
       z:=[Valuation(g):g in gs];d:=-&+z;
       assert d ge 0 and &and[D[i]+z[i] ge 0:i in [1..13]];
       printf "section %o %o residual %o cusps %o\n",ix,id,d,[D[i]+z[i]:i in [1..13]];
       if d gt 0 then
          u:=[gs[i]*t^(-z[i]):i in [1..13]];
          assert &and[AbsolutePrecision(a) gt d:a in u];rr:=&*u;
          fprintf out,"%o %o",id,d;
          for j in [0..d] do fprintf out," %o",Integers()!Coefficient(rr,j);end for;
          fprintf out,"\n";
       end if;
       id+:=1;
   end for;
end for;
quit;
