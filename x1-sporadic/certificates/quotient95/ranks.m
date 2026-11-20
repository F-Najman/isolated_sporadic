SetQuitOnError(true);SetColumns(0);k:=GF(3);
f:=Open(file,"r");rows:=[];
for i in [0..12] do
    a:=[];
    for n in [1..100] do Append(~a,[StringToInteger(x):x in Eltseq(Gets(f))]);end for;Append(~rows,Matrix(k,a));
end for;
polars:=Split(Read(patterns),"\n");
for p in polars do
    v:=[StringToInteger(s):s in Split(p," ")];b:=v[2..14];
    M:=ZeroMatrix(k,0,131);
    for i in [1..13] do if b[i] gt 0 then M:=VerticalJoin(M,Submatrix(rows[i],1,1,b[i],131));end if;end for;
    r:=Rank(Submatrix(M,1,1,Nrows(M),130));printf "fibre %o degree %o trace-zero rank %o full rank %o\n",v[1],&+b,r,Rank(M);
end for;
quit;
