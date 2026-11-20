load "forms.m";
// High-order characteristic-two Tate slopes at the split cusps. The constant
// terms in y and a4 cancel since s_1(q)=s_3(q) in characteristic two.
L<t>:=LaurentSeriesRing(k,450); sl:=[L|];
for a in [1..65] do
    wseq:=[t^(131*n+a):n in [0..3]];
    vseq:=[t^(131*n-a):n in [1..3]];
    x:=&+[z/(1+z)^2:z in wseq cat vseq];
    ym:=&+[z^2/(1+z)^3:z in wseq]+&+[z/(1+z)^3:z in vseq];
    s:=x+ym/x;
    assert AbsolutePrecision(s) gt 129;Append(~sl,s);
end for;
sl cat:=[1+sl[131-a]:a in [66..130]];
weights0:=[L!1] cat [sl[a]+sl[1]:a in [2..65]];
products0:=products(weights0,pp);
A:=Matrix(k,[[Coefficient(z,n):z in products0]:n in [1..129]]);
assert Submatrix(A,1,1,22,Ncols(A)) eq Submatrix(B,2,1,22,Ncols(B));
anti:=(IdentityMatrix(k,#cols)-delta)*canon;
J:=[];T0:=A;
for i in [0..64] do Append(~J,T0*anti); T0:=T0*delta;end for;
function first(B,n)
    return Submatrix(B,1,1,n,Ncols(B));
end function;
// Largest, second and third multiplicities among at most three split cusps.
J129:=J[1]; J64:=[first(M,64):M in J]; J43:=[first(M,43):M in J];
bad:=0;count:=0;
for i in [2..65] do
    T:=VerticalJoin(J129,J64[i]);
    for j in [2..65] do
        if i eq j then continue;end if;
        count+:=1;
        if Rank(VerticalJoin(T,J43[j])) lt 129+64+43 then bad+:=1;end if;
    end for;
end for;
printf "three split cusps: rectangles %o failures %o\n",count,bad;
// The first differential rows at the unique degree-65 cusp above c_infty.
L:=GF(2,130);z:=PrimitiveElement(L)^((#L-1) div 131);
sl:=[L|];ds:=[L|];
for a in [1..130] do
    w:=z^a; x:=w/(1+w)^2;y:=w^2/(1+w)^3;
    x1:=w+1/w;y1:=1+1/w;
    Append(~sl,x+y/x);Append(~ds,x1+(y1+1)/x+y*x1/x^2);
end for;
u:=weights(sl,1,1);du:=weights(ds,1,0);
I:=absoluteRows([u[v[1]]*du[v[2]]+du[v[1]]*u[v[2]]:v in pp])*anti;
I:=BasisMatrix(RowSpace(I));assert Rank(I) eq 64;
J64:=first(J[1],64); J32:=[first(M,32):M in J];J21:=[first(M,21):M in J];
badI:=0;countI:=0;
for i in [2..65] do
    T:=VerticalJoin(VerticalJoin(I,J64),J32[i]);
    for j in [2..65] do
        if i eq j then continue;end if;
        countI+:=1;
        if Rank(VerticalJoin(T,J21[j])) lt 64+64+32+21 then badI+:=1;end if;
    end for;
end for;
printf "inert orbit and three split cusps: rectangles %o failures %o\n",countI,badI;
assert count eq 4032 and bad eq 0 and countI eq 4032 and badI eq 0;
function initialBlock(J,n)
    M:=first(J[1],n);
    for i in [2..#J] do M:=VerticalJoin(M,first(J[i],n));end for;
    return M;
end function;
function extraRectangles(J,initial,top)
    m:=6;
    N:=Transpose(BasisMatrix(Nullspace(Transpose(initial))));
    U:=[Submatrix(M,m+1,1,top[1]-m,Ncols(M))*N:M in J];
    V:=[first(M,top[2]-m):M in U]; W:=[first(M,top[3]-m):M in U];
    count:=0;failed:=0;target:=&+top-3*m;
    for i in [2..65] do
        T:=VerticalJoin(U[1],V[i]);
        for j in [2..65] do
            if i eq j then continue;end if;
            count+:=1;
            if Rank(VerticalJoin(T,W[j])) ne target then failed+:=1;end if;
        end for;
    end for;
    return count,failed;
end function;
A6:=initialBlock(J,6); r6:=Rank(A6);
A6I:=VerticalJoin(A6,I);r6I:=Rank(A6I);
printf "six jets on all split cusps: rank %o with inert differences %o\n",r6,r6I;
count,bad:=extraRectangles(J,A6,[129,64,43]);
countI,badI:=extraRectangles(J,A6I,[64,32,21]);
printf "six jets and three higher multiplicities: rectangles %o failures %o\n",count,bad;
printf "six jets, inert orbit and three higher multiplicities: rectangles %o failures %o\n",countI,badI;
assert r6 eq 384 and r6I eq 448 and bad eq 0 and badI eq 0;
for m in [7..9] do
    A0:=initialBlock(J,m);
    assert Rank(VerticalJoin(A0,I)) eq 64*(m+1);
    printf "all split jets through %o and inert differences: rank %o\n",m,Rank(VerticalJoin(A0,I));
end for;
quit;
