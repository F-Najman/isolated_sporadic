// Product evaluations on X_1(131) in characteristic two.
SetQuitOnError(true); SetColumns(0); SetSeed(131);
k:=GF(2); R<T>:=PolynomialRing(k); k131:=GF(131); U131<W>:=PolynomialRing(k131);
p:=131; pairs:=[<a,b>:a,b in [1..65]|a le b];
function pivots(B)
    A:=BasisMatrix(RowSpace(B));
    return [Min([j:j in [1..Ncols(A)]|A[i,j] ne 0]):i in [1..Nrows(A)]];
end function;
function choose(B,cols)
    return Submatrix(B,[1..Nrows(B)],cols);
end function;
function products(w,pp)
    return [w[z[1]]*w[z[2]]:z in pp];
end function;
function absoluteRows(vals)
    L:=Universe(vals);V,rep:=VectorSpace(L,k);
    return BasisMatrix(RowSpace(Transpose(Matrix(k,[Eltseq(rep(v)):v in vals]))));
end function;
function weights(m,c,unit)
    return [Universe(m)!unit] cat [m[c*a mod 131]+m[c]:a in [2..65]];
end function;
function splitCuspRows(length)
    L<t>:=LaurentSeriesRing(k,130); sl:=[L|];
    for a in [1..65] do
        w:=t^a; v:=t^(131-a);
        x:=w/(1+w)^2+v/(1+v)^2;
        y:=w^2/(1+w)^3+v/(1+v)^3;
        m:=x+y/x;
        assert AbsolutePrecision(m) gt length;Append(~sl,m);
    end for;
    sl cat:=[1+sl[131-a]:a in [66..130]];
    jets:=[[Coefficient(m,i):m in sl]:i in [0..length]];
    rows:=[];
    for i in [0..64] do
        c:=3^i mod 131;
        U:=[Matrix(k,1,65,weights(jets[j+1],c,j eq 0 select 1 else 0)):j in [0..length]];
        for n in [0..length] do
            C:=&+[Transpose(U[j+1])*U[n-j+1]:j in [0..n]];
            Append(~rows,Vector(k,[C[v[1],v[2]]:v in pairs]));
        end for;
    end for;
    return Matrix(k,[Eltseq(v):v in rows]);
end function;
precision:=23; assert 65*precision gt 2*651-2+130;
assert Order((k131!3)^2) eq 65;
stack:=splitCuspRows(precision-1);
cols:=pivots(stack); pp:=[pairs[i]:i in cols]; B:=choose(stack,cols);
base:=[i*precision+1:i in [0..64]];
L:=GF(2,130);z:=PrimitiveElement(L)^((#L-1) div 131);
m:=[(z^a/(1+z^a))^2:a in [1..130]];
cusps:=VerticalJoin(Submatrix(B,base,[1..#cols]),absoluteRows(products(weights(m,1,1),pp)));
canon:=Transpose(BasisMatrix(Nullspace(Transpose(cusps))));
assert #cols eq 779 and Rank(cusps) eq 129 and Ncols(canon) eq 650;
idx:=pivots(Transpose(B));next:=[(i-1+precision) mod Nrows(B)+1:i in [1..Nrows(B)]];
delta:=Submatrix(B,idx,[1..#cols])^-1*Submatrix(B,[next[j]:j in idx],[1..#cols]);
assert B*delta eq Submatrix(B,next,[1..#cols]) and delta^65 eq IdentityMatrix(k,#cols);
rot:=[delta^i*canon:i in [0..64]];
printf "forms %o cusp rank %o regular %o\n",#cols,Rank(cusps),Ncols(canon);
function genPoint(f,r,jcode,dual)
    K:=GF(2,f); pol:=R!Intseq(jcode,2);
    assert IsIrreducible(pol) and Degree(pol) eq f;
    j:=Roots(ChangeRing(pol,K))[1][1]; E:=EllipticCurve([K|1,0,0,0,1/j]);
    q:=#K; tr:=TraceOfFrobenius(E); eig:=[v[1]:v in Roots(W^2-tr*W+q)];
    eligible:=[z:z in eig|Order(z^2) eq r];
    assert #eligible eq 1 and #eig eq 2;
    lambda:=eligible[1]; other:=[v:v in eig|v ne lambda][1];
    KK:=GF(2,f*r); Embed(K,KK);
    three:=Factorization(DivisionPolynomial(BaseChange(E,KK),3));
    n:=Lcm([Order(lambda)] cat [r*Degree(v[1]):v in three]);
    LL:=GF(2,f*n); Embed(KK,LL); EE:=BaseChange(E,LL);
    t0:=2;t1:=tr;for i in [2..n] do t0,t1:=t1,tr*t1-q*t0;end for;
    cofactor:=q^n+1-t1; cofactor div:=131^Valuation(cofactor,131);
    repeat
        P:=cofactor*Random(EE);
        while P ne EE!0 and 131*P ne EE!0 do P:=131*P;end while;
        if P ne EE!0 then P:=EE![P[1]^q,P[2]^q]-(Integers()!other)*P;end if;
    until P ne EE!0;
    assert EE![P[1]^q,P[2]^q] eq (Integers()!lambda)*P;
    UU<X>:=PolynomialRing(LL); img:=[* *];
    for term in three do
        root:=Roots(UU!term[1])[1][1]; target,phi:=IsogenyFromKernel(EE,X-root);
        Append(~img,<phi(P),Degree(term[1])>);
    end for;
    assert &+[t[2]:t in img] eq 4;
    return P,img;
end function;
function ordinaryRows(P,c,jet)
    E:=Curve(P);a1,a2,a3,a4,a6:=Explode(aInvariants(E));
    assert a1 eq 1 and a3 eq 0;
    Q:=P; m:=[]; der:=[];
    for a in [1..65] do
        x:=Q[1];Append(~m,(x^2+a4+Q[2])/x);
        if jet then Append(~der,1/x^2);end if;
        Q+:=P;
    end for;
    m cat:=[1+m[131-a]:a in [66..130]];
    w:=weights(m,c,1); A:=absoluteRows(products(w,pp));
    if jet then
        der cat:=[der[131-a]:a in [66..130]];
        wd:=weights(der,c,0);
        D:=absoluteRows([w[v[1]]*wd[v[2]]+wd[v[1]]*w[v[2]]:v in pp]);
        A:=VerticalJoin(A,D);
    end if;
    return A;
end function;
function primaryBlock(f,r,j,jw)
    P,img:=genPoint(f,r,j,jw);
    A:=VerticalJoin(ordinaryRows(P,1,false),ordinaryRows(P,3,true));
    for entry in img do A:=VerticalJoin(A,ordinaryRows(entry[1],3,false));end for;
    return BasisMatrix(RowSpace(A));
end function;
