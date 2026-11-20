// Modular forms, ordinary marked curves, and split cusp expansions on Y.
SetQuitOnError(true);
SetColumns(0);
SetSeed(131);
k := GF(2);
R<T> := PolynomialRing(k);
H := [Integers()!((GF(131)!2)^(26*i)):i in [0..4]];
assert #Seqset(H) eq 5 and Order((GF(131)!3)^10) eq 13;
width := 65*66 div 2;

// The columns are traces of the products u_a*u_b, 1 <= a <= b <= 65.
// Each input row gives the coefficient of a successive power of the local
// parameter in the tangent slopes m_a. Expand field elements over F_2.
function conditions(jets)
    L := Universe(jets[1]);
    weights := [Matrix(L,[[d eq 1 select L!1 else L!0] cat
                  [jets[d][a*h mod 131]+jets[d][h] : a in [2..65]]
                  : h in H]) : d in [1..#jets]];
    blocks := ZeroMatrix(k,0,width);
    for n in [0..#weights-1] do
        coeff := &+[Transpose(weights[d+1])*weights[n-d+1]:d in [0..n]];
        vals := [coeff[a,b]:a,b in [1..65] | a le b];
        mat := Matrix(k,[Eltseq(v):v in vals]);
        blocks := VerticalJoin(blocks,BasisMatrix(RowSpace(Transpose(mat))));
    end for;
    return blocks;
end function;

function cuspConditions()
    Hpm := Seqset(H cat [131-h:h in H]);
    todo := {1..130};
    blocks := ZeroMatrix(k,0,width);
    while #todo ne 0 do
        c := Min(todo);
        todo diff:= {c*h mod 131:h in Hpm};
        slopes := [k | (a*c mod 131 gt 65) select 1 else 0:a in [1..130]];
        blocks := VerticalJoin(blocks,conditions([slopes]));
    end while;
    L := GF(2,130);
    zeta := PrimitiveElement(L)^((#L-1) div 131);
    assert Order(zeta) eq 131 and Order((GF(131)!2)^10) eq 13;
    slopes := [(zeta^a/(1+zeta^a))^2:a in [1..130]];
    return BasisMatrix(RowSpace(VerticalJoin(blocks,conditions([slopes]))));
end function;

function curveConditions(p,firstJet)
    e := Curve(p);
    a1,a2,a3,a4,a6 := Explode(aInvariants(e));
    assert a1 eq 1 and a2 eq 0 and a3 eq 0 and p ne e!0 and 131*p eq e!0;
    mults := [a*p:a in [1..130]];
    slopes := [(pt[1]^2+a4+pt[2])/pt[1]:pt in mults];
    assert &and[slopes[a]+slopes[131-a] eq 1:a in [1..65]];
    if firstJet then
        assert a4 eq 0 and a6 eq 1/jInvariant(e);
        first := [1/pt[1]^2:pt in mults];
        return conditions([slopes,first]);
    end if;
    return conditions([slopes]);
end function;

R131<W> := PolynomialRing(GF(131));
function cuspSeries(length)
L<t> := LaurentSeriesRing(k,130);
slopes := [L | ];
for r in [1..65] do
    w := t^r;
    v := t^(131-r);
    x := w/(1+w)^2 + v/(1+v)^2;
    y := w^2/(1+w)^3 + v/(1+v)^3;
    m := x+y/x;
    assert AbsolutePrecision(m) gt length;
    Append(~slopes,m);
end for;
slopes cat:= [1+slopes[131-r]:r in [66..130]];
rows := [];
for j in [0..12] do
    c := 3^j mod 131;
    weights := [Matrix(k,[[n eq 0 select 1 else 0] cat
          [Coefficient(slopes[a*c*h mod 131]+slopes[c*h mod 131],n):a in [2..65]]
             :h in H]) : n in [0..length]];
    jets := [];
    for n in [1..length] do
        C := &+[Transpose(weights[i+1])*weights[n-i+1] : i in [0..n]];
        Append(~jets,Vector(k,[C[a,b]: a,b in [1..65] | a le b]));
    end for;
    Append(~rows,jets);
end for;

    return rows;
end function;
function markedPoints(b,code,dual)
    K := GF(2,b);
    poly := R!Intseq(code,2);
    assert IsIrreducible(poly) and Degree(poly) eq b;
    q := #K;
    j := Roots(ChangeRing(poly,K))[1][1];
    curve := EllipticCurve([K|1,0,0,0,1/j]);
    tr := TraceOfFrobenius(curve);
    eig := Roots(W^2-tr*W+q);
    assert #eig ge 1;
    good := [r[1]:r in eig | r[1]^10 eq 1];
    assert #good ge 1;
    result:=[* *];
    for r in good do
    three := Factorization(DivisionPolynomial(curve,3));
    assert &and[f[2] eq 1:f in three] and &+[Degree(f[1]):f in three] eq 4;
    n := Lcm([Order(r)] cat [Degree(f[1]):f in three]);
    L := GF(2,b*n);
    Embed(K,L);
    E := BaseChange(curve,L);
    U<X> := PolynomialRing(L);
    t0,t1 := Explode([2,tr]);
    for m in [2..n] do t0,t1 := t1,tr*t1-q*t0; end for;
    size := q^n+1-t1;
    assert Valuation(size,131) ge 1;
    repeat
        P := (size div 131^Valuation(size,131))*Random(E);
        while P ne E!0 and 131*P ne E!0 do P := 131*P; end while;
        if #eig eq 2 and P ne E!0 then
            other:=[v[1]:v in eig|v[1] ne r][1];
            P:=E![P[1]^q,P[2]^q]-(Integers()!other)*P;
        end if;
    until P ne E!0;
    assert E![P[1]^q,P[2]^q] eq (Integers()!r)*P;
    if #eig eq 1 then assert (tr^2-4*q) mod 131^2 ne 0;end if;
    if #good gt 1 then
        ker:=&*[X-(a*P)[1]:a in [1..65]];
        quotient:=IsogenyFromKernel(E,ker);
        if Evaluate(R!Intseq(dual,2),jInvariant(quotient)) ne 0 then continue;end if;
    end if;
    images := [* *];
    for f in three do
        roots := Roots(U!f[1]);
        target,phi := IsogenyFromKernel(E,X-roots[1][1]);
        Append(~images,<phi(P),Degree(f[1])>);
    end for;
    Append(~result,<P,images>);
    end for;
    assert #result gt 0;return result;
end function;
function chooseCols(A,cols)
    return Matrix(k,[[A[i,j]:j in cols]:i in [1..Nrows(A)]]);
end function;
function pivots(B)
    row:=BasisMatrix(RowSpace(B));
    return [Min([j:j in [1..Ncols(row)]|row[i,j] ne 0]):i in [1..Nrows(row)]];
end function;
raw:=cuspSeries(26); allc:=cuspConditions();
stack:=VerticalJoin(allc,Matrix(k,[Eltseq(v):v in &cat raw]));
cols:=pivots(stack);
cc:=chooseCols(allc,cols);
canon:=Transpose(BasisMatrix(Nullspace(Transpose(cc))));
assert #cols eq 155 and Rank(cc) eq 25;
pairs:=[<a,b>:a,b in [1..65]|a le b];
selectPairs:=[pairs[j]:j in cols];
function valsAt(slopes,v)
    L:=Universe(slopes); ww:=[[L!1] cat [slopes[a*h*v mod 131]+slopes[h*v mod 131]:a in [2..65]]:h in H];
    return [&+[w[P[1]]*w[P[2]]:w in ww]:P in selectPairs];
end function;
function slopesAt(P)
    E:=Curve(P); a1,a2,a3,a4,a6:=Explode(aInvariants(E)); assert a1 eq 1 and a2 eq 0 and a3 eq 0;
    Q:=P; m:=[];
    for i in [1..65] do Append(~m,(Q[1]^2+a4+Q[2])/Q[1]); Q+:=P;end for;
    return m cat [1+m[131-i]:i in [66..130]];
end function;
function absRows(values)
    L:=Universe(values);V,f:=VectorSpace(L,k);
    return BasisMatrix(RowSpace(Transpose(Matrix(k,[Eltseq(f(v)):v in values]))));
end function;
function evalAt(P,v)
    return absRows(valsAt(slopesAt(P),v));
end function;
function cuspCoords()
    at0:=chooseCols(Matrix(k,[Eltseq(c[1]):c in raw]),cols)*canon;
    L:=GF(2,130); z:=PrimitiveElement(L)^((#L-1) div 131);
    m:=[L|]; jets:=[L|];
    for a in [1..130] do
        w:=z^a; x0:=w/(1+w)^2;y0:=w^2/(1+w)^3;
        x1:=w+1/w;y1:=1+1/w;
        Append(~m,x0+y0/x0);
        Append(~jets,x1+(y1+1)/x0+y0*x1/x0^2);
    end for;
    m0:=conditions([m]);m01:=conditions([m,jets]);
    // Only the first coefficient matters on cusp forms. The constants
    // are already zero after multiplication by canon.
    atI:=BasisMatrix(RowSpace(chooseCols(m01,cols)*canon));
    assert Rank(at0) eq 13 and Rank(atI) eq 13;
    assert Rank(VerticalJoin(at0,atI)) eq 25;
    return VerticalJoin(at0,atI),atI;
end function;
fibreO,atI:=cuspCoords();
// Passing from P to 3P permutes split cusp expansions. All coefficients
// beyond the constants transform by this cyclic permutation (the Tate q is fixed).
// At the nonsplit cusp the effect is the same field automorphism on all
// coefficient values, so compute translated raw columns algebraically.
function translatedProducts(full,v)
    // u_a(vP)=u_{av}(P)+u_v(P) in characteristic two; unfolding a
    // negative residue contributes u_1 (the Hasse invariant).
    if full eq 1 then return [1];end if;
    r:=full*v mod 131; vs:=v mod 131; terms:=[];
    if (r gt 65) xor (vs gt 65) then Append(~terms,1);end if;
    ar:=Min(r,131-r); av:=Min(vs,131-vs);
    if ar ne 1 then Append(~terms,ar);end if;
    if av ne 1 then Append(~terms,av);end if;
    return terms;
end function;
tr:=ZeroMatrix(k,Nrows(stack),#cols);
for l in [1..#selectPairs] do
    P:=selectPairs[l];
    for a in translatedProducts(P[1],3), b in translatedProducts(P[2],3) do
        pos:=Index(pairs,<Min(a,b),Max(a,b)>);
        for i in [1..Nrows(stack)] do tr[i,l]+:=stack[i,pos];end for;
    end for;
end for;
basis:=chooseCols(stack,cols);
idx:=pivots(Transpose(basis)); delta:=Matrix(k,[Eltseq(basis[i]):i in idx])^-1*Matrix(k,[Eltseq(tr[i]):i in idx]);
assert basis*delta eq tr and delta^13 eq IdentityMatrix(k,#cols);
rot:=[delta^i*canon:i in [0..12]];
function block3(data)
    P,imgs:=Explode(data);
    B:=ZeroMatrix(k,0,#cols);
    for v in imgs do B:=VerticalJoin(B,evalAt(v[1],1));end for;
    return BasisMatrix(RowSpace(B));
end function;
