// Degree-13 residue classes on Y = X_1(131)/H_5 in characteristic two.
// Check l(s*T_3*D + D + 2*s*D) = 1, where s acts on markings by 3.
// Run with magma -b hecke13_131.m. No input files are needed.
SetQuitOnError(true);
SetColumns(0);
SetSeed(131);
k := GF(2);
R<T> := PolynomialRing(k);
polys := [T^13+T^11+T^10+T^9+T^8+T^2+1,
          T^13+T^12+T^9+T^6+T^5+T^2+1];
K := GF(2,13);
q := #K;
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

// Over F_(q^20) all 131-torsion is rational and the 131-primary part
// has order 131^2. Project onto the desired Frobenius eigenspace.
function eigenPoint(e,r,other,cofactor)
    repeat
        t := cofactor*Random(e);
        if t ne e!0 then
            assert 131*t eq e!0;
            t := e![t[1]^q,t[2]^q]-other*t;
        end if;
    until t ne e!0;
    assert 131*t eq e!0 and e![t[1]^q,t[2]^q] eq r*t;
    return t;
end function;

cusps := cuspConditions();
assert Rank(cusps) eq 25;
L := GF(2,260);
Embed(K,L);
U<X> := PolynomialRing(L);
for i in [1,2] do
    assert IsIrreducible(polys[i]);
    j := Roots(ChangeRing(polys[i],K))[1][1];
    e0 := EllipticCurve([K|1,0,0,0,1/j]);
    tr := TraceOfFrobenius(e0);
    assert tr eq -71;
    R131<W> := PolynomialRing(GF(131));
    eig := Sort([Integers()!v[1]:v in Roots(W^2-tr*W+q)]);
    assert eig eq [61,130] and &and[(GF(131)!r)^10 eq 1:r in eig];
    t0 := 2;
    t1 := tr;
    for n in [2..20] do t0,t1 := t1,tr*t1-q*t0; end for;
    order := q^20+1-t1;
    assert Valuation(order,131) eq 2;
    e := BaseChange(e0,L);
    three := DivisionPolynomial(e0,3);
    assert IsIrreducible(three) and Degree(three) eq 4;
    roots3 := Roots(U!three);
    assert #roots3 eq 4;
    target, phi := IsogenyFromKernel(e,X-roots3[1][1]);
    assert Degree(MinimalPolynomial(jInvariant(target))) eq 52;
    for n in [1,2] do
        p := eigenPoint(e,eig[n],eig[3-n],order div 131^2);
        kernel := &*[X-(a*p)[1]:a in [1..65]];
        quotient := IsogenyFromKernel(e,kernel);
        assert Evaluate(polys[3-i],jInvariant(quotient)) eq 0;
        simple := curveConditions(p,false);
        double := curveConditions(3*p,true);
        neighbours := curveConditions(phi(3*p),false);
        stacked := VerticalJoin(VerticalJoin(cusps,simple),VerticalJoin(double,neighbours));
        ranks := [Rank(cusps),Rank(simple),Rank(double),Rank(neighbours),Rank(stacked)];
        assert ranks eq [25,13,26,52,116];
        printf "j-polynomial %o, eigenvalue %o: ranks %o, canonical rank increase %o\n",
               i,eig[n],ranks,Rank(stacked)-Rank(cusps);
    end for;
end for;
quit;
