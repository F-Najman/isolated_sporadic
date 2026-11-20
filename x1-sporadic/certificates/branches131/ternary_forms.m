SetQuitOnError(true); SetColumns(0); SetSeed(131);
k:=GF(3); R<T>:=PolynomialRing(k);
// The original quotient is the default. The full-cover Hecke test uses H=1.
fullMarking:=assigned useFullCurve and useFullCurve;
H:=fullMarking select [1] else [Integers()!((GF(131)!2)^(26*i)):i in [0..4]];
deckSize:=65 div #H;assert #Seqset(H)*deckSize eq 65;
pairs:=[<a,b>:a,b in [1..65]|a le b];
multIndex:=func<i|Min(i mod 131,(-i) mod 131)>;
function symCoefficients(mat)
    return [mat[p[1],p[2]]:p in pairs];
end function;
function allProducts(slopes)
    L:=Universe(slopes);
    weights:=Matrix(L,[[slopes[a*h mod 131]:a in [1..65]]:h in H]);
    return symCoefficients(Transpose(weights)*weights);
end function;
function absoluteRows(values)
    L:=Universe(values); V,ev:=VectorSpace(L,k);
    return BasisMatrix(RowSpace(Transpose(Matrix(k,[Eltseq(ev(v)):v in values]))));
end function;
function splitSeries(length)
    S<u>:=LaurentSeriesRing(k,150);
    m:=[S|];
    for r in [1..65] do
        w:=u^r; v:=u^(131-r);
        x:=w/(1-w)^2 + v/(1-v)^2;
        y:=w^2/(1-w)^3 - v/(1-v)^3;
        z:=(x/2)/(2*y+x);
        assert AbsolutePrecision(z) gt length and Coefficient(z,0) eq 1/k!2;
        Append(~m,z);
    end for;
    m cat:=[-m[131-r]:r in [66..130]];
    rows:=[];
    for i in [0..deckSize-1] do
        c:=2^i mod 131;
        ws:=[Matrix(k,[[Coefficient(m[a*c*h mod 131],n):a in [1..65]]:h in H]):n in [0..length]];
        coeffs:=[];
        for n in [0..length] do
            mat:=&+[Transpose(ws[j+1])*ws[n-j+1]:j in [0..n]];
            Append(~coeffs,Vector(k,symCoefficients(mat)));
        end for;
        Append(~rows,coeffs);
    end for;
    return rows;
end function;
function fullConstants(raw)
    split:=Matrix(k,[Eltseq(v[1]):v in raw]);
    L:=GF(3,65); zet:=PrimitiveElement(L)^((#L-1) div 131); assert Order(zet) eq 131;
    m:=[(1-zet^a)/(2*(1+zet^a)):a in [1..130]];
    return BasisMatrix(RowSpace(VerticalJoin(split,absoluteRows(allProducts(m)))));
end function;
function selectCols(A,cols)
    return Matrix(k,[[A[i,j]:j in cols]:i in [1..Nrows(A)]]);
end function;
raw:=splitSeries(22);
constants:=fullConstants(raw);
stack:=VerticalJoin(constants,Matrix(k,[Eltseq(v):v in &cat raw]));
ech:=BasisMatrix(RowSpace(stack));
cols:=[Min([j:j in [1..Ncols(ech)]|ech[i,j] ne 0]):i in [1..Nrows(ech)]];
assert #Seqset(cols) eq #cols;
cc:=selectCols(constants,cols);
selected:=[pairs[j]:j in cols];
assert #cols eq 12*deckSize-1 and Rank(cc) eq 2*deckSize-1;
print "modular forms",#cols,"cusps",Rank(cc),"available cusp forms",#cols-Rank(cc);
function productValues(slopes,shift)
    L:=Universe(slopes);
    v:=[ [slopes[a*h*shift mod 131]:a in [1..65]]:h in H ];
    return [ &+[row[p[1]]*row[p[2]]:row in v] :p in selected];
end function;
function tangentSlopes(P)
    E:=Curve(P);
    a1,a2,a3,a4,a6:=Explode(aInvariants(E));
    assert a1 eq 0 and a2 eq 1 and a3 eq 0;
    out:=[]; Q:=P;
    for n in [1..65] do
        Append(~out,(2*Q[1]+a4)/(2*Q[2])); Q+:=P;
    end for;
    return out cat [-out[131-r]:r in [66..130]];
end function;
Q131<W>:=PolynomialRing(GF(131));
function marked(f,jcode,dualcode)
    K:=GF(3,f); q:=#K;
    jf:=R!Intseq(jcode,3); assert IsIrreducible(jf) and Degree(jf) eq f;
    j:=Roots(ChangeRing(jf,K))[1][1]; assert j ne 0;
    E:=EllipticCurve([K|0,1,0,0,-1/j]);
    tr:=TraceOfFrobenius(E);
    eig:=[v[1]:v in Roots(W^2-tr*W+q)];
    eligible:=[r:r in eig|r^(2*#H) eq 1]; assert #eligible gt 0;
    two:=Factorization(DivisionPolynomial(E,2));
    assert &+[Degree(z[1]):z in two] eq 3 and &and[z[2] eq 1:z in two];
    data:=[* *];
    for r in eligible do
        // Factor the two-kernel polynomial over its splitting extension.
        // Each irreducible kernel orbit contributes one evaluation block.
        n:=Lcm([Order(r)] cat [Degree(z[1]):z in two]);
        L:=GF(3,f*n); Embed(K,L); EL:=BaseChange(E,L);
        RR<X>:=PolynomialRing(L);
        t0:=2;t1:=tr;
        for m in [2..n] do t0,t1:=t1,tr*t1-q*t0;end for;
        card:=q^n+1-t1;
        repeat
            P:=(card div 131^Valuation(card,131))*Random(EL);
            while P ne EL!0 and 131*P ne EL!0 do P:=131*P;end while;
            if #eig eq 2 and P ne EL!0 then
              other:=[v:v in eig|v ne r][1];
              P:=EL![P[1]^q,P[2]^q]-(Integers()!other)*P;
            end if;
        until P ne EL!0;
        if #eig eq 1 then
            assert (tr^2-4*q) mod 131^2 ne 0;
        end if;
        assert EL![P[1]^q,P[2]^q] eq (Integers()!r)*P;
        if #eligible gt 1 then
            ker:=&*[X-(a*P)[1]:a in [1..65]];
            iso:=IsogenyFromKernel(EL,ker);
            if Evaluate(R!Intseq(dualcode,3),jInvariant(iso)) ne 0 then continue; end if;
        end if;
        neighbours:=[* *]; orig:=tangentSlopes(P);
        for fac in two do
            root:=Roots(RR!fac[1])[1][1];
            target,phi:=IsogenyFromKernel(EL,X-root);
            Append(~neighbours,<tangentSlopes(phi(P)),Degree(fac[1]),phi(P)>);
        end for;
        Append(~data,<orig,neighbours,P>);
    end for;
    assert #data gt 0; return data;
end function;
function mblock(comp,shift)
    slopes:=comp[1]; ngh:=comp[2];
    B:=VerticalJoin(absoluteRows(productValues(slopes,shift)),
                    absoluteRows(productValues(slopes,2*shift)));
    for v in ngh do B:=VerticalJoin(B,absoluteRows(productValues(v[1],2*shift)));end for;
    return BasisMatrix(RowSpace(B));
end function;
