// Run in certificates/ with magma -b fibres131.m.
SetQuitOnError(true);
SetColumns(0);
load "fricke131-forms.m";
// Put z=s+2u in the integral equation z^2=R. Coefficients are reduced
// in the basis 1,y using y^2+y=x^3-x^2+x before division by 2 or 4.
IZ<t> := PolynomialRing(Integers());
v := t^3-t^2+t;
s0 := t^5+t^4+t+1;
s1 := t^2+1;
function divide(P,n)
    assert &and[c mod n eq 0 : c in Coefficients(P)];
    return Parent(P)![c div n : c in Coefficients(P)];
end function;
H0 := divide(R(t,0)-s0^2-s1^2*v,4);
H1 := divide(R(t,1)-R(t,0)-2*s0*s1+s1^2,4);
jm0 := Jm(t,0);
jm1 := Jm(t,1)-jm0;
jj0 := divide(Jp(t,0)+s0*jm0+s1*jm1*v,2);
jj1 := divide(Jp(t,1)-Jp(t,0)+s0*jm1+s1*jm0-s1*jm1,2);

k := GF(2);
K<X> := FunctionField(k);
S<U> := PolynomialRing(K);
A := U^2+Evaluate(s0,X)*U-Evaluate(H0,X);
B := Evaluate(H1,X)-Evaluate(s1,X)*U;
F<u> := ext<K | A^2+A*B-B^2*(X^3-X^2+X)>;
x := F!X;
y := Evaluate(A,u)/Evaluate(B,u);
assert y^2+y eq x^3-x^2+x and Genus(F) eq 11;
jfun := Evaluate(jj0,x)+y*Evaluate(jj1,x)+u*Jm(x,y);
cusps := Poles(x);
assert #cusps eq 2 and Seqset(Poles(jfun)) eq Seqset(cusps);
c0 := [P:P in cusps | Valuation(jfun,P) eq -131][1];
ci := [P:P in cusps | Valuation(jfun,P) eq -1][1];
assert Degree(c0) eq 1 and Degree(ci) eq 1;
E := EllipticCurve([k|0,-1,1,1,0]);
G := E![0,0];
assert Order(G) eq 3;
P131<W> := PolynomialRing(GF(131));
function minimumResidueDegree(P)
    if P eq c0 then return 1; end if;
    if P eq ci then return Order((GF(131)!2)^10); end if;
    L,ev := ResidueClassField(P);
    j0 := ev(jfun);
    curve := EllipticCurveWithjInvariant(L!j0);
    curves := j0 eq 0 select Twists(curve) else [curve];
    orders := {};
    for e in curves do
        roots := Roots(W^2-TraceOfFrobenius(e)*W+2^Degree(P));
        orders join:= {Order(r[1]^10):r in roots};
    end for;
    assert #orders gt 0 and orders subset {1,13};
    return Min(orders);
end function;
cache := AssociativeArray();
for b in [11..14] do
    total := 0;
    cuspOnly := 0;
    irreducible := 0;
    withCusp := 0;
    degree13classes := [];
    for h in [-26,-13,0,13,26], n in [0..2] do
        Q := n*G;
        if n eq 0 then
            fibre := c0+ci;
        else
            places := [P:P in Zeros(x-Q[1]) | Evaluate(y,P) eq Q[2]];
            fibre := &+[Valuation(x-Q[1],P)*P:P in places];
        end if;
        assert Degree(fibre) eq 2;
        D := b*c0+h*(c0-ci)+fibre-(c0+ci);
        V, section := RiemannRochSpace(D);
        if b eq 13 and h in {-13,0} and n eq 2 then
            assert Dimension(V) eq 3;
        end if;
        count := 0;
        for pivot in [1..Dimension(V)] do
            for tail in VectorSpace(k,Dimension(V)-pivot) do
                vec := V!([k!0:i in [1..pivot-1]] cat [k!1] cat Eltseq(tail));
                eff := D+Divisor(section(vec));
                count +:= 1;
                success := true;
                for P in Support(eff) do
                    if not IsDefined(cache,P) then cache[P] := minimumResidueDegree(P); end if;
                    if Valuation(eff,P) mod cache[P] ne 0 then success := false; break; end if;
                end for;
                if success then
                    support := Support(eff);
                    if Seqset(support) subset Seqset(cusps) then
                        assert eff eq b*c0 or (b ge 13 and eff eq (b-13)*c0+13*ci);
                        cuspOnly +:= 1;
                    elif #support eq 1 and Degree(support[1]) eq b then
                        irreducible +:= 1;
                        if b eq 13 then
                            Append(~degree13classes,<h,n>);
                            residue,ev := ResidueClassField(support[1]);
                            Rj<T> := PolynomialRing(k);
                            first := T^13+T^11+T^10+T^9+T^8+T^2+1;
                            second := T^13+T^12+T^9+T^6+T^5+T^2+1;
                            assert Rj!MinimalPolynomial(ev(jfun)) eq
                                (h eq -13 select first else second);
                            assert Rj!MinimalPolynomial(ev(Jp(x,y)-jfun)) eq
                                (h eq -13 select second else first);
                        end if;
                    else
                        assert b eq 14 and Valuation(eff,c0) eq 1 and
                            #support eq 2 and exists{Q:Q in support | Degree(Q) eq 13};
                        withCusp +:= 1;
                    end if;
                end if;
            end for;
        end for;
        assert count eq 2^Dimension(V)-1;
        total +:= count;
    end for;
    assert [total,cuspOnly,irreducible,withCusp] eq
           [[23,1,0,0],[57,1,0,0],[129,2,2,0],[241,2,5,2]][b-10];
    if b eq 13 then assert degree13classes eq [<-13,2>,<0,2>]; end if;
    printf "degree %o: %o tested, %o cuspidal, %o noncuspidal places, %o mixed\n",
           b,total,cuspOnly,irreducible,withCusp;
end for;
quit;
