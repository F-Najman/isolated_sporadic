// Gonality of X_1(36) over F_5, used for the arithmetic bound at level 180.
// The plane equation is Sutherland's X1opt36new.txt, reduced modulo 5.
SetQuitOnError(true);
k := GF(5);
K<u> := FunctionField(k);
R<v> := PolynomialRing(K);
model := u^4*(u-1)*(u^2-u+1)*v^8
    + 4*u^4*(u-1)*(u+1)*(u^2-u+1)*v^7
    + u^2*(6*u^7+2*u^6-10*u^5+10*u^4+u^3-13*u^2+4*u-1)*v^6
    + 2*u^2*(2*u^8+5*u^7-5*u^6+12*u^4-18*u^3+u^2+u-1)*v^5
    + u^3*(u^8+8*u^7+u^6-9*u^5+27*u^4-9*u^3-55*u^2+34*u-13)*v^4
    + 2*u*(u^10+2*u^9-2*u^8+4*u^7+13*u^6-27*u^5-12*u^4+25*u^3-19*u^2+6*u-1)*v^3
    + (u^11+15*u^8-3*u^7-43*u^6+9*u^5+29*u^4-39*u^3+22*u^2-7*u+1)*v^2
    + 2*u^3*(u^3-3*u+1)*(u^3+3*u^2-2*u+1)*v + u^6*(u^3-3*u+1);
assert IsIrreducible(model);
F := FunctionField(model);
g := Genus(F);
assert g eq 17 and Degree(F) eq 8;

// Identify the reduced model over F_5 by its marked curve and j-degree.
U := F!u;
V := F.1;
G := (U-1)/U;
X := (U+V)/(U-1);
Y := (G^5*X^2-G^4*X^4-4*G^4*X^2-G^4+2*G^3*X^4+8*G^3*X^2+2*G^3*X+2*G^3
      -12*G^2*X^2-2*G^2*X-2*G^2+6*G*X^2+2*G*X+G-2*X^2)
  / (G^5*X-G^4*X^3-G^4*X^2-3*G^4*X-G^4-G^3*X^4+5*G^3*X^3+G^3*X^2
      +7*G^3*X+2*G^3+G^2*X^4-3*G^2*X^3-5*G^2*X^2-7*G^2*X-2*G^2
      -G*X^4+3*G*X^3+G*X^2+5*G*X+G-X^2-X);
r := (X^2*Y-X*Y+Y-1)/(X^2*Y-X);
s := (X*Y-Y+1)/(X*Y);
E := EllipticCurve([1+s-r*s,r*s-r^2*s,r*s-r^2*s,0,0]);
P := E![0,0];
assert 36*P eq E!0 and 18*P ne E!0 and 12*P ne E!0;
placesJ, ordersJ := Support(Divisor(jInvariant(E)));
assert -&+[Degree(placesJ[i])*Min(0,ordersJ[i]) : i in [1..#ordersJ]] eq 432;

places := [Places(F,d) : d in [1..4]];
assert [#s : s in places] eq [15,5,49,144];
canonical := CanonicalDivisor(F);
space, inclusion := RiemannRochSpace(canonical);
basis := [inclusion(v) : v in Basis(space)];
assert #basis eq g;

function Coordinates(value, degree)
    coordinates := ChangeUniverse(Eltseq(value,k),k);
    return coordinates cat [k!0 : i in [#coordinates+1..degree]];
end function;

// Row spaces of the first m jet conditions at a closed point. In the
// representation L(K), a section has local expansion starting at -v_P(K).
function Flags(P, m)
    d := Degree(P);
    vK := Valuation(canonical,P);
    if vK eq 0 and m eq 1 then
        entries := &cat[Coordinates(Evaluate(f,P),d) : f in basis];
        return [Transpose(Matrix(k,g,d,entries))];
    end if;
    series := [Expand(f,P : AbsPrec:=Max(2,m-vK+1)) : f in basis];
    assert &and[AbsolutePrecision(s) ge m-vK : s in series];
    result := [* *];
    rows := ZeroMatrix(k,0,g);
    for j in [0..m-1] do
        entries := &cat[Coordinates(Coefficient(s,j-vK),d) : s in series];
        rows := VerticalJoin(rows,Transpose(Matrix(k,g,d,entries)));
        Append(~result,rows);
    end for;
    return result;
end function;

rationalFlags := [Flags(P,5) : P in places[1]];
quadraticFlags := [Flags(P,2) : P in places[2]];
extras := AssociativeArray();
empty := ZeroMatrix(k,0,g);
extras[0] := [empty];
extras[2] := [s[1] : s in quadraticFlags];
extras[3] := [Flags(P,1)[1] : P in places[3]];
extras[4] := [Flags(P,1)[1] : P in places[4]]
           cat [s[2] : s in quadraticFlags]
           cat [VerticalJoin(quadraticFlags[i][1],quadraticFlags[j][1])
                : i,j in [1..#quadraticFlags] | i lt j];

procedure TestRationalPart(start, remaining, used, rows, extra, ~number)
    if remaining eq 0 then
        if used ge 3 then
            for B in extra do
                assert Rank(VerticalJoin(rows,B)) eq 7;
            end for;
            number +:= #extra;
        end if;
        return;
    end if;
    if used+Min(remaining,#rationalFlags-start+1) lt 3 then
        return;
    end if;
    for i in [start..#rationalFlags] do
        // At least two other distinct places occur in the full rational part.
        for m in [1..Min(5,remaining)] do
            TestRationalPart(i+1,remaining-m,used+1,
                VerticalJoin(rows,rationalFlags[i][m]),extra,~number);
        end for;
    end for;
end procedure;

total := 0;
for h in [0,2,3,4] do
    count := 0;
    TestRationalPart(1,7-h,0,empty,extras[h],~count);
    expected := #extras[h]*&+[Binomial(15,j)*Binomial(6-h,j-1) : j in [3..7-h]];
    assert count eq expected;
    total +:= count;
    print "Degree of nonrational part and number of divisors tested:", h,count;
end for;
print "All degree-seven divisors tested have l(D)=1:", total;
print "Gonality over F_5:", Degree(F);
