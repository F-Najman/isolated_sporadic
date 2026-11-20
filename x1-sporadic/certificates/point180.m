// Three fibres over census points on X_1(60) and X_1(45).
SetQuitOnError(true);
Q<T> := PolynomialRing(Rationals());
function CensusCurve(polynomial, B, bDen, C, cDen, n)
    assert IsIrreducible(Q!polynomial);
    K<z> := NumberField(Q!polynomial);
    b := Evaluate(Q!B,z)/bDen;
    c := Evaluate(Q!C,z)/cDen;
    E := EllipticCurve([1-c,-b,-b,0,0]);
    P := E![0,0];
    assert n*P eq E!0 and &and[(n div p)*P ne E!0 : p in PrimeDivisors(n)];
    // j,b,c are rational functions of the marked isomorphism class.
    assert Max([Degree(MinimalPolynomial(t)) : t in [jInvariant(E),b,c]]) eq Degree(K);
    return E;
end function;

function DivisionPolynomials(E)
    R<x> := PolynomialRing(BaseRing(E));
    b2,b4,b6,b8 := Explode(bInvariants(E));
    psi2Squared := 4*x^3 + b2*x^2 + 2*b4*x + b6;
    psi3 := 3*x^4 + b2*x^3 + 3*b4*x^2 + 3*b6*x + b8;
    psi4OverPsi2 := 2*x^6 + b2*x^5 + 5*b4*x^4 + 10*b6*x^3
                  + 10*b8*x^2 + (b2*b8-b4*b6)*x + b4*b8-b6^2;
    return x,psi2Squared,psi3,psi4OverPsi2;
end function;

procedure Tripling(E, expected)
    K := BaseRing(E);
    x,D,psi3,psi4OverPsi2 := DivisionPolynomials(E);
    // The x-coordinate of [3] is numerator/psi3^2.
    numerator := x*psi3^2 - D*psi4OverPsi2;
    assert Degree(numerator) eq 9 and IsIrreducible(numerator);
    assert GCD(numerator,psi3) eq 1;
    L<u> := ext<K | numerator>;
    EL := BaseChange(E,L);
    // Pulling back dx/(2y+a1*x+a3) recovers the required y-coordinate.
    a1,a2,a3,a4,a6 := Explode(aInvariants(EL));
    y := (3*a3*Evaluate(psi3,u)^2/Evaluate(Derivative(numerator),u) - a1*u-a3)/2;
    P180 := EL![u,y];
    assert 3*P180 eq EL![0,0];
    assert 180*P180 eq EL!0 and &and[(180 div p)*P180 ne EL!0 : p in [2,3,5]];
    assert Degree(K)*Degree(L) eq expected;
    print "Tripling fibre: base degree, relative degree, residue degree",
        Degree(K),Degree(L),expected;
end procedure;

E := CensusCurve([-1,-1,0,6,5,-5,-4,-7,6,6,-5,3,-3,1],
    [-3185,-4960,12292,12874,-5257,-6986,-22178,12721,15921,-12764,8660,-7785,2407],59,
    [-66,-916,509,2689,-1079,345,-3051,-488,4370,-3223,2156,-1437,379],59,60);
assert Degree(MinimalPolynomial(jInvariant(E))) eq 13;
Tripling(E,117);

E := CensusCurve([-2,8,-20,32,-5,-8,-56,47,2,32,-22,5,-14,2,-3,1],
    [-5118765286984,10347780302453,-27196557622320,20291953540247,22966163991731,
      32851392757372,-61627321002265,-10000031592676,-24837679795463,26864802025524,
      1135919756101,16197709622111,-1412669275541,3093224947121,-1174456886897],414372187720,
    [-951130654218,1392074773512,-4134345886833,3823717754347,4331802399786,
      5108324552951,-11113569696615,-2501017835012,-3764147948422,4826772355481,
      543426226373,2757947596612,-120292659741,483271182730,-200166951401],207186093860,60);
Tripling(E,135);

E := CensusCurve([-1,3,-1,1,-12,10,6,-14,11,-5,1],
    [-694,3205,-6203,5358,1098,-6782,6995,-4072,1443,-236],27,
    [86,-368,1006,-1782,732,1255,-1912,1325,-549,103],27,45);
K := BaseRing(E);
x,D,psi3,primitive4 := DivisionPolynomials(E);
assert Degree(primitive4) eq 6 and IsIrreducible(primitive4);
L<u> := ext<K | primitive4/LeadingCoefficient(primitive4)>;
S<t> := PolynomialRing(L);
quadratic := t^2-Evaluate(D,u);
assert IsIrreducible(quadratic);
M<w> := ext<L | quadratic>;
EM := BaseChange(E,M);
a1,a2,a3,a4,a6 := Explode(aInvariants(EM));
T4 := EM![u,(w-a1*u-a3)/2];
assert 4*T4 eq EM!0 and 2*T4 ne EM!0;
P180 := 34*(EM![0,0]) + T4;
assert 4*P180 eq EM![0,0];
assert 180*P180 eq EM!0 and &and[(180 div p)*P180 ne EM!0 : p in [2,3,5]];
assert Degree(K)*Degree(L)*Degree(M) eq 120;
print "Primitive order-four fibre: base degree, relative degree, residue degree",10,12,120;
