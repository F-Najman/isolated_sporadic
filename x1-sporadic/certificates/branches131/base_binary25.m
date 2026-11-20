// Run in certificates/ with magma -b fibres131.m.
SetQuitOnError(true);
SetColumns(0);
load "../fricke131-forms.m";
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
kX<X> := FunctionField(k);
S<U> := PolynomialRing(kX);
A := U^2+Evaluate(s0,X)*U-Evaluate(H0,X);
B := Evaluate(H1,X)-Evaluate(s1,X)*U;
K<u> := ext<kX | A^2+A*B-B^2*(X^3-X^2+X)>;
x := K!X;
y := Evaluate(A,u)/Evaluate(B,u);
assert y^2+y eq x^3-x^2+x and Genus(K) eq 11;
jfun := Evaluate(jj0,x)+y*Evaluate(jj1,x)+u*Jm(x,y);
cusps := Poles(x);
assert #cusps eq 2 and Seqset(Poles(jfun)) eq Seqset(cusps);
c0 := [P:P in cusps | Valuation(jfun,P) eq -131][1];
ci := [P:P in cusps | Valuation(jfun,P) eq -1][1];
assert Degree(c0) eq 1 and Degree(ci) eq 1;
E := EllipticCurve([k|0,-1,1,1,0]);
G := E![0,0];
assert Order(G) eq 3;

function classdiv(b,h,n)
    Q:=n*G;
    if Q eq E!0 then fibre:=c0+ci;
    else
        places := [P:P in Zeros(x-Q[1]) | Evaluate(y,P) eq Q[2]];
        fibre:=&+[Valuation(x-Q[1],P)*P:P in places];
    end if;
    assert Degree(fibre) eq 2;
    return b*c0+h*(c0-ci)+fibre-(c0+ci);
end function;
