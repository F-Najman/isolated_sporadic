SetQuitOnError(true);SetColumns(0);load "../fricke131-forms.m";
v:=[StringToInteger(a):a in Split(Split(Read(file),"\n")[1]," ")];
assert v[1] eq 17 and v[2] eq 95;k:=GF(17);B<X,Y>:=PolynomialRing(k,2);
mon:=[B!1,0] cat [IsEven(i) select X^(i div 2) else X^((i-3) div 2)*Y:i in [2..95]];
t:=&+[v[2+i]*mon[i]:i in [1..96]];w:=&+[v[98+i]*mon[i]:i in [1..96]];
x:=k!14;y:=k!13;z:=k!1;xy:=[x,y];R0:=R(X,Y);
assert y^2+y eq x^3-x^2+x and z^2 eq Evaluate(R0,xy);
assert Evaluate(t,xy)+z*Evaluate(w,xy) eq 0;
yp:=(3*x^2-2*x+1)/(2*y+1);
zp:=(Evaluate(Derivative(R0,1),xy)+Evaluate(Derivative(R0,2),xy)*yp)/(2*z);
dh:=Evaluate(Derivative(t,1),xy)+Evaluate(Derivative(t,2),xy)*yp+
    zp*Evaluate(w,xy)+z*(Evaluate(Derivative(w,1),xy)+Evaluate(Derivative(w,2),xy)*yp);
assert dh ne 0;
j:=(Jp(x,y)+z*Jm(x,y))/2;
assert j eq 6 and Jp(x,y)-j eq 1;
for v in [j,Jp(x,y)-j] do
    tr:=TraceOfFrobenius(EllipticCurveWithjInvariant(v));
    assert Abs(tr) eq 6;
    F<U>:=PolynomialRing(GF(131));assert Sort([Order(r[1]^2):r in Roots(U^2-tr*U+17)]) eq [5,65];
end for;
print "simple zero at 17, j-invariants 6 and 1, absolute trace 6, projective marking orders 5 and 65";
quit;
