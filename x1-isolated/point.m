// Check one explicit closed point on Sutherland's alternate X_1(n) model.
SetColumns(0);
SetQuitOnError(true);
if not assigned n then n:=17; else n:=StringToInteger(n); end if;
if not assigned modeldir then modeldir:="models"; end if;
R<t>:=PolynomialRing(Rationals());
input:=eval Read("points/point" cat IntegerToString(n) cat ".txt");
f:=R!input[1]; assert IsIrreducible(f);
K<a>:=NumberField(f);
x0:=Evaluate(R!input[2],a); y0:=Evaluate(R!input[3],a);
A<x,y>:=AffineSpace(Rationals(),2);
equation:=eval Read(modeldir cat "/FFFc" cat IntegerToString(n) cat ".txt");
C:=ProjectiveClosure(Curve(A,equation));
pt0:=C(K)![x0,y0,1]; assert not IsSingular(pt0);
assert Max(Degree(MinimalPolynomial(x0)),Degree(MinimalPolynomial(y0))) eq Degree(f);
r:=(x0^2*y0-x0*y0+y0-1)/(x0*(x0*y0-1)); s:=(x0*y0-y0+1)/(x0*y0);
b:=r*s*(r-1); c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,0,0]); Q:=E![0,0,1];
assert Discriminant(E) ne 0 and n*Q eq E!0;
assert &and[(n div l)*Q ne E!0:l in PrimeDivisors(n)];
if assigned prime then
    p:=StringToInteger(prime); assert IsPrime(p) and p gt 3 and GCD(p,n) eq 1;
    k:=GF(p); S:=PolynomialRing(k); reduced:=ChangeRing(C,k);
    fbar:=S!f; assert Degree(fbar) eq Degree(f) and IsSquarefree(fbar);
    divisors:=[];
    for z in Factorization(fbar) do
        h:=z[1];
        if Degree(h) eq 1 then L:=k; a0:=-Coefficient(h,0); else L<a0>:=ext<k|h>; end if;
        X:=Evaluate(S!input[2],a0); Y:=Evaluate(S!input[3],a0);
        pt:=reduced(L)![X,Y,1]; assert not IsSingular(pt);
        rr:=(X^2*Y-X*Y+Y-1)/(X*(X*Y-1)); ss:=(X*Y-Y+1)/(X*Y);
        bb:=rr*ss*(rr-1); cc:=ss*(rr-1);
        assert Discriminant(EllipticCurve([1-cc,-bb,-bb,0,0])) ne 0;
        assert ss^2-ss-rr+1 ne 0 and rr*ss-2*rr+1 ne 0;
        q:=Place(pt); assert Degree(q) eq Degree(h);
        Append(~divisors,Divisor(q));
    end for;
    D:=&+divisors;
    assert Degree(D) eq Degree(f);
    ell:=Dimension(RiemannRochSpace(D));
    print "good prime",p,"reduced dimension",ell;
else
    P:=Place(pt0); assert Degree(P) eq Degree(f);
    ell:=Dimension(RiemannRochSpace(Divisor(P)));
end if;
assert ell eq 1;
print "level",n,"degree",Degree(f),"l(P)",ell,"marked order",n;
quit;
