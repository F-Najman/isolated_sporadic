// A degree-45 point above y = -1/(x^3+x^2) on Sutherland's plane model
// https://math.mit.edu/~drew/X1/X1opt87.txt
// Checking the pointed elliptic curve needs no external files.
Q<t>:=PolynomialRing(Rationals());
g:=t^45 + 33*t^44 + 510*t^43 + 4884*t^42 + 32280*t^41 + 154893*t^40 + 
    550411*t^39 + 1437567*t^38 + 2615876*t^37 + 2671355*t^36 - 801777*t^35 - 
    7999663*t^34 - 12594901*t^33 - 4736202*t^32 + 14509197*t^31 + 24451029*t^30 
    + 6060501*t^29 - 24148133*t^28 - 25443805*t^27 + 7407090*t^26 + 
    27433066*t^25 + 5827851*t^24 - 20027124*t^23 - 9830840*t^22 + 12963444*t^21 
    + 9218439*t^20 - 7884260*t^19 - 7075846*t^18 + 4209681*t^17 + 4432894*t^16 -
    1925940*t^15 - 2183163*t^14 + 793816*t^13 + 823181*t^12 - 306834*t^11 - 
    226414*t^10 + 104770*t^9 + 39297*t^8 - 27319*t^7 - 1898*t^6 + 4389*t^5 - 
    752*t^4 - 243*t^3 + 117*t^2 - 18*t + 1;
assert Degree(g) eq 45 and IsIrreducible(g);
K<x>:=NumberField(g);
y:=-1/(x^3+x^2);
q:=(1+x)/x;
v:=1/((y+1)*(x+1));
E:=EllipticCurve([(1+q)*v+2-q,0,q*v^2+(1-q)*v,0,0]);
P:=E![-v,v^2];
assert Discriminant(E) ne 0;
assert 87*P eq E!0 and 29*P ne E!0 and 3*P ne E!0;
assert Degree(MinimalPolynomial(jInvariant(E))) eq 45;
print "VERIFIED", Degree(K), 87;
quit;
