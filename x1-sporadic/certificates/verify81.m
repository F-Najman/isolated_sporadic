// Intersect the earlier Sutherland model with y=-x.
// https://math.mit.edu/~drew/X1/X1opt81.txt
Q<t>:=PolynomialRing(Rationals());
g:=t^48 + t^47 + 19*t^46 + 39*t^45 + 200*t^44 + 955*t^43 + 2832*t^42 + 
    9793*t^41 + 36045*t^40 + 104129*t^39 + 259458*t^38 + 623313*t^37 + 
    1317455*t^36 + 2214155*t^35 + 2822259*t^34 + 1685926*t^33 - 5921221*t^32 - 
    29696585*t^31 - 79899709*t^30 - 157075513*t^29 - 247325490*t^28 - 
    335096985*t^27 - 425607576*t^26 - 549450408*t^25 - 738072737*t^24 - 
    994316643*t^23 - 1279938432*t^22 - 1514396686*t^21 - 1589555818*t^20 - 
    1424228830*t^19 - 1036561543*t^18 - 556782590*t^17 - 148308631*t^16 + 
    92917543*t^15 + 172336247*t^14 + 155003511*t^13 + 106082449*t^12 + 
    61614075*t^11 + 31724577*t^10 + 14580801*t^9 + 5864665*t^8 + 1994487*t^7 + 
    550455*t^6 + 117292*t^5 + 17870*t^4 + 1630*t^3 + 25*t^2 - 12*t - 1;
assert Degree(g) eq 48 and IsIrreducible(g);
K<x>:=NumberField(g);
q:=-x;
v:=-(x+1);
E:=EllipticCurve([(1+q)*v+2-q,0,q*v^2+(1-q)*v,0,0]);
P:=E![-v,v^2];
assert Discriminant(E) ne 0 and 81*P eq E!0 and 27*P ne E!0;
assert Degree(MinimalPolynomial(jInvariant(E))) eq 48;
print "VERIFIED", Degree(K), 81;
quit;
