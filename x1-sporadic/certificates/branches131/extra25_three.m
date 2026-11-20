load "base_ternary.m";
Cl,rep,dl:=ClassGroup(K); T:=dl(c0-ci); G1:=dl(classdiv(0,0,1));
assert Invariants(Cl) eq [5,125515,0] and Order(G) eq 5;
extra:=[Q:Q in Places(K,1)|Q notin cusps];assert #extra eq 2;
for Q in extra do
    assert &and[h*T+n*G1 ne 13*dl(Q-c0):h in [-26,-13,0,13,26],n in [0..4]];
end for;
ss:=Zeros(jfun); assert #ss eq 11;
for Q in ss do
    assert Degree(Q) eq 2 and (Integers()!Eltseq(dl(Q)-2*dl(c0))[2]) mod 13 ne 0;
end for;
eta:=x^5*Differential(x)/((2*y+1)*z);
assert IsEffective(Divisor(eta)) and &and[Valuation(Divisor(eta),Q) eq 0:Q in cusps];
print "empty shifted classes excluded, eleven nonsplit supersingular places, nonvanishing anti-Fricke differential at both cusps";
quit;
