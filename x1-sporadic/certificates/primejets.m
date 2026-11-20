SetQuitOnError(true);
p:=StringToInteger(N); mult:=p eq 59 select 3 else 10; basic:=p eq 59 select 2 else 3; r:=(p-1) div 2;
S:=CuspidalSubspace(ModularSymbolsH(p,[1],2,-1)); print "genus",Dimension(S); assert Dimension(S) eq ((p-5)*(p-7)) div 24;
ell:=1;repeat ell:=NextPrime(ell); until Order(GF(p)!ell) eq p-1;
T:=func<a| IntegralHeckeOperator(S,a)>; Q:=Rationals();
G:=ChangeRing((ChangeRing(T(ell)^2-T(ell^2),Q)/ell),Integers());
k:=GF(1009); g:=ChangeRing(G,k); assert g^r eq Parent(g)!1;
V:=VectorSpace(k,Dimension(S)); SetSeed(20261118); v:=Random(V);
H:=[ChangeRing(T(i),k):i in [1..mult]];
R:=[[v*g^j*Tn:Tn in H]:j in [0..r-1]];
A:=&cat [row[1..basic]: row in R]; assert Rank(Matrix(A)) eq #A;
if p eq 59 then
 w:=A cat R[1][3..3]; assert Rank(Matrix(w)) eq #w;
else
 for i in [2..r-1] do for j in [i+1..r] do
 w:=A cat R[1][4..mult] cat R[i][4..mult] cat R[j][4..mult];
 assert Rank(Matrix(w)) eq #w;
 end for;end for;
end if;
print "JET_PASS",p,r,basic,mult;quit;
