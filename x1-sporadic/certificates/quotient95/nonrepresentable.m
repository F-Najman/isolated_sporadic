load "artin.m";
D:=5*c0+35*(ci-c0);V,rr:=RiemannRochSpace(D);assert Dimension(V) eq 1;
Z:=D+Divisor(rr(V.1));W0:=Support(Z)[1];
assert Degree(W0) eq 5 and Valuation(Z,W0) eq 1 and Order(im(dl(W0))) eq 13;
jc:=MinimalPolynomial(Evaluate(jfun,W0),k);assert Degree(jc) eq 5;
jcode:=Seqint([Integers()!a:a in Coefficients(jc)],3);
load "../branches131/ternary_forms.m";
load "../branches131/ternary_cols.m";
precision:=100;
load "tate_jets.m";
K:=GF(3,5);q:=#K; j:=Roots(ChangeRing(R!Intseq(jcode,3),K))[1][1];
E:=EllipticCurve([K|0,1,0,0,-1/j]);tr:=TraceOfFrobenius(E);
eig:=[z[1]:z in Roots(W^2-tr*W+q)];good:=[v:v in eig|Order(v^2) eq 13];assert #good gt 0;
for lambda in good do
    n:=Order(lambda);L:=GF(3,5*n);Embed(K,L); EL:=BaseChange(E,L);
    a:=2;b:=tr;for i in [2..n] do a,b:=b,tr*b-q*a;end for;
    card:=q^n+1-b;cof:=card div 131^Valuation(card,131);
    repeat
        P:=cof*Random(EL);while P ne EL!0 and 131*P ne EL!0 do P:=131*P;end while;
        if P ne EL!0 and #eig eq 2 then
            other:=[v:v in eig|v ne lambda][1];P:=EL![P[1]^q,P[2]^q]-(Integers()!other)*P;
        end if;
    until P ne EL!0 and EL![P[1]^q,P[2]^q] eq (Integers()!lambda)*P;
    if #eig eq 1 then assert (tr^2-4*q) mod 131^2 ne 0;end if;
    A:=absoluteRows(productValues(tangentSlopes(P),1))*canon;
    C:=Submatrix(J,1,1,30,Ncols(J))*canon; rr:=Rank(VerticalJoin(A,C));
    printf "degree-five base place, trace %o, eigenvalue %o, ranks %o %o %o\n",tr,lambda,Rank(A),Rank(C),rr;
    assert Rank(A) eq 65 and rr eq 95;
end for;
quit;