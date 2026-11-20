/* Evaluate all 26 points of an ordinary Fricke fibre over F_(2^13).
   Both eigenlines are computed on one curve so the two kernels are dual. */
function ordinaryFibre(code)
    f:=13; K:=GF(2,f); pol:=R!Intseq(code,2);
    j:=Roots(ChangeRing(pol,K))[1][1]; C:=EllipticCurve([K|1,0,0,0,1/j]);
    q:=#K; tr:=TraceOfFrobenius(C);
    rr:=Roots(W^2-tr*W+q); assert #rr eq 2;
    roots:=[Integers()!v[1]:v in rr];
    assert code eq 12707 and tr eq 5 and Seqset(roots) eq {58,78};
    roots:=[58,78];
    n:=Lcm([Order(v[1]):v in rr]); L:=GF(2,f*n);Embed(K,L);
    EC:=BaseChange(C,L); t0:=2;t1:=tr;
    for i in [2..n] do t0,t1:=t1,tr*t1-q*t0;end for;
    ord:=q^n+1-t1;
    pp:=[];
    for i in [1,2] do
        repeat
            P:=(ord div 131^Valuation(ord,131))*Random(EC);
            while P ne EC!0 and 131*P ne EC!0 do P:=131*P;end while;
            if P ne EC!0 then P:=EC![P[1]^q,P[2]^q]-roots[3-i]*P;end if;
        until P ne EC!0;
        assert EC![P[1]^q,P[2]^q] eq roots[i]*P;
        Append(~pp,P);
    end for;
    U<X>:=PolynomialRing(L);
    target,phi:=IsogenyFromKernel(EC,&*[X-(a*pp[1])[1]:a in [1..65]]);
    assert jInvariant(target) ne j;
    assert MinimalPolynomial(jInvariant(target),k) eq R!Intseq(14243,2);
    markings:=[*pp[1],phi(pp[2]) *];
    rows:=[];
    for P in markings do
        m:=slopesAt(P);
        for i in [0..12] do
            vals:=valsAt(m,3^i mod 131);
            assert &and[v^q eq v:v in vals];
            Append(~rows,[K!v:v in vals]);
        end for;
    end for;
    F:=Matrix(K,rows)*ChangeRing(canon,K);
    assert Nrows(F) eq 26 and Rank(F) eq 25;
    printf "ordinary fibre %o rank %o\n",code,Rank(F);
    return F;
end function;
