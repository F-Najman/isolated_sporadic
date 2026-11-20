// Enumerate the divisor classes on B compatible with a full lift.
SetQuitOnError(true); SetColumns(0);
load "../branches131/base_binary25.m";
Cl,rep,dl:=ClassGroup(K); assert Invariants(Cl) eq [2,4290,0]; c:=dl(c0);
function residueDegree(Q)
    a:=Integers()!Eltseq(dl(Q)-Degree(Q)*c)[2];
    return 65 div GCD(65,a);
end function;
function polyCode(j)
    return Seqint([Integers()!a:a in Coefficients(MinimalPolynomial(j,k))],2);
end function;
assert &and[Degree(P)*residueDegree(P) ge 65:P in Zeros(jfun)];
assert residueDegree(c0) eq 1 and residueDegree(ci) eq 65;
supersingular:=assigned ss and ss eq "1";
shiftSS:=0*(c0-ci);
if supersingular then
    Qss:=[P:P in Places(K,1)|P notin cusps][1];
    assert Evaluate(jfun,Qss) eq 0;
    ssplaces:=Zeros(jfun);
    assert &and[residueDegree(Q) eq 65:Q in ssplaces];
    assert {Q:Q in ssplaces|Degree(Q)*residueDegree(Q) lt 130} eq {Qss};
    shiftSS:=65*(Qss-c0);
    G1:=dl(classdiv(0,0,1)); Vss:=dl(shiftSS);
    assert Order(Vss) eq 6 and Vss notin {i*G1:i in [0..2]};
end if;
upper:=StringToInteger(hi); lower:=StringToInteger(lo);
assert upper le 38 and lower ge 1;
n:=StringToInteger(cl); a:=StringToInteger(shift);
shard:=StringToInteger(part); bits:=StringToInteger(cut);
cache:=AssociativeArray(); bad:=[]; exceptional:=[];
for f in [1..8] do
    for Q in Places(K,f) do
        r:=residueDegree(Q);cache[Q]:=r;
        if r ne 1 or Q in cusps then Append(~bad,Q);end if;
        if Q notin cusps and r ne 1 and r*f le 38 then Append(~exceptional,Q);end if;
    end for;
end for;
assert #exceptional eq 1;
Q0:=exceptional[1]; assert Degree(Q0) eq 3 and cache[Q0] eq 5;
code0:=polyCode(Evaluate(jfun,Q0)); dual0:=polyCode(Evaluate(Jp(x,y)-jfun,Q0));
assert code0 eq 13;
Sort(~bad,func<u,v|Degree(u)-Degree(v)>);
for e in [lower..upper] do
    b:=e-15*a; if b lt 0 then continue;end if;
    D:=classdiv(b,0,n)-5*a*(Q0-3*c0)-shiftSS;
    assert (Integers()!Eltseq(dl(D)-b*c)[2]) mod 65 eq 0;
    if b eq 0 then
        if shard eq 0 and IsPrincipal(D) then
            printf "case %o %o %o 3 %o 5 %o %o\n",e,n,a,a,code0,dual0;
        end if;
        continue;
    end if;
    V,m:=RiemannRochSpace(D); F:=[m(v):v in Basis(V)]; dim:=#F;
    if dim eq 0 then continue;end if;
    one:=ZeroMatrix(k,dim,0); others:=[* *]; skip:=false;
    for Q in bad do
        vals:=[Valuation(g,Q):g in F]; v:=Min(vals);
        if v+Valuation(D,Q) gt 0 then skip:=true; break;end if;
        assert v+Valuation(D,Q) eq 0;
        idx:=Index(vals,v); L,ev:=ResidueClassField(Q); W,coords:=VectorSpace(L,k);
        col:=Matrix(k,[Eltseq(coords(L!ev(g/F[idx]))):g in F]);
        if Degree(Q) eq 1 then one:=HorizontalJoin(one,col);else Append(~others,col);end if;
    end for;
    if skip then continue;end if;
    ok,v0:=IsConsistent(one,Vector(k,[1:i in [1..Ncols(one)]]));
    if not ok then continue;end if;
    B:=BasisMatrix(Nullspace(one)); t:=Min(bits,Nrows(B));
    if shard ge 2^t then continue;end if;
    if Nrows(B) gt 0 then
        head:=Intseq(shard,2,t) cat [0:i in [t+1..Nrows(B)]];
        v0+:=Vector(k,head)*B;
    end if;
    rest:=Nrows(B)-t;
    if rest gt 0 then B:=Submatrix(B,t+1,1,rest,dim);end if;
    survivors:=0; processed:=0;
    for vec in VectorSpace(k,rest) do
        v:=v0; if rest gt 0 then v+:=vec*B;end if;
        if exists{A:A in others|IsZero(v*A)} then continue;end if;
        processed+:=1;
        eff:=D+Divisor(m(V!Eltseq(v))); supp:=Support(eff);
        Sort(~supp,func<P,Q|Degree(P)-Degree(Q)>);
        good:=true;
        for i in [1..#supp-1] do
            Q:=supp[i];
            if not IsDefined(cache,Q) then cache[Q]:=residueDegree(Q);end if;
            if cache[Q] ne 1 then good:=false;break;end if;
        end for;
        if not good then continue;end if;
        assert GCD(Valuation(eff,supp[#supp]),65) eq 1;
        printf "case %o %o %o",e,n,a;
        if a gt 0 then printf " 3 %o 5 %o %o",a,code0,dual0;end if;
        for Q in supp do
            printf " %o %o 1 %o %o",Degree(Q),Valuation(eff,Q),polyCode(Evaluate(jfun,Q)),polyCode(Evaluate(Jp(x,y)-jfun,Q));
        end for;
        printf "\n"; survivors+:=1;
    end for;
    printf "summary %o %o %o %o %o %o\n",e,n,a,shard,processed,survivors;
end for;
printf "complete %o %o %o %o %o %o\n",lo,hi,cl,shift,part,cut;
quit;
