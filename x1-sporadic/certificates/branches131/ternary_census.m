SetQuitOnError(true); SetColumns(0);
load "base_ternary.m";
if not assigned lo then lo:="19"; end if;
if not assigned hi then hi:=lo; end if;
if not assigned cls then cls:="-1"; end if;
if not assigned shift then shift:="0"; end if;
shiftindex:=StringToInteger(shift);
extra:=[Q:Q in Places(K,1)|Q notin cusps];
assert #extra eq 2;
Cl,cr,dl:=ClassGroup(K);
assert Invariants(Cl) eq [5,125515,0];
cclass:=dl(c0);
function splits(P)
    cl:=dl(P)-Degree(P)*cclass;
    return (Integers()!Eltseq(cl)[2]) mod 13 eq 0;
end function;
function code(j)
    pol:=MinimalPolynomial(j,k);
    return Seqint([Integers()!v:v in Coefficients(pol)],3);
end function;
assert splits(c0) and not splits(ci);
cache:=AssociativeArray();
for f in [1..5] do
    for P in Places(K,f) do
        cache[P]:=splits(P);
    end for;
end for;
bad:=[P:P in Keys(cache) | not cache[P] or P eq c0];
Sort(~bad,func<a,b|Degree(a)-Degree(b)>);
for b in [StringToInteger(lo)..StringToInteger(hi)], h in [-26,-13,0,13,26], n in [0..4] do
    id:=(h div 13+2)*5+n;
    if StringToInteger(cls) ge 0 and id ne StringToInteger(cls) then continue; end if;
    D:=classdiv(b,h,n);
    if shiftindex gt 0 then D-:=13*(extra[shiftindex]-c0);end if;
    coord:=Eltseq(dl(D)-b*cclass); assert (Integers()!coord[2]) mod 13 eq 0;
    V,m:=RiemannRochSpace(D);
    F:=[m(v):v in Basis(V)]; dim:=#F;
    if dim eq 0 then continue; end if;
    blocks:=[* *];
    for Q in bad do
        vals:=[Valuation(a,Q):a in F]; e:=Min(vals);
        if e+Valuation(D,Q) gt 0 then blocks:=[* ZeroMatrix(k,dim,1) *];break; end if;
        assert e+Valuation(D,Q) eq 0;
        pivot:=Index(vals,e);
        L,ev:=ResidueClassField(Q);
        W,coords:=VectorSpace(L,k); mat:=Matrix(k,[Eltseq(coords(L!ev(a/F[pivot]))):a in F]);
        Append(~blocks,mat);
    end for;
    if exists{M:M in blocks|IsZero(M)} then continue; end if;
    processed:=0; survived:=0;
    for pivot in [1..dim] do
        for tail in VectorSpace(k,dim-pivot) do
            v:=Vector(k,[0:i in [1..pivot-1]] cat [1] cat Eltseq(tail));
            if exists{M:M in blocks|IsZero(v*M)} then continue; end if;
            processed+:=1;
            eff:=D+Divisor(m(V!Eltseq(v)));
            supp:=Support(eff);
            assert c0 notin supp and ci notin supp;
            // The class of eff has trivial Artin symbol. Once all smaller
            // places split, so does the last if its multiplicity is prime to 13.
            Sort(~supp,func<a,b|Degree(a)-Degree(b)>);
            good:=true;
            for i in [1..#supp-1] do
                Q:=supp[i]; if not IsDefined(cache,Q) then cache[Q]:=splits(Q); end if;
                if not cache[Q] then good:=false;break;end if;
            end for;
            if not good then continue;end if;
            assert Valuation(eff,supp[#supp]) mod 13 ne 0;
            survived+:=1;
            printf "case %o %o %o",b,h,n;
            for Q in supp do
                L,ev:=ResidueClassField(Q);
                printf " %o %o %o %o",Degree(Q),Valuation(eff,Q),code(ev(jfun)),code(ev(Jp(x,y)-jfun));
            end for;
            printf "\n";
        end for;
    end for;
    printf "summary %o %o %o %o %o %o\n",b,h,n,dim,processed,survived;
end for;
printf "complete %o %o %o %o\n",lo,hi,cls,shift;
quit;
