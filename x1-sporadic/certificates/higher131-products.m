SetQuitOnError(true);
load "higher131-forms.m";
// Output a key for the closed point of Y and its absolute residue degree.
// On an ordinary normalized curve the five x-coordinates of H.P determine
// the marking modulo +/-H. The chosen Frobenius representative is immaterial.
function placeKey(P)
    E := Curve(P); L := BaseRing(E); U<X> := PolynomialRing(L);
    ker := &*[X-(h*P)[1]:h in H];
    tuple := [jInvariant(E)] cat [Coefficient(ker,i):i in [0..4]];
    f := Lcm([Degree(MinimalPolynomial(v)):v in tuple]);
    F := GF(2,f);
    Embed(F,L);
    v := [F!x:x in tuple];
    reps := [];
    for i in [1..f] do
        Append(~reps,[Seqint([Integers()!a:a in Eltseq(x)],2):x in v]);
        v := [x^2:x in v];
    end for;
    assert #Seqset(reps) eq f;
    return Sprint([f] cat Min(reps)),f;
end function;

// In the ordinary family, X_a denotes x(aP) as a formal series in u=eta-eta0.
// Doubling gives X_a=X_{a/2}^2+(eta0+u)/X_{a/2}^2.
// Squaring doubles the precision at each iteration.
function jetBlock(P,m)
    E := Curve(P); L := BaseRing(E);
    a1,a2,a3,a4,a6 := Explode(aInvariants(E));
    assert a1 eq 1 and a2 eq 0 and a3 eq 0;
    eta0 := a6+a4^2;
    assert eta0*jInvariant(E) eq 1;
    mults := [i*P:i in [1..130]];
    slope := [(Q[1]^2+Q[2]+a4)/Q[1]:Q in mults];
    if m eq 1 then return conditions([slope]); end if;
    S<u> := PowerSeriesRing(L,m+1);
    xx := [S!Q[1]:Q in mults];
    prec := 1;
    while prec lt m do
        half := [xx[(66*a) mod 131]^2:a in [1..130]];
        xx := [v+(eta0+u)/v:v in half];
        prec *:= 2;
    end while;
    jets := [slope];
    assert &and[slope[a]^2+slope[a] eq Coefficient(xx[(2*a) mod 131],0):a in [1..130]];
    for n in [1..m-1] do
        next := [Coefficient(xx[(2*a) mod 131],n):a in [1..130]];
        if IsEven(n) then next := [next[a]+jets[1+n div 2][a]^2:a in [1..130]]; end if;
        Append(~jets,next);
    end for;
    assert jets[2] eq [1/Q[1]^2:Q in mults];
    return conditions(jets);
end function;

procedure record(P,weight,~points,~degrees,~coeffs)
    name,f := placeKey(P);
    assert weight mod f eq 0;
    if not IsDefined(points,name) then points[name]:=P; degrees[name]:=f; end if;
    coeffs[name] := (IsDefined(coeffs,name) select coeffs[name] else 0) + weight div f;
end procedure;

procedure components(b,code,shifts,~points,~degrees,~memo)
    P,images := markedPoint(b,code);
    for n in shifts do
        coeffs := AssociativeArray();
        v := 3^n mod 131;
        record(v*P,-b,~points,~degrees,~coeffs);
        record(3*v*P,-2*b,~points,~degrees,~coeffs);
        record(9*v*P,3*b,~points,~degrees,~coeffs);
        for pair in images do
            Q,t := Explode(pair);
            record(v*Q,b*t,~points,~degrees,~coeffs);
            record(3*v*Q,-b*t,~points,~degrees,~coeffs);
        end for;
        memo[<code,n>]:=coeffs;
    end for;
end procedure;


function selectedColumns(A,cols)
    return Matrix(k,[[A[i,j]:j in cols]:i in [1..Nrows(A)]]);
end function;
function productRank(bases,qeval)
    left,right:=Explode(bases);
    a:=[left*Transpose(rows):rows in qeval];
    b:=[right*Transpose(rows):rows in qeval];
    length:=Ncols(a[1]);
    products:=[];
    for i in [1..Nrows(left)], j in [1..Nrows(right)] do
        row:=[];
        for c in [1..13] do
            for n in [2..length+1] do
                Append(~row,&+[a[c][i,v]*b[c][j,n-v]:v in [1..n-1]]);
            end for;
        end for;
        Append(~products,row);
    end for;
    return Rank(Matrix(k,products));
end function;

procedure obstruct(cases,label,cusps,memo,points,degrees,qcols,qeval,~jets,~results)
    D:=AssociativeArray();
    for pair in cases do
        one:=memo[pair];
        for key in Keys(one) do
            D[key]:=(IsDefined(D,key) select D[key] else 0)+one[key];
        end for;
    end for;
    degree:= &+[Integers()|D[key]*degrees[key]:key in Keys(D)|D[key] gt 0];
    assert &+[Integers()|D[key]*degrees[key]:key in Keys(D)] eq 0;
    assert degree gt 0;
    success:=false; ranks:=[]; both:=[* *];
    for sign in [-1,1] do
        B:=cusps;
        for key in Keys(D) do
            m:=D[key]*sign;
            if m le 0 then continue; end if;
            if not IsDefined(jets,<key,m>) then
                jets[<key,m>]:=BasisMatrix(RowSpace(jetBlock(points[key],m)));
            end if;
            B:=VerticalJoin(B,jets[<key,m>]);
        end for;
        rank:=Rank(B)-Rank(cusps);
        assert rank le degree;
        Append(~ranks,rank); Append(~both,B);
        if rank eq degree then success:=true; break; end if;
    end for;
    if not success then
        assert ranks eq [degree-1,degree-1];
        bases:=[BasisMatrix(Nullspace(Transpose(selectedColumns(B,qcols)))):B in both];
        upper:=131-(degree-1);
        rank:=productRank(bases,qeval);
        printf "canonical product degree %o dimensions %o product rank %o symmetric bound %o\n",
            degree,[Nrows(B):B in bases],rank,upper*(upper+1) div 2;
        success:=rank gt upper*(upper+1) div 2;
    end if;
    if not success then printf "remaining %o %o degree %o ranks %o\n",label,cases,degree,ranks; end if;
    Append(~results,success);
end procedure;

cusps:=cuspConditions(); assert Rank(cusps) eq 25;
qs:=cuspSeries(23);
stack:=VerticalJoin(cusps,Matrix(k,[Eltseq(v):v in &cat qs]));
ech:=BasisMatrix(RowSpace(stack));
qcols:=[Min([j:j in [1..Ncols(ech)] | ech[i,j] ne 0]):i in [1..Nrows(ech)]];
print "weight-two evaluation rank",#qcols; assert #Seqset(qcols) eq #qcols;
qeval:=[selectedColumns(Matrix(k,[Eltseq(v):v in a]),qcols):a in qs];
points:=AssociativeArray(); degrees:=AssociativeArray(); memo:=AssociativeArray(); jets:=AssociativeArray();
single:=[378895,524145];
pairs:=[<3,13,15,35825>,<9,617,9,865>,<8,419,10,1591>];
repeated:=[5451,4411];
for code in single do components(18,code,[0],~points,~degrees,~memo); end for;
todo:={ <a[1],a[2]>:a in pairs } join { <a[3],a[4]>:a in pairs } join {<12,code>:code in repeated};
for item in todo do
    b,code:=Explode(item);
    components(b,code,[0..12],~points,~degrees,~memo);
end for;
printf "hard singles %o pairs %o repeated %o\n",single,pairs,repeated;
out:=[];
for code in single do
    obstruct([<code,0>],"single",cusps,memo,points,degrees,qcols,qeval,~jets,~out);
end for;
for row in pairs do
    a,j,b,l:=Explode(row);
    for k0 in [0..12] do
        obstruct([<j,0>,<l,k0>],"pair",cusps,memo,points,degrees,qcols,qeval,~jets,~out);
    end for;
end for;
for l in repeated do
    for i in [0..12], j in [i..12] do
        obstruct([<13,i>,<13,j>,<l,0>],"repeated",cusps,memo,points,degrees,qcols,qeval,~jets,~out);
    end for;
end for;
printf "confirmed %o of %o\n",#[v:v in out|v],#out;
assert #out eq 223 and &and out;
quit;
