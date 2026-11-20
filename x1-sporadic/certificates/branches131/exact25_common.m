// Common ordinary-place and divisor routines for the hard Hecke configurations.
function placeKey(P)
    E:=Curve(P); L:=BaseRing(E); U<X>:=PolynomialRing(L);
    a1,a2,a3,a4,a6:=Explode(aInvariants(E));
    r0:=Characteristic(k) eq 3 select -a4/2 else L!0;
    ker:=&*[X-((h*P)[1]-r0):h in H];
    val:=[jInvariant(E)] cat [Coefficient(ker,i):i in [0..4]];
    f:=Lcm([Degree(MinimalPolynomial(v,k)):v in val]);
    F:=GF(Characteristic(k),f); Embed(F,L);
    u:=[F!v:v in val];
    seqs:=[];
    for i in [1..f] do
        Append(~seqs,[Seqint([Integers()!a:a in Eltseq(x)],Characteristic(k)):x in u]);
        u:=[x^Characteristic(k):x in u];
    end for;
    assert #Seqset(seqs) eq f;
    return Sprint([f] cat Min(seqs)),f;
end function;
function jetRows(P,m)
    E:=Curve(P); L:=BaseRing(E); a1,a2,a3,a4,a6:=Explode(aInvariants(E));
    if m eq 1 then return simpleEval(P);end if;
    pts:=[];Q:=P;
    for i in [1..130] do Append(~pts,Q);Q+:=P;end for;
    S<u>:=PowerSeriesRing(L,m+1);
    eta0:=1/jInvariant(E); eta:=eta0+u;
    if Characteristic(k) eq 2 then
        assert a1 eq 1 and a2 eq 0 and a3 eq 0 and eta0 eq a6+a4^2;
        sl:=[(Q[1]^2+a4+Q[2])/Q[1]:Q in pts];
        xx:=[S!Q[1]:Q in pts]; n:=1;
        while n lt m do
            half:=[xx[66*a mod 131]^2:a in [1..130]];
            xx:=[v+eta/v:v in half];n*:=2;
        end while;
        jets:=[sl];
        for i in [1..m-1] do
            col:=[Coefficient(xx[2*a mod 131],i):a in [1..130]];
            if IsEven(i) then col:=[col[a]+jets[i div 2+1][a]^2:a in [1..130]];end if;
            Append(~jets,col);
        end for;
        B:=chooseCols(conditions(jets),cols)*canon;
    else
        assert a1 eq 0 and a2 eq 1 and a3 eq 0;
        r0:=-a4/2;
        xx:=[S!(Q[1]-r0):Q in pts]; yy:=[S!Q[2]:Q in pts]; n:=1;
        while n lt m do
            third:=[xx[44*a mod 131]^3:a in [1..130]];
            xx:=[(v^3+eta*v-eta^3)/(v-eta)^2:v in third];n*:=3;
        end while;
        assert &and[Coefficient(xx[a],0) eq pts[a][1]-r0:a in [1..130]];
        rhs:=[v^3+v^2-eta:v in xx];n:=1;
        while n lt m do yy:=[(yy[a]+rhs[a]/yy[a])/2:a in [1..130]];n*:=2;end while;
        sl:=[xx[a]/yy[a]:a in [1..130]];
        B:=ZeroMatrix(k,0,Ncols(canon));
        for n in [0..m-1] do
            values:=[&+[Coefficient(sl[h*v[1] mod 131]*sl[h*v[2] mod 131],n):h in H]:v in selected];
            B:=VerticalJoin(B,absRows(values)*canon);
        end for;
    end if;
    return BasisMatrix(RowSpace(B));
end function;
procedure addPoint(P,weight,~places,~degrees,~D)
    key,f:=placeKey(P);assert weight mod f eq 0;
    if not IsDefined(places,key) then places[key]:=P; degrees[key]:=f;end if;
    D[key]:=(IsDefined(D,key) select D[key] else 0)+(weight div f);
end procedure;
procedure recComponents(f,j,dual,shifts,~places,~degrees,~ans)
    data:=allMarked(f,j,dual); opts:=[* *];
    for details in data do
        P,imgs:=Explode(details);
        for n in shifts do
            d:=diamond^n mod 131; D:=AssociativeArray();
            if Characteristic(k) eq 3 then
                addPoint(d*P,-f,~places,~degrees,~D);
                addPoint(2*d*P,-f,~places,~degrees,~D);
                addPoint(4*d*P,2*f,~places,~degrees,~D);
                for entry in imgs do
                    Q,t:=Explode(entry);
                    addPoint(d*Q,t*f,~places,~degrees,~D);
                    addPoint(2*d*Q,-t*f,~places,~degrees,~D);
                end for;
            else
                addPoint(d*P,-f,~places,~degrees,~D);
                addPoint(3*d*P,-3*f,~places,~degrees,~D);
                for entry in imgs do
                    Q,t:=Explode(entry);
                    addPoint(d*Q,t*f,~places,~degrees,~D);
                end for;
            end if;
            Append(~opts,D);
        end for;
    end for;
    ans:=opts;
end procedure;
function fullChoices(opts,number)
    D:=AssociativeArray();
    for comp in opts do
        v:=comp[(number mod #comp)+1];number:=number div #comp;
        for key in Keys(v) do D[key]:=(IsDefined(D,key) select D[key] else 0)+v[key];end for;
    end for;
    return D;
end function;
procedure sideRank(D,sign,fixed,places,degrees,~jets,~ans)
    B:=fixed;
    for key in Keys(D) do
        m:=sign*D[key]; if m le 0 then continue;end if;
        if not IsDefined(jets,<key,m>) then jets[<key,m>]:=jetRows(places[key],m);end if;
        B:=VerticalJoin(B,jets[<key,m>]);
    end for;
    ans:=Rank(B)-Rank(fixed);
end procedure;
input:=Open(file,"r");
places:=AssociativeArray(); degrees:=AssociativeArray(); jets:=AssociativeArray();
cache:=AssociativeArray();
while true do
    line:=Gets(input); if IsEof(line) then break;end if;
    ints:=[StringToInteger(w):w in Split(line," ")]; use:=ints[1];Remove(~ints,1);
    factors:=[<ints[4*i-3],ints[4*i-2],ints[4*i-1],ints[4*i]>:i in [1..#ints div 4]];
    Sort(~factors,func<a,b|b[1]-a[1]>);
    opts:=[* *];
    for t in factors do
        f,m,j,jw:=Explode(t);
        for rep in [1..m] do
            moves:=#opts eq 0 select [0] else [0..12]; tag:=<j,jw,#moves>;
            if not IsDefined(cache,tag) then
                res:=[* *];recComponents(f,j,jw,moves,~places,~degrees,~res);cache[tag]:=res;
            end if;
            Append(~opts,cache[tag]);
        end for;
    end for;
    count:=&*[#v:v in opts]; missed:=0;
    fixed:=use eq 0 select ZeroMatrix(k,0,Ncols(canon)) else fibreO;
    for choice in [0..count-1] do
        D:=fullChoices(opts,choice);
        e:=&+[Integers()|D[key]*degrees[key]:key in Keys(D)|D[key] gt 0];
        assert &+[Integers()|D[key]*degrees[key]:key in Keys(D)] eq 0;
        ranks:=[];
        for sign in [1,-1] do
            rank:=0;sideRank(D,sign,fixed,places,degrees,~jets,~rank); Append(~ranks,rank);
            if e gt 0 and rank eq e then break;end if;
        end for;
        if e gt 26 and ranks[#ranks] ne e and assigned move then
            assert use eq 1;
            for sign in [1,-1] do
                B:=ZeroMatrix(k,0,Ncols(canon));
                for key in Keys(D) do
                    m:=sign*D[key]; if m le 0 then continue;end if;
                    B:=VerticalJoin(B,jets[<key,m>]);
                end for;
                A:=ChangeRing(B,BaseRing(moved));
                rk:=Rank(VerticalJoin(A,moved));
                printf "moved %o sign %o degree %o rank %o\n",choice,sign,e,rk;
                if rk eq e+Rank(moved) then Append(~ranks,e);break;end if;
            end for;
        end if;
        if e eq 0 or ranks[#ranks] ne e then
            missed+:=1;printf "remaining %o degree %o ranks %o row %o\n",choice,e,ranks,line;
        end if;
    end for;
    printf "exact %o %o row %o\n",count,missed,line;
end while;
quit;
