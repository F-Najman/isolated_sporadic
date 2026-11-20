load "forms.m";
// x(P) distinguishes the marking modulo sign on the normalized ordinary model.
function placeKey(P)
    L:=BaseRing(Curve(P)); vals:=[jInvariant(Curve(P)),P[1]];
    f:=Lcm([Degree(MinimalPolynomial(v,k)):v in vals]);F:=GF(2,f);Embed(F,L);
    a:=[F!v:v in vals];vv:=[];
    for i in [1..f] do
        Append(~vv,[Seqint([Integers()!z:z in Eltseq(v)],2):v in a]);a:=[v^2:v in a];
    end for;
    assert #Seqset(vv) eq f;return Sprint([f] cat Min(vv)),f;
end function;
function jetRows(P,m)
    if m le 2 then return BasisMatrix(RowSpace(ordinaryRows(P,1,m eq 2)*canon));end if;
    E:=Curve(P);L:=BaseRing(E);a1,a2,a3,a4,a6:=Explode(aInvariants(E));
    assert a1 eq 1 and a2 eq 0 and a3 eq 0 and a6+a4^2 eq 1/jInvariant(E);
    Q:=P;pts:=[];for a in [1..130] do Append(~pts,Q);Q+:=P;end for;
    S<u>:=PowerSeriesRing(L,m+1);eta:=a6+a4^2+u;
    xx:=[S!Q[1]:Q in pts];n:=1;
    while n lt m do
        half:=[xx[66*a mod 131]^2:a in [1..130]];
        xx:=[v+eta/v:v in half];n*:=2;
    end while;
    slope:=[(Q[1]^2+a4+Q[2])/Q[1]:Q in pts];jets:=[slope];
    for i in [1..m-1] do
        res:=[Coefficient(xx[2*a mod 131],i):a in [1..130]];
        if IsEven(i) then res:=[res[a]+jets[i div 2+1][a]^2:a in [1..130]];end if;
        Append(~jets,res);
    end for;
    w:=[weights(jets[n+1],1,n eq 0 select 1 else 0):n in [0..m-1]];
    M:=ZeroMatrix(k,0,Ncols(canon));
    for n in [0..m-1] do
        val:=[&+[w[i+1][t[1]]*w[n-i+1][t[2]]:i in [0..n]]:t in pp];
        M:=VerticalJoin(M,absoluteRows(val)*canon);
    end for;
    return BasisMatrix(RowSpace(M));
end function;
procedure addPoint(P,weight,~places,~degs,~D)
    key,f:=placeKey(P);assert weight mod f eq 0;
    if not IsDefined(places,key) then places[key]:=P;degs[key]:=f;end if;
    D[key]:=(IsDefined(D,key) select D[key] else 0)+weight div f;
end procedure;
procedure component(P,img,h,moves,~places,~degs,~opts)
    opts:=[* *];
    for i in moves do
        b:=3^i mod 131;S:=AssociativeArray();
        addPoint(b*P,-h,~places,~degs,~S);
        addPoint(3*b*P,-2*h,~places,~degs,~S);
        addPoint(9*b*P,3*h,~places,~degs,~S);
        for q in img do
            Q,m:=Explode(q);
            addPoint(b*Q,h*m,~places,~degs,~S);
            addPoint(3*b*Q,-h*m,~places,~degs,~S);
        end for;
        Append(~opts,S);
    end for;
end procedure;
function combined(opts,number)
    D:=AssociativeArray();
    for group in opts do
        C:=group[number mod #group+1];number div:=#group;
        for key in Keys(C) do D[key]:=(IsDefined(D,key) select D[key] else 0)+C[key];end for;
    end for;
    return D;
end function;
input:=Open(file,"r"); places:=AssociativeArray();degs:=AssociativeArray();
jets:=AssociativeArray(); opts_cache:=AssociativeArray(); pts_cache:=AssociativeArray();
while true do
    line:=Gets(input);if IsEof(line) then break;end if;
    v:=[StringToInteger(s):s in Split(line," ")];
    terms:=[<v[5*i-4],v[5*i-3],v[5*i-2],v[5*i-1],v[5*i]>:i in [1..#v div 5]];
    Sort(~terms,func<a,b|b[1]*b[3]-a[1]*a[3]>); opts:=[* *];
    for z in terms do
        f,m,r,j,jw:=Explode(z);tag:=<f,r,j,jw>;
        if not IsDefined(pts_cache,tag) then P,img:=genPoint(f,r,j,jw);pts_cache[tag]:=<P,img>;end if;
        P,img:=Explode(pts_cache[tag]);
        for rep in [1..m] do
            moves:=#opts eq 0 select [0] else [0..65 div r-1];id:=<tag,#moves>;
            if not IsDefined(opts_cache,id) then
                A:=[* *];component(P,img,f*r,moves,~places,~degs,~A);opts_cache[id]:=A;
            end if;
            Append(~opts,opts_cache[id]);
        end for;
    end for;
    count:=&*[#v:v in opts];miss:=0;
    for choice in [0..count-1] do
        D:=combined(opts,choice);
        deg:=&+[Integers()|D[key]*degs[key]:key in Keys(D)|D[key] gt 0];
        assert deg gt 0 and &+[Integers()|D[key]*degs[key]:key in Keys(D)] eq 0;
        ranks:=[];
        for sign in [-1,1] do
            B:=ZeroMatrix(k,0,Ncols(canon));
            for key in Keys(D) do
                m:=sign*D[key];if m le 0 then continue;end if;
                if not IsDefined(jets,<key,m>) then jets[<key,m>]:=jetRows(places[key],m);end if;
                B:=VerticalJoin(B,jets[<key,m>]);
            end for;
            Append(~ranks,Rank(B)); if ranks[#ranks] eq deg then break;end if;
        end for;
        if ranks[#ranks] ne deg then
            miss+:=1;printf "survivor %o degree %o ranks %o row %o\n",choice,deg,ranks,line;
        end if;
    end for;
    printf "exact %o %o row %o\n",count,miss,line;
end while;
quit;
