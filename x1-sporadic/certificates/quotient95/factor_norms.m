load "artin.m";
mon:=[K!1,0] cat [IsEven(i) select x^(i div 2) else x^((i-3) div 2)*y:i in [2..180]];
fd:=Open(file,"r");
while true do
    l:=Gets(fd);if IsEof(l) then break;end if;
    vv:=[StringToInteger(a):a in Split(l," ")]; idx:=vv[1];d:=vv[2];
    aa:=&+[vv[2+i]*mon[i]:i in [1..d+1]];bb:=&+[vv[d+3+i]*mon[i]:i in [1..d+1]];
    h:=aa+z*bb;
    if Valuation(h,c0) ne -d then h:=aa-z*bb;end if;
    assert Valuation(h,c0) eq -d and Poles(h) eq [c0];
    zeros:=Zeros(h);pat:=[];ok:=true;unamb:=true;
    for W in zeros do
        r:=Order(im(dl(W)));m:=Valuation(h,W);
        Append(~pat,<Degree(W),m,r>);if m mod r ne 0 then ok:=false;end if;
        if r mod 5 eq 0 and r mod 13 ne 0 then unamb:=false;end if;
    end for;
    if ok then
        assert unamb;
        for W in zeros do
            if W in cusps then
                assert W eq ci and Valuation(h,W) eq 65;
                printf "inert %o 65\n",idx;continue;
            end if;
            jj:=MinimalPolynomial(Evaluate(jfun,W),k);jd:=MinimalPolynomial(Evaluate(Jp(x,y)-jfun,W),k);
            assert Degree(jj) eq Degree(W) and Valuation(h,W) eq 1 and im(dl(W)) eq Art!0;
            printf "data %o %o %o %o\n",idx,Degree(W),Seqint([Integers()!a:a in Coefficients(jj)],3),Seqint([Integers()!a:a in Coefficients(jd)],3);
        end for;
    end if;
    Sort(~pat);
    printf "norm %o %o pass %o unamb %o places",idx,d,ok,unamb;
    for t in pat do printf " %o %o %o",t[1],t[2],t[3]; end for;printf "\n";
end while;
quit;
