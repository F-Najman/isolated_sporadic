load "forms.m";
input:=Open(file,"r");
cache:=AssociativeArray();
while true do
    line:=Gets(input);if IsEof(line) then break;end if;
    v:=[StringToInteger(s):s in Split(line," ")];
    factors:=[<v[5*i-4],v[5*i-3],v[5*i-2],v[5*i-1],v[5*i]>:i in [1..#v div 5]];
    Sort(~factors,func<a,b|b[1]*b[3]-a[1]*a[3]>);
    blocks:=[* *];degree:=0;
    for term in factors do
        f,m,r,j,jw:=Explode(term);degree+:=f*m*r;tag:=<f,r,j,jw>;
        if IsDefined(cache,tag) then A:=cache[tag];
        else A:=primaryBlock(f,r,j,jw);
            if f*r le 21 then cache[tag]:=A;end if;
        end if;
        for copy in [1..m] do
            choices:=#blocks eq 0 select [0] else [0..(65 div r)-1];
            Append(~blocks,[*BasisMatrix(RowSpace(A*rot[a+1])):a in choices*]);
        end for;
    end for;
    num:=&*[#a:a in blocks];fail:=0;
    for c in [0..num-1] do
        ix:=c;M:=ZeroMatrix(k,0,Ncols(canon));
        for b in blocks do M:=VerticalJoin(M,b[1+ix mod #b]);ix div:=#b;end for;
        if Rank(M) ne 7*degree then fail+:=1;end if;
    end for;
    printf "checked %o %o %o row %o\n",degree,num,fail,line;
end while;
quit;
