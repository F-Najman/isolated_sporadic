load "binary25_forms.m";
cache:=AssociativeArray();
input:=Open(file,"r");
while true do
    line:=Gets(input); if IsEof(line) then break;end if;
    ints:=[StringToInteger(t):t in Split(line," ")];
    use:=ints[1];Remove(~ints,1);
    types:=[<ints[4*i-3],ints[4*i-2],ints[4*i-1],ints[4*i]>:i in [1..#ints div 4]];
    Sort(~types,func<a,b|b[1]-a[1]>);
    b:=&+[v[1]*v[2]:v in types];
    options:=[* *]; labels:=[* *];miss:=false;
    for typ in types do
        f,m,j,jw:=Explode(typ); key:=<f,j,jw>;
        if not IsDefined(cache,key) then
            vals:=markedPoints(f,j,jw);
            mats:=[block3(a):a in vals];
            if f le 13 then cache[key]:=mats;end if;
        else mats:=cache[key]; end if;
        for rep in [1..m] do
            choices:=#options eq 0 select [0] else [0..12];
            blocks:=[* *]; tags:=[];
            for z in [1..#mats],i in choices do
                Append(~blocks,BasisMatrix(RowSpace(mats[z]*rot[i+1])));
                Append(~tags,<f,j,jw,z,i>);
            end for;
            Append(~options,blocks);Append(~labels,tags);
        end for;
    end for;
    combinations:=&*[#v:v in options]; failed:=0;
    for number in [0..combinations-1] do
        idx:=number; stacked:=use eq 1 select fibreO else ZeroMatrix(k,0,Ncols(canon)); tags:=[];
        for col in [1..#options] do
            z:=1+idx mod #options[col];idx:=idx div #options[col];
            stacked:=VerticalJoin(stacked,options[col][z]);
            Append(~tags,labels[col][z]);
        end for;
        rank:=Rank(stacked);
        if rank ne 4*b+Rank(fibreO)*use then
            failed+:=1;
            printf "hard %o %o %o", b,use,rank;
            for t in tags do printf " %o %o %o %o %o",t[1],t[2],t[3],t[4],t[5];end for;
            printf "\n";
        end if;
    end for;
    printf "checked %o %o %o %o\n",b,combinations,failed,line;
end while;
quit;
