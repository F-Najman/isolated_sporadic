load "ternary_forms.m";
// The selected weight-two products are determined by these cusp expansions.
// Translate columns by twice the marking and solve for the diamond matrix.
function pivots(B)
    e:=BasisMatrix(RowSpace(B));
    return [Min([j:j in [1..Ncols(e)]|e[i,j] ne 0]):i in [1..Nrows(e)]];
end function;
cl:=[]; signs:=[];
for P in selected do
    a:=2*P[1] mod 131; b:=2*P[2] mod 131;
    Append(~signs,((a gt 65) xor (b gt 65)) select -1 else 1);
    u:=multIndex(a);v:=multIndex(b);
    Append(~cl,Index(pairs,<Min(u,v),Max(u,v)>));
end for;
base:=selectCols(stack,cols);
translated:=selectCols(stack,cl)*DiagonalMatrix(k,signs);
pos:=pivots(Transpose(base));
delta:=Matrix(k,[Eltseq(base[i]):i in pos])^-1*Matrix(k,[Eltseq(translated[i]):i in pos]);
assert base*delta eq translated and delta^13 eq IdentityMatrix(k,#cols);
canon:=Transpose(BasisMatrix(Nullspace(Transpose(cc))));
rot:=[delta^i*canon:i in [0..12]];
cache:=AssociativeArray();
input:=Open(file,"r");
while true do
    line:=Gets(input); if IsEof(line) then break;end if;
    ints:=[StringToInteger(t):t in Split(line," ")];
    types:=[<ints[4*i-3],ints[4*i-2],ints[4*i-1],ints[4*i]>:i in [1..#ints div 4]];
    Sort(~types,func<a,b|b[1]-a[1]>);
    b:=&+[v[1]*v[2]:v in types];
    options:=[* *]; labels:=[* *];miss:=false;
    for typ in types do
        f,m,j,jw:=Explode(typ); key:=<f,j,jw>;
        if not IsDefined(cache,key) then
            vals:=marked(f,j,jw);
            mats:=[mblock(a,1):a in vals];
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
        idx:=number; stacked:=ZeroMatrix(k,0,Ncols(canon)); tags:=[];
        for col in [1..#options] do
            z:=1+idx mod #options[col];idx:=idx div #options[col];
            stacked:=VerticalJoin(stacked,options[col][z]);
            Append(~tags,labels[col][z]);
        end for;
        rank:=Rank(stacked);
        if rank ne 5*b then
            failed+:=1;
            printf "hard %o %o", b,rank;
            for t in tags do printf " %o %o %o %o %o",t[1],t[2],t[3],t[4],t[5];end for;
            printf "\n";
        end if;
    end for;
    printf "checked %o %o %o %o\n",b,combinations,failed,line;
end while;
quit;
