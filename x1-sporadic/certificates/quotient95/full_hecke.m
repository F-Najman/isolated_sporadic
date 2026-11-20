useFullCurve:=true;
load "../branches131/ternary_forms.m";
load "../branches131/ternary_cols.m";

fromFile:=Open(file,"r");groups:=AssociativeArray();degrees:=AssociativeArray();
while true do
   line:=Gets(fromFile);if IsEof(line) then break;end if;
   v:=Split(line," ");if #v eq 0 or v[1] ne "data" then continue;end if;
   label,f,j,jw:=Explode([StringToInteger(s):s in v[2..5]]);
   comps:=marked(f,j,jw);assert #comps eq 1;
   A:=mblock(comps[1],1);
   printf "right-side conditions %o %o rank %o\n",label,f,Rank(A*canon);
   if not IsDefined(groups,label) then groups[label]:=[* *];degrees[label]:=0;end if;
   degrees[label]+:=f;
   g:=groups[label];Append(~g,A);groups[label]:=g;
end while;
for label in Sort(Setseq(Keys(groups))) do
    g:=groups[label];rows:=g[1]*canon;ranks:=[];
    if #g eq 1 then Append(~ranks,Rank(rows));
    else
        assert #g eq 2;
        for i in [0..64] do Append(~ranks,Rank(VerticalJoin(rows,g[2]*rot[i+1])));end for;
    end if;
    assert Seqset(ranks) eq {5*degrees[label]};
    printf "full Hecke fibre %o relative choices %o ranks %o\n",label,#ranks,Sort(Setseq(Seqset(ranks)));
end for;
quit;
