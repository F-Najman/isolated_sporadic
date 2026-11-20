SetQuitOnError(true);
load "higher131-forms.m";
cusps := cuspConditions(); assert Rank(cusps) eq 25;
function component(b,code,shifts)
    P,images := markedPoint(b,code);
    blocks := [];
    for k0 in shifts do
        mult := 3^k0 mod 131;
        simple := curveConditions(mult*P,false);
        double := curveConditions(3*mult*P,true);
        neigh := ZeroMatrix(k,0,width);
        for pair in images do
            neigh := VerticalJoin(neigh,curveConditions(3*mult*pair[1],false));
        end for;
        Append(~blocks,BasisMatrix(RowSpace(VerticalJoin(VerticalJoin(simple,double),neigh))));
    end for;
    return blocks;
end function;
pairs := [
    <3,13,14,19963>,
    <3,13,14,27717>,
    <8,451,9,971>,
    <8,357,9,865>,
    <3,13,14,30405>,
    <3,13,15,34591>,
    <9,967,9,865>,
    <3,13,15,56373>,
    <3,13,15,35825>,
    <3,13,15,50915>,
    <9,617,9,865>,
    <8,419,10,1591>
];
cache := AssociativeArray();
for input in pairs do
    d1,j1,d2,j2 := Explode(input);
    if not IsDefined(cache,j1) then cache[j1]:=component(d1,j1,[0..12]); end if;
    if not IsDefined(cache,j2) then cache[j2]:=component(d2,j2,[0..12]); end if;
    bs := VerticalJoin(cusps,cache[j1][1]);
    ranks := [Rank(VerticalJoin(bs,block))-Rank(cusps):block in cache[j2]];
    printf "pair %o %o %o %o ranks %o expected %o\n",d1,j1,d2,j2,ranks,7*(d1+d2);
end for;
quit;
