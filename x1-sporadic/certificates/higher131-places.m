SetQuitOnError(true);
load "higher131-forms.m";
cusps := cuspConditions(); assert Rank(cusps) eq 25;
rows := [
    <14,24673>,
    <14,22515>,
    <14,24139>,
    <14,17101>,
    <14,23007>,
    <15,51661>,
    <15,49573>,
    <15,59245>,
    <15,64223>,
    <15,45749>,
    <16,67897>,
    <16,84625>,
    <16,113819>,
    <16,103403>,
    <16,97659>,
    <16,99355>,
    <16,123443>,
    <16,114907>,
    <16,115945>,
    <16,114227>,
    <16,96751>,
    <16,114191>,
    <17,209695>,
    <17,137381>,
    <17,215717>,
    <17,259501>,
    <17,249319>,
    <17,154719>,
    <17,222251>,
    <17,257387>,
    <17,203013>,
    <17,141083>,
    <17,216013>,
    <17,134573>,
    <17,244835>,
    <17,161223>,
    <17,181255>,
    <18,504165>,
    <18,277663>,
    <18,309265>,
    <18,377083>,
    <18,303161>,
    <18,356695>,
    <18,272335>,
    <18,378895>,
    <18,443511>,
    <18,524145>,
    <18,325911>,
    <18,447913>,
    <18,355697>,
    <18,332395>,
    <18,355897>,
    <18,473083>,
    <18,375415>,
    <18,410283>,
    <18,503009>,
    <18,468549>
];
for entry in rows do
    b, code := Explode(entry);
    P,images,tr,r := markedPoint(b,code);
    poly := R!Intseq(code,2);
    simple := curveConditions(P,false);
    double := curveConditions(3*P,true);
    neighbours := ZeroMatrix(k,0,width);
    same := false;
    for pair in images do
        Q,t := Explode(pair);
        same := same or Evaluate(poly,jInvariant(Curve(Q))) eq 0;
        neighbours := VerticalJoin(neighbours,curveConditions(3*Q,false));
    end for;
    all := VerticalJoin(VerticalJoin(cusps,simple),VerticalJoin(double,neighbours));
    ranks := [Rank(cusps),Rank(simple),Rank(double),Rank(neighbours),Rank(all)];
    printf "case %o %o trace %o eigenvalue %o self %o ranks %o %o %o %o %o expected %o\n",
       b,code,tr,Integers()!r,same,ranks[1],ranks[2],ranks[3],ranks[4],ranks[5],25+7*b;
end for;
quit;
