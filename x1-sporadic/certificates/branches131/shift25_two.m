SetQuitOnError(true);
load "../higher131-forms.m";
cusps:=cuspConditions(); P,imgs,tr,r:=markedPoint(12,7963);
B:=VerticalJoin(cusps,curveConditions(P,false));
B:=VerticalJoin(B,curveConditions(3*P,true));
for im in imgs do B:=VerticalJoin(B,curveConditions(3*im[1],false));end for;
printf "shifted supersingular orbit, residual degree 12, rank %o\n",Rank(B)-Rank(cusps);
assert Rank(B)-Rank(cusps) eq 84;
quit;
