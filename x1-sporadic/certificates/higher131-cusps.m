SetQuitOnError(true);
load "higher131-forms.m";
cusps := cuspConditions(); assert Rank(cusps) eq 25;
length:=18; rows:=cuspSeries(length);
prym := [[rows[c][n]+rows[1+(c mod 13)][n] : n in [1..length]] : c in [1..13]];
function test(orders)
    extras := &cat[[prym[c][n]:n in [1..orders[c]]]:c in [1..13]];
    return Rank(VerticalJoin(cusps,Matrix(k,[Eltseq(v):v in extras])))-Rank(cusps);
end function;
assert test([0] cat [10:i in [1..12]]) eq 120;
for b in [14..18] do
    for m in [11..b], empty in [2..13] do
        bounds := [i eq 1 select m else (i eq empty select 0 else b-m):i in [1..13]];
        assert test(bounds) eq &+bounds;
    end for;
    printf "cuspidal cotangent test %o verified\n",b;
end for;
quit;
