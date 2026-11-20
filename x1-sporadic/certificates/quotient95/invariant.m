load "artin.m";
assert Invariants(Art) eq [65] and Order(dl(5*(ci-c0))) eq 13;
for b in [1..7] do
    for k in [0..12] do
        D:=b*c0+5*k*(ci-c0);V,rr:=RiemannRochSpace(D);if Dimension(V) eq 0 then continue;end if;
        for a in Points(ProjectiveSpace(GF(3),Dimension(V)-1)) do
            eff:=D+Divisor(rr(V!Eltseq(a)));supp:=Support(eff);if c0 in supp then continue;end if;
            lift:=&and[(Order(im(dl(W))) mod 5 ne 0) or (Valuation(eff,W) mod 5 eq 0):W in supp];
            printf "invariant quotient %o parameter %o lifts %o places",b,k,lift;
            for W in supp do printf " %o %o %o %o", Degree(W),Valuation(eff,W),Order(im(dl(W))),W in cusps;end for;
            printf "\n";
        end for;
    end for;
end for;
quit;
