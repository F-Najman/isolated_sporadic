// Run after the Fricke model and FrobeniusWitness have been defined.
c0 := 1*[P : P in cusps | Evaluate(z/x^5,P) eq -1][1];
cinf := 1*[P : P in cusps | Evaluate(z/x^5,P) eq 1][1];
fibreO := c0+cinf;
fricke := hom<F -> F | -z>;
assert R(0,0) eq -19 and Jp(0,0) eq -1769472 and Jm(0,0) eq 0;
tenPoints := {@ @};
for h in [-32..32] do
    D := 10*c0+h*(c0-cinf);
    V, section := RiemannRochSpace(D+fibreO);
    if Dimension(V) ne 2 then continue; end if;
    u := section(V.1);
    v := section(V.2);
    W := (u*fricke(v)-fricke(u)*v)/z;
    if W eq 0 then continue; end if;
    fibres := [];
    if W in k then
        Append(~fibres,fibreO);
    else
        places := Zeros(W);
        rationalValues := {};
        for P in places do
            L, ev := ResidueClassField(P);
            if ev(x) in k and ev(y) in k then
                Include(~rationalValues,<k!ev(x),k!ev(y)>);
            end if;
        end for;
        for coords in rationalValues do
            // No rational two-torsion on E, so x-x(A) is a parameter at A.
            assert 2*coords[2]+1 ne 0;
            fibre := 0*c0;
            for P in places do
                L, ev := ResidueClassField(P);
                if ev(x) eq coords[1] and ev(y) eq coords[2] then
                    fibre +:= Valuation(x-coords[1],P)*P;
                end if;
            end for;
            assert Degree(fibre) eq 2;
            Append(~fibres,fibre);
        end for;
    end if;
    for fibre in fibres do
        V0, s0 := RiemannRochSpace(D+fibreO-fibre);
        assert Dimension(V0) eq 1;
        effective := D+fibreO-fibre+Divisor(s0(V0.1));
        assert IsEffective(effective);
        for P in Support(effective) do
            if Degree(P) ne 10 or P in tenPoints or not SameResidueField(P) then continue; end if;
            Include(~tenPoints,P);
            print "degree 10, class, Frobenius witness", h, FrobeniusWitness(P);
        end for;
    end for;
end for;
assert #tenPoints eq 22;
print "degree 10 points with full relative degree", #tenPoints;
