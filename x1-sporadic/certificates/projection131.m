// Run from the certificates directory, with fricke131-forms.m in that directory.
SetQuitOnError(true);
SetColumns(0);
load "fricke131-forms.m";
k := Rationals();
K<X> := FunctionField(k);
S<Z> := PolynomialRing(K);
A := R(X,0);
B := R(X,1)-A;
F<z> := ext<K | (Z^2-A)^2+B*(Z^2-A)-B^2*(X^3-X^2+X)>;
x := F!X;
y := (z^2-F!A)/(F!B);
assert y^2+y eq x^3-x^2+x and z^2 eq R(x,y);
assert Genus(F) eq 11;
cusps := Poles(x);
assert #cusps eq 2 and &and[Degree(P) eq 1 : P in cusps];
c := 1*[P : P in cusps | Evaluate(z/x^5,P) eq 1][1];
c1 := 1*[P : P in cusps | Evaluate(z/x^5,P) eq -1][1];
T := c-c1;
branch := Zeros(z);
assert Sort([Degree(P) : P in branch]) eq [5,15];

function SameResidueField(P)
    L, ev := ResidueClassField(P);
    u := ev(x);
    v := ev(y);
    b := Degree(P);
    if b eq 1 then return true; end if;
    // A primitive linear combination exists among these b+1 choices.
    return &or[Degree(MinimalPolynomial(u+i*v)) eq b : i in [0..b]];
end function;

function FrobeniusWitness(P)
    L, ev := ResidueClassField(P);
    j := (Jp(ev(x),ev(y)) + ev(z)*Jm(ev(x),ev(y)))/2;
    f := MinimalPolynomial(j);
    assert Degree(f) eq Degree(P);
    T131<t> := PolynomialRing(GF(131));
    for ell in [263,787,1049,2621] do
        assert IsPrime(ell) and ell mod 131 eq 1;
        if &or[Denominator(a) mod ell eq 0 : a in Coefficients(f)] then continue; end if;
        residue := GF(ell);
        for root in Roots(ChangeRing(f,residue)) do
            jbar := root[1];
            if root[2] ne 1 or jbar in {residue!0,residue!1728} then continue; end if;
            curve := EllipticCurveWithjInvariant(jbar);
            a := TraceOfFrobenius(curve);
            eigenvalues := Roots(t^2-a*t+ell);
            if #eigenvalues eq 0 then continue; end if;
            orders := {Order(r[1]) div GCD(Order(r[1]),2) : r in eigenvalues};
            if orders eq {65} then return <ell, Integers()!jbar, a>; end if;
        end for;
    end for;
    error "Missing Frobenius witness.";
end function;

found := {@ @};
degreeCounts := [];
twists := [];
for m in [10,11] do
    // For the exceptional degree-11 pencils, also remove the fibre at O.
    for a in [0..m-10] do
        exponents := a eq 0 select [-32..32] else twists;
        for h in exponents do
            D := m*c+h*T-a*(c+c1);
            space, section := RiemannRochSpace(D);
            if m eq 11 and a eq 0 and Dimension(space) gt 1 then
                Append(~twists,h);
            end if;
            if Dimension(space) eq 0 then continue; end if;
            effective := D+Divisor(section(space.1));
            assert IsEffective(effective);
            for P in Support(effective) do
                if Degree(P) gt 9 or P in cusps or P in found or not SameResidueField(P) then continue; end if;
                Include(~found,P);
                Append(~degreeCounts,Degree(P));
                if P in branch then
                    assert Degree(P) eq 5;
                    print "Fricke fixed", Degree(P);
                else
                    witness := FrobeniusWitness(P);
                    print "degree, representative, Frobenius witness", Degree(P), <m,a,h>, witness;
                end if;
            end for;
        end for;
    end for;
end for;
print "additional degree-11 representatives", twists;
assert Sort(degreeCounts) eq [5,8,8,8,8,9,9,9,9];
print "degrees", Sort(degreeCounts);
load "ten131.m";
