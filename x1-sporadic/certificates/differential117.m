// Run with Magma V2.30. The base field is F_2 throughout.
// x=(eta(9*tau)/eta(tau))^3  and  y=(eta(13*tau)/eta(tau))^2
// are Hauptmoduln of X_0(9) and X_0(13). Both have a simple zero at infinity.
SetColumns(0);
k := GF(2);
precision := 420;
L<q> := LaurentSeriesRing(k, precision);
q := q + O(q^precision);
euler := &*[1-q^i : i in [1..precision-1]];
function Hauptmodul(n, power)
    euler_n := &*[1-q^(n*i) : i in [1..precision div n]];
    return q*(euler_n/euler)^power;
end function;
xq := Hauptmodul(9,3);
yq := Hauptmodul(13,2);
assert Valuation(xq) eq 1 and Valuation(yq) eq 1;

// Since E_4 = 1 mod 2, 1/j = Delta = q prod (1-q^n)^24.
// Each rational function below has the known degree of the j-map.
// The difference of two degree-d functions has pole degree at most 2d.
R<z> := PolynomialRing(k);
A := z^3+z^2+z;
B := z^12+z^8+z^4+1;
C := z;
D := z^14+z^12+z^2+1;
assert GCD(A,B) eq 1 and Max(Degree(A),Degree(B)) eq 12;
assert GCD(C,D) eq 1 and Max(Degree(C),Degree(D)) eq 14;
Jq := q*euler^24;
assert Valuation(Jq-Evaluate(A,xq)/Evaluate(B,xq)) gt 24;
assert Valuation(Jq-Evaluate(C,yq)/Evaluate(D,yq)) gt 28;

// As 9 and 13 are coprime, their cyclic subgroups determine the cyclic
// subgroup of order 117. The normalization of this fiber product is X_0(117).
S<x> := FunctionField(k);
P<Y> := PolynomialRing(S);
relation := Evaluate(A,x)*Evaluate(D,Y)-Evaluate(B,x)*Evaluate(C,Y);
assert Degree(relation) eq 14 and IsIrreducible(relation);
F<v> := ext<S | relation>;
u := F!x;
assert Genus(F) eq 11;
j := Evaluate(B,u)/Evaluate(A,u);
assert j eq Evaluate(D,v)/Evaluate(C,v);

// Evaluate the coefficient of a differential, relative to du, at infinity.
function EvaluateBase(a, xq)
    a := S!a;
    return Evaluate(Numerator(a),xq)/Evaluate(Denominator(a),xq);
end function;
function AtCusp(a, xq, yq)
    coefficients := Eltseq(a);
    return &+[EvaluateBase(coefficients[i],xq)*yq^(i-1)
              : i in [1..#coefficients] | coefficients[i] ne 0];
end function;
basis := BasisOfHolomorphicDifferentials(F);
assert #basis eq 11;
du := Differential(u);
images := [AtCusp(w/du,xq,yq)*Derivative(xq) : w in basis];
terms := 70;
assert terms gt 2*Genus(F)-2;
assert AbsolutePrecision(Jq) gt 28;
assert forall{h : h in images | AbsolutePrecision(h) ge terms};
M := Matrix(k, [[Coefficient(h,i) : i in [0..terms-1]] : h in images]);
assert Rank(M) eq 11;

// The optimal curve has Manin constant 1 (Cremona, Elliptic curve data,
// Theorem 1(a)). Thus the pullback of its Neron differential is f(q)dq/q.
E := EllipticCurve("117a1");
assert Conductor(E) eq 117;
fq := qExpansion(ModularForm(E),terms+2);
fvector := Vector(k,[Coefficient(fq,i) : i in [1..terms]]);
coefficients := Solution(M,fvector);
assert coefficients*M eq fvector;
omega := &+[basis[i] : i in [1..#basis] | coefficients[i] ne 0];

// In characteristic two the supersingular locus is precisely j=0.
supersingular := Zeros(j);
assert #supersingular eq 7;
assert forall{p : p in supersingular | Degree(p) eq 2};
zero_divisor := Divisor(omega);
places := Support(zero_divisor);
assert #places eq 3 and Degree(zero_divisor) eq 20;
assert forall{p : p in places | p in supersingular and Degree(p) eq 2};
orders := Sort([Valuation(omega,p) : p in places]);
assert orders eq [2,2,6];
assert &+[Degree(p) : p in places] eq 6;
print "Supersingular places: seven of degree 2.";
print "Differential zeros: three of degree 2, of orders", orders;
quit;
