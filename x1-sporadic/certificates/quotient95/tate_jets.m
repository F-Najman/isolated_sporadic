// Weight-two products at one split cusp, with common parameter t (Tate q=t^131).
// In characteristic three s_1(q)=s_3(q), so the constants in x/2+a4 cancel.
L<v>:=LaurentSeriesRing(k,720);
slope:=[L|];
for a in [1..65] do
    pos:=[v^(131*n+a):n in [0..4]]; neg:=[v^(131*n-a):n in [1..5]];
    sx:=&+[w/(1-w)^2:w in pos cat neg];
    sy:=&+[w^2/(1-w)^3:w in pos]-&+[w/(1-w)^3:w in neg];
    m:=(sx/2)/(2*sy+sx);
    assert AbsolutePrecision(m) gt precision;
    Append(~slope,m);
end for;
slope cat:=[-slope[131-a]:a in [66..130]];
s0:=[&+[slope[h*a mod 131]*slope[h*b mod 131]:h in H] where a is pair[1] where b is pair[2]:pair in selected];
J:=Matrix(k,[[Coefficient(f,n):f in s0]:n in [1..precision]]);
assert Submatrix(J,1,1,22,Ncols(J)) eq Submatrix(selectCols(Matrix(k,[Eltseq(r):r in raw[1]]),cols),2,1,22,Ncols(J));
