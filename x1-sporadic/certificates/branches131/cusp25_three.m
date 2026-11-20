load "ternary_forms.m";load "ternary_cols.m";
anti:=(IdentityMatrix(k,#cols)-delta)*canon;
sraw:=splitSeries(25);
jj:=[[Vector(k,[v[l]:l in cols])*anti:v in a[2..#a]]:a in sraw];
L:=GF(3,65);zeta:=PrimitiveElement(L)^((#L-1) div 131);
m:=[L|];first:=[L|];
for a in [1..130] do
    w:=zeta^a; x0:=w/(1-w)^2; y0:=w^2/(1-w)^3;
    x1:=w+1/w-2;y1:=1-1/w;
    f0:=x0/2; f1:=x1/2+1; d0:=2*y0+x0;d1:=2*y1+x1;
    Append(~m,f0/d0);Append(~first,f1/d0-f0*d1/d0^2);
end for;
v:=[&+[m[p[1]*h mod 131]*first[p[2]*h mod 131]+first[p[1]*h mod 131]*m[p[2]*h mod 131]:h in H]:p in selected];
inert:=BasisMatrix(RowSpace(absoluteRows(v)*anti));
assert Rank(inert) eq 12;
load "cusp25_common.m";
