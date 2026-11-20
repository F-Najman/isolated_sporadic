load "../branches131/ternary_forms.m";
load "../branches131/ternary_cols.m";
precision:=260;
load "tate_jets.m";
anti:=(IdentityMatrix(k,#cols)-delta)*canon;
pick:=pivots(anti);assert #pick eq 120;
U:=Submatrix(anti,[1..Nrows(anti)],pick);
tseries:=[[StringToInteger(a):a in Split(l," ")]:l in Split(Read("fricke-series.txt"),"\n")];
x:=&+[L!(k!tseries[2][i])*v^(i-3):i in [1..340]]+O(v^338);
y:=&+[L!(k!tseries[3][i])*v^(i-4):i in [1..340]]+O(v^337);
z:=&+[L!(k!tseries[4][i])*v^(i-11):i in [1..340]]+O(v^330);
f:=v*Derivative(x)/(2*y+1); // dx/(2y+1) as a multiple of dv/v.
assert Valuation(f) eq 1;
basisB:=[x^i*f/z:i in [0..5]] cat [x^i*y*f/z:i in [0..3]];
B0:=Matrix(k,[[Coefficient(b,n):b in basisB]:n in [1..260]]);
assert Ncols(B0) eq 10;
trace:=Matrix(k,260,1,[Coefficient(f,n):n in [1..260]]);
A:=J;matrices:=[];out:=Open(file,"w");
for i in [0..12] do
    M:=HorizontalJoin(A*U,B0); Append(~matrices,M);
    for n in [1..260] do
        fprintf out,"%o\n",&cat[IntegerToString(Integers()!t):t in Eltseq(M[n]) cat [trace[n,1]]];
    end for;
    A:=A*delta;
end for;
stack:=VerticalJoin(matrices);
assert Rank(stack) eq 130 and Rank(HorizontalJoin(stack,VerticalJoin([trace:i in [0..12]]))) eq 131;
print "trace-zero space 130, full space 131, jets 13 by 260";
quit;
