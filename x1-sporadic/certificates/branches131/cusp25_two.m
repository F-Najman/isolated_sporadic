load "binary25_forms.m";
anti:=(IdentityMatrix(k,#cols)-delta)*canon;
jj:=[[ (Vector(k,[v[l]:l in cols]))*anti :v in a]:a in raw];
piv:=pivots(Transpose(canon));
D:=Matrix(k,[Eltseq(canon[i]):i in piv])^-1*Matrix(k,[Eltseq((delta*canon)[i]):i in piv]);
assert canon*D eq delta*canon;
inert:=BasisMatrix(RowSpace(atI*(IdentityMatrix(k,Ncols(canon))-D)));
assert Rank(inert) eq 12;
load "cusp25_common.m";
