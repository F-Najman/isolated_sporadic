function pivots(B)
    e:=BasisMatrix(RowSpace(B));
    return [Min([j:j in [1..Ncols(e)]|e[i,j] ne 0]):i in [1..Nrows(e)]];
end function;
cl:=[]; signs:=[];
for P in selected do
    a:=2*P[1] mod 131; b:=2*P[2] mod 131;
    Append(~signs,((a gt 65) xor (b gt 65)) select -1 else 1);
    u:=multIndex(a);v:=multIndex(b);
    Append(~cl,Index(pairs,<Min(u,v),Max(u,v)>));
end for;
base:=selectCols(stack,cols);
translated:=selectCols(stack,cl)*DiagonalMatrix(k,signs);
pos:=pivots(Transpose(base));
delta:=Matrix(k,[Eltseq(base[i]):i in pos])^-1*Matrix(k,[Eltseq(translated[i]):i in pos]);
assert base*delta eq translated and delta^deckSize eq IdentityMatrix(k,#cols);
canon:=Transpose(BasisMatrix(Nullspace(Transpose(cc))));
rot:=[delta^i*canon:i in [0..deckSize-1]];
