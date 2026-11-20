load "ternary_forms.m"; load "ternary_cols.m";
diamond:=2;fibreO:=ZeroMatrix(k,0,Ncols(canon)); chooseCols:=selectCols;absRows:=absoluteRows;
function simpleEval(P) return absoluteRows(productValues(tangentSlopes(P),1))*canon;end function;
function allMarked(f,j,jw)
    return [* <d[3],[* <e[3],e[2]>:e in d[2] *]>:d in marked(f,j,jw) *];
end function;
// The binary branch of the shared routine is never called in characteristic three.
function conditions(v) error "binary only";end function;
moved:=ZeroMatrix(k,0,Ncols(canon));
load "exact25_common.m";
