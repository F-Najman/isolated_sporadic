useFullCurve:=true;
load "../branches131/ternary_forms.m";
load "../branches131/ternary_cols.m";
precision:=100;
load "tate_jets.m";
anti:=(IdentityMatrix(k,#cols)-delta^13)*canon;
pick:=pivots(anti);assert #pick eq 520;
U:=Submatrix(anti,[1..Nrows(anti)],pick);
M:=J; jj:=[];
for i in [0..51] do
    Append(~jj,M*U); M:=M*delta;
end for;
fd:=Open(patterns,"r");
while true do
    line:=Gets(fd);if IsEof(line) then break;end if;
    v:=[StringToInteger(a):a in Split(line," ")];D:=v[2..14];W:=ZeroMatrix(k,0,520);
    for i in [0..12] do
        if D[i+1] eq 0 then continue;end if;
        for a in [0..3] do
            W:=VerticalJoin(W,Submatrix(jj[i+13*a+1],1,1,D[i+1],520));
        end for;
    end for;
    r:=Rank(W); printf "anti %o %o %o\n",v[1],&+D,r;
end while;
quit;