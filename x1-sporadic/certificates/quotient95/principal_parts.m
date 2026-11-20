SetQuitOnError(true);SetColumns(0);k:=GF(3);
input:=Open(file,"r");J:=[];
for i in [1..13] do
    A:=[];
    for n in [1..100] do Append(~A,[StringToInteger(t):t in Eltseq(Gets(input))]);end for;
    Append(~J,Matrix(k,A));
end for;
for line in Split(Read(patterns),"\n") do
    v:=[StringToInteger(x):x in Split(line," ")];D:=v[2..14];
    A:=ZeroMatrix(k,0,131);
    for i in [1..13] do
        if D[i] gt 0 then A:=VerticalJoin(A,Submatrix(J[i],1,1,D[i],131));end if;
    end for;
    r:=Rank(Submatrix(A,1,1,Nrows(A),130));N:=Nullspace(A);dim:=Dimension(N);
    printf "fibre %o degree %o trace-zero rank %o full rank %o\n",v[1],&+D,r,&+D-dim;
    if dim eq 0 then continue;end if;
    // Each line is the principal part of a function, up to a scalar.
    // Residue duality identifies this space with L(D) modulo constants.
    for a in Points(ProjectiveSpace(k,dim-1)) do
        pp:=Eltseq(Vector(k,Eltseq(a))*BasisMatrix(N));offset:=0;
        printf "poles %o",v[1];
        for i in [1..13] do
            ord:=Max([0] cat [j:j in [1..D[i]]|pp[offset+j] ne 0]);
            printf " %o",ord;offset+:=D[i];
        end for;
        printf "\n";
    end for;
end for;
quit;
