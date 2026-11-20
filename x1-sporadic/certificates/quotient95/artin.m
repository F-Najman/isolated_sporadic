load "../branches131/base_ternary.m";
Cl,rep,dl:=ClassGroup(K);
Art,im:=quo<Cl|[dl(c0),dl(classdiv(0,0,1)),65*Cl.1,65*Cl.2]>;
assert Invariants(Art) eq [65] and Order(im(dl(ci))) eq 65;
assert &and[Degree(P) eq 2 and Order(im(dl(P))) eq 65:P in Zeros(jfun)];
