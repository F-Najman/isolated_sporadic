// Enumerate cuspidal reduced divisors on the smooth projective curve over Q.
SetColumns(0);
SetQuitOnError(true);
if not assigned n then n:=13; else n:=StringToInteger(n); end if;
if not assigned modeldir then modeldir:="models"; end if;
A<x,y>:=AffineSpace(Rationals(),2);
equation:=eval Read(modeldir cat "/FFFc" cat IntegerToString(n) cat ".txt");
C:=ProjectiveClosure(Curve(A,equation));
g:=Genus(C);
F:=FunctionField(C); x:=F.1; y:=F.2;
r:=(x^2*y-x*y+y-1)/(x*(x*y-1)); s:=(x*y-y+1)/(x*y);
b:=r*s*(r-1); c:=s*(r-1);
cusps:=Poles(jInvariant(EllipticCurve([1-c,-b,-b,0,0])));
base:=Rep({P:P in cusps|Degree(P) eq 1});
allowed:=[P:P in cusps|P ne base and Degree(P) le g];
procedure Count(i,D,~counts)
    if i gt #allowed then
        if Dimension(RiemannRochSpace(D)) eq 1 then counts[1+Degree(D)]+:=1; end if;
        return;
    end if;
    Q:=Divisor(allowed[i]);
    for k in [0..(g-Degree(D)) div Degree(Q)] do Count(i+1,D+k*Q,~counts); end for;
end procedure;
counts:=[0:d in [0..g]];
Count(1,0*Divisor(base),~counts);
expected:=case<n | 13:19, 16:20, 18:21, 20:60, default:0>;
assert &+counts eq expected;
printf "negative level %o genus %o counts %o total %o\n",n,g,counts,&+counts;
quit;
