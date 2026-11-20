load "../branches131/base_binary25.m";
S0:=x^5+x^4+x+1+(x^2+1)*y;
h:=K!0;
for z in [x^5,x^3*y,x^4,x^2*y,x^3,x*y,x^2,y,x] do
    if Valuation(u+h,c0) eq Valuation(z,c0) then h+:=z;end if;
end for;
if Valuation(u+h,c0) eq 0 then h+:=1;end if;
beta:=(u+h)/S0*Differential(x);
assert IsEffective(Divisor(beta));
assert &and[Valuation(Divisor(Differential(x)),P) eq 0:P in cusps];
assert Valuation(Divisor(beta),c0) eq 11 and Valuation(Divisor(beta),ci) eq 0;
printf "trace differential split vanishing %o nonsplit vanishing %o\n",Valuation(Divisor(beta),c0),Valuation(Divisor(beta),ci);
quit;
