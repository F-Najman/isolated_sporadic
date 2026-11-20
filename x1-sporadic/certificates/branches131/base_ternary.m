SetQuitOnError(true);
SetColumns(0);
load "../fricke131-forms.m";
k:=GF(3); kX<X>:=FunctionField(k); kXS<S>:=PolynomialRing(kX);
A0:=R(X,0); A1:=R(X,1)-A0;
a:=S^2-A0;
K<z>:=ext<kX|a^2+a*A1-A1^2*(X^3-X^2+X)>;
x:=K!X; y:=(z^2-K!A0)/(K!A1);
assert y^2+y eq x^3-x^2+x and z^2 eq R(x,y) and Genus(K) eq 11;
jfun:=(Jp(x,y)+z*Jm(x,y))/2;
cusps:=Poles(x); assert #cusps eq 2;
c0:=[P:P in cusps|Valuation(jfun,P) eq -131][1];
ci:=[P:P in cusps|Valuation(jfun,P) eq -1][1];
E:=EllipticCurve([k|0,-1,1,1,0]); G:=E![0,0];
function classdiv(b,h,n)
    Q:=n*G;
    if Q eq E!0 then fibre:=c0+ci;
    else
        places := [P:P in Zeros(x-Q[1]) | Evaluate(y,P) eq Q[2]];
        fibre:=&+[Valuation(x-Q[1],P)*P:P in places];
    end if;
    assert Degree(fibre) eq 2;
    return b*c0+h*(c0-ci)+fibre-(c0+ci);
end function;
