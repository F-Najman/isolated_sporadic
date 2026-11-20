// Check all support patterns by containing them in divisors with full cotangent rank.
function rankOf(orders,extra)
    vec:=&cat[[jj[i][r]:r in [1..orders[i]]]:i in [1..13]];
    if #vec eq 0 then return Rank(extra);end if;
    M:=Matrix(k,[Eltseq(v):v in vec]);
    return Rank(VerticalJoin(extra,M));
end function;
for ss in [0,1] do
    limit:=ss eq 0 select 25 else 12;
    extra:=ss eq 0 select ZeroMatrix(k,0,Ncols(canon)) else inert;
    goal:=Rank(extra);
    cap:=ss eq 0 select 10 else 6;
    while rankOf([0] cat [cap:i in [1..12]],extra) lt 12*cap+goal do cap-:=1;end while;
    assert cap ge 5;
    failures:=0; tests:=1;
    for m in [cap+1..limit] do
        for n in [0..Min(m,limit-m)] do
            r:=Min(n,limit-m-n);
            for empty in [2..13], peak in [2..13] do
                if peak eq empty then continue;end if;
                // If peak and the other ten slots have the same cap,
                // its index does not distinguish the resulting array.
                if r eq n and peak ne (empty eq 2 select 3 else 2) then continue;end if;
                dims:=[i eq 1 select m else (i eq empty select 0 else (i eq peak select n else r)):i in [1..13]];
                tests+:=1; rank:=rankOf(dims,extra);
                if rank ne &+dims+goal then
                    failures+:=1;
                    if failures le 6 then print "rectangle failed",ss,m,n,r,empty,peak,rank,&+dims+goal;end if;
                end if;
            end for;
        end for;
    end for;
    assert failures eq 0;
    printf "cuspidal test at %o inert orbit %o cap %o rectangles %o failures %o\n",Characteristic(k),ss,cap,tests,failures;
end for;
quit;
