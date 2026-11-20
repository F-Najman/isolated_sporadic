"""Principal parts and formal-immersion certificates for the exceptional cusp divisors."""

from collections import defaultdict


def records(text,prefix):
    return [line for line in text.splitlines() if line.startswith(prefix)]


def pencil_layout(divisors,rank_text):
    dimensions={}
    directions=defaultdict(list)
    for line in records(rank_text,"fibre "):
        v=line.split()
        i,d,r,s=map(int,[v[1],v[3],v[6],v[9]])
        assert d==95 and r==divisors[i][0] and r<=s<=95 and i not in dimensions
        dimensions[i]=96-s
    assert set(dimensions)==set(range(len(divisors)))
    for line in records(rank_text,"poles "):
        i,*m=map(int,line[6:].split())
        assert len(m)==13 and 0<sum(m)<=95 and all(x<=y for x,y in zip(m,divisors[i][1:]))
        directions[i].append(tuple(m))
    for i,e in dimensions.items():
        assert len(directions[i])==(3**(e-1)-1)//2
    primitive=sorted({m for parts in directions.values() for m in parts})
    return primitive,directions


def cusp_and_ordinary(divisors,primitive,directions,pencil_text):
    dimensions={}
    for line in records(pencil_text,"cuspidal representative "):
        v=line.split()
        i,d,e=map(int,[v[2],v[4],v[6]])
        assert d==sum(primitive[i]) and e>=2
        dimensions[i]=e
    assert set(dimensions)==set(range(len(primitive)))
    sections=defaultdict(list)
    ordinary={}
    endpoints=defaultdict(set)
    for line in records(pencil_text,"section "):
        lhs,rhs=line.split(" cusps ")
        _,i,idx,_,d=lhs.split()
        i,idx,d=map(int,(i,idx,d))
        caps=tuple(map(int,rhs.strip("[] ").split(", ")))
        assert len(caps)==13 and d+sum(caps)==sum(primitive[i]) and min(caps)>=0
        sections[i].append(idx)
        if d:
            ordinary[idx]=d
        else:
            endpoints[primitive[i]].add(caps)
    for i,e in dimensions.items():
        assert len(sections[i])==(3**e-1)//2 and primitive[i] in endpoints[primitive[i]]
    pure={tuple(v[1:]) for v in divisors}
    for i,v in enumerate(divisors):
        for m in set(directions[i]):
            for endpoint in endpoints[m]:
                pure.add(tuple(a-b+c for a,b,c in zip(v[1:],m,endpoint)))
    return sorted(pure),ordinary


def partition(group,pure):
    widths=[max(pure[i][j] for i in group)-min(pure[i][j] for i in group) for j in range(13)]
    axis=max(range(13),key=widths.__getitem__)
    assert widths[axis]>0
    order=sorted(group,key=lambda i:pure[i][axis])
    mid=len(order)//2
    return order[:mid],order[mid:]


def certify_boxes(pure,dest,magma,write_rows):
    def corner(group):
        return tuple(max(pure[i][j] for i in group) for j in range(13))

    groups=[]
    pending=[list(range(len(pure)))]
    while pending:
        group=pending.pop()
        if sum(corner(group))<=113:
            groups.append(group)
        else:
            pending.extend(partition(group,pure))

    boxes=[]
    refs=[None]*len(pure)
    anti=[]
    floors=[]
    while groups:
        corners=[corner(g) for g in groups]
        write_rows(dest/"corners.dat",[(i,*v) for i,v in enumerate(corners)])
        write_rows(dest/"floor-input.dat",[(i,*(m//5 for m in v)) for i,v in enumerate(corners)])
        atr=records(magma("full_cusps","corner-anti",patterns="corners.dat"),"anti ")
        fl=records(magma("ranks","corner-floors",file="jets100.dat",patterns="floor-input.dat"),"fibre ")
        assert len(atr)==len(fl)==len(groups)
        retry=[]
        for group,v,a,b in zip(groups,corners,atr,fl):
            av=a.split()
            bv=b.split()
            degree=sum(v)
            if int(av[3])!=4*degree or bv[6]!=bv[9]:
                assert len(group)>1,(v,a,b)
                retry.extend(partition(group,pure))
                continue
            j=len(boxes)
            boxes.append(v)
            for i in group:
                refs[i]=j
            anti.append(f"anti {j} {degree} {av[3]}")
            bv[1]=str(j)
            floors.append(" ".join(bv))
        groups=retry
    return boxes,refs,"\n".join(anti),"\n".join(floors)

