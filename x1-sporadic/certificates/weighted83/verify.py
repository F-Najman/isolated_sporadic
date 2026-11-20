"""Replay or audit the modular factorizations for the weighted level-83 search."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

GROUPS=[
    (5,7,11,13,17,19,23,29),
    (31,37,41,43,47,53,59,61),
    (67,71,73,79,89,97,101),
    (103,107,109,113,127,131,137),
    (139,149,151,157,163,167,173),
    (179,181,191,193,197,199,211),
]
BOUND=56


def replay(root,inputs,primes,log,series,threads):
    args=[str(root/'search'),str(root/'units.txt'),str(inputs),str(series),
          ','.join(map(str,primes)),str(BOUND),str(threads)]
    temporary=log.with_suffix('.partial')
    with gzip.open(temporary,'wb',compresslevel=1) as zipped, subprocess.Popen(args,stdout=subprocess.PIPE) as proc:
        for block in iter(lambda:proc.stdout.read(1<<17),b''):
            zipped.write(block)
        if proc.wait():
            raise RuntimeError('norm computation failed')
    temporary.replace(log)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdir',type=Path)
    parser.add_argument('series',type=Path)
    parser.add_argument('--threads',type=int,default=12)
    parser.add_argument('--replay',action='store_true')
    parser.add_argument('--make-record',action='store_true')
    args=parser.parse_args()
    root=args.workdir.resolve()
    ids={}
    entries={}
    counts=Counter()
    with (root/'fibres.txt').open() as f:
        for line in f:
            i,a,b,p,d=map(int,line.split())
            key=i,a,b
            assert key not in entries
            entries[key]=p,d
            if d:
                ids[key]=(1<<(min(BOUND,d//2)+1))-2
            counts['positive' if d else 'all_cuspidal']+=1
    def digest(file):
        with file.open('rb') as f:
            return hashlib.file_digest(f,'sha256').hexdigest()
    record={'counts':counts, 'units_sha256':digest(root/'units.txt'),
            'fibres_sha256':digest(root/'fibres.txt'),'stages':[]}
    record['minimum_positive_degree']=min(d for _,d in entries.values() if d)
    assert record['minimum_positive_degree']==57
    for group,primes in enumerate(GROUPS):
        if not ids: break
        inputs=root/'fibres.txt'
        if group:
            inputs=root/f'stage{group}.in'
            with inputs.open('w') as out:
                for key in sorted(ids):
                    print(*key,*entries[key],file=out)
        log=root/f'stage{group}.gz'
        if args.replay or not log.exists():
            replay(root,inputs,primes,log,args.series,args.threads)
        left={}
        hashes={}
        with gzip.open(log,'rt') as f:
            for line in f:
                head,tail=line.split(' left')
                data=head.split()
                i,a,b,p,d=map(int,data[:5])
                key=i,a,b
                assert key in ids and key not in hashes and entries[key]==(p,d)
                mask=(1<<(min(BOUND,d//2)+1))-2
                normalized=[]
                previous=-1
                for step in data[5:]:
                    prime,pattern=step.split(':')
                    index=primes.index(int(prime))
                    assert index>previous
                    previous=index
                    fac=Counter()
                    for s in pattern.rstrip(',').split(','):
                        j,m=map(int,s.split('^'))
                        assert j>0 and m>0
                        fac[j]+=m
                    assert sum(j*m for j,m in fac.items())==d
                    possible=1
                    for j,m in fac.items():
                        for _ in range(m): possible|=possible<<j
                    mask&=possible
                    normalized.append((int(prime),sorted(fac.items())))
                assert [j for j in range(1,BOUND+1) if mask>>j&1]==list(map(int,tail.split()))
                if mask & ids[key]:
                    left[key]=mask & ids[key]
                hashes[key]=hashlib.sha256(repr((key,p,d,normalized)).encode('ascii')).digest()
        assert hashes.keys()==ids.keys()
        h=hashlib.sha256()
        for key in sorted(hashes):h.update(hashes[key])
        record['stages'].append(dict(tested=len(ids),remaining=len(left),sha256=h.hexdigest()))
        print('stage',group,'tested',len(ids),'remaining',len(left),flush=True)
        ids=left
    record['remaining']=[[*key,*entries[key],[j for j in range(1,BOUND+1) if m>>j&1]]
                         for key,m in sorted(ids.items())]
    if args.make_record:
        (root/'record.json').write_text(json.dumps(record,indent=1)+'\n')
    else:
        assert record==json.loads((root/'record.json').read_text())
    print('total signed and nonsigned fibres',sum(counts.values()),'pending',len(ids),flush=True)


if __name__=='__main__':
    main()
