// Exact cusp contacts of linear combinations of three eta products at level 83.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

using namespace std;
using I = int64_t;
constexpr int N=83, H=41, max_factors=5;
constexpr I prime=83000831, zeta=1495617;
using Indices=array<int,max_factors>;
using Cyclo=array<I,N>;
int factors, height, cutoff, limit;
I xi;

int fold(int x) { x%=N; return min(x,N-x); }
I power(I a,int n) {
    I value=1;
    for(;n;n/=2,a=a*a%prime) if(n%2) value=value*a%prime;
    return value;
}

struct Monomial {
    Indices idx{};
    uint64_t mask=0;
    array<int,H> val{},sign{};
    array<vector<int>,H> jets;
    I sine=1;
};
struct Coeff { int a,b,c; };
struct Case { int i,j,k; Coeff w; int raw,degree; };

vector<int> qpart(const Indices& m,int c,int length) {
    vector<int> series(length);
    series[0]=1;
    for(int j=0;j<factors;++j) {
        int r=c*m[j]%N;
        for(int s:{r,N-r}) for(int n=s;n<length;n+=N)
            for(int t=length-1;t>=n;--t) series[t]-=series[t-n];
    }
    return series;
}

Monomial monomial(Indices m,bool build) {
    Monomial x;
    x.idx=m;
    for(int j=0;j<factors;++j) {
        int i=m[j];
        x.mask|=uint64_t(1)<<(i-1);
        I v=(power(xi,i)+prime-power(xi,2*N-i))%prime;
        x.sine=x.sine*v%prime;
    }
    for(int j=0;j<H;++j) {
        int c=j+1; if(c%2==0) c+=N;
        int sum=0,par=0;
        for(int t=0;t<factors;++t) {
            int r=c*m[t]%N;
            sum+=r*(N-r);
            par+=c*m[t]/N;
        }
        x.val[j]=-sum/(2*N);
        x.sign[j]=par%2 ? -1:1;
        if(build) x.jets[j]=qpart(m,c,limit);
    }
    return x;
}

bool is_zero(const Cyclo& z) {
    return all_of(z.begin()+1,z.end(),[&](I a){return a==z[0];});
}
Cyclo shifted(const Cyclo& v,int e) {
    Cyclo y;
    e=(e%N+N)%N;
    for(int i=0;i<N;++i) y[(i+e)%N]=v[i];
    return y;
}
// Leading factors of E_i at the nonsplit cusp are proportional to
// xi^i-xi^-i, with xi=-zeta^42. The proportionality cancels in the ratio.
vector<Cyclo> other_jets(const Indices& m,int n) {
    vector<Cyclo> s(n+1); s[0][0]=1;
    for(int j=0;j<factors;++j) {
        int i=m[j];
        Cyclo a=shifted(s[0],42*i), b=shifted(s[0],-42*i);
        for(int u=0;u<N;++u) s[0][u]=(i%2?-1:1)*(a[u]-b[u]);
    }
    for(int j=0;j<factors;++j) for(int e:{m[j],-m[j]})
        for(int t=1;t<=n;++t) for(int u=n;u>=t;--u) {
            Cyclo v=shifted(s[u-t],e);
            for(int h=0;h<N;++h) s[u][h]-=v[h];
        }
    return s;
}
int other_contact(const Monomial& a,const Monomial& b,const Monomial& c,
                  Coeff w,int raw) {
    if((w.a*a.sine+w.b*b.sine+w.c*c.sine)%prime) return 0;
    int n=raw/H;
    auto A=other_jets(a.idx,n), B=other_jets(b.idx,n), C=other_jets(c.idx,n);
    for(int i=0;i<=n;++i) {
        Cyclo z;
        for(int k=0;k<N;++k) z[k]=w.a*A[i][k]+w.b*B[i][k]+w.c*C[i][k];
        if(!is_zero(z)) return i;
    }
    // A nonzero combination cannot vanish at this orbit to order >raw/41.
    abort();
}

template<size_t Terms>
bool least_translate(const array<Indices,Terms>& original) {
    if(original[0][0]!=1) return false;
    // The lexicographically least translate must include 1.
    for(const auto& side:original) for(int j=0;j<factors;++j) {
        int m=1;
        while(fold(m*side[j])!=1) ++m;
        auto v=original;
        for(auto& row:v) {
            for(int k=0;k<factors;++k) row[k]=fold(m*row[k]);
            sort(row.begin(),row.begin()+factors);
        }
        sort(v.begin(),v.end());
        if(v<original) return false;
    }
    return true;
}
bool representative(const array<Indices,3>& original) { return least_translate(original); }

struct Summary {
    I triples=0,identities=0,nonzero=0,nonsplit=0;
    map<int,I> histogram;
    vector<Case> selected;
};

void examine(const vector<Monomial>& g,int i,int j,int k,Coeff w,Summary& s) {
    const auto& A=g[i]; const auto& B=g[j]; const auto& C=g[k];
    int raw=0, cancel=0;
    array<int,H> lo;
    for(int h=0;h<H;++h) {
        lo[h]=min({A.val[h],B.val[h],C.val[h]});
        raw+=A.val[h]-lo[h];
    }
    for(int h=0;h<H;++h) {
        int x=A.val[h]-lo[h], y=B.val[h]-lo[h], z=C.val[h]-lo[h];
        int a=w.a*A.sign[h], b=w.b*B.sign[h], c=w.c*C.sign[h];
        if((x==0?a:0)+(y==0?b:0)+(z==0?c:0)!=0) continue;
        auto leading=[&](const vector<int>& u,const vector<int>& v,const vector<int>& t,int start,int stop) {
            int e=start;
            for(;e<stop;++e)
                if((x<=e?a*u[e-x]:0)+(y<=e?b*v[e-y]:0)+(z<=e?c*t[e-z]:0)) break;
            return e;
        };
        int bound=raw-cancel+1;
        int e=leading(A.jets[h],B.jets[h],C.jets[h],1,min(bound,limit));
        if(e==limit && bound>limit) {
            int cusp=h+1;
            e=leading(qpart(A.idx,cusp,bound),qpart(B.idx,cusp,bound),
                      qpart(C.idx,cusp,bound),limit,bound);
        }
        cancel+=e;
        if(cancel>raw) { ++s.identities; return; }
    }
    raw-=cancel;
    int nonsplit=other_contact(A,B,C,w,raw), degree=raw-H*nonsplit;
    ++s.nonzero;
    ++s.histogram[degree];
    if(nonsplit) ++s.nonsplit;
    if(degree && degree<=cutoff) s.selected.push_back({i,j,k,w,raw,degree});
}

void fill(int pos,int least,int square,Indices indices,array<vector<Monomial>,3>& groups) {
    if(pos==factors) {
        if(square<=2) groups[square].push_back(monomial(indices,true));
        return;
    }
    for(int i=least;i<=H;++i) {
        indices[pos]=i;
        fill(pos+1,i,(square+i*i)%N,indices,groups);
    }
}

#ifndef ETA_CONTACT_LIBRARY
int main(int argc,char** argv) {
    if(argc!=5) return 2;
    factors=stoi(argv[1]); height=stoi(argv[2]); cutoff=stoi(argv[3]); limit=48;
    assert(2<=factors && factors<=3 && 1<=height && height<=8);
    xi=prime-power(zeta,42);
    assert(power(xi,83)==prime-1 && power(xi,2)==zeta);
    array<vector<Monomial>,3> groups;
    fill(0,1,0,Indices{},groups);
    vector<Coeff> weights;
    for(int a=1;a<=height;++a) for(int b=-height;b<=height;++b)
        for(int c=-height;c<=height;++c)
            if(b && c && gcd(a,gcd(abs(b),abs(c)))==1) weights.push_back({a,b,c});
    ofstream out(argv[4]); assert(out);
    Summary total;
    for(int bucket=0;bucket<=2;++bucket) {
        const auto& g=groups[bucket];
        vector<Summary> summaries(g.size());
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<int(g.size());++i) {
            auto& s=summaries[i];
            if(!bucket && g[i].idx[0]!=1) continue;
            for(int j=i+1;j<int(g.size());++j) for(int k=j+1;k<int(g.size());++k) {
                if(g[i].mask&g[j].mask&g[k].mask) continue;
                if(!bucket && !representative({g[i].idx,g[j].idx,g[k].idx})) continue;
                ++s.triples;
                for(auto w:weights) examine(g,i,j,k,w,s);
            }
        }
        for(auto& s:summaries) {
            total.triples+=s.triples; total.identities+=s.identities;
            total.nonzero+=s.nonzero; total.nonsplit+=s.nonsplit;
            for(auto [d,n]:s.histogram) total.histogram[d]+=n;
            for(auto ca:s.selected) {
                for(int i:{ca.i,ca.j,ca.k}) for(int n=0;n<factors;++n) out<<g[i].idx[n]<<' ';
                out<<ca.w.a<<' '<<ca.w.b<<' '<<ca.w.c<<' '<<ca.raw<<' '<<ca.degree<<'\n';
            }
        }
    }
    assert(total.triples*I(weights.size())==total.identities+total.nonzero);
    cout<<"triples "<<total.triples<<" coefficient_choices "<<weights.size()
        <<" zero_functions "<<total.identities<<" nonzero_functions "<<total.nonzero
        <<" nonsplit_contact "<<total.nonsplit<<'\n';
    for(auto [d,n]:total.histogram) cout<<d<<' '<<n<<'\n';
}
#endif
