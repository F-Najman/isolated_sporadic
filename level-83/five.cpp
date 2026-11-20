// Rational cusp orders and signed contact orders for five-over-five units on X_1(83).
// Compile with c++ -O3 -std=c++17. Run with a residual-degree cutoff and an output file.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

constexpr int p = 83, cusps = 41;
using Quint = std::array<int,5>;
using Pair = std::pair<Quint,Quint>;
using Integer = std::int64_t;
constexpr Integer split_prime = 83000831, zeta = 1495617;
int sine_value[p], inverse[p];

Integer power(Integer a, unsigned n, Integer mod) {
    Integer result = 1;
    while (n) {
        if (n&1) result = result*a%mod;
        a = a*a%mod;
        n /= 2;
    }
    return result;
}
int fold(int a) {
    a %= p;
    return std::min(a,p-a);
}
Integer magnitude(const Quint& q, int k) {
    Integer result = 1;
    for (int a : q) result = result*sine_value[(k*a)%p]%split_prime;
    return result;
}
bool zero_representative(const Quint& a, const Quint& b) {
    auto translate = [](const Quint& q, int k) {
        Quint t;
        for (int j=0;j<5;++j) t[j]=fold(q[j]*k);
        std::sort(t.begin(),t.end());
        return t;
    };
    const Pair original{a,b};
    for (const Quint* side : {&a,&b})
        for (int k : *side) {
            auto left=translate(a,inverse[k]), right=translate(b,inverse[k]);
            if (right<left) std::swap(left,right);
            if (Pair{left,right}<original) return false;
        }
    return true;
}
struct Monomial {
    Quint indices;
    std::array<std::int16_t,cusps> order, parity, first;
    std::uint64_t mask=0;
    Integer mag;
    explicit Monomial(const Quint& q) : indices(q), mag(magnitude(q,1)) {
        for (int i:q) {
            mask |= std::uint64_t(1)<<(i-1);
        }
        for (int c=1;c<=cusps;++c) {
            int u=c%2 ? c : c+p;
            int v=0, s=0, m=cusps;
            for (int i:q) {
                int r=u*i%p;
                v+=r*(p-r);
                s+=u*i/p;
                m=std::min(m,fold(r));
            }
            // In a square-sum bucket v modulo 2p is the same on each side.
            order[c-1]=v/(2*p);
            parity[c-1]=s%2;
            first[c-1]=m;
        }
    }
};
int main(int argc,char** argv) {
    if (argc!=3) return 2;
    int cutoff=std::stoi(argv[1]);
    std::ofstream selected(argv[2]);
    assert(selected);
    assert(power(zeta,p,split_prime)==1 && zeta!=1);
    assert(power(2,41,p)==82);
    for (int i=1;i<p;++i) {
        inverse[i]=power(i,p-2,p);
        sine_value[i]=(2+2*split_prime-power(zeta,i,split_prime)
                       -power(zeta,p-i,split_prime))%split_prime;
    }
    std::array<Integer,p> sizes{};
    std::array<std::vector<Monomial>,3> groups;
    for (int a=1;a<=41;++a) for (int b=a;b<=41;++b)
      for (int c=b;c<=41;++c) for (int d=c;d<=41;++d) for(int e=d;e<=41;++e) {
        int s=(a*a+b*b+c*c+d*d+e*e)%p;
        ++sizes[s];
        if (s<=2) groups[s].emplace_back(Quint{a,b,c,d,e});
      }
    Integer raw=0,units=0,collisions=0;
    for (Integer count:sizes) raw+=count*(count-1)/2;
    std::array<Integer,2500> hist{};
    for (int s=0;s<=2;++s) {
        const auto& g=groups[s];
        for (std::size_t i=0;i<g.size();++i) {
            const auto& a=g[i];
            if (!s && a.indices[0]!=1) break;
            for (std::size_t j=i+1;j<g.size();++j) {
                const auto& b=g[j];
                if (a.mask&b.mask) continue;
                if (!s && !zero_representative(a.indices,b.indices)) continue;
                if (a.mag==b.mag) {
                    bool different=false;
                    for (int k=2;k<=41 && !different;++k)
                        different=magnitude(a.indices,k)!=magnitude(b.indices,k);
                    assert(different);
                    ++collisions;
                }
                ++units;
                int pole=0,plus=0,minus=0;
                for (int k=0;k<cusps;++k) {
                    int diff=a.order[k]-b.order[k];
                    pole+=std::max(diff,0);
                    int contact=(diff==0) ? std::min(a.first[k],b.first[k]) : 0;
                    if (a.parity[k]==b.parity[k]) plus+=contact; else minus+=contact;
                }
                for (int sign:{1,-1}) {
                    int degree=pole-(sign==1 ? plus : minus);
                    assert(degree>=0 && degree<static_cast<int>(hist.size()));
                    ++hist[degree];
                    if (degree && degree<=cutoff) {
                        for (int v:a.indices) selected<<v<<' ';
                        for (int v:b.indices) selected<<v<<' ';
                        selected<<sign<<' '<<degree<<' '<<pole<<'\n';
                    }
                }
            }
        }
    }
    // Separate the raw pairs by the size of their exact common multiset.
    assert(raw==41*(units+41*1876505LL+861*18132LL+12341*100LL));
    
    std::cout<<"units "<<units<<" signed "<<2*units<<" prime-value collisions resolved "<<collisions
             <<" raw "<<raw<<"\n";
    for (size_t d=0;d<hist.size();++d) if(hist[d]) std::cout<<d<<" "<<hist[d]<<"\n";
}
