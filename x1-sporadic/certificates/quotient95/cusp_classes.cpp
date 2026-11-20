// Exhaust degree-95 divisors on specified subsets of the 13 split cusps.
// A full-rank upper bound on multiplicities certifies every smaller divisor.
#include <algorithm>
#include <array>
#include <cassert>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
constexpr int target=95, nforms=130, ncusps=13;
struct Row {
    array<uint64_t,3> a{},b{};
    int pivot() const {
        for(int w=0;w<3;++w) if(auto z=a[w]|b[w]) return 64*w+__builtin_ctzll(z);
        return -1;
    }
    bool negative(int p) const {return b[p/64] & (uint64_t(1)<<(p%64));}
    void add(const Row& v,bool negate) {
        for(int w=0;w<3;++w) {
            auto c=negate?v.b[w]:v.a[w],d=negate?v.a[w]:v.b[w];
            auto u=a[w],v0=b[w],zero=~(u|v0),otherzero=~(c|d);
            a[w]=(u&otherzero)|(zero&c)|(v0&d);
            b[w]=(v0&otherzero)|(zero&d)|(u&c);
        }
    }
};
array<array<Row,100>,ncusps> jets;
struct Counter {uint64_t passed=0,pencils=0,failed=0,checks=0;};
using Unit=array<int,13>;
vector<Unit> units;
ofstream unresolved;
using Shape=vector<int>;
int rankCaps(const Shape& positions,const Shape& caps) {
    array<Row,nforms> basis{}; array<bool,nforms> seen{};int rank=0;
    for(unsigned i=0;i<caps.size();++i) for(int j=0;j<caps[i];++j){
        Row v=jets[positions[i]][j];
        for(int p=v.pivot();p>=0;p=v.pivot()){
            if(!seen[p]){
                if(v.negative(p))swap(v.a,v.b);
                basis[p]=v;seen[p]=true;++rank;break;
            }
            v.add(basis[p],!v.negative(p));
        }
    }
    return rank;
}
uint64_t ways(const Shape& lo,const Shape& hi,int sumlo){
    // Coefficients of the product (1+...+X^{hi[i]-lo[i]}).
    array<uint64_t,target+1> prev{},next{};prev[0]=1;int room=target-sumlo;
    for(unsigned i=0;i<lo.size();++i){
        uint64_t across=0;
        for(int j=0;j<=room;++j){
            across+=prev[j];
            if(j>hi[i]-lo[i])across-=prev[j-hi[i]+lo[i]-1];
            next[j]=across;
        }
        prev=next;
    }
    return prev[room];
}
bool pencilTest(const Shape& positions,const Shape& caps,const vector<Unit>& localUnits){
    array<int,13> D{};for(unsigned i=0;i<positions.size();++i)D[positions[i]]=caps[i];
    for(const Unit& u:localUnits){
        Shape pos,high;bool effective=true;
        for(int i=0;i<13;++i){
            int a=D[i]+u[i];if(a<0){effective=false;break;}
            if(a>0){pos.push_back(i);high.push_back(a);}
        }
        if(effective&&rankCaps(pos,high)==94)return true;
    }
    return false;
}
void cover(const Shape& positions,Shape lo,Shape hi,const vector<Unit>& localUnits,Counter& count){
    int l=0,h=0;for(int n:lo)l+=n;for(int n:hi)h+=n;
    if(l>target||h<target)return;
    int sl=0,sh=0,widest=0;
    for(unsigned i=0;i<lo.size();++i){
        int a=lo[i],b=hi[i];
        lo[i]=max(a,target-h+b);hi[i]=min(b,target-l+a);
        sl+=lo[i];sh+=hi[i];
        if(hi[i]-lo[i]>hi[widest]-lo[widest])widest=i;
    }
    if(sh<=(positions.size()>=8?126:119)){
        ++count.checks;
        int rank=rankCaps(positions,hi);
        if(rank==sh){count.passed+=ways(lo,hi,sl);return;}
        if(sh==target){
            if(rank==94&&pencilTest(positions,hi,localUnits)){++count.pencils;return;}
            ++count.failed;
            array<int,13> d{};
            for(unsigned j=0;j<positions.size();++j)d[positions[j]]=hi[j];
            #pragma omp critical(unsolved)
            {unresolved<<rank;for(int x:d)unresolved<<" "<<x;unresolved<<endl;}
            return;
        }
    }
    assert(hi[widest]>lo[widest]);
    int mid=(lo[widest]+hi[widest])/2;
    Shape small=hi;small[widest]=mid;
    cover(positions,lo,small,localUnits,count);lo[widest]=mid+1;cover(positions,lo,hi,localUnits,count);
}
uint64_t binom(int n,int k){uint64_t v=1;for(int i=1;i<=k;++i)v=v*(n-i+1)/i;return v;}
int main(int argc,char**argv){
    assert(argc==6);unresolved.open(argv[5]);ifstream in(argv[1]);string s;
    for(auto& jet:jets)for(auto& r:jet){
        assert(bool(in>>s)&&s.size()==131);
        for(int p=0;p<nforms;++p){
            if(s[p]=='1')r.a[p/64]|=uint64_t(1)<<(p%64);
            else if(s[p]=='2')r.b[p/64]|=uint64_t(1)<<(p%64);
            else assert(s[p]=='0');
        }
    }
    ifstream unitFile(argv[2]);int id;
    while(unitFile>>id){Unit v{};v[0]=id;for(int i=1;i<13;++i)unitFile>>v[i];units.push_back(v);}
    assert(units.size()==468);
    int first=atoi(argv[3]),upto=atoi(argv[4]);assert(upto<=9);
    for(int k=first;k<=upto;++k){
        int all=(1<<ncusps)-1,num=0;Counter count;
        vector<Counter> byMask(all+1);
        #pragma omp parallel for schedule(dynamic) reduction(+:num)
        for(int mask=1;mask<=all;++mask){
            if(__builtin_popcount(mask)!=k)continue;
            bool rep=true;int shift=mask;
            for(int i=1;i<ncusps;++i){shift=((shift<<1)&all)|(shift>>(ncusps-1));if(shift<mask)rep=false;}
            if(!rep)continue;
            ++num;Shape pos;
            for(int j=0;j<ncusps;++j)if((1<<j)&mask)pos.push_back(j);
            vector<Unit> localUnits;
            for(const Unit& u:units){
                bool fits=true;for(int j=0;j<13;++j)if(u[j]<0&&!(mask&(1<<j)))fits=false;
                if(fits)localUnits.push_back(u);
            }
            cover(pos,Shape(k,1),Shape(k,target),localUnits,byMask[mask]);
        }
        for(const auto& c:byMask){count.passed+=c.passed;count.pencils+=c.pencils;count.failed+=c.failed;count.checks+=c.checks;}
        assert(13*num==int(binom(13,k)));
        assert(count.passed+count.pencils+count.failed==num*binom(94,k-1));
        cout<<"support "<<k<<" diamond representatives "<<num<<" rank boxes "<<count.checks
            <<" divisors passed "<<13*count.passed<<" pencils "<<13*count.pencils<<" unresolved "<<13*count.failed<<endl;
    }
}
