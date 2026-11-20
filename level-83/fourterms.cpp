// Signed sums of four products of three generalized eta functions.
#define ETA_CONTACT_LIBRARY
#include "three.cpp"

using Four=array<int,4>;
using Support=array<const Monomial*,4>;
struct Selected { Four terms; int mask,raw,degree; };
struct Tally {
    I supports=0,zero=0,other=0;
    array<I,500> hist{};
    vector<Selected> cases;
};
array<Four,8> weights;

int fourth_contact(const Support& s,const Four& w,int raw) {
    I initial=0;
    for(int i=0;i<4;++i) initial+=w[i]*s[i]->sine;
    if(initial%prime) return 0;
    const int n=raw/H;
    vector<Cyclo> value(n+1);
    for(int i=0;i<4;++i) {
        auto v=other_jets(s[i]->idx,n);
        for(int k=0;k<=n;++k) for(int j=0;j<N;++j) value[k][j]+=w[i]*v[k][j];
    }
    for(int k=0;k<=n;++k) if(!is_zero(value[k])) return k;
    abort();
}

void examine_four(const vector<Monomial>& g,Four ids,Tally& sums) {
    const Support s={&g[ids[0]],&g[ids[1]],&g[ids[2]],&g[ids[3]]};
    ++sums.supports;
    array<int,H> lows;
    int base=0;
    for(int h=0;h<H;++h) {
        lows[h]=min({s[0]->val[h],s[1]->val[h],s[2]->val[h],s[3]->val[h]});
        base+=s[0]->val[h]-lows[h];
    }
    array<int,8> contacts{};
    for(int h=0;h<H;++h) {
        Four offsets,lead;
        int nr=0;
        for(int i=0;i<4;++i) {
            offsets[i]=s[i]->val[h]-lows[h];
            lead[i]=offsets[i]?0:s[i]->sign[h];
            nr+=!offsets[i];
        }
        if(nr%2) continue;
        array<vector<int>,4> extra;
        for(int m=0;m<8;++m) {
            const auto& w=weights[m];
            if(lead[0]+w[1]*lead[1]+w[2]*lead[2]+w[3]*lead[3]) continue;
            if(contacts[m]>base) continue;
            Four v;
            for(int i=0;i<4;++i) v[i]=s[i]->sign[h]*w[i];
            int e=1, cap=base-contacts[m]+1;
            for(;e<cap;++e) {
                if(e==limit && extra[0].empty())
                    for(int i=0;i<4;++i) extra[i]=qpart(s[i]->idx,h+1,base+1);
                int z=0;
                for(int i=0;i<4;++i) if(offsets[i]<=e) {
                    const auto& q=(e<limit?s[i]->jets[h]:extra[i]);
                    z+=v[i]*q[e-offsets[i]];
                }
                if(z) break;
            }
            contacts[m]+=e;
        }
    }
    for(int m=0;m<8;++m) {
        int raw=base-contacts[m];
        if(raw<0) { ++sums.zero; continue; }
        int add=fourth_contact(s,weights[m],raw);
        int degree=raw-H*add;
        if(add) ++sums.other;
        assert(0<=degree && degree<int(sums.hist.size()));
        ++sums.hist[degree];
        if(degree && degree<=cutoff) sums.cases.push_back({ids,m,raw,degree});
    }
}

int main(int argc,char** argv) {
    if(argc!=3) return 2;
    factors=3; cutoff=stoi(argv[1]); limit=48;
    xi=prime-power(zeta,42);
    for(int m=0;m<8;++m) weights[m]={1,(m&1)?-1:1,(m&2)?-1:1,(m&4)?-1:1};
    array<vector<Monomial>,3> groups;
    fill(0,1,0,Indices{},groups);
    ofstream out(argv[2]); assert(out);
    Tally total;
    for(int bucket=0;bucket<=2;++bucket) {
        const auto& g=groups[bucket]; const int n=int(g.size());
        vector<Tally> sums(n);
        #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<n;++i) {
            if(!bucket && g[i].idx[0]!=1) continue;
            for(int j=i+1;j<n;++j) for(int k=j+1;k<n;++k) for(int l=k+1;l<n;++l) {
                if(g[i].mask&g[j].mask&g[k].mask&g[l].mask) continue;
                if(!bucket && !least_translate<4>({g[i].idx,g[j].idx,g[k].idx,g[l].idx})) continue;
                examine_four(g,{i,j,k,l},sums[i]);
            }
        }
        for(auto& s:sums) {
            total.supports+=s.supports; total.zero+=s.zero; total.other+=s.other;
            for(int j=0;j<int(s.hist.size());++j) total.hist[j]+=s.hist[j];
            for(auto a:s.cases) {
                for(int i:a.terms) for(int j=0;j<3;++j) out<<g[i].idx[j]<<' ';
                for(auto w:weights[a.mask]) out<<w<<' ';
                out<<a.raw<<' '<<a.degree<<'\n';
            }
        }
    }
    assert(total.supports*8==total.zero+accumulate(total.hist.begin(),total.hist.end(),I(0)));
    cout<<"supports "<<total.supports<<" signed_sections "<<total.supports*8
        <<" zero_functions "<<total.zero<<" nonsplit_contact "<<total.other<<'\n';
    for(int i=0;i<int(total.hist.size());++i) if(total.hist[i]) cout<<i<<' '<<total.hist[i]<<'\n';
}
