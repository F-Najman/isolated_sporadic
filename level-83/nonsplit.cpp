// Select trinomials that meet the cusp of degree 41 before computing the other contacts.
#define ETA_CONTACT_LIBRARY
#include "three.cpp"
#include <unordered_map>

int main(int argc,char** argv) {
    if(argc!=5) return 2;
    factors=stoi(argv[1]); height=stoi(argv[2]); cutoff=stoi(argv[3]); limit=48;
    assert(4<=factors && factors<=5 && 1<=height && height<=4);
    xi=prime-power(zeta,42);
    assert(power(xi,83)==prime-1 && power(xi,2)==zeta);
    array<vector<Monomial>,3> groups;
    fill(0,1,0,Indices{},groups);
    ofstream out(argv[4]); assert(out);
    Summary total;
    I probes=0,modhits=0,genuine=0;
    for(int bucket=0;bucket<=2;++bucket) {
        const auto& g=groups[bucket];
        // All conjugates satisfy any true relation. The separate integer
        // cyclotomic check below certifies the surviving leading relations.
        unordered_map<I,vector<int>> index;
        for(int k=0;k<int(g.size());++k) index[g[k].sine].push_back(k);
        vector<Coeff> weights;
        vector<I> neg_c_inverse;
        for(int a=1;a<=height;++a) for(int b=-height;b<=height;++b)
            for(int c=-height;c<=height;++c) if(b && c && gcd(a,gcd(abs(b),abs(c)))==1) {
                weights.push_back({a,b,c});
                neg_c_inverse.push_back(power((prime-c)%prime,prime-2));
            }
        vector<Cyclo> at0;
        for(const auto& m:g) at0.push_back(other_jets(m.idx,0)[0]);
        vector<Summary> summaries(g.size());
        #pragma omp parallel for schedule(dynamic) reduction(+:probes,modhits,genuine)
        for(int i=0;i<int(g.size());++i) {
            auto& s=summaries[i];
            if(!bucket && g[i].idx[0]!=1) continue;
            for(int j=i+1;j<int(g.size());++j) for(size_t t=0;t<weights.size();++t) {
                const auto w=weights[t];
                I v=(w.a*g[i].sine+w.b*g[j].sine)%prime;
                I wanted=(v+prime)*neg_c_inverse[t]%prime;
                ++probes;
                auto it=index.find(wanted);
                if(it==index.end()) continue;
                for(int k:it->second) {
                    if(k<=j || (g[i].mask&g[j].mask&g[k].mask)) continue;
                    if(!bucket && !representative({g[i].idx,g[j].idx,g[k].idx})) continue;
                    ++modhits;
                    Cyclo z;
                    for(int u=0;u<N;++u) z[u]=w.a*at0[i][u]+w.b*at0[j][u]+w.c*at0[k][u];
                    if(!is_zero(z)) continue;
                    ++genuine;
                    examine(g,i,j,k,w,s);
                }
            }
        }
        for(auto& s:summaries) {
            total.identities+=s.identities; total.nonzero+=s.nonzero; total.nonsplit+=s.nonsplit;
            for(auto [d,n]:s.histogram) total.histogram[d]+=n;
            for(auto ca:s.selected) {
                for(int i:{ca.i,ca.j,ca.k}) for(int n=0;n<factors;++n) out<<g[i].idx[n]<<' ';
                out<<ca.w.a<<' '<<ca.w.b<<' '<<ca.w.c<<' '<<ca.raw<<' '<<ca.degree<<'\n';
            }
        }
    }
    assert(genuine==total.identities+total.nonzero && total.nonsplit==total.nonzero);
    cout<<"lookups "<<probes<<" matches_mod_p "<<modhits<<" exact_cusp_relations "<<genuine
        <<" identically_zero "<<total.identities<<" nonconstant_or_all_cusps "<<total.nonzero<<'\n';
    for(auto [d,n]:total.histogram) cout<<d<<' '<<n<<'\n';
}
