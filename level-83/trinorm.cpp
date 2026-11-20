// Norms of nonzero linear combinations of balanced eta monomials.
#include "norm_core.hpp"
#define ETA_CONTACT_LIBRARY
#include "three.cpp"

struct NormRow {
    vector<int> terms, weights;
    int raw, degree, base;
};

int main(int argc,char** argv) {
    // CASES SERIES START STOP THREADS PRIMES DEGREE_LIMIT FACTORS [TERMS]
    if(argc<9) return 2;
    const int start=stoi(argv[3]), stop=stoi(argv[4]), threshold=stoi(argv[7]);
    factors=stoi(argv[8]);
    const int count=argc>9?stoi(argv[9]):3;
    assert(factors>=2 && factors<=5 && count>=2 && threshold<143);
    omp_set_num_threads(stoi(argv[5]));
    stringstream inputs(argv[6]);
    string token;
    vector<U> primes;
    ::modulus=1;
    while(getline(inputs,token,',')) {
        U p=stoul(token);
        assert(p>3 && p!=N && n_is_prime(p) && ::modulus%p && ::modulus<(1UL<<61)/p);
        primes.push_back(p); ::modulus*=p;
    }
    for(Poly* f:{&cubic,&branch.a,&branch.b}) nmod_init(&f->data.mod,::modulus);
    xi=prime-power(zeta,42);
    ifstream in(argv[1]),series(argv[2]); assert(in && series);
    map<Indices,int> ids;
    vector<Monomial> mons;
    vector<NormRow> rows;
    int depth=0, max_raw=0;
    for(int n=0;n<stop;++n) {
        string line;
        if(!getline(in,line)) break;
        if(n<start) continue;
        istringstream s(line);
        NormRow row;
        for(int i=0;i<count;++i) {
            Indices a{};
            for(int j=0;j<factors;++j) s>>a[j];
            assert(s && a[0]>=1 && a[factors-1]<=H);
            auto [it,new_key]=ids.emplace(a,mons.size());
            if(new_key) mons.push_back(monomial(a,false));
            row.terms.push_back(it->second);
        }
        row.weights.resize(count);
        for(int& w:row.weights) s>>w;
        s>>row.raw>>row.degree;
        assert(s && row.degree>0 && row.raw>=row.degree && (row.raw-row.degree)%H==0);
        row.base=0;
        for(int h=0;h<H;++h) {
            int ref=mons[row.terms[0]].val[h], low=ref;
            for(auto i:row.terms) low=min(low,mons[i].val[h]);
            row.base+=ref-low;
        }
        assert(row.base>=row.raw);
        // The contact at any cusp is at most the sum of all contacts.
        depth=max(depth,row.raw+(row.base-row.raw)+1);
        max_raw=max(max_raw,row.raw);
        rows.push_back(row);
    }
    if(rows.empty()) return 0;
    make_basis(series,max_raw);
    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0;i<mons.size();++i)
        for(int h=0;h<H;++h) mons[i].jets[h]=qpart(mons[i].idx,h+1,depth);
    Poly euler,inverse;
    vector<U> den(max_raw+1);
    den[0]=1;
    for(int n=1;n<=max_raw;++n) if(n%N)
        for(int t=0;t<factors;++t) for(int k=max_raw;k>=n;--k)
            den[k]=(den[k]+::modulus-den[k-n]) % ::modulus;
    for(int k=0;k<=max_raw;++k) euler.set(k,den[k]);
    nmod_poly_inv_series(inverse,euler,max_raw+1);
    #pragma omp parallel for schedule(dynamic,2)
    for(size_t n=0;n<rows.size();++n) {
        const auto& row=rows[n];
        const int d=row.raw;
        Poly prod,local;
        nmod_poly_one(prod);
        int contacts=0;
        for(int h=0;h<H;++h) {
            int low=mons[row.terms[0]].val[h];
            for(auto i:row.terms) low=min(low,mons[i].val[h]);
            vector<int> coeff(row.base+1);
            for(int j=0;j<count;++j) {
                const auto& a=mons[row.terms[j]];
                int offset=a.val[h]-low, w=a.sign[h]*row.weights[j];
                for(int t=offset;t<=row.base;++t) coeff[t]+=w*a.jets[h][t-offset];
            }
            int g=0;
            while(g<=row.base && !coeff[g]) ++g;
            assert(g+d<=row.base);
            contacts+=g;
            nmod_poly_fit_length(local,d+1);
            for(int k=0;k<=d;++k) local.data.coeffs[k]=residue(coeff[k+g]);
            local.data.length=d+1; _nmod_poly_normalise(local);
            nmod_poly_mullow(prod,prod,local,d+1);
        }
        assert(row.base-contacts==d);
        nmod_poly_mullow(prod,prod,inverse,d+1);
        vector<U> f(d+1);
        for(int k=0;k<=d;++k) f[k]=prod.at(k);
        Poly norm;
        norm_from_series(norm,f,row.degree<d,0);
        assert(nmod_poly_degree(norm)<=row.degree);
        ostringstream out;
        for(auto i:row.terms) for(int j=0;j<factors;++j) out<<mons[i].idx[j]<<' ';
        for(int w:row.weights) out<<w<<' ';
        out<<row.degree;
        bitset<143> possible;
        for(int e=1;e<=min(threshold,row.degree/2);++e) possible.set(e);
        for(U p:primes) {
            if(!possible.any()) break;
            if(norm.at(row.degree)%p) possible &= factor_degrees(norm,row.degree,p,out);
        }
        out<<" left";
        for(int e=1;e<=threshold;++e) if(possible[e]) out<<' '<<e;
        #pragma omp critical
        { cout<<out.str()<<'\n'; }
    }
}
