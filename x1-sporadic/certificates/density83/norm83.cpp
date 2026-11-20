// Norm of a signed Siegel-unit fibre to the Fricke elliptic curve
// and then to the x-line. Compute simultaneously modulo several primes.
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>
#include <gmpxx.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <bitset>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <vector>

using U = unsigned long;
// Set -DUNIT_FACTORS to the number of indices on each side.
#ifndef UNIT_FACTORS
#define UNIT_FACTORS 3
#endif
constexpr int factors = UNIT_FACTORS;
static_assert(factors >= 2 && factors <= 4);
using Indices = std::array<int,factors>;
using std::vector;
constexpr int level=83, nc=41;
U modulus=5;
int precision;

U residue(int value) {
    if (value >= 0 && U(value) < modulus) return U(value);
    if (value < 0 && U(-value) <= modulus) return modulus-U(-value);
    long r = value % static_cast<long>(modulus);
    return r < 0 ? modulus-U(-r) : U(r);
}

struct Poly {
    nmod_poly_struct data;
    explicit Poly(U m=modulus) { nmod_poly_init(&data,m); }
    Poly(const Poly& other) : Poly(other.data.mod.n) { nmod_poly_set(&data,&other.data); }
    Poly& operator=(const Poly& other) { nmod_poly_set(&data,&other.data); return *this; }
    ~Poly() { nmod_poly_clear(&data); }
    operator nmod_poly_struct*() { return &data; }
    operator const nmod_poly_struct*() const { return &data; }
    U at(int n) const { return nmod_poly_get_coeff_ui(&data,n); }
    void set(int n,U x) { nmod_poly_set_coeff_ui(&data,n,x); }
};

struct EllFunction { Poly a,b; }; // a(x)+b(x)*y
Poly cubic;
EllFunction branch;

void multiply(EllFunction& result,const EllFunction& f,const EllFunction& g) {
    // y^2 = x^3+x^2+x-(x+1)y.
    Poly ac,bd,mixed,t;
    nmod_poly_mul(ac,f.a,g.a);
    nmod_poly_mul(bd,f.b,g.b);
    nmod_poly_add(mixed,f.a,f.b);
    nmod_poly_add(t,g.a,g.b);
    nmod_poly_mul(mixed,mixed,t);
    nmod_poly_mul(t,cubic,bd);
    nmod_poly_add(result.a,ac,t);
    nmod_poly_sub(mixed,mixed,ac);
    nmod_poly_scalar_mul_nmod(t,bd,2);
    nmod_poly_sub(result.b,mixed,t);
    nmod_poly_shift_left(t,bd,1); nmod_poly_sub(result.b,result.b,t);
}

struct Expansion {
    int valuation,sign,contact;
    vector<short> coeff;
};
struct Monomial { Indices indices; std::array<Expansion,nc> at; };
struct Fibre { int a,b,sign,degree,poles; };
vector<Monomial> monomials;
vector<Poly> basis, anti;
Poly denominator_inverse;

void build_monomial(Monomial& mon) {
    const int length=precision+nc+1;
    for(int j=0;j<nc;++j) {
        int c=j+1;
        if(c%2==0) c+=level;
        auto& v=mon.at[j];
        v.valuation=0; v.sign=1; v.contact=nc;
        vector<int> coeff(length);
        coeff[0]=1;
        for(int i:mon.indices) {
            int r=c*i%level;
            v.valuation+=r*(r-level);
            if((c*i/level)%2) v.sign=-v.sign;
            v.contact=std::min({v.contact,r,level-r});
            for(int first:{r,level-r}) for(int n=first;n<length;n+=level)
                for(int k=length-1;k>=n;--k) coeff[k]-=coeff[k-n];
        }
        v.coeff.resize(length);
        for(int k=0;k<length;++k) {
            assert(std::abs(coeff[k])<32768);
            v.coeff[k]=coeff[k];
        }
    }
}

void make_basis(std::istream& in,int max_degree) {
    int length; in>>length;
    assert(length>=max_degree+12);
    Poly data[3];
    for(auto& series:data) for(int i=0;i<length;++i) {
        mpz_class v; in>>v;
        series.set(i,mpz_fdiv_ui(v.get_mpz_t(),modulus));
    }
    assert(in && data[0].at(0)==1 && data[2].at(0)==1);
    basis.resize(max_degree+1);
    anti.resize(max_degree+1);
    nmod_poly_one(basis[0]);
    basis[2]=data[0]; basis[3]=data[1];
    for(int d=4;d<=max_degree;++d)
        nmod_poly_mullow(basis[d],basis[d-2],data[0],max_degree+11);
    for(int i=0;i<=max_degree-6;++i)
        nmod_poly_mullow(anti[i],basis[i],data[2],i+11);
    // Every monomial has the same unsigned product of its denominator
    // unit parts across the 41 rational cusps.
    vector<U> den(max_degree+1); den[0]=1;
    for(int n=1;n<=max_degree;++n) if(n%level)
        for(int t=0;t<factors;++t) for(int k=max_degree;k>=n;--k)
            den[k]=(den[k]+modulus-den[k-n])%modulus;
    Poly euler;
    for(int k=0;k<=max_degree;++k) euler.set(k,den[k]);
    nmod_poly_inv_series(denominator_inverse,euler,max_degree+1);
    cubic.set(1,1); cubic.set(2,1); cubic.set(3,1);
    const int aa[]={-19,-39,-31,3,16,8,1};
    const int bb[]={-11,-13,4,10,4};
    for(int i=0;i<=6;++i) branch.a.set(i,residue(aa[i]));
    for(int i=0;i<=4;++i) branch.b.set(i,residue(bb[i]));
}

void subtract_row(vector<U>& f,int pos,const Poly& row,U value) {
    if(value==0) return;
    int len=int(f.size())-pos;
    assert(row.at(0)==1 || row.at(0)==modulus-1);
    int n=std::min(len,int(row.data.length));
    _nmod_vec_scalar_addmul_nmod(f.data()+pos,row.data.coeffs,n,modulus-value,row.data.mod);
}

void put_basis_term(EllFunction& out,int order,U coefficient) {
    if(order%2) out.b.set((order-3)/2,coefficient);
    else out.a.set(order/2,coefficient);
}

void fricke_polynomial(Poly& out,const Fibre& F) {
    const int d=F.degree;
    const auto& A=monomials[F.a];
    const auto& B=monomials[F.b];
    Poly product,local;
    nmod_poly_one(product);
    int degree_check=0,pole_check=0;
    for(int j=0;j<nc;++j) {
        const auto& a=A.at[j]; const auto& b=B.at[j];
        int diff=a.valuation-b.valuation;
        assert(diff%(2*level)==0);
        int v=diff/(2*level);
        int cancel= v==0 && a.sign==F.sign*b.sign ? std::min(a.contact,b.contact):0;
        pole_check+=std::max(-v,0); degree_check+=std::max(-v,0)-cancel;
        int sa=std::max(v,0), sb=std::max(-v,0);
        nmod_poly_fit_length(local,d+1);
        for(int k=0;k<=d;++k) {
            int ta=k+cancel-sa,tb=k+cancel-sb;
            int value=(ta<0 ? 0:a.coeff[ta]*a.sign)
                      -(tb<0 ? 0:b.coeff[tb]*b.sign*F.sign);
            local.data.coeffs[k]=residue(value);
        }
        local.data.length=d+1; _nmod_poly_normalise(local);
        assert(std::gcd(local.at(0),modulus)==1);
        nmod_poly_mullow(product,product,local,d+1);
    }
    assert(d==degree_check && F.poles==pole_check);
    nmod_poly_mullow(product,product,denominator_inverse,d+1);
    vector<U> f(d+1);
    for(int k=0;k<=d;++k) f[k]=product.at(k);
    vector<U> rem=f;
    EllFunction W,T;
    // F-zW is Fricke-conjugate to F and is regular at the split cusp.
    for(int j=d;j>=1;--j) {
        U lead=rem[d-j];
        if(j<6 || j==7) { assert(!lead); continue; }
        const Poly& row=anti[j-6];
        U coeff=row.at(0)==1 ? lead : (modulus-lead)%modulus;
        put_basis_term(W,j-6,coeff);
        subtract_row(rem,d-j,row,coeff);
    }
    // T=F+wF. Only the constant correction is required at this cusp.
    f[d]=(f[d]+rem[d])%modulus;
    for(int j=d;j>=0;--j) {
        if(j==1) { assert(f[d-j]==0); continue; }
        U lead=f[d-j];
        U coeff=basis[j].at(0)==1 ? lead : (modulus-lead)%modulus;
        put_basis_term(T,j,coeff);
        subtract_row(f,d-j,basis[j],coeff);
    }
    EllFunction T2,W2,RW;
    multiply(T2,T,T); multiply(W2,W,W); multiply(RW,branch,W2);
    nmod_poly_sub(T2.a,T2.a,RW.a);
    nmod_poly_sub(T2.b,T2.b,RW.b);
    // Norm of a+b*y = a*(a-(x+1)b) - (x^3+x^2+x)*b^2.
    Poly c,tail;
    nmod_poly_sub(c,T2.a,T2.b);
    nmod_poly_shift_left(tail,T2.b,1); nmod_poly_sub(c,c,tail);
    nmod_poly_mul(out,c,T2.a);
    nmod_poly_mul(tail,T2.b,T2.b);
    nmod_poly_mul(tail,tail,cubic);
    nmod_poly_sub(out,out,tail);
    assert(nmod_poly_degree(out)<=d);
}

std::bitset<143> factor_degrees(const Poly& integral,int d,U p,std::ostream& log) {
    Poly f(p);
    for(int i=0;i<=d;++i) f.set(i,integral.at(i)%p);
    assert(nmod_poly_degree(f)==d);
    nmod_poly_make_monic(f,f);
    nmod_poly_factor_t sf,dd;
    nmod_poly_factor_init(sf); nmod_poly_factor_init(dd);
    nmod_poly_factor_squarefree(sf,f);
    std::bitset<143> mask; mask.set(0);
    int total=0;
    log<<' '<<p<<':';
    for(int i=0;i<sf->num;++i) {
        int n=nmod_poly_degree(sf->p+i);
        vector<slong> degrees(n+1); slong* ptr=degrees.data();
        dd->num=0;
        nmod_poly_factor_distinct_deg(dd,sf->p+i,&ptr);
        for(int j=0;j<dd->num;++j) {
            int degree=degrees[j]; assert(degree>0 && nmod_poly_degree(dd->p+j)%degree==0);
            int count=(nmod_poly_degree(dd->p+j)/degree)*sf->exp[i];
            for(int k=0;k<count;++k) mask|=mask<<degree;
            log<<degree<<'^'<<count<<',';
            total+=degree*count;
        }
    }
    nmod_poly_factor_clear(sf); nmod_poly_factor_clear(dd);
    assert(total==d);
    return mask;
}

int main(int argc,char** argv) {
    // CASES SERIES START STOP THREADS prime,prime,... [degree limit [polynomial]]
    if(argc<7) return 2;
    int start=std::stoi(argv[3]),stop=std::stoi(argv[4]);
    int limit=argc>7 ? std::stoi(argv[7]):129;
    assert(limit>0 && limit<=142);
    omp_set_num_threads(std::stoi(argv[5]));
    std::stringstream choices(argv[6]);
    vector<U> primes;
    std::string token;
    modulus=1;
    while(std::getline(choices,token,',')) {
        U p=std::stoul(token);
        assert(p>3 && p!=level && n_is_prime(p) && modulus%p && modulus < (1UL<<61)/p);
        primes.push_back(p); modulus*=p;
    }
    // These empty polynomial objects were constructed before the modulus
    // was read. Initialize their coefficient rings here.
    for(Poly* f:{&cubic,&branch.a,&branch.b,&denominator_inverse}) nmod_init(&f->data.mod,modulus);
    std::ifstream in(argv[1]),series(argv[2]); assert(in && series);
    std::map<Indices,int> lookup;
    vector<Fibre> cases;
    Indices a,b; int sign,d,pole;
    auto index=[&](const Indices& t) {
        auto [it,added]=lookup.emplace(t,monomials.size());
        if(added) { monomials.emplace_back(); monomials.back().indices=t; }
        return it->second;
    };
    for(int i=0; i<stop && in>>a[0]; ++i) {
        for(int j=1;j<factors;++j) in>>a[j];
        for(int& t:b) in>>t;
        in>>sign>>d>>pole;
        assert(in && (sign==1 || sign==-1));
        if(i<start || !d) continue;
        cases.push_back({index(a),index(b),sign,d,pole});
        precision=std::max(precision,d);
    }
    if(cases.empty()) return 0;
    make_basis(series,precision);
    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0;i<monomials.size();++i) build_monomial(monomials[i]);
    #pragma omp parallel for schedule(dynamic,2)
    for(size_t i=0;i<cases.size();++i) {
        const auto& c=cases[i];
        Poly f; fricke_polynomial(f,c);
        auto a=monomials[c.a].indices, b=monomials[c.b].indices;
        std::ostringstream line;
        for(int t:a) line<<t<<' ';
        for(int t:b) line<<t<<' ';
        line<<c.sign<<' '<<c.degree;
        std::bitset<143> possible;
        for(int j=1;j<=std::min(limit,c.degree/2);++j) possible.set(j);
        for(U p:primes) {
            if(!possible.any()) break;
            if(f.at(c.degree)%p==0) continue;
            possible &= factor_degrees(f,c.degree,p,line);
        }
        line<<" left";
        for(int k=1;k<=limit;++k) if(possible[k]) line<<' '<<k;
        if(argc>8) {
            line<<" poly";
            for(int j=0;j<=c.degree;++j) line<<' '<<f.at(j);
        }
        #pragma omp critical
        { std::cout<<line.str()<<'\n'; }
    }
}
