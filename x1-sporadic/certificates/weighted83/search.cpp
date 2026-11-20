// Norms of rational fibres of level-83 modular units with arbitrary exponents.
#include "norm_core.hpp"
#include <string>

struct Fibre {
    int numerator, denominator, poles, degree;
};
struct Unit {
    std::array<int,nc> exponents{};
    std::vector<Fibre> fibres;
    bool zero_at_nonsplit_cusp=false;
};
struct Cusp {
    int order, first, sign;
    Poly unit;
};

inline U add_mod(U a,U b) {
    U z=a+b;
    return z>=modulus ? z-modulus : z;
}
inline U subtract_mod(U a,U b) {
    return a<b ? modulus-(b-a) : a-b;
}

std::vector<Cusp> cusp_data(const Unit& unit,int limit) {
    std::vector<Cusp> result(nc);
    int total=std::accumulate(unit.exponents.begin(),unit.exponents.end(),0);
    assert(total>=0 && total%12==0);
    for(int c=0;c<nc;++c) {
        int a=c+1;
        if(a%2==0) a+=level;
        auto& point=result[c];
        int numerator=level*level*total;
        int parity=0;
        point.first=level;
        std::array<int,level> exponents{};
        for(int j=0;j<nc;++j) {
            int e=unit.exponents[j];
            if(!e) continue;
            int m=a*(j+1), r=m%level;
            numerator+=6*e*r*(r-level);
            parity+=e*(m/level);
            point.first=std::min({point.first,r,level-r});
            exponents[r]+=e;
            exponents[level-r]+=e;
        }
        assert(numerator%(12*level)==0);
        point.order=numerator/(12*level);
        point.sign=parity%2==0 ? 1:-1;
        int length=limit+point.first+1;
        std::vector<U> v(length);
        v[0]=1;
        for(int i=1;i<length;++i) {
            int e=exponents[i%level];
            for(int j=0;j<e;++j)
                for(int n=length-1;n>=i;--n) v[n]=subtract_mod(v[n],v[n-i]);
            for(int j=0;j< -e;++j)
                for(int n=i;n<length;++n) v[n]=add_mod(v[n],v[n-i]);
        }
        for(int n=0;n<length;++n) point.unit.set(n,v[n]);
    }
    return result;
}

void fibre_polynomial(Poly& norm,const std::vector<Cusp>& cusps,const Fibre& test,
                      bool zero_at_nonsplit_cusp) {
    int d=test.degree;
    Poly product,local;
    nmod_poly_one(product);
    int poles=0,contacts=0;
    for(const auto& p:cusps) {
        int v=p.order;
        int cancel=(v==0 && test.denominator*p.sign==test.numerator) ? p.first:0;
        poles+=std::max(0,-v);
        contacts+=cancel;
        // Expand b*U-a after removing its signed order at this cusp.
        int offset=cancel-std::max(v,0);
        int constant=std::max(-v,0)-cancel;
        nmod_poly_fit_length(local,d+1);
        U scale=residue(test.denominator*p.sign);
        for(int i=0;i<=d;++i) {
            U z=0;
            if(i+offset>=0 && i+offset<p.unit.data.length)
                z=nmod_mul(scale,p.unit.data.coeffs[i+offset],p.unit.data.mod);
            if(i==constant) z=subtract_mod(z,residue(test.numerator));
            local.data.coeffs[i]=z;
        }
        local.data.length=d+1;
        _nmod_poly_normalise(local);
        nmod_poly_mullow(product,product,local,d+1);
    }
    assert(poles==test.poles && poles-contacts==d);
    std::vector<U> coefficients(d+1);
    for(int i=0;i<=d;++i) coefficients[i]=product.at(i);
    // All 41 conjugate cusp values of b*U-a equal -a when sum(e_i)>0.
    U expected=nmod_pow_ui(residue(-test.numerator),nc,norm.data.mod);
    norm_from_series(norm,std::move(coefficients),zero_at_nonsplit_cusp,expected);
}

int main(int argc,char** argv) {
    // Exponents, selected fibres, elliptic series, primes, degree bound, threads.
    if(argc!=7) return 2;
    std::ifstream exponents(argv[1]),targets(argv[2]),series(argv[3]);
    assert(exponents && targets && series);
    std::stringstream list(argv[4]);
    std::string token;
    std::vector<U> primes;
    modulus=1;
    while(std::getline(list,token,',')) {
        U p=std::stoul(token);
        assert(p>3 && p!=level && n_is_prime(p) && modulus%p && modulus<(1UL<<61)/p);
        modulus*=p;
        primes.push_back(p);
    }
    for(Poly* f:{&cubic,&branch.a,&branch.b}) nmod_init(&f->data.mod,modulus);
    int bound=std::stoi(argv[5]);
    assert(bound>0 && bound<=142);
    omp_set_num_threads(std::stoi(argv[6]));
    std::vector<Unit> units;
    int label;
    while(exponents>>label) {
        assert(label==static_cast<int>(units.size()));
        units.emplace_back();
        for(int& e:units.back().exponents) exponents>>e;
        assert(exponents);
        const auto& e=units.back().exponents;
        units.back().zero_at_nonsplit_cusp=std::accumulate(e.begin(),e.end(),0)>0;
    }
    int max_degree=0;
    Fibre t;
    while(targets>>label>>t.numerator>>t.denominator>>t.poles>>t.degree) {
        assert(label>=0 && label<static_cast<int>(units.size()));
        if(t.degree) {
            units[label].fibres.push_back(t);
            max_degree=std::max(max_degree,t.degree);
        }
    }
    if(!max_degree) return 0;
    make_basis(series,max_degree);
    #pragma omp parallel for schedule(dynamic,1)
    for(size_t i=0;i<units.size();++i) {
        const auto& u=units[i];
        if(u.fibres.empty()) continue;
        int precision=0;
        for(const auto& f:u.fibres) precision=std::max(precision,f.degree);
        auto cusps=cusp_data(u,precision);
        for(const auto& f:u.fibres) {
            std::ostringstream log;
            log<<i<<' '<<f.numerator<<' '<<f.denominator<<' '<<f.poles<<' '<<f.degree;
            Poly norm;
            fibre_polynomial(norm,cusps,f,u.zero_at_nonsplit_cusp);
            std::bitset<143> allowed;
            for(int j=1;j<=std::min(bound,f.degree/2);++j) allowed.set(j);
            for(U p:primes) {
                if(allowed.none()) break;
                if(norm.at(f.degree)%p) allowed &= factor_degrees(norm,f.degree,p,log);
            }
            log<<" left";
            for(int j=1;j<=bound;++j) if(allowed[j]) log<<' '<<j;
            #pragma omp critical
            { std::cout<<log.str()<<'\n'; }
        }
    }
}
