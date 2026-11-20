// Elliptic algebra from the level-83 signed-fibre certificate.
#pragma once
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
using std::vector;
constexpr int level=83, nc=41;
U modulus=5;

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

vector<Poly> basis, anti;
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

// From the first norm's Laurent expansion at the pole cusp, compute
// its norm on the Fricke quotient and then on the elliptic x-line.
void norm_from_series(Poly& out, vector<U> f, bool check_other_cusp, U other_cusp_value) {
    const int d=static_cast<int>(f.size())-1;
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
    if(check_other_cusp) assert(rem[d]==other_cusp_value);
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

