// Push the two degree-57 divisors through the ell-isogeny correspondence.
#include "norm_core.hpp"
using Pair=std::array<int,2>;
const Pair nums[]={{35,40},{2,9}}, dens[]={{3,34},{18,33}};
const int signs[]={-1,1};
int ell,D,length;
vector<Poly> theta;
Poly inverse_den;
std::map<Pair,Poly> products;

int folded(int x) {
    x%=83; return std::min(x,83-x);
}
void build_theta() {
    theta.resize(42);
    Poly denominator;
    nmod_poly_one(denominator);
    for(int r=1;r<=41;++r) {
        auto& q=theta[r]; q.set(0,1);
        for(int k=1;;++k) {
            int a=83*k*(k-1)/2+r*k, b=83*k*(k+1)/2-r*k;
            if(std::min(a,b)>length+42) break;
            U value=k%2?modulus-1:1;
            if(a<=length+42) q.set(a,value);
            if(b<=length+42) q.set(b,value);
        }
        nmod_poly_mullow(denominator,denominator,q,length+1);
    }
    nmod_poly_mullow(denominator,denominator,denominator,length+1);
    nmod_poly_inv_series(inverse_den,denominator,length+1);
}
struct Local { Pair indices; int val,sgn; };
Local location(Pair a,int h) {
    if(h%2==0) h+=83;
    int sum=0,par=0;
    for(int& i:a) {
        int r=h*i%83;
        sum+=r*(83-r); par+=h*i/83;
        i=folded(r);
    }
    std::sort(a.begin(),a.end());
    if(!products.count(a)) nmod_poly_mullow(products[a],theta[a[0]],theta[a[1]],length+43);
    return {a,-(sum/166),par%2?-1:1};
}
void first_norm(Poly& out,int id) {
    nmod_poly_one(out);
    int order=0;
    for(int h=1;h<=41;++h) {
        auto a=location(nums[id],h), b=location(dens[id],h);
        int lo=std::min(a.val,b.val);
        int g=(a.val==b.val && a.sgn+signs[id]*b.sgn==0)?
              std::min(a.indices[0],b.indices[0]):0;
        order+=b.val-lo-g;
        Poly f,t;
        nmod_poly_shift_left(f,products[a.indices],a.val-lo);
        nmod_poly_shift_left(t,products[b.indices],b.val-lo);
        if(a.sgn==-1) nmod_poly_neg(f,f);
        if(signs[id]*b.sgn==1) nmod_poly_add(f,f,t); else nmod_poly_sub(f,f,t);
        nmod_poly_shift_right(f,f,g);
        nmod_poly_mullow(out,out,f,length+1);
    }
    assert(order==57);
    nmod_poly_mullow(out,out,inverse_den,length+1);
}
void hecke_series(Poly& out,const Poly& a) {
    if(ell==1) { out=a; return; }
    U root=n_powmod(n_primitive_root_prime(modulus),(modulus-1)/ell,modulus);
    assert(root!=1 && n_powmod(root,ell,modulus)==1);
    Poly prod,local,compressed,scaled;
    nmod_poly_one(prod);
    vector<U> roots(ell); roots[0]=1;
    for(int k=1;k<ell;++k) roots[k]=roots[k-1]*root%modulus;
    for(int k=0;k<ell;++k) {
        nmod_poly_fit_length(local,length+1);
        for(int i=0;i<=length;++i)
            local.data.coeffs[i]=a.at(i)*roots[(i%ell)*k%ell]%modulus;
        local.data.length=length+1; _nmod_poly_normalise(local);
        nmod_poly_mullow(prod,prod,local,length+1);
    }
    for(int i=0;i<=length;++i) {
        if(i%ell) assert(!prod.at(i)); else compressed.set(i/ell,prod.at(i));
    }
    for(int i=0;ell*i<=D;++i) scaled.set(ell*i,a.at(i));
    nmod_poly_mullow(out,compressed,scaled,D+1);
}
int main(int argc,char** argv) {
    if(argc!=4) return 2;
    ell=std::stoi(argv[1]); modulus=std::stoul(argv[2]);
    assert((ell==1 || (ell!=83 && n_is_prime(ell))) && n_is_prime(modulus)
           && modulus>3 && modulus!=83 && (modulus-1)%ell==0 && modulus<1000000);
    D=ell==1?57:57*(ell+1);
    length=D*ell;
    for(Poly* f:{&cubic,&branch.a,&branch.b,&inverse_den}) nmod_init(&f->data.mod,modulus);
    std::ifstream in(argv[3]); assert(in);
    make_basis(in,D);
    build_theta();
    for(int id=0;id<2;++id) {
        Poly a,b,finalnorm;
        first_norm(a,id);
        hecke_series(b,a);
        vector<U> z(D+1);
        for(int k=0;k<=D;++k) z[k]=b.at(k);
        norm_from_series(finalnorm,z,false,0);
        std::cout<<id+1<<' '<<D;
        if(nmod_poly_degree(finalnorm)==D) {
            if(ell==1) {
                nmod_poly_make_monic(finalnorm,finalnorm);
                for(int i=0;i<=D;++i) std::cout<<' '<<finalnorm.at(i);
            } else factor_degrees(finalnorm,D,modulus,std::cout);
        }
        std::cout<<'\n';
    }
}
