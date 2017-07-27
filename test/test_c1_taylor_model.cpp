#include <iostream>

#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "function/taylor_model.hpp"
#include "function/c1_taylor_model.hpp"
#include "function/chebyshev_model.hpp"

using namespace Ariadne;

using UnivariateC1TaylorModelType = Ariadne::UnivariateC1TaylorModel<FloatDP,FloatDP>;
using C1TaylorModelType = Ariadne::C1TaylorModel<FloatDP,FloatDP>;

Int test_taylor_series() {
    UnivariateC1TaylorModelType f1=UnivariateC1TaylorModelType::constant({2,dp});
    f1._polynomial=UnivariatePolynomial<FloatDPValue>({2,3,5,7},dp);
//    f1._errors._uniform=Positive<Decimal>(0.1_dec);
    f1._errors._uniform=Error<FloatDP>(0.1_dec,dp);
    std::cout << "f1="<<f1<<"\n";

    UnivariateC1TaylorModelType f2=UnivariateC1TaylorModelType::coordinate(dp);
    f2._polynomial=UnivariatePolynomial<FloatDPValue>({1,2,3,4,5},dp);
    f2._errors._uniform=0.01;
    std::cout << "f2="<<f2<<"\n";

    const FloatDPValue c1={3.0_exact,dp};
    const FloatDPBounds c2=1/c1;
    UnivariateC1TaylorModelType f1pc1=f1; f1pc1+=c1;
    std::cout << "f1pc1="<<f1pc1<<"\n";
    UnivariateC1TaylorModelType f1tc2=f1; f1tc2*=c2;
    std::cout << "f1tc2="<<f1tc2<<"\n";


    UnivariateC1TaylorModelType f1pf2=f1+f2;
    std::cout << "f1pf2="<<f1pf2<<"\n";
    UnivariateC1TaylorModelType f1tf2=f1*f2;
    std::cout << "f1tf2="<<f1tf2<<"\n";
    std::cout << "\n\n";
    //std::cout << "f2*f2="<<f2*f2<<"\n";

    FloatDPBounds x(2,dp);
    std::cout << "f1(x)="<<evaluate(f1,x)<<"\n";
}

Polynomial<FloatDPValue> make_polynomial(Expansion<FloatDPValue> e) {
    return Polynomial<FloatDPValue>(e); }
//Expansion<FloatDPValue> make_expansion(Expansion<Int> e, DoublePrecision pr) {
//    return Expansion<FloatDPValue>(e,pr); }

Int test_taylor_function() {
    DoublePrecision pr;

    C1TaylorModelType f1=C1TaylorModelType(2,dp);
//    f1._polynomial=make_polynomial({{{1,1},7},{{1,0},5},{{0,1},3},{{0,0},2}},dp);
    std::cout << "f1="<<f1<<"\n";
    std::cout << "(f1+=13)="<<(f1+=13)<<"\n";
    std::cout << "(f1+=1/5)="<<(f1+=FloatDPValue(0.2))<<"\n";
    std::cout << "(f1*=1/5)="<<(f1*=FloatDPValue(0.2))<<"\n";

    C1TaylorModelType f2=C1TaylorModelType::coordinate(2,1,pr);
    f2._errors._uniform=0.01;
    std::cout << "f2="<<f2<<"\n";

    C1TaylorModelType f3=C1TaylorModelType::constant(2,3.0_exact);
    std::cout << "f3="<<f3<<"\n\n";

    f1._polynomial=make_polynomial({{{{1,1},7},{{1,0},5},{{0,1},3},{{0,0},2}},pr});
//    f2._polynomial=make_polynomial({{{2,0},1.3},{{1,1},1.1},{{1,0},0.7},{{0,2},0.5},{{0,1},0.3},{{0,0},0.2}},pr);
//    f1._polynomial._expansion=make_expansion({{{1,1},7},{{1,0},5},{{0,1},3},{{0,0},2}},pr);
//    f2._polynomial._expansion=make_expansion({{{2,0},1.3},{{1,1},1.1},{{1,0},0.7},{{0,2},0.5},{{0,1},0.3},{{0,0},0.2}},pr);
//    f1._polynomial._expansion.sort(ReverseLexicographicKeyLess());
//    f2._polynomial._expansion.sort(ReverseLexicographicKeyLess());
    std::cout << "f1="<<f1<<"\n";
    std::cout << "f2="<<f2<<"\n";
    C1TaylorModelType f1pf2=f1+f2;
    std::cout << "f1pf2="<<f1pf2<<"\n";
    C1TaylorModelType f1tf2=f1*f2;
    std::cout << "f1tf2="<<f1tf2<<"\n";
    std::cout << "\n\n";
    //std::cout << "f2*f2="<<f2*f2<<"\n";

    Vector<FloatDPValue> x={{2,3},pr};
    std::cout << "x="<<x<<"\n";
    std::cout << "f1(x)="<<evaluate(f1,x)<<"\n";

    UnivariateC1TaylorModelType s(pr);
    s._polynomial=UnivariatePolynomial<FloatDPValue>{{1,2,3,4},pr};
    s/=FloatDPValue(8,pr);
    std::cout << "s="<<s<<"\n";
    std::cout << "f="<<f1<<"\n";
    std::cout << "compose(s,f)="<<compose(s,f1)<<"\n";
    compose(f1,{f1,f2});
}

void test_chebyshev_model() {
//    UnivariateChebyshevModel<FloatDP> cm({2,3,5,7,11},dp);
    UnivariateChebyshevModel<FloatDP> cm({0,0,0,1,1},dp);
    cm*=rec(FloatDPBounds(5,dp));
    std::cout << "cm=" << cm << "\n";
    UnivariateTaylorModel<FloatDP> tm(cm);
    std::cout << "tm=" << tm << "\n";
    FloatDPBounds x(3,dp);
    Vector<FloatDPBounds> v({x});
    std::cout << "v="<<v<<" tm(v)=" << tm(v) << "\n";
    std::cout << "x="<<x<<" cm(x)=" << cm(x) << "\n";
    std::cerr<<"norm(tm)="<<norm(tm)<<", norm(cm)=" << norm(cm) << "\n";
    UnivariateChebyshevModel<FloatDP> cmp(tm);
    std::cout << "cmp=" << cmp << "\n";
    std::cout << "cmp=" << cmp.error().raw()-cm.error().raw() << "\n";
}

Int main() {
    test_taylor_series();
    test_taylor_function();
    test_chebyshev_model();
    std::cout << "Done\n";
}
