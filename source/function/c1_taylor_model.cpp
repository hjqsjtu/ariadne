/***************************************************************************
 *            c1_taylor_model.cpp
 *
 *  Copyright 2008--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "function/functional.hpp"
#include "config.h"

#include <iostream>
#include <iomanip>

#include "numeric/float.hpp"
#include "numeric/rational.hpp"

#include "utility/macros.hpp"
#include "utility/exceptions.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/expansion.hpp"
#include "function/polynomial.tpl.hpp"

#include "function/c1_taylor_model.hpp"

#define VOLATILE ;

namespace Ariadne {

template<class PRE, class PR> PRE cast_error_precision(PR) = delete;
template<> inline DP cast_error_precision<DP,DP>(DP pr) { return pr; }
template<> inline DP cast_error_precision<DP,MP>(MP pr) { return DP(); }
template<> inline MP cast_error_precision<MP,MP>(MP pr) { return pr; }

template<class FE, class F> Error<FE> cast_error(Error<F>) = delete;
template<> inline Error<FloatDP> cast_error<FloatDP,FloatDP>(Error<FloatDP> x) { return x; }
template<> inline Error<FloatDP> cast_error<FloatDP,FloatMP>(Error<FloatMP> x) { return Error<FloatDP>(Dyadic(x.raw()),DoublePrecision()); }
template<> inline Error<FloatMP> cast_error<FloatMP,FloatMP>(Error<FloatMP> x) { return x; }

template<class FE, class F> Error<FE> cast_error(PositiveUpperBound<F> x) { return cast_error<FE>(Error<F>(x)); }

template<class FE, class F> FE cast_raw_float(F) = delete;
template<> inline FloatDP cast_raw_float<FloatDP,FloatDP>(FloatDP x) { return x; }
template<> inline FloatDP cast_raw_float<FloatDP,FloatMP>(FloatMP x) { return FloatDP(Dyadic(x.raw()),DoublePrecision()); }
template<> inline FloatMP cast_raw_float<FloatMP,FloatMP>(FloatMP x) { return x; }

static const char* plusminus = u8"\u00B1";

typedef UnivariateC1TaylorModel<FloatDP,FloatDP> UnivariateC1TaylorModelType;

// FIXME: Change FloatDPError so that operations below are not ambiguous
inline FloatDPError operator*(Nat i, FloatDPError e) { return FloatDPError(i,e.precision())*e; }
inline FloatDPError operator*(DegreeType i, FloatDPError e) { return FloatDPError(i,e.precision())*e; }
inline FloatDPError operator*(SizeType i, FloatDPError e) { return FloatDPError(i,e.precision())*e; }
inline FloatMPError operator*(Nat i, FloatMPError e) { return FloatMPError(i,e.precision())*e; }
inline FloatMPError operator*(DegreeType i, FloatMPError e) { return FloatMPError(i,e.precision())*e; }
inline FloatMPError operator*(SizeType i, FloatMPError e) { return FloatMPError(i,e.precision())*e; }

template<class X> struct MultiIndexProduct : CovectorExpression<MultiIndexProduct<X>> {
    MultiIndex const& _a; X const& _x;
    typedef decltype(_a[0]*_x) ScalarType;
    MultiIndexProduct(MultiIndex const& a, X const& x) : _a(a), _x(x) { }
    SizeType size() const {return _a.size(); }
    ScalarType operator[](SizeType i) const { return _a[i]*_x; }
};
template<class FE> MultiIndexProduct<Error<FE>> operator*(MultiIndex const& a, Error<FE> const& e) {
    return MultiIndexProduct<Error<FE>>{a,e};
}


template<class FE> inline C1Errors<FE> operator+(C1Errors<FE> const& e1, C1Errors<FE> const& e2) {
    C1Errors<FE> r(e1._gradient.size(),min(e1.precision(),e2.precision()));
    r._at_zero=e1._at_zero+e2._at_zero, r._uniform=e1._uniform+e2._uniform, r._gradient=e1._gradient+e2._gradient; return r; }
template<class FE> inline C1Errors<FE> operator*(C1Errors<FE> const& e1, C1Errors<FE> const& e2) {
    C1Errors<FE> r(e1._gradient.size(),min(e1.precision(),e2.precision()));
    r._gradient=e1._gradient*e2._uniform+e1._uniform*e2._gradient; r._at_zero=e1._at_zero*e2._at_zero; r._uniform=e1._uniform*=e2._uniform; return r; }

template<class FE> inline C1Errors<FE>& operator+=(C1Errors<FE>& e1, C1Errors<FE> const& e2) {
    e1._at_zero +=e1._at_zero ; e1._uniform += e2._uniform; e1._gradient += e2._gradient; return e1; }
template<class FE> inline C1Errors<FE>& operator+=(C1Errors<FE>& e, Error<SelfType<FE>> const& s) {
    e._at_zero +=s ; e._uniform += s; return e; }
template<class FE> inline C1Errors<FE>& iaddmul(C1Errors<FE>& e, MultiIndex const& a, Error<SelfType<FE>> const& s) {
    e._uniform += s; if(a.degree()==0) { e._at_zero+=s; } else { e._gradient += a*s; } return e; }
template<class FE> inline C1Errors<FE>& operator+=(C1Errors<FE>& e, Pair<MultiIndex const&, Error<SelfType<FE>>const&> ac) {
    return iaddmul(e,ac.first,ac.second); }
template<class FE> inline C1Errors<FE>& operator*=(C1Errors<FE>& e, Error<SelfType<FE>> const& s) {
    e._at_zero *=s ; e._uniform *= s; e._gradient *=s; return e; }
template<class FE> inline C1Errors<FE>& operator*=(C1Errors<FE>& e1, C1Errors<SelfType<FE>> const& e2) {
    e1._gradient*=e2._uniform; e1._gradient+=e1._uniform*e2._gradient; e1._at_zero*=e2._at_zero; e1._uniform*=e2._uniform; return e1; }

template<class FE> inline UnivariateC1Errors<FE>& iaddmul(UnivariateC1Errors<FE>& e, DegreeType const& a, Error<SelfType<FE>> const& s) {
    e._uniform += s; if(a==0) { e._at_zero+=s; } else { e._derivative += a*s; } return e; }

template<class FE> inline UnivariateC1Errors<FE> operator+(UnivariateC1Errors<FE> const& e1, UnivariateC1Errors<FE> const& e2) {
    UnivariateC1Errors<FE> r(min(e1.precision(),e2.precision()));
    r._derivative=e1._derivative+e2._derivative; r._at_zero=e1._at_zero+e2._at_zero; r._uniform=e1._uniform+e2._uniform; return r; }
template<class FE> inline UnivariateC1Errors<FE> operator*(UnivariateC1Errors<FE> const& e1, UnivariateC1Errors<FE> const& e2) {
    UnivariateC1Errors<FE> r(min(e1.precision(),e2.precision()));
    r._derivative=e1._derivative*e2._uniform+e1._uniform*e2._derivative; r._at_zero=e1._at_zero*e2._at_zero; r._uniform=e1._uniform*e2._uniform; return r; }
template<class FE> inline UnivariateC1Errors<FE>& operator+=(UnivariateC1Errors<FE>& e1, UnivariateC1Errors<FE> const& e2) {
    e1._derivative+=e2._derivative; e1._at_zero+=e2._at_zero; e1._uniform+=e2._uniform; return e1; }
template<class FE> inline UnivariateC1Errors<FE>& operator+=(UnivariateC1Errors<FE>& e, Pair<DegreeType const&, Error<SelfType<FE>>const&> ac) {
    return iaddmul(e,ac.first,ac.second); }

template<class FE> OutputStream& operator<<(OutputStream& os, UnivariateC1Errors<FE> const& e) {
    return os << "{@0:" << e._at_zero << ",C0=" << e._uniform << ",C1=" << e._derivative << "}"; }

template<class FE> OutputStream& operator<<(OutputStream& os, C1Errors<FE> const& e) {
    return os << "{@0:" << e._at_zero << ",C0=" << e._uniform << ",C1=" << e._gradient << "}"; }



template<class FLT> FLT add_rnd(FLT const& x1, FLT const& x2) { return add(x1,x2); }
template<class FLT> FLT sub_rnd(FLT const& x1, FLT const& x2) { return sub(x1,x2); }
template<class FLT> FLT mul_rnd(FLT const& x1, FLT const& x2) { return mul(x1,x2); }
template<class FLT> FLT mul_rnd(Nat n1, FLT const& x2) { return mul(n1,x2); }
template<class FLT> FLT fma_rnd(FLT const& x1, FLT const& x2, FLT y) { return fma(x1,x2,std::move(y)); }

double fma_rnd(double x1, double x2, double y) { return fma(x1,x2,y); }

template<class F, class FE>
Value<F> add_err(Value<SelfType<F>> const& v1, Value<F> const& v2, Error<FE>& e) {
    F::set_rounding_downward();
    F l = add_rnd(v1.raw(),v2.raw());
    F::set_rounding_upward();
    F u = add_rnd(v1.raw(),v2.raw());
    F d = hlf(sub_rnd(u,l));
    e.raw() = add_rnd(e.raw(),cast_raw_float<FE>(d));
    F::set_rounding_to_nearest();
    return Value<F>(add_rnd(v1.raw(),v2.raw()));
}

template<class F, class FE>
Value<F> mul_err(Value<SelfType<F>> const& v1, Value<F> const& v2, Error<FE>& e) {
    F::set_rounding_downward();
    F l = mul_rnd(v1.raw(),v2.raw());
    F::set_rounding_upward();
    F u = mul_rnd(v1.raw(),v2.raw());
    F d = hlf(sub_rnd(u,l));
    e.raw() = add_rnd(e.raw(),cast_raw_float<FE>(d));
    F::set_rounding_to_nearest();
    return Value<F>(mul_rnd(v1.raw(),v2.raw()));
}

template<class F, class FE> inline
Value<F> fma_err(Value<F> const& v1, Value<F> const& v2, Value<F> const& w, Error<FE>& e) {
    F::set_rounding_downward();
    F l = fma_rnd(v1.raw(),v2.raw(),w.raw());
    F::set_rounding_upward();
    F u = fma_rnd(v1.raw(),v2.raw(),w.raw());
    F d = hlf(sub_rnd(u,l));
    e.raw() = add(e.raw(),cast_raw_float<FE>(d));
    F::set_rounding_to_nearest();
    return Value<F>(fma_rnd(v1.raw(),v2.raw(),w.raw()));
}

//FloatDP fma(FloatDP x1, FloatDP x2, FloatDP y) { return add(mul(x1,x2),y); }

template<class PR, class PRE>
FloatValue<PR> add_errs(MultiIndex const& a, FloatValue<PR> const& v1, FloatValue<PR> const& v2, FloatError<PRE>& e0, Covector<FloatError<PRE>>& e1) {
    typedef RawFloat<PR> FLT;
    FLT::set_rounding_downward();
    FLT l = add_rnd(v1.raw(),v2.raw());
    FLT::set_rounding_upward();
    FLT u = add_rnd(v1.raw(),v2.raw());
    FLT e = hlf(sub_rnd(u,l));
    e0.raw() = add_rnd(e0.raw(),e);
    for(SizeType i=0; i!=a.size(); ++i) {
        e1[i].raw() = add_rnd(e1[i].raw(),mul_rnd(a[i],e));
    }
    FLT::set_rounding_to_nearest();
    return FloatValue<PR>(add_rnd(v1.raw(),v2.raw()));
}

template<class F> inline Error<F> abssum(UnivariatePolynomial<Value<F>> const& a) {
    ARIADNE_DEBUG_ASSERT(a.size()>=1);
    Error<F> e=PositiveValue<F>(abs(a[0]));
    for(SizeType i=1; i!=a.size(); ++i) {
        e+=abs(a[i]);
    }
    return e;
}

template<class F> inline Error<F> indabssum(UnivariatePolynomial<Value<F>> const& a) {
    ARIADNE_DEBUG_ASSERT(a.size()>=1);
    Error<F> e=0u*Error<F>(abs(a[0]));
    for(SizeType i=1; i!=a.size(); ++i) {
        e+=i*abs(a[i]);
    }
    return e;
}

template<class I, class X> Pair<I const&, X const&> make_monomial(I const& a, X const& c) { return Pair<I const&, X const&>(a,c); }

template<class F, class PRE> UnivariateC1Errors<RawFloatType<PRE>> c1_seminorms(UnivariatePolynomial<Value<F>> const& p, PRE pre) {
    using FE=RawFloatType<PRE>;
    UnivariateC1Errors<FE> r(pre);
//    for (auto iter=p.expansion().begin(); iter!=p.expansion().end(); ++iter) {
    for (DegreeType a=0; a!=p.degree(); ++a) {
        Value<F> const& c=p[a];
        Error<FE> e=cast_error<FE>(mag(c));
        r+=make_monomial(a,e);
    }
    return r;
}

template<class F, class PRE> C1Errors<RawFloatType<PRE>> c1_seminorms(Polynomial<Value<F>> const& p, PRE pre) {
    using FE=RawFloatType<PRE>;
    C1Errors<FE> r(p.argument_size(),pre);
    for (auto iter=p.expansion().begin(); iter!=p.expansion().end(); ++iter) {
        Error<FE> e=cast_positive(abs(iter->coefficient()));
        r+=make_monomial(iter->index(),e);
    }
    return r;
}

template<class F> inline Error<F> c0_norm(UnivariatePolynomial<Value<F>> const& p) {
    return abssum(p);
}
template<class F> Error<F> c0_seminorm(UnivariatePolynomial<Value<F>> const& p) {
    return abssum(p);
}
template<class F> Error<F> c1_seminorm(UnivariatePolynomial<Value<F>> const& p) {
    return indabssum(p);
}

template<class F, class FE> Error<FE> UnivariateC1TaylorModel<F,FE>::_c1_seminorm() const{
    return cast_error<FE>(c1_seminorm(this->_polynomial))+this->_errors._derivative;
}

template<class X, class Y>
static Y polynomial_evaluate(const std::vector<X>& f, const Y& x)
{
    if(f.size()>=2) {
        DegreeType i=f.size()-2;
        Y r=x*f[i+1]+f[i];
        while(i!=0) {
            i=i-1;
            r=r*x;
            r+=f[i];
        }
        return r;
    } else if(f.size()==1) {
        return x*X(0)+f[0];
    } else {
        return x*X(0);
    }
}

template<class F, class FE> UnivariateC1TaylorModel<F,FE>::UnivariateC1TaylorModel(PR pr)
    : _polynomial(1u,Value<F>(0,pr))
    , _errors(cast_error_precision<PRE>(pr))
{
}


template<class F, class FE> UnivariateC1TaylorModel<F,FE>::UnivariateC1TaylorModel(DegreeType d, PR pr)
    : _polynomial(d+1u,FloatValue<PR>(pr))
    , _errors(cast_error_precision<PRE>(pr))
{
}


template<class F, class FE> UnivariateC1TaylorModel<F,FE> UnivariateC1TaylorModel<F,FE>::constant(NumericType c) {
    UnivariateC1TaylorModel<F,FE> result(c.precision());
    result._polynomial[0]=c.value();
    result._errors._at_zero=cast_error<FE>(c.error());
    result._errors._uniform=cast_error<FE>(c.error());
    return result;
}

template<class F, class FE> UnivariateC1TaylorModel<F,FE> UnivariateC1TaylorModel<F,FE>::coordinate(PR pr) {
    UnivariateC1TaylorModel<F,FE> result(1u,pr);
    result._polynomial[1]=1;
    return result;
}

template<class F, class FE> UnivariateC1TaylorModel<F,FE> UnivariateC1TaylorModel<F,FE>::uniform_ball(PR pr) {
    UnivariateC1TaylorModel<F,FE> result(0u,pr);
    result._errors._uniform=1u;
    return result;
}

template<class F, class FE> UnivariateC1TaylorModel<F,FE> UnivariateC1TaylorModel<F,FE>::derivative_ball(PR pr) {
    UnivariateC1TaylorModel<F,FE> result(0u,pr);
    result._errors._derivative=1u;
    return result;
}

template<class F, class FE> PrecisionType<F> UnivariateC1TaylorModel<F,FE>::precision() const {
    assert(not this->_polynomial.empty());
    return this->_polynomial[0].precision();
}

template<class F, class FE> PrecisionType<FE> UnivariateC1TaylorModel<F,FE>::error_precision() const {
    return this->_errors.precision();
}

template<class F, class FE> UnitInterval UnivariateC1TaylorModel<F,FE>::domain() const {
    return UnitInterval();
}

template<class F, class FE> DegreeType UnivariateC1TaylorModel<F,FE>::degree() const {
    return this->_polynomial.size()-1;
}

template<class F, class FE> UnivariateC1Errors<FE> UnivariateC1TaylorModel<F,FE>::_c1_seminorms(UnivariateC1TaylorModel<F,FE> const& f) {
    return c1_seminorms(f.polynomial(),f.error_precision())+f.errors();
}

#define ARIADNE_BOUNDS_INTERVAL_SUM
#if defined ARIADNE_MIDPOINT_INTERVAL_SUM
template<class F, class FE> UnivariateC1TaylorModel<F,FE>& UnivariateC1TaylorModel<F,FE>::_iadd(UnivariateC1TaylorModel<F,FE>& f, Bounds<F> const& ic) {
    F::set_rounding_upward();
    F& fv=f._polynomial[0];
    F c=ic.midpoint();
    F::set_rounding_upward();
    VOLATILE F fvu=fv+c;
    VOLATILE F mfvl=(-fv)-c;
    F e=(fvu+mfvl)/2;
    e+=max(ic.upper()-c,c-ic.lower());
    f._errors._at_zero+=(fvu+mfvl)/2;
    f._errors._uniform+=(fvu+mfvl)/2;
    F::set_rounding_to_nearest();
    fv+=c;
    ARIADNE_ASSERT_MSG(f._errors._at_zero>=0,"f="<<f<<" c="<<c);
    return f;
}

#elif defined ARIADNE_BOUNDS_INTERVAL_SUM
template<class F, class FE> UnivariateC1TaylorModel<F,FE>& UnivariateC1TaylorModel<F,FE>::_iadd(UnivariateC1TaylorModel<F,FE>& f, Bounds<F> const& ic) {
    using PR=typename UnivariateC1TaylorModel<F,FE>::PR; using PRE=typename UnivariateC1TaylorModel<F,FE>::PRE;
    FloatValue<PR>& fv=f._polynomial[0];
    FloatError<PRE>& fe=f._errors._uniform;
    FloatError<PRE> e(fe.precision());
    FloatBall<PR,PRE> c(ic);
    fv=add_err(fv,c.value(),e);
    e+=c.error();
    f._errors._at_zero+=e;
    f._errors._uniform+=e;
    return f;
}
#endif

template<class F, class FE> UnivariateC1TaylorModel<F,FE>& UnivariateC1TaylorModel<F,FE>::_imul(UnivariateC1TaylorModel<F,FE>& f, Bounds<F> const& ic) {
    using PR=typename UnivariateC1TaylorModel<F,FE>::PR; using PRE=typename UnivariateC1TaylorModel<F,FE>::PRE;
    FloatBall<PR,PRE> c(ic);
    const FloatValue<PR> cv=c.value();
    const FloatError<PRE> ce=c.error();

    const PositiveFloatUpperBound<PRE> ac=cast_error<FE>(Error<F>(mag(c)));
    f._errors._at_zero*=ac;
    f._errors._uniform*=ac;
    f._errors._derivative*=ac;

    FloatError<PRE> e(f.error_precision());
    for(DegreeType i=0; i!=f._polynomial.size(); ++i) {
        e=0;
        f._polynomial[i]=mul_err(f._polynomial[i],cv,e);
        f._errors._uniform+=e;
        f._errors._derivative+=i*e;
        if(i==0) { f._errors._at_zero+=e; }
    }
    return f;
}


template<class F, class FE> UnivariateC1TaylorModel<F,FE> UnivariateC1TaylorModel<F,FE>::_add(UnivariateC1TaylorModel<F,FE> const& f1, UnivariateC1TaylorModel<F,FE> const& f2) {
    using PR=typename UnivariateC1TaylorModel<F,FE>::PR; using PRE=typename UnivariateC1TaylorModel<F,FE>::PRE;
    PrecisionType pr=min(f1.precision(),f2.precision());
    UnivariateC1TaylorModel<F,FE> fr(std::max(f1.degree(),f2.degree()),pr);

    FloatError<PRE> e(fr._errors._uniform.precision());
    for(DegreeType i=0u; i!=std::min(f1._polynomial.size(),f2._polynomial.size()); ++i) {
        fr._polynomial[i]=add_err(f1._polynomial[i],f2._polynomial[i],e);
        fr._errors._uniform+=e;
        fr._errors._derivative+=i*e;
        if(i==0) { fr._errors._at_zero=e; }
        e=0u;
    }
    for(DegreeType i=std::min(f1._polynomial.size(),f2._polynomial.size()); i!=f1._polynomial.size(); ++i) {
        fr._polynomial[i]=f1._polynomial[i];
    }
    for(DegreeType i=std::min(f1._polynomial.size(),f2._polynomial.size()); i!=f2._polynomial.size(); ++i) {
        fr._polynomial[i]=f2._polynomial[i];
    }
    fr._errors._at_zero+=(f1._errors._at_zero+f2._errors._at_zero);
    fr._errors._uniform+=f1._errors._uniform+f2._errors._uniform;
    fr._errors._derivative+=f1._errors._derivative+f2._errors._derivative;
    return fr;
}

// (f*g)' - (p*q)' = f'g+fg'-p'q-pq'
//   = f'g-f'q+f'q-p'q + fg'-pg'+pg'-pq'
//   = f'(g-q)+(f'-p')q + (f-p)g'+p(g'-q')
//   = (f'-p')(g-q)+p'(g-q)+(f'-p')q + (f-p)(g'-q')+(f-p)q'+p(g'-q')
//   = (f'-p')(g-q)+(f-p)(g'-q') +p'(g-q)+(f-p)q' + (f'-p')q+p(g'-q')
template<class F, class FE> UnivariateC1TaylorModel<F,FE> UnivariateC1TaylorModel<F,FE>::_mul(UnivariateC1TaylorModel<F,FE> const& f1, UnivariateC1TaylorModel<F,FE> const& f2) {
    using PR=typename UnivariateC1TaylorModel<F,FE>::PR; using PRE=typename UnivariateC1TaylorModel<F,FE>::PRE;
    PrecisionType pr=min(f1.precision(),f2.precision());
    ErrorPrecisionType pre=min(f1.error_precision(),f2.error_precision());
    UnivariateC1TaylorModel<F,FE> fr(f1.degree()+f2.degree(),pr);
    // std::cerr<<"d0="<<fr.degree()<<", d1="<<f1.degree()<<", d2="<<f2.degree()<<"\n";

    FloatError<PRE> e(fr.error_precision());
    for(DegreeType k=0u; k<=fr.degree(); ++k) {
        e=0u;
        for(DegreeType i=std::max(k,f2.degree())-f2.degree(); i<=std::min(k,f1.degree()); ++i) {
            DegreeType j=k-i;
            fr._polynomial[k]=fma_err(f1._polynomial[i],f2._polynomial[j],fr._polynomial[k],e);
        }
        fr._errors += make_monomial(k,e);
    }

    fr._errors += f1.errors()*f2.errors();
    fr._errors += c1_seminorms(f1.polynomial(),pre)*f2.errors()+f1.errors()*c1_seminorms(f2.polynomial(),pre);

    return fr;
}


template<class F, class FE> UnivariateC1TaylorModel<F,FE> UnivariateC1TaylorModel<F,FE>::_compose(UnivariateC1TaylorModel<F,FE> const& f, UnivariateC1TaylorModel<F,FE> const& g) {
    using PR=typename UnivariateC1TaylorModel<F,FE>::PR; using PRE=typename UnivariateC1TaylorModel<F,FE>::PRE;
    DegreeType i=f.degree();

    UnivariateC1TaylorModel<F,FE> h=polynomial_evaluate(f.polynomial(),g);

    UnivariateC1Errors<FE> g_seminorms=c1_seminorms(g);

    // Use the formulae: |fog(x)-pog(x)|<=|f-p|, |fog0-pog0|<=min(|f-p|,|f0-p0|+|f'-p'|*|g0|)
    // |fog'(x)-pog'(x)|=|f'(g(x))g'(x)-p'(g(x))g'(x)|<=|f'(g(x))-p'(g(x))|*|g'(x)| <= |f'-p'|*|g'|
    h._errors._at_zero+=min(f._errors._at_zero+g_seminorms._at_zero*f._errors._derivative, f._errors._uniform);
    h._errors._uniform+=f._errors._uniform;
    h._errors._derivative+=f._errors._derivative*g_seminorms._derivative;

    return h;
}

template<class F> Error<F> operator*(PositiveValue<F> x1, Error<F> x2) { return Error<F>(mul(up,x1._v,x2._e)); }

template<class F, class FE> Bounds<F> UnivariateC1TaylorModel<F,FE>::_evaluate(UnivariateC1TaylorModel<F,FE> const& f, Bounds<F> const& x) {
    DegreeType i=f.degree();
    Bounds<F> r = polynomial_evaluate(f.polynomial(),x);

    // Find the best error value to use at x
    Error<FE> xe = f._errors._at_zero + cast_error<FE>(mag(x))*f._errors._derivative;
    Error<FE> const& ue = f._errors._uniform;
    if(refines(xe,ue)) {
        r+=Ball<F,FE>(F(x.precision()),xe.raw());
    } else {
        r+=Ball<F,FE>(F(x.precision()),ue.raw());
    }
    return r;
}


template<class F, class FE> OutputStream& UnivariateC1TaylorModel<F,FE>::_write(OutputStream& os, const UnivariateC1TaylorModel<F,FE>& f) {
//    os << "(" << f._polynomial[0] << "±" << f._errors._at_zero << "/" << f._errors._uniform << ")+(" << f._polynomial[1] << plusminus << f._errors._derivative <<")*x";
    return os << "(" << f._polynomial << "±" << f._errors << ")";
    os << "(" << f._polynomial[0] << "±" << f._errors._uniform << ")+(" << f._polynomial[1] << plusminus << f._errors._derivative <<")*x";
    for(DegreeType i=2; i<=f._polynomial.size(); ++i) {
        if (f._polynomial[i]>=0) { os << "+"; }
        os << f._polynomial[i] << "*x^" << i;
    }
    return os;
}


template<class F, class FE> C1TaylorModel<F,FE>::C1TaylorModel()
    : C1TaylorModel(0u,PR()) { }

template<class F, class FE> C1TaylorModel<F,FE>::C1TaylorModel(SizeType as, PrecisionType pr)
    : _polynomial(as)
    , _errors(as,PRE())
{
    FloatValue<PR> z(dp);
    MultiIndex ind(as);
    for(SizeType i=0; i!=as; ++i) {
        ind[i]=1;
        _polynomial.expansion().append(ind,z);
        ind[i]=0;
    }
    _polynomial.expansion().append(ind,z);
    _polynomial.expansion().reverse_lexicographic_sort();
}

template<class F, class FE> PrecisionType<F> C1TaylorModel<F,FE>::precision() const {
    return this->_polynomial.begin()->coefficient().precision();
}

template<class F, class FE> PrecisionType<FE> C1TaylorModel<F,FE>::error_precision() const {
    return this->_errors._uniform.precision();
}


template<class F, class FE> C1TaylorModel<F,FE> C1TaylorModel<F,FE>::constant(SizeType as, NumericType c) {
    C1TaylorModel<F,FE> result(as,c.precision());
    result = c;
    return result;
}

template<class F, class FE> C1TaylorModel<F,FE> C1TaylorModel<F,FE>::constant(SizeType as, NumberType c, PrecisionType pr) {
    C1TaylorModel<F,FE> result(as,pr);
    result = c;
    return result;
}

template<class F, class FE> C1TaylorModel<F,FE> C1TaylorModel<F,FE>::coordinate(SizeType as, SizeType j, PrecisionType pr) {
    C1TaylorModel<F,FE> result(as,pr);
    MultiIndex ind=MultiIndex::unit(as,j);
    //result._expansion[ind]=1;
    result._polynomial[ind]=1;
    return result;
}

template<class F, class FE> SizeType C1TaylorModel<F,FE>::argument_size() const {
    return this->_polynomial.argument_size();
}

template<class F, class FE> Void C1TaylorModel<F,FE>::clear() {
    this->_polynomial.expansion().clear();
    this->_errors._at_zero=0u;
    this->_errors._uniform=0u;
    for(SizeType j=0; j!=this->_errors._gradient.size(); ++j) {
        this->_errors._gradient[j]=0u;
    }
}

template<class F, class FE> C1TaylorModel<F,FE>& C1TaylorModel<F,FE>::operator=(NumberType c) {
    return this->operator=(NumericType(c,this->precision()));
}

template<class F, class FE> C1TaylorModel<F,FE>& C1TaylorModel<F,FE>::operator=(NumericType c) {
    this->clear();
    MultiIndex ind=MultiIndex::zero(this->argument_size());
    this->_polynomial[ind]=c.value();
    this->_errors._at_zero=c.error();
    this->_errors._uniform=c.error();
    return *this;
}

template<class F, class FE> C1Errors<FE> C1TaylorModel<F,FE>::_c1_seminorms(C1TaylorModel<F,FE> const& f) {
    return c1_seminorms(f.polynomial(),f.error_precision())+f.errors();
}

template<class F, class FE> C1TaylorModel<F,FE>& C1TaylorModel<F,FE>::_iadd(C1TaylorModel<F,FE>& f, Bounds<F> const& c) {
    Ball<F,FE> nv=Ball<F,FE>(static_cast<Value<F>const&>(f._polynomial.value())+c,f.error_precision());
    f._polynomial.value()=nv.value();
    f._errors._at_zero += nv.error();
    f._errors._uniform += nv.error();
    return f;
}

template<class F, class FE> C1TaylorModel<F,FE>& C1TaylorModel<F,FE>::_imul(C1TaylorModel<F,FE>& f, Bounds<F> const& c) {
    using PR=typename C1TaylorModel<F,FE>::PR; using PRE=typename C1TaylorModel<F,FE>::PRE;
    FloatValue<PR> cv=c.value();
    FloatError<PRE> ce=c.error();
    PositiveFloatUpperBound<PRE> ac=cast_positive(mag(c));
    f._errors *= ac;

    FloatError<PRE> e(f.error_precision());
    for(auto data : f._polynomial) {
        e=0u;
        data.coefficient()=mul_err(data.coefficient(),cv,e);
        MultiIndexReference a=data.index();
        f._errors+=make_monomial(a,e);
    }
    return f;
}



template<class F, class FE> C1TaylorModel<F,FE> C1TaylorModel<F,FE>::_add(C1TaylorModel<F,FE> const& f1, C1TaylorModel<F,FE> const& f2) {
    ARIADNE_PRECONDITION(f1.argument_size()==f2.argument_size());
    const SizeType n=f1.argument_size();
    const PrecisionType pr=min(f1.precision(),f2.precision());
    C1TaylorModel<F,FE> f0(n,pr);
    f0._polynomial.reserve(f1._polynomial.number_of_terms()+f2._polynomial.number_of_terms());

    Expansion<FloatDPValue>::ConstIterator i1=f1._polynomial.begin();
    Expansion<FloatDPValue>::ConstIterator i2=f2._polynomial.begin();
    FloatDPError e(f0.error_precision());
    while(i1!=f1._polynomial.end() && i2!=f2._polynomial.end()) {
        if(i1->key()==i2->key()) {
            e=0u;
            const MultiIndex& a = i1->key();
            f0._polynomial.expansion().append(a,add_err(i1->coefficient(),i2->coefficient(),e));
            f0._errors += std::make_pair(a,e);
        } else if(reverse_lexicographic_less(i1->key(),i2->key())) {
            f0._polynomial.expansion().append(i1->key(),i1->data());
            ++i1;
        } else {
            f0._polynomial.expansion().append(i2->key(),i2->data());
            ++i2;
        }
    }
    while(i1!=f1._polynomial.end()) {
        f0._polynomial.expansion().append(i1->key(),i1->data());
        ++i1;
    }
    while(i2!=f2._polynomial.end()) {
        f0._polynomial.expansion().append(i2->key(),i2->data());
        ++i2;
    }
    f0._errors+=(f1._errors+f2._errors);

    return f0;
}

template<class F, class FE> Void fma(C1TaylorModel<F,FE>& f0, const C1TaylorModel<F,FE>& f1, const C1TaylorModel<F,FE>& f2, const FloatDPValue& c3, const MultiIndex& a3) {
    ARIADNE_PRECONDITION(f0.argument_size()==f1.argument_size());
    ARIADNE_PRECONDITION(f1.argument_size()==f2.argument_size());
    const SizeType n=f1.argument_size();
    MultiIndex a(n);
    f0.clear();
    PositiveFloatDPValue ac3=cast_positive(abs(c3));

    auto iter1=f1._polynomial.begin();
    auto iter2=f2._polynomial.begin();
    while(iter1!=f1._polynomial.end() && iter2!=f2._polynomial.end()) {
        if(iter1->key()==iter2->key()+a3) {
            const MultiIndex& a = iter1->key();
            FloatDP::set_rounding_upward();
            FloatDP fvu=iter1->data().raw()+iter2->data().raw()*c3.raw();
            FloatDP mfvl=(-iter1->data().raw())+iter2->data().raw()*(-c3.raw());
            const FloatDP e=(fvu+mfvl)/2;
            if(a.degree()==0) {
                f0._errors._at_zero.raw()+=e;
            }
            f0._errors._uniform.raw()+=e;
            for(SizeType j=0; j!=n; ++j) {
                f0._errors._gradient[j].raw()+=a[j]*e;
            }
            FloatDP::set_rounding_to_nearest();
            f0._polynomial.expansion().append(a,cast_exact(iter1->data().raw()+iter2->data().raw()*c3.raw()));
            ++iter1;
            ++iter2;
        } else if(reverse_lexicographic_less(iter1->key(),iter2->key()+a3)) {
            f0._polynomial.expansion().append(iter1->key(),iter1->data());
            ++iter1;
        } else {
            a=iter2->key()+a3;
            FloatDP::set_rounding_upward();
            FloatDP fvu=iter2->data().raw()*c3.raw();
            FloatDP mfvl=iter2->data().raw()*(-c3.raw());
            const FloatDP e=(fvu+mfvl)/2;
            if(a.degree()==0) {
                f0._errors._at_zero.raw()+=e;
            }
            f0._errors._uniform.raw()+=e;
            for(SizeType j=0; j!=n; ++j) {
                f0._errors._gradient[j].raw()+=static_cast<DegreeType>(a[j])*e;
            }
            FloatDP::set_rounding_to_nearest();
            f0._polynomial.expansion().append(a,cast_exact(iter2->data().raw()*c3.raw()));
            ++iter2;
        }
    }
    while(iter1!=f1._polynomial.end()) {
        f0._polynomial.expansion().append(iter1->key(),iter1->data());
        ++iter1;
    }
    while(iter2!=f2._polynomial.end()) {
        FloatDP::set_rounding_upward();
        FloatDP fvu=iter2->data().raw()*c3.raw();
        FloatDP mfvl=iter2->data().raw()*(-c3.raw());
        const FloatDP e=(fvu+mfvl)/2;
        if(a.degree()==0) {
            f0._errors._at_zero.raw()+=e;
        }
        f0._errors._uniform.raw()+=e;
        for(SizeType j=0; j!=n; ++j) {
            f0._errors._gradient[j].raw()+=static_cast<DegreeType>(a[j])*e;
        }
        FloatDP::set_rounding_to_nearest();
        f0._polynomial.expansion().append(iter2->key()+a3,cast_exact(iter2->data().raw()*c3.raw()));
        ++iter2;
    }

    FloatDP::set_rounding_upward();
    f0.zero_error().raw()+=(f1.zero_error().raw()+f2.zero_error().raw()*ac3.raw());
    f0.uniform_error().raw()+=(f1.uniform_error().raw()+f2.uniform_error().raw()*ac3.raw());
    for(SizeType j=0; j!=n; ++j) {
        f0.derivative_error(j).raw()+=(f1.derivative_error(j).raw()+f2.derivative_error(j).raw()*ac3.raw());
    }
    FloatDP::set_rounding_to_nearest();

}

template<class F, class FE> C1TaylorModel<F,FE> C1TaylorModel<F,FE>::_mul(C1TaylorModel<F,FE> const& f1, C1TaylorModel<F,FE> const& f2) {
    ARIADNE_PRECONDITION(f1.argument_size()==f2.argument_size());
    const SizeType n=f1.argument_size();
    const PrecisionType pr=min(f1.precision(),f2.precision());
    C1TaylorModel<F,FE> f0a(n,pr);
    f0a._polynomial.expansion().clear();
    C1TaylorModel<F,FE> f0b(f0a);
    f0b.clear();
    C1TaylorModel<F,FE>* ftp=&f0a;
    C1TaylorModel<F,FE>* frp=&f0b;
    for(Expansion<FloatDPValue>::ConstIterator i2=f2._polynomial.begin();
        i2!=f2._polynomial.end(); ++i2)
    {
        fma(*frp,*ftp,f1,i2->data(),i2->key());
        std::swap(ftp,frp);
    }
    return *frp;
}

template<class F, class FE> Bounds<F> C1TaylorModel<F,FE>::_evaluate(C1TaylorModel<F,FE> const& f, Vector<Bounds<F>> const& x) {
    Bounds<F> r=horner_evaluate(f._polynomial.expansion(),x);
    r += Bounds<F>(-f.uniform_error(),+f.uniform_error());
    return r;
}


template<class F, class FE> C1TaylorModel<F,FE> C1TaylorModel<F,FE>::_compose(UnivariateC1TaylorModel<F,FE> const& f, C1TaylorModel<F,FE> const& g) {
    C1TaylorModel<F,FE> r = polynomial_evaluate(f.polynomial(),g);
    // Use the formula |fog,j(x)-pog,j(x)|=|f'(g(x))g,j(x)-p'(g(x))g,j(x)|<=|f'(g(x))-p'(g(x))|*|g,j'(x)| <= |f'-p'|*|g'|
    C1Errors<FE> g_seminorms=c1_seminorms(g);
    r._errors._at_zero+=min(f._errors._at_zero+g_seminorms._at_zero*f._errors._derivative, f._errors._uniform);
    r._errors._uniform+=f._errors._uniform;
    r._errors._gradient+=f._errors._derivative*g_seminorms._gradient;
}

template<class F, class FE> C1TaylorModel<F,FE> C1TaylorModel<F,FE>::_compose(C1TaylorModel<F,FE> const& f, Vector<C1TaylorModel<F,FE>> const& g) {
    // Use the formula |fog,j(x)-pog,j(x)|=|f,i(g(x))gi,j(x)-p,i(g(x))gi,j(x)|<=|f,i(g(x))-p,i(g(x))|*|g,ij(x)| <= |f'-p'|*||
    ARIADNE_NOT_IMPLEMENTED;
    C1TaylorModel<F,FE> gz=f;//g.zero_element();
    C1TaylorModel<F,FE> r(gz.argument_size(),gz.precision());
        //r=horner_evaluate(f._polynomial.expansion(),g);
    std::cerr<<"intermediate="<<r<<"\n";
    r._errors += f._errors;
    // TODO: How do first derivatives change?
    return r;
}

template<class F, class FE> OutputStream& C1TaylorModel<F,FE>::_write(OutputStream& os, C1TaylorModel<F,FE> const& f) {
    return os << f._polynomial << "±" << f._errors; }


template<class T> struct ListForm {
    ListForm(const T& t) : value(t) { } const T& value;
};
template<class T> ListForm<T> list_form(const T& t) { return ListForm<T>(t); }

template<class X> OutputStream& operator<<(OutputStream& os, const ListForm<Expansion<X>>& lfe) {
    const Expansion<X>& e=lfe.value;
    os << "{ ";
    for(auto iter=e.begin(); iter!=e.end(); ++iter)
    {
        if(iter!=e.begin()) { os << ", "; }
        for(SizeType i=0; i!=iter->key().size(); ++i) {
            os << DegreeType(iter->key()[i]);
            if(i+1!=iter->key().size()) { os << ","; }
        }
        os << ":" << iter->data();
    }
    return os << " }";
}

template<class F, class FE> OutputStream& operator<<(OutputStream& os, const C1TaylorModel<F,FE>& f) {
    os << "C1TaylorModel<F,FE>( coefficients=" << list_form(f._polynomial.expansion())
       << ", errors=" << f.errors()
       << ")";
    return os;

    for(auto iter=f._polynomial.expansion().begin();
        iter!=f._polynomial.expansion().end(); ++iter)
    {
        os << *iter;
    }
}


template class UnivariateC1TaylorModel<FloatDP,FloatDP>;
template class UnivariateC1TaylorModel<FloatMP,FloatDP>;
template class C1TaylorModel<FloatDP,FloatDP>;

} // namespace Ariadne
