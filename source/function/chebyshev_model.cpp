/***************************************************************************
 *            chebyshev_model.cpp
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

#include "function/taylor_series.hpp"
#include "function/taylor_model.hpp"
#include "function/chebyshev_model.hpp"

#define VOLATILE ;

namespace Ariadne {

static const char* plusminus = u8"\u00B1";

typedef UnivariateChebyshevModel<FloatDP,FloatDP> UnivariateChebyshevModelType;

FloatDP make_float(FloatDP x, DoublePrecision pr) { return x; }

template<class F, class PRE> Ball<F,RawFloatType<PRE>> mul(Value<F> const& v1, Value<F> const& v2, PRE const& pre) {
    typedef RawFloatType<PRE> FE;
    F::set_rounding_downward();
    F l = mul_rnd(v1.raw(),v2.raw());
    F::set_rounding_upward();
    F u = mul_rnd(v1.raw(),v2.raw());
    F e = hlf(sub_rnd(u,l));
    F::set_rounding_to_nearest();
    F v = mul_rnd(v1.raw(),v2.raw());
    return Ball<F,FE>(v,make_float(e,pre));
}
template<class F, class FE> Value<F> fma_err(Value<F> const& v1, Value<F> const& v2, Value<F> const& w, Error<FE>& e);

template<class F, class FE>
Value<F> add_err(Value<SelfType<F>> const& v1, Value<F> const& v2, Error<FE>& e) {
    F::set_rounding_downward();
    F l = add_rnd(v1.raw(),v2.raw());
    F::set_rounding_upward();
    F u = add_rnd(v1.raw(),v2.raw());
    F d = hlf(sub_rnd(u,l));
    e.raw() = add_rnd(e.raw(),d);
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
    e.raw() = add_rnd(e.raw(),d);
    F::set_rounding_to_nearest();
    return Value<F>(mul_rnd(v1.raw(),v2.raw()));
}


template<class FLT> FLT add_rnd(FLT const& x1, FLT const& x2) { return add(x1,x2); }
template<class FLT> FLT sub_rnd(FLT const& x1, FLT const& x2) { return sub(x1,x2); }
template<class FLT> FLT mul_rnd(FLT const& x1, FLT const& x2) { return mul(x1,x2); }
template<class FLT> FLT mul_rnd(Nat n1, FLT const& x2) { return mul(n1,x2); }
template<class FLT> FLT fma_rnd(FLT const& x1, FLT const& x2, FLT y) { return fma(x1,x2,std::move(y)); }

inline double fma_rnd(double x1, double x2, double y) { return fma(x1,x2,y); }

template<class F, class FE> inline
Value<F> fma_err(Value<F> const& v1, Value<F> const& v2, Value<F> const& w, Error<FE>& e) {
    F::set_rounding_downward();
    F l = fma_rnd(v1.raw(),v2.raw(),w.raw());
    F::set_rounding_upward();
    F u = fma_rnd(v1.raw(),v2.raw(),w.raw());
    e.raw() = add(e.raw(),hlf(sub(u,l)));
    F::set_rounding_to_nearest();
    return Value<F>(fma_rnd(v1.raw(),v2.raw(),w.raw()));
}

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

template<class F> inline Error<F> abssum(UnivariateExpansion<Value<F>> const& a) {
    ARIADNE_DEBUG_ASSERT(a.size()>=1);
    Error<F> e=PositiveValue<F>(abs(a[0]));
    for(SizeType i=1; i!=a.size(); ++i) {
        e+=abs(a[i]);
    }
    return e;
}

template<class F> inline Error<F> norm(UnivariateExpansion<Value<F>> const& p) {
    return abssum(p);
}

template<class X, class Y>
static Y horner_evaluate(const Array<X>& f, const Y& x)
{
    Y z=nul(x);
    Y r=z+f[0];
    Y t=z+1;
    for(SizeType i=1; i!=f.size(); ++i) {
        t=x*t;
        r+=f[i]*t;
    }
    return r;
}

template<class X, class Y>
static Y chebyshev_polynomial_evaluate(const List<X>& f, const Y& x)
{
    Y z=nul(x);
    Y r=z+f[0];
    if(f.size()<=1) { return r; }
    Y tp=z+1;
    Y tc=x;
    Y& tn=tp;
    r+=f[1]*tc;
    for(SizeType i=2; i!=f.size(); ++i) {
        tn=2*x*tc-tp;
        r+=f[i]*tn;
        std::swap(tp,tc);
    }
    return r;
}


template<class F, class FE> UnivariateChebyshevModel<F,FE>::UnivariateChebyshevModel()
    : UnivariateChebyshevModel<F,FE>(PR())
{
}

template<class F, class FE> UnivariateChebyshevModel<F,FE>::UnivariateChebyshevModel(Nat c)
    : UnivariateChebyshevModel<F,FE>(PR())
{
    *this = NumericType(c,this->precision());
}

template<class F, class FE> UnivariateChebyshevModel<F,FE>::UnivariateChebyshevModel(PR pr)
    : _expansion(0u,Value<F>(0,pr))
    , _error(Ariadne::PrecisionType<FE>())
{
}

template<class F, class FE> UnivariateChebyshevModel<F,FE>::UnivariateChebyshevModel(InitializerList<Dyadic> cs, PR pr)
    : _expansion(cs,pr)
    , _error(Ariadne::PrecisionType<FE>())
{
}

template<class F, class FE> UnivariateChebyshevModel<F,FE>::UnivariateChebyshevModel(UnivariateTaylorModel<F> const& tm)
    : UnivariateChebyshevModel<F,FE>(tm.precision())
{
    auto t=UnivariateChebyshevModel<F,FE>::coordinate(tm.precision());
    Vector<UnivariateChebyshevModel<F,FE>> vt(1u,t);
    *this=horner_evaluate(tm.expansion(),vt);
    this->_error+=tm.error();
}


template<class F, class FE> UnivariateChebyshevModel<F,FE>::operator UnivariateTaylorModel<F> () const
{
    ThresholdSweeper<FloatDP> swp(this->precision(),1e-12);
    UnivariateTaylorModel<F> r=chebyshev_polynomial_evaluate(this->_expansion,UnivariateTaylorModel<F>::coordinate(1u,0u,swp));
    r.error()+=this->_error;
    return r;
}

template<class F, class FE> UnivariateChebyshevModel<F,FE>::UnivariateChebyshevModel(DegreeType d, PR pr)
    : _expansion(d,FloatValue<PR>(pr))
    , _error(Ariadne::PrecisionType<FE>())
{
}


template<class F, class FE> UnivariateChebyshevModel<F,FE>& UnivariateChebyshevModel<F,FE>::operator=(ValidatedNumericType c) {
    this->_expansion.resize(1);
    this->_expansion[0]=c.value();
    this->_error=c.error();
    return *this;
}

template<class F, class FE> UnivariateChebyshevModel<F,FE>& UnivariateChebyshevModel<F,FE>::operator=(ValidatedNumber c) {
    return *this = ValidatedNumericType(c,this->precision());
}

template<class F, class FE> UnivariateChebyshevModel<F,FE> UnivariateChebyshevModel<F,FE>::constant(ValidatedNumericType c) {
    UnivariateChebyshevModel<F,FE> result(c.precision());
    result._expansion[0]=c.value();
    result._error=c.error();
    return result;
}

template<class F, class FE> UnivariateChebyshevModel<F,FE> UnivariateChebyshevModel<F,FE>::coordinate(PR pr) {
    UnivariateChebyshevModel<F,FE> result(1u,pr);
    result._expansion[1]=1;
    return result;
}

template<class F, class FE> UnivariateChebyshevModel<F,FE> UnivariateChebyshevModel<F,FE>::ball(PR pr) {
    UnivariateChebyshevModel<F,FE> result(0u,pr);
    result._error=1u;
    return result;
}

template<class F, class FE> PrecisionType<F> UnivariateChebyshevModel<F,FE>::precision() const {
    assert(not this->_expansion.empty());
    return this->_expansion[0].precision();
}

template<class F, class FE> PrecisionType<FE> UnivariateChebyshevModel<F,FE>::error_precision() const {
    return this->_error.precision();
}

template<class F, class FE> ExactIntervalType UnivariateChebyshevModel<F,FE>::domain() const {
    return ExactIntervalType(-1,+1);
}

template<class F, class FE> DegreeType UnivariateChebyshevModel<F,FE>::degree() const {
    assert(this->_expansion.size()>=1u);
    return this->_expansion.size()-1u;
}

template<class F, class FE> Error<FE> UnivariateChebyshevModel<F,FE>::error() const {
    return this->_error;
}

template<class F, class FE> Error<FE> UnivariateChebyshevModel<F,FE>::_norm(UnivariateChebyshevModel<F,FE> const& f) {
    return abssum(f._expansion);
}

#define ARIADNE_BOUNDS_INTERVAL_SUM
#if defined ARIADNE_MIDPOINT_INTERVAL_SUM
template<class F, class FE> UnivariateChebyshevModel<F,FE>& UnivariateChebyshevModel<F,FE>::_iadd(UnivariateChebyshevModel<F,FE>& f, Bounds<F> const& ic) {
    F::set_rounding_upward();
    F& fv=f._expansion[0];
    F c=ic.midpoint();
    F::set_rounding_upward();
    VOLATILE F fvu=fv+c;
    VOLATILE F mfvl=(-fv)-c;
    F e=(fvu+mfvl)/2;
    e+=max(ic.upper()-c,c-ic.lower());
    f._error+=(fvu+mfvl)/2;
    F::set_rounding_to_nearest();
    fv+=c;
    ARIADNE_ASSERT_MSG(f._error>=0,"f="<<f<<" c="<<c);
    return f;
}

#elif defined ARIADNE_BOUNDS_INTERVAL_SUM
template<class F, class FE> UnivariateChebyshevModel<F,FE>& UnivariateChebyshevModel<F,FE>::_iadd(UnivariateChebyshevModel<F,FE>& f, Bounds<F> const& ic) {
    using PR=typename UnivariateChebyshevModel<F,FE>::PR; using PRE=typename UnivariateChebyshevModel<F,FE>::PRE;
    FloatValue<PR>& fv=f._expansion[0];
    FloatError<PRE>& fe=f._error;
    FloatError<PRE> e(fe.precision());
    FloatBall<PR> c(ic);
    fv=add_err(fv,c.value(),e);
    e+=c.error();
    return f;
}
#endif

template<class F, class FE> UnivariateChebyshevModel<F,FE>& UnivariateChebyshevModel<F,FE>::_imul(UnivariateChebyshevModel<F,FE>& f, Bounds<F> const& ic) {
    using PR=typename UnivariateChebyshevModel<F,FE>::PR; using PRE=typename UnivariateChebyshevModel<F,FE>::PRE;
    FloatBall<PR> c(ic);
    const FloatValue<PR> cv=c.value();
    const FloatError<PRE> ce=c.error();
    FloatError<PRE>& fe=f._error;

    const PositiveFloatUpperBound<PRE> ac=mag(c);
    const PositiveFloatUpperBound<PRE> af=norm(f);
    fe*=ac;
    fe+=af*ce;

    FloatError<PRE> e;
    for(DegreeType i=0; i!=f._expansion.size(); ++i) {
        f._expansion[i]=mul_err(f._expansion[i],cv,e);
    }
    fe+=e;
    return f;
}

template<class F, class FE> UnivariateChebyshevModel<F,FE> UnivariateChebyshevModel<F,FE>::_add(UnivariateChebyshevModel<F,FE> const& f1, UnivariateChebyshevModel<F,FE> const& f2) {
    using PR=typename UnivariateChebyshevModel<F,FE>::PR; using PRE=typename UnivariateChebyshevModel<F,FE>::PRE;
    PrecisionType pr=min(f1.precision(),f2.precision());
    UnivariateChebyshevModel<F,FE> fr(std::max(f1.degree(),f2.degree()),pr);

    for(DegreeType i=0u; i!=std::min(f1._expansion.size(),f2._expansion.size()); ++i) {
        fr._expansion[i]=add_err(f1._expansion[i],f2._expansion[i],fr._error);
    }
    for(DegreeType i=std::min(f1._expansion.size(),f2._expansion.size()); i!=f1._expansion.size(); ++i) {
        fr._expansion[i]=f1._expansion[i];
    }
    for(DegreeType i=std::min(f1._expansion.size(),f2._expansion.size()); i!=f2._expansion.size(); ++i) {
        fr._expansion[i]=f2._expansion[i];
    }
    fr._error+=(f1._error+f2._error);
    return fr;
}

template<class F, class FE> UnivariateChebyshevModel<F,FE> UnivariateChebyshevModel<F,FE>::_mul(UnivariateChebyshevModel<F,FE> const& f1, UnivariateChebyshevModel<F,FE> const& f2) {
    using PR=typename UnivariateChebyshevModel<F,FE>::PR; using PRE=typename UnivariateChebyshevModel<F,FE>::PRE;
    PrecisionType pr=min(f1.precision(),f2.precision());
    ErrorPrecisionType pre=min(f1.error_precision(),f2.error_precision());
    UnivariateChebyshevModel<F,FE> fr(f1.degree()+f2.degree(),pr);
//    std::cerr<<"dr="<<fr.degree()<<", d1="<<f1.degree()<<", d2="<<f2.degree()<<"\n";
    Ball<F,FE> tb(pr,pre);
    UnivariateExpansion<Ball<F,FE>> rc(f1.degree()+f2.degree(), tb);
    for(DegreeType i1=0u; i1<=f1.degree(); ++i1) {
        for(DegreeType i2=0u; i2<=f2.degree(); ++i2) {
            tb=mul(f1._expansion[i1],f2._expansion[i2],pre);
            rc[abs(i1-i2)]+=tb;
            rc[i1+i2]+=tb;
        }
    }
    for(DegreeType k=0u; k<=fr.degree(); ++k) {
        fr._expansion[k]=hlf(rc[k].value());
        fr._error+=rc[k].error();
    }
    fr._error=hlf(fr._error);
    return fr;
}

template<class F, class FE> UnivariateChebyshevModel<F,FE> UnivariateChebyshevModel<F,FE>::_compose(UnivariateChebyshevModel<F,FE> const& f, UnivariateChebyshevModel<F,FE> const& g) {
    using PR=typename UnivariateChebyshevModel<F,FE>::PR; using PRE=typename UnivariateChebyshevModel<F,FE>::PRE;
    return chebyshev_polynomial_evaluate(f._expansion,g)+Bounds<FE>(-f._error,+f._error);
}

template<class F> Error<F> operator*(PositiveValue<F> x1, Error<F> x2) { return Error<F>(mul(up,x1._v,x2._e)); }

template<class F, class FE> Bounds<F> UnivariateChebyshevModel<F,FE>::_evaluate(UnivariateChebyshevModel<F,FE> const& f, Bounds<F> const& x) {
    return chebyshev_polynomial_evaluate(f._expansion,x)+Bounds<F>(-f._error,+f._error);
}


template<class F, class FE> OutputStream& UnivariateChebyshevModel<F,FE>::_write(OutputStream& os, const UnivariateChebyshevModel<F,FE>& f) {
//    os << "(" << f._expansion[0] << "±" << f._error << "/" << f._error << ")+(" << f._expansion[1] << plusminus << f._error <<")*x";
    return os << "(" << f._expansion << "±" << f._error << ")";
    os << "(" << f._expansion[0] << "±" << f._error << ")+(" << f._expansion[1] << plusminus << f._error <<")*x";
    for(DegreeType i=2; i<=f._expansion.size(); ++i) {
        if (f._expansion[i]>=0) { os << "+"; }
        os << f._expansion[i] << "*x^" << i;
    }
    return os;
}


template<class F, class FE> ChebyshevModel<F,FE>::ChebyshevModel()
    : ChebyshevModel(0u,PR()) { }

template<class F, class FE> ChebyshevModel<F,FE>::ChebyshevModel(SizeType as, PrecisionType pr)
    : _expansion(as)
    , _error(as,PRE())
{
    FloatValue<PR> z(dp);
    MultiIndex ind(as);
    for(SizeType i=0; i!=as; ++i) {
        ind[i]=1;
        _expansion.expansion().append(ind,z);
        ind[i]=0;
    }
    _expansion.expansion().append(ind,z);
    _expansion.expansion().reverse_lexicographic_sort();
}

template<class F, class FE> PrecisionType<F> ChebyshevModel<F,FE>::precision() const {
    return this->_expansion.begin()->coefficient().precision();
}

template<class F, class FE> PrecisionType<FE> ChebyshevModel<F,FE>::error_precision() const {
    return this->_error.precision();
}


template<class F, class FE> ChebyshevModel<F,FE> ChebyshevModel<F,FE>::constant(SizeType as, NumericType c) {
    ChebyshevModel<F,FE> result(as,c.precision());
    result = c;
    return result;
}

template<class F, class FE> ChebyshevModel<F,FE> ChebyshevModel<F,FE>::constant(SizeType as, GenericNumericType c, PrecisionType pr) {
    ChebyshevModel<F,FE> result(as,pr);
    result = c;
    return result;
}

template<class F, class FE> ChebyshevModel<F,FE> ChebyshevModel<F,FE>::coordinate(SizeType as, SizeType j, PrecisionType pr) {
    ChebyshevModel<F,FE> result(as,pr);
    MultiIndex ind=MultiIndex::unit(as,j);
    //result._expansion[ind]=1;
    result._expansion[ind]=1;
    return result;
}

template<class F, class FE> SizeType ChebyshevModel<F,FE>::argument_size() const {
    return this->_expansion.argument_size();
}

template<class F, class FE> Void ChebyshevModel<F,FE>::clear() {
    this->_expansion.clear();
    this->_error=0u;
}

template<class F, class FE> ChebyshevModel<F,FE>& ChebyshevModel<F,FE>::operator=(GenericNumericType c) {
    return this->operator=(NumericType(c,this->precision()));
}

template<class F, class FE> ChebyshevModel<F,FE>& ChebyshevModel<F,FE>::operator=(NumericType c) {
    this->clear();
    MultiIndex ind=MultiIndex::zero(this->argument_size());
    this->_expansion[ind]=c.value();
    this->_error=c.error();
    return *this;
}

template<class F, class FE> Error<FE> ChebyshevModel<F,FE>::_norm(ChebyshevModel<F,FE> const& f) {
    ARIADNE_NOT_IMPLEMENTED;
    //return abssum(f._expansion);
}

template<class F, class FE> ChebyshevModel<F,FE>& ChebyshevModel<F,FE>::_iadd(ChebyshevModel<F,FE>& f, Bounds<F> const& c) {
    Ball<F,FE> nv=Ball<F,FE>(static_cast<Value<F>const&>(f._expansion.value())+c,f.error_precision());
    f._expansion.value()=nv.value();
    f._error += nv.error();
    return f;
}

template<class F, class FE> ChebyshevModel<F,FE>& ChebyshevModel<F,FE>::_imul(ChebyshevModel<F,FE>& f, Bounds<F> const& c) {
    using PR=typename ChebyshevModel<F,FE>::PR; using PRE=typename ChebyshevModel<F,FE>::PRE;
    FloatValue<PR> cv=c.value();
    FloatError<PRE> ce=c.error();
    PositiveFloatUpperBound<PRE> ac=cast_positive(mag(c));
    PositiveFloatUpperBound<PRE> af=cast_positive(norm(f));
    f._error *= ac;
    f._error += af*ce;

    FloatError<PRE> e(f.error_precision());
    for(auto data : f._expansion) {
        e=0u;
        data.coefficient()=mul_err(data.coefficient(),cv,e);
        MultiIndexReference a=data.index();
    }
    f._error+=e;
    return f;
}



template<class F, class FE> ChebyshevModel<F,FE> ChebyshevModel<F,FE>::_add(ChebyshevModel<F,FE> const& f1, ChebyshevModel<F,FE> const& f2) {
    ARIADNE_PRECONDITION(f1.argument_size()==f2.argument_size());
    const SizeType n=f1.argument_size();
    const PrecisionType pr=min(f1.precision(),f2.precision());
    ChebyshevModel<F,FE> f0(n,pr);
    f0._expansion.reserve(f1._expansion.number_of_terms()+f2._expansion.number_of_terms());

    Expansion<FloatDPValue>::ConstIterator i1=f1._expansion.begin();
    Expansion<FloatDPValue>::ConstIterator i2=f2._expansion.begin();
    FloatDPError e(f0.error_precision());
    while(i1!=f1._expansion.end() && i2!=f2._expansion.end()) {
        if(i1->key()==i2->key()) {
            e=0u;
            const MultiIndex& a = i1->key();
            f0._expansion.expansion().append(a,add_err(i1->coefficient(),i2->coefficient(),e));
            f0._error += e;
        } else if(reverse_lexicographic_less(i1->key(),i2->key())) {
            f0._expansion.expansion().append(i1->key(),i1->data());
            ++i1;
        } else {
            f0._expansion.expansion().append(i2->key(),i2->data());
            ++i2;
        }
    }
    while(i1!=f1._expansion.end()) {
        f0._expansion.expansion().append(i1->key(),i1->data());
        ++i1;
    }
    while(i2!=f2._expansion.end()) {
        f0._expansion.expansion().append(i2->key(),i2->data());
        ++i2;
    }
    f0._error+=(f1._error+f2._error);

    return f0;
}

template<class F, class FE> Void fma(ChebyshevModel<F,FE>& f0, const ChebyshevModel<F,FE>& f1, const ChebyshevModel<F,FE>& f2, const FloatDPValue& c3, const MultiIndex& a3) {
    ARIADNE_NOT_IMPLEMENTED;
}
/*
template<class F, class FE> Void fma(ChebyshevModel<F,FE>& f0, const ChebyshevModel<F,FE>& f1, const ChebyshevModel<F,FE>& f2, const FloatDPValue& c3, const MultiIndex& a3) {
    ARIADNE_PRECONDITION(f0.argument_size()==f1.argument_size());
    ARIADNE_PRECONDITION(f1.argument_size()==f2.argument_size());
    const SizeType n=f1.argument_size();
    MultiIndex a(n);
    f0.clear();
    PositiveFloatDPValue ac3=cast_positive(abs(c3));

    auto iter1=f1._expansion.begin();
    auto iter2=f2._expansion.begin();
    while(iter1!=f1._expansion.end() && iter2!=f2._expansion.end()) {
        if(iter1->key()==iter2->key()+a3) {
            const MultiIndex& a = iter1->key();
            FloatDP::set_rounding_upward();
            FloatDP fvu=iter1->data().raw()+iter2->data().raw()*c3.raw();
            FloatDP mfvl=(-iter1->data().raw())+iter2->data().raw()*(-c3.raw());
            const FloatDP e=(fvu+mfvl)/2;
            if(a.degree()==0) {
                f0._error.raw()+=e;
            }
            f0._error.raw()+=e;
            for(SizeType j=0; j!=n; ++j) {
                f0._error[j].raw()+=a[j]*e;
            }
            FloatDP::set_rounding_to_nearest();
            f0._expansion.expansion().append(a,cast_exact(iter1->data().raw()+iter2->data().raw()*c3.raw()));
            ++iter1;
            ++iter2;
        } else if(reverse_lexicographic_less(iter1->key(),iter2->key()+a3)) {
            f0._expansion.expansion().append(iter1->key(),iter1->data());
            ++iter1;
        } else {
            a=iter2->key()+a3;
            FloatDP::set_rounding_upward();
            FloatDP fvu=iter2->data().raw()*c3.raw();
            FloatDP mfvl=iter2->data().raw()*(-c3.raw());
            const FloatDP e=(fvu+mfvl)/2;
            if(a.degree()==0) {
                f0._error.raw()+=e;
            }
            f0._error.raw()+=e;
            for(SizeType j=0; j!=n; ++j) {
                f0._error[j].raw()+=static_cast<DegreeType>(a[j])*e;
            }
            FloatDP::set_rounding_to_nearest();
            f0._expansion.expansion().append(a,cast_exact(iter2->data().raw()*c3.raw()));
            ++iter2;
        }
    }
    while(iter1!=f1._expansion.end()) {
        f0._expansion.expansion().append(iter1->key(),iter1->data());
        ++iter1;
    }
    while(iter2!=f2._expansion.end()) {
        FloatDP::set_rounding_upward();
        FloatDP fvu=iter2->data().raw()*c3.raw();
        FloatDP mfvl=iter2->data().raw()*(-c3.raw());
        const FloatDP e=(fvu+mfvl)/2;
        if(a.degree()==0) {
            f0._error.raw()+=e;
        }
        f0._error.raw()+=e;
        for(SizeType j=0; j!=n; ++j) {
            f0._error[j].raw()+=static_cast<DegreeType>(a[j])*e;
        }
        FloatDP::set_rounding_to_nearest();
        f0._expansion.expansion().append(iter2->key()+a3,cast_exact(iter2->data().raw()*c3.raw()));
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
*/

template<class F, class FE> ChebyshevModel<F,FE> ChebyshevModel<F,FE>::_mul(ChebyshevModel<F,FE> const& f1, ChebyshevModel<F,FE> const& f2) {
    ARIADNE_PRECONDITION(f1.argument_size()==f2.argument_size());
    const SizeType n=f1.argument_size();
    const PrecisionType pr=min(f1.precision(),f2.precision());
    ChebyshevModel<F,FE> f0a(n,pr);
    f0a._expansion.expansion().clear();
    ChebyshevModel<F,FE> f0b(f0a);
    f0b.clear();
    ChebyshevModel<F,FE>* ftp=&f0a;
    ChebyshevModel<F,FE>* frp=&f0b;
    for(Expansion<FloatDPValue>::ConstIterator i2=f2._expansion.begin();
        i2!=f2._expansion.end(); ++i2)
    {
        fma(*frp,*ftp,f1,i2->data(),i2->key());
        std::swap(ftp,frp);
    }
    return *frp;
}

/*
template<class F, class FE> Bounds<F> ChebyshevModel<F,FE>::_evaluate(ChebyshevModel<F,FE> const& f, Vector<Bounds<F>> const& x) {
    Bounds<F> r=horner_evaluate(f._expansion.expansion(),x);
    r += Bounds<F>(-f.uniform_error(),+f.uniform_error());
    return r;
}
*/

/*
template<class F, class FE> ChebyshevModel<F,FE> ChebyshevModel<F,FE>::_compose(UnivariateChebyshevModel<F,FE> const& f, ChebyshevModel<F,FE> const& g) {
    ChebyshevModel<F,FE> r = g;
    r.clear();

    DegreeType i=f.degree();
    r+=f._expansion[i];

    while(i!=0) {
        r=r*g;
        --i;
        r+=f._expansion[i];
    }

    // TODO: How do first derivatives change?
    ARIADNE_NOT_IMPLEMENTED;
    return r;
}
*/

template<class F, class FE> OutputStream& ChebyshevModel<F,FE>::_write(OutputStream& os, ChebyshevModel<F,FE> const& f) {
    return os << f._expansion << "±" << f._error; }

template<class F, class FE> ChebyshevModel<F,FE> ChebyshevModel<F,FE>::_compose(ChebyshevModel<F,FE> const& f, Vector<ChebyshevModel<F,FE>> const& g) {
    ARIADNE_NOT_IMPLEMENTED;
}

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

template<class F, class FE> OutputStream& operator<<(OutputStream& os, const ChebyshevModel<F,FE>& f) {
    os << "ChebyshevModel<F,FE>( coefficients=" << list_form(f._expansion.expansion())
       << ", zero_error=" << f.zero_error()
       << ", uniform_error=" << f.uniform_error()
       << ", derivative_error=" << f.derivative_error()
       << ")";
    return os;

    for(auto iter=f._expansion.expansion().begin();
        iter!=f._expansion.expansion().end(); ++iter)
    {
        os << *iter;
    }
}


template class UnivariateChebyshevModel<FloatDP,FloatDP>;
//template class UnivariateChebyshevModel<FloatMP,FloatDP>;
template class ChebyshevModel<FloatDP,FloatDP>;

} // namespace Ariadne
