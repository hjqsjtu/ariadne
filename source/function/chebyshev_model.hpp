/***************************************************************************
 *            chebyshev_model.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file chebyshev_model.hpp
 *  \brief Over-approximations of continuously-differentiable functions based on Chebyshev polynomials.
 */

#ifndef ARIADNE_CHEBYSHEV_FUNCTION_HPP
#define ARIADNE_CHEBYSHEV_FUNCTION_HPP

#include <iosfwd>
#include "utility/declarations.hpp"
#include "utility/container.hpp"
#include "numeric/numeric.hpp"
#include "numeric/float.hpp"
#include "algebra/expansion.hpp"
#include "function/domain.hpp"

namespace Ariadne {

class MultiIndex;

template<class F, class FE=F> class UnivariateChebyshevModel;
template<class F, class FE=F> class ChebyshevModel;

template<class P, class F> class TaylorModel;
template<class F> using UnivariateTaylorModel = TaylorModel<ValidatedTag,F>;


template<class X> struct UnivariateExpansion : public List<X> {
  public:
    template<class... PRS, EnableIf<IsConstructible<X,PRS...>> =dummy> UnivariateExpansion(DegreeType deg, PRS... prs) : List<X>(deg+1,X(prs...)) { }
    template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>> =dummy> UnivariateExpansion(InitializerList<Y> vals, PRS... prs) {
        for(auto val : vals) { this->append(X(val,prs...)); } }
    template<class Y, EnableIf<IsAssignable<Y,X>> =dummy> UnivariateExpansion<X>& operator=(Array<Y> vals) {
        X x=nul(this->at(0u)); this->clear(); for(auto val : vals) { x=val; this->append(x); } }
    DegreeType degree() const { return this->size()-1u; }
    X const& operator[](SizeType i) const { if(i>=this->size()) { std::cerr<<"cs="<<*this<<" i="<<i<<"\n"; assert(i<=this->degree()); } return this->List<X>::operator[](i); }
    X& operator[](SizeType i) { if(i>=this->size()) { std::cerr<<"cs="<<*this<<" i="<<i<<"\n"; assert(i<=this->degree()); } return this->List<X>::operator[](i); }
};

template<class F, class FE=F> class UnivariateTaylorMod
{
    typedef Ariadne::PrecisionType<F> PR;
    typedef Ariadne::PrecisionType<FE> PRE;
    UnivariateExpansion<FloatValue<PR>> _expansion;
    FloatError<PRE> _error;
  public:

};

/*! \ingroup FunctionModelSubModule
 *  \brief A UnivariateChebyshevModel is a univariate function \f$f:\R\rightarrow\R\f$ on an interval \f$[a,b]\f$ is approximated by polynomial \f$p\f$ with error bounds \f$e_0 \geq |f(c)-p(c)|\f$ and \f$e_1\geq sup_{\xi\in D} |f'(\xi)-p'(\xi)|\f$.
 */
template<class F, class FE> class UnivariateChebyshevModel
{
    typedef Ariadne::PrecisionType<F> PR;
    typedef Ariadne::PrecisionType<FE> PRE;
  public:
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef Bounds<F> ValidatedNumericType;
    typedef Bounds<F> NumericType;
    typedef UnivariateChebyshevModel<F,FE> UnivariateChebyshevModelType;
  public:
    UnivariateExpansion<FloatValue<PR>> _expansion;
    FloatError<PRE> _error;
  private:
    explicit UnivariateChebyshevModel<F,FE>(DegreeType d, PrecisionType pr);
  public:
    UnivariateChebyshevModel<F,FE>(); //< Deprecated
    explicit UnivariateChebyshevModel<F,FE>(Nat); //< Deprecated
    UnivariateChebyshevModel<F,FE>(PR pr);
    UnivariateChebyshevModel<F,FE>(InitializerList<Dyadic> c, PR pr);
    UnivariateChebyshevModel<F,FE>(UnivariateTaylorModel<F> const& tm);
    operator UnivariateTaylorModel<F>() const;
    UnivariateChebyshevModel<F,FE>& operator=(ValidatedNumber c);
    UnivariateChebyshevModel<F,FE>& operator=(ValidatedNumericType c);
    static UnivariateChebyshevModelType constant(ValidatedNumericType c);
    static UnivariateChebyshevModelType constant(ValidatedNumber, PrecisionType pr);
    static UnivariateChebyshevModelType coordinate(PrecisionType pr);
    static UnivariateChebyshevModelType ball(PrecisionType pr);
  public:
    IntervalDomainType domain() const;
    DegreeType degree() const;
    Void sweep(RawFloat<PR> threshold);
    Value<F> coefficient(SizeType i) const;
    Error<FE> error() const;
    PrecisionType precision() const;
    ErrorPrecisionType error_precision() const;
  public:
    NumericType operator()(NumericType const& x) const { return evaluate(*this,x); }
  public:
    friend FloatError<PRE> norm(UnivariateChebyshevModelType const& f) { return _norm(f); }
    friend UnivariateChebyshevModelType nul(UnivariateChebyshevModelType f) { f=0; return std::move(f); }
    friend UnivariateChebyshevModelType neg(UnivariateChebyshevModelType f) { _ineg(f); return std::move(f); }
    friend UnivariateChebyshevModelType operator-(UnivariateChebyshevModelType const& f) { _ineg(f); return std::move(f); }
    friend UnivariateChebyshevModelType operator+(UnivariateChebyshevModelType const& f1, UnivariateChebyshevModelType const& f2) { return _add(f1,f2); }
    friend UnivariateChebyshevModelType operator-(UnivariateChebyshevModelType const& f1, UnivariateChebyshevModelType const& f2) { return f1+(-1)*f2; }
    friend UnivariateChebyshevModelType operator*(UnivariateChebyshevModelType const& f1, UnivariateChebyshevModelType const& f2) { return _mul(f1,f2); }
    friend UnivariateChebyshevModelType& operator+=(UnivariateChebyshevModelType& f1, UnivariateChebyshevModelType const& f2) { return f1=f1+f2; }
    friend UnivariateChebyshevModelType& operator-=(UnivariateChebyshevModelType& f1, UnivariateChebyshevModelType const& f2) { return f1=f1-f2; }
    friend UnivariateChebyshevModelType& operator+=(UnivariateChebyshevModelType& f, ValidatedNumericType c) { return _iadd(f,c); }
    friend UnivariateChebyshevModelType& operator-=(UnivariateChebyshevModelType& f, ValidatedNumericType c) { return _iadd(f,neg(c)); }
    friend UnivariateChebyshevModelType& operator*=(UnivariateChebyshevModelType& f, ValidatedNumericType c) { return _imul(f,c); }
    friend UnivariateChebyshevModelType& operator/=(UnivariateChebyshevModelType& f, ValidatedNumericType c) { return _imul(f,rec(c)); }
    friend UnivariateChebyshevModelType operator+(ValidatedNumericType c, UnivariateChebyshevModelType f) { _iadd(f,c); return std::move(f); }
    friend UnivariateChebyshevModelType operator*(ValidatedNumericType c, UnivariateChebyshevModelType f) { _imul(f,c); return std::move(f); }
    friend UnivariateChebyshevModelType operator+(UnivariateChebyshevModelType f, ValidatedNumericType c) { _iadd(f,c); return std::move(f); }
    friend UnivariateChebyshevModelType operator-(UnivariateChebyshevModelType f, ValidatedNumericType c) { _iadd(f,neg(c)); return std::move(f); }
    friend UnivariateChebyshevModelType operator*(UnivariateChebyshevModelType f, ValidatedNumericType c) { _imul(f,c); return std::move(f); }
    friend UnivariateChebyshevModelType operator/(UnivariateChebyshevModelType f, ValidatedNumericType c) { _imul(f,rec(c)); return std::move(f); }
    friend UnivariateChebyshevModelType& operator+=(UnivariateChebyshevModelType& f, ValidatedNumber c) { return _iadd(f,f.make_concrete(c)); }
    friend UnivariateChebyshevModelType& operator-=(UnivariateChebyshevModelType& f, ValidatedNumber c) { return _iadd(f,neg(f.make_concrete(c))); }
    friend UnivariateChebyshevModelType& operator*=(UnivariateChebyshevModelType& f, ValidatedNumber c) { return _imul(f,f.make_concrete(c)); }
    friend UnivariateChebyshevModelType& operator/=(UnivariateChebyshevModelType& f, ValidatedNumber c) { return _imul(f,rec(f.make_concrete(c))); }
    friend UnivariateChebyshevModelType operator+(ValidatedNumber c, UnivariateChebyshevModelType f) { _iadd(f,f.make_concrete(c)); return std::move(f); }
    friend UnivariateChebyshevModelType operator*(ValidatedNumber c, UnivariateChebyshevModelType f) { _imul(f,f.make_concrete(c)); return std::move(f); }
    friend UnivariateChebyshevModelType operator+(UnivariateChebyshevModelType f, ValidatedNumber c) { _iadd(f,f.make_concrete(c)); return std::move(f); }
    friend UnivariateChebyshevModelType operator-(UnivariateChebyshevModelType f, ValidatedNumber c) { _iadd(f,f.make_concrete(neg(c))); return std::move(f); }
    friend UnivariateChebyshevModelType operator*(UnivariateChebyshevModelType f, ValidatedNumber c) { _imul(f,f.make_concrete(c)); return std::move(f); }
    friend UnivariateChebyshevModelType operator/(UnivariateChebyshevModelType f, ValidatedNumber c) { _imul(f,f.make_concrete(rec(c))); return std::move(f); }
    friend ValidatedNumericType evaluate(UnivariateChebyshevModelType const& f, ValidatedNumericType const& c) { return _evaluate(f,c); }
    friend UnivariateChebyshevModelType compose(UnivariateChebyshevModelType const& f, UnivariateChebyshevModelType const& g) { return _compose(f,g); }
    friend OutputStream& operator<< (OutputStream& os, const UnivariateChebyshevModelType& f) { return _write(os,f); }
  public:
    ValidatedNumericType make_concrete(ValidatedNumber c) { return ValidatedNumericType(c,this->precision()); }
    static Error<FE> _norm(UnivariateChebyshevModel<F,FE> const&);
    static UnivariateChebyshevModel<F,FE>& _iadd(UnivariateChebyshevModel<F,FE>&, Bounds<F> const&);
    static UnivariateChebyshevModel<F,FE>& _imul(UnivariateChebyshevModel<F,FE>&, Bounds<F> const&);
    static UnivariateChebyshevModel<F,FE>& _ineg(UnivariateChebyshevModel<F,FE>&);
    static UnivariateChebyshevModel<F,FE> _add(UnivariateChebyshevModel<F,FE> const&, UnivariateChebyshevModel<F,FE> const&);
    static UnivariateChebyshevModel<F,FE> _mul(UnivariateChebyshevModel<F,FE> const&, UnivariateChebyshevModel<F,FE> const&);
    static Bounds<F> _evaluate(UnivariateChebyshevModel<F,FE> const&, Bounds<F> const&);
    static UnivariateChebyshevModel<F,FE> _compose(UnivariateChebyshevModel<F,FE> const&, UnivariateChebyshevModel<F,FE> const&);
    static OutputStream& _write(OutputStream& os, const UnivariateChebyshevModel<F,FE>& f);
};


template<class F, class FE> class ChebyshevModel
{
    typedef Ariadne::PrecisionType<F> PR; typedef Ariadne::PrecisionType<FE> PRE;
  public:
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef FloatBounds<PR> NumericType;
    typedef ValidatedNumber GenericNumericType;
    typedef UnivariateChebyshevModel<F,FE> UnivariateChebyshevModelType;
    typedef ChebyshevModel<F,FE> ChebyshevModelType;
  public:
    Polynomial<FloatValue<PR>> _expansion;
    FloatError<PRE> _error;
  public:
    ChebyshevModel();
    ChebyshevModel(SizeType as, PR pr);
  public:
    static ChebyshevModel constant(SizeType as, NumericType c);
    static ChebyshevModel constant(SizeType as, GenericNumericType c, PrecisionType pr);
    static ChebyshevModel coordinate(SizeType as, SizeType ind, PrecisionType pr);
  public:
    BoxDomainType domain() const;
    SizeType argument_size() const;
    PrecisionType precision() const;
    ErrorPrecisionType error_precision() const;
    Void sweep(RawFloat<PR> threshold);
    ChebyshevModel& operator=(NumericType c);
    ChebyshevModel& operator=(GenericNumericType c);
    Void clear();
    FloatError<PRE> const& error() const { return _error; }
    FloatError<PRE>& error() { return _error; }
  public:
    friend PositiveUpperBound<FE> norm(ChebyshevModelType const& f) { return _norm(f); }
    friend ChebyshevModelType& operator+=(ChebyshevModelType& f, ValidatedNumber c) { return _iadd(f,ValidatedNumericType(c,f.precision())); }
    friend ChebyshevModelType& operator+=(ChebyshevModelType& f, NumericType c) { return _iadd(f,c); }
    friend ChebyshevModelType& operator*=(ChebyshevModelType& f, NumericType c) { return _imul(f,c); }
    friend ChebyshevModelType operator+(ChebyshevModelType f1, ChebyshevModelType f2) { return _add(f1,f2); }
    friend ChebyshevModelType operator*(ChebyshevModelType f1, ChebyshevModelType f2) { return _mul(f1,f2); }
    friend NumericType evaluate(ChebyshevModelType f, Vector<NumericType> x) { return _evaluate(f,x); }
    friend ChebyshevModelType compose(UnivariateChebyshevModelType const& f, ChebyshevModelType g) { return _compose(f,g); }
    friend ChebyshevModelType compose(ChebyshevModelType f, Vector<ChebyshevModelType> g) { return _compose(f,g); }
    friend OutputStream& operator<< (OutputStream& os, const ChebyshevModelType& f) { return _write(os,f); }
  public:
    static Error<FE> _norm(ChebyshevModel<F,FE> const&);
    static ChebyshevModel<F,FE>& _iadd(ChebyshevModel<F,FE>&, Bounds<F> const&);
    static ChebyshevModel<F,FE>& _imul(ChebyshevModel<F,FE>&, Bounds<F> const&);
    static ChebyshevModel<F,FE> _add(ChebyshevModel<F,FE> const&, ChebyshevModel<F,FE> const&);
    static ChebyshevModel<F,FE> _mul(ChebyshevModel<F,FE> const&, ChebyshevModel<F,FE> const&);
    static Bounds<F> _evaluate(ChebyshevModel<F,FE> const&, Vector<Bounds<F>> const&);
    static ChebyshevModel<F,FE> _compose(UnivariateChebyshevModel<F,FE> const&, ChebyshevModel<F,FE> const&);
    static ChebyshevModel<F,FE> _compose(ChebyshevModel<F,FE> const&, Vector<ChebyshevModel<F,FE>> const&);
    static OutputStream& _write(OutputStream& os, const ChebyshevModel<F,FE>& f);
};

} // namespace Ariadne

#endif // ARIADNE_CHEBYSHEV_FUNCTION_HPP
