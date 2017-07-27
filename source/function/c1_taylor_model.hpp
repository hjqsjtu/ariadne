/***************************************************************************
 *            c1_taylor_model.hpp
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

/*! \file c1_taylor_model.hpp
 *  \brief Over-approximations of continuously-differentiable functions based on Taylor expansions.
 */

#ifndef ARIADNE_C1_TAYLOR_MODEL_HPP
#define ARIADNE_C1_TAYLOR_MODEL_HPP

#include <iosfwd>
#include "utility/declarations.hpp"
#include "utility/container.hpp"
#include "numeric/numeric.hpp"
#include "numeric/float.hpp"
#include "algebra/expansion.hpp"
#include "function/domain.hpp"
#include "function/affine.hpp"
#include "function/polynomial.hpp"

namespace Ariadne {


class MultiIndex;

template<class F, class FE> class UnivariateC1TaylorModel;
template<class F, class FE> class C1TaylorModel;

template<class FE> struct C0Error;
template<class FE> struct UnivariateC1Errors;
template<class FE> struct C1Errors;

template<class X> struct UnivariatePolynomial : public List<X> {
  public:
    template<class... PRS, EnableIf<IsConstructible<X,PRS...>> =dummy> UnivariatePolynomial(DegreeType deg, PRS... prs) : List<X>(deg+1,X(prs...)) { }
    template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>> =dummy> UnivariatePolynomial(InitializerList<Y> vals, PRS... prs) {
        for(auto val : vals) { this->append(X(val,prs...)); } }
    template<class Y, EnableIf<IsAssignable<Y,X>> =dummy> UnivariatePolynomial<X>& operator=(Array<Y> vals) {
        X x=nul(this->at(0u)); this->clear(); for(auto val : vals) { x=val; this->append(x); } }
    DegreeType degree() const { return this->size()-1u; }
};

template<class FE> struct UnivariateC1Errors {
    typedef PrecisionType<FE> PRE;
    UnivariateC1Errors(PRE pr) : _at_zero(pr), _uniform(pr), _derivative(pr) { }
    PRE precision() const { return _uniform.precision(); }
    Error<FE> _at_zero;
    Error<FE> _uniform;
    Error<FE> _derivative;
};

template<class FE> struct C1Errors {
    typedef PrecisionType<FE> PRE;
    C1Errors(SizeType n, PRE pr) : _at_zero(pr), _uniform(pr), _gradient(n,pr) { }
    PRE precision() const { return _uniform.precision(); }
    Error<FE> _at_zero;
    Error<FE> _uniform;
    Covector<Error<FE>> _gradient;
};

template<class F, class PRE> UnivariateC1Errors<RawFloatType<PRE>> c1_seminorms(UnivariatePolynomial<Value<F>> const& p, PRE pre);
template<class F, class PRE> C1Errors<RawFloatType<PRE>> c1_seminorms(Polynomial<Value<F>> const& p, PRE pre);


/*! \ingroup FunctionModelSubModule
 *  \brief A UnivariateC1TaylorModel is a univariate function \f$f:\R\rightarrow\R\f$ on an interval \f$[a,b]\f$ is approximated by polynomial \f$p\f$ with error bounds \f$e_0 \geq |f(c)-p(c)|\f$ and \f$e_1\geq sup_{\xi\in D} |f'(\xi)-p'(\xi)|\f$.
 */
template<class F, class FE> class UnivariateC1TaylorModel
{
    typedef Ariadne::PrecisionType<F> PR;
    typedef Ariadne::PrecisionType<FE> PRE;
  public:
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef UnivariateC1TaylorModel<F,FE> UnivariateC1TaylorModelType;
    typedef Bounds<F> NumericType;
    typedef ValidatedNumber NumberType;
  public:
    UnivariatePolynomial<Value<F>> _polynomial;
    UnivariateC1Errors<FE> _errors;
  private:
    explicit UnivariateC1TaylorModel<F,FE>(DegreeType d, PrecisionType pr);
  public:
    UnivariateC1TaylorModel<F,FE>(PR pr);
    static UnivariateC1TaylorModelType constant(NumericType c);
    static UnivariateC1TaylorModelType constant(NumberType, PrecisionType pr);
    static UnivariateC1TaylorModelType coordinate(PrecisionType pr);
    static UnivariateC1TaylorModelType uniform_ball(PrecisionType pr);
    static UnivariateC1TaylorModelType derivative_ball(PrecisionType pr);
  public:
    UnitInterval domain() const;
    DegreeType degree() const;
    Void sweep(RawFloat<PR> threshold);
    FloatValue<PR> coefficient(SizeType i) const;
    UnivariatePolynomial<Value<F>> const& polynomial() const { return _polynomial; }
    UnivariateC1Errors<FE> const& errors() const { return _errors; }
    PrecisionType precision() const;
    ErrorPrecisionType error_precision() const;
  public:
    friend Error<FE> c0_norm(UnivariateC1TaylorModelType const& f) { return f._c0_norm(); }
    friend Error<FE> c1_seminorm(UnivariateC1TaylorModelType const& f) { return f._c1_seminorm(); }
    friend UnivariateC1TaylorModelType& operator+=(UnivariateC1TaylorModelType& f, NumericType c) { return _iadd(f,c); }
    friend UnivariateC1TaylorModelType& operator-=(UnivariateC1TaylorModelType& f, NumericType c) { return _iadd(f,neg(c)); }
    friend UnivariateC1TaylorModelType& operator*=(UnivariateC1TaylorModelType& f, NumericType c) { return _imul(f,c); }
    friend UnivariateC1TaylorModelType& operator/=(UnivariateC1TaylorModelType& f, NumericType c) { return _imul(f,rec(c)); }
    friend UnivariateC1TaylorModelType operator+(NumericType c, UnivariateC1TaylorModelType f) { _iadd(f,c); return std::move(f); }
    friend UnivariateC1TaylorModelType operator*(NumericType c, UnivariateC1TaylorModelType f) { _imul(f,c); return std::move(f); }
    friend UnivariateC1TaylorModelType operator+(UnivariateC1TaylorModelType f, NumericType c) { _iadd(f,c); return std::move(f); }
    friend UnivariateC1TaylorModelType operator-(UnivariateC1TaylorModelType f, NumericType c) { _iadd(f,neg(c)); return std::move(f); }
    friend UnivariateC1TaylorModelType operator*(UnivariateC1TaylorModelType f, NumericType c) { _imul(f,c); return std::move(f); }
    friend UnivariateC1TaylorModelType operator/(UnivariateC1TaylorModelType f, NumericType c) { _imul(f,rec(c)); return std::move(f); }
    friend UnivariateC1TaylorModelType operator+(NumberType c, UnivariateC1TaylorModelType f) { _iadd(f,c); return std::move(f); }
    friend UnivariateC1TaylorModelType operator*(NumberType c, UnivariateC1TaylorModelType f) { _imul(f,c); return std::move(f); }
    friend UnivariateC1TaylorModelType operator+(UnivariateC1TaylorModelType f, NumberType c) { _iadd(f,c); return std::move(f); }
    friend UnivariateC1TaylorModelType operator-(UnivariateC1TaylorModelType f, NumberType c) { _iadd(f,neg(c)); return std::move(f); }
    friend UnivariateC1TaylorModelType operator*(UnivariateC1TaylorModelType f, NumberType c) { _imul(f,c); return std::move(f); }
    friend UnivariateC1TaylorModelType operator/(UnivariateC1TaylorModelType f, NumberType c) { _imul(f,rec(c)); return std::move(f); }
    friend UnivariateC1TaylorModelType operator+(UnivariateC1TaylorModelType const& f1, UnivariateC1TaylorModelType const& f2) { return _add(f1,f2); }
    friend UnivariateC1TaylorModelType operator*(UnivariateC1TaylorModelType const& f1, UnivariateC1TaylorModelType const& f2) { return _add(f1,f2); }

    friend UnivariateC1Errors<FE> c1_seminorms(UnivariateC1TaylorModelType const& f) { return _c1_seminorms(f); }
    friend ValidatedNumericType evaluate(UnivariateC1TaylorModelType const& f, ValidatedNumericType const& c) { return _evaluate(f,c); }
    friend UnivariateC1TaylorModelType compose(UnivariateC1TaylorModelType const& f, UnivariateC1TaylorModelType const& g) { return _compose(f,g); }
    friend OutputStream& operator<< (OutputStream& os, const UnivariateC1TaylorModelType& f) { return _write(os,f); }
  public:
    Error<FE> _c0_norm() const;
    Error<FE> _c1_seminorm() const;
    static UnivariateC1Errors<FE> _c1_seminorms(UnivariateC1TaylorModel<F,FE> const&);
    static UnivariateC1TaylorModel<F,FE>& _iadd(UnivariateC1TaylorModel<F,FE>&, Bounds<F> const&);
    static UnivariateC1TaylorModel<F,FE>& _imul(UnivariateC1TaylorModel<F,FE>&, Bounds<F> const&);
    static UnivariateC1TaylorModel<F,FE> _add(UnivariateC1TaylorModel<F,FE> const&, UnivariateC1TaylorModel<F,FE> const&);
    static UnivariateC1TaylorModel<F,FE> _mul(UnivariateC1TaylorModel<F,FE> const&, UnivariateC1TaylorModel<F,FE> const&);
    static Bounds<F> _evaluate(UnivariateC1TaylorModel<F,FE> const&, Bounds<F> const&);
    static UnivariateC1TaylorModel<F,FE> _compose(UnivariateC1TaylorModel<F,FE> const&, UnivariateC1TaylorModel<F,FE> const&);
    static OutputStream& _write(OutputStream& os, const UnivariateC1TaylorModel<F,FE>& f);
};

template<class F, class FE> class C1TaylorModel
{
    typedef Ariadne::PrecisionType<F> PR; typedef Ariadne::PrecisionType<FE> PRE;
  public:
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef Bounds<F> NumericType;
    typedef ValidatedNumber NumberType;
    typedef UnivariateC1TaylorModel<F,FE> UnivariateC1TaylorModelType;
    typedef C1TaylorModel<F,FE> C1TaylorModelType;
  public:
    Polynomial<Value<F>> _polynomial;
    C1Errors<FE> _errors;
  public:
    C1TaylorModel();
    C1TaylorModel(SizeType as, PR pr);
  public:
    static C1TaylorModel constant(SizeType as, NumericType c);
    static C1TaylorModel constant(SizeType as, NumberType c, PrecisionType pr);
    static C1TaylorModel coordinate(SizeType as, SizeType ind, PrecisionType pr);
  public:
    UnitBox domain() const;
    SizeType argument_size() const;
    PrecisionType precision() const;
    ErrorPrecisionType error_precision() const;
    Void sweep(RawFloat<PR> threshold);
    C1TaylorModel& operator=(NumericType c);
    C1TaylorModel& operator=(NumberType c);
    Void clear();
    Polynomial<Value<F>> const& polynomial() const { return _polynomial; }
    C1Errors<FE> const& errors() const { return _errors; }
  public:
    Error<FE> const& zero_error() const { return _errors._at_zero; }
    Error<FE> const& uniform_error() const { return _errors._uniform; }
    Error<FE> const& derivative_error(SizeType j) const { return _errors._gradient[j]; }
    Covector<Error<FE>> const& derivative_errors() const { return _errors._gradient; }
    Error<FE>& zero_error() { return _errors._at_zero; }
    Error<FE>& uniform_error() { return _errors._uniform; }
    Error<FE>& derivative_error(SizeType j) { return _errors._gradient[j]; }
  public:
    friend C1TaylorModelType& operator+=(C1TaylorModelType& f, ValidatedNumber c) { return _iadd(f,ValidatedNumericType(c,f.precision())); }
    friend C1TaylorModelType& operator*=(C1TaylorModelType& f, ValidatedNumber c) { return _imul(f,ValidatedNumericType(c,f.precision())); }
    friend C1TaylorModelType& operator+=(C1TaylorModelType& f, NumericType c) { return _iadd(f,c); }
    friend C1TaylorModelType& operator*=(C1TaylorModelType& f, NumericType c) { return _imul(f,c); }
    friend C1TaylorModelType operator+(NumericType const& c, C1TaylorModelType f) { _iadd(f,c); return std::move(f); }
    friend C1TaylorModelType operator*(NumericType const& c, C1TaylorModelType f) { _imul(f,c); return std::move(f); }
    friend C1TaylorModelType operator+(C1TaylorModelType f, NumericType c) { _iadd(f,c); return std::move(f); }
    friend C1TaylorModelType operator-(C1TaylorModelType f, NumericType c) { _iadd(f,neg(c)); return std::move(f); }
    friend C1TaylorModelType operator*(C1TaylorModelType f, NumericType c) { _imul(f,c); return std::move(f); }
    friend C1TaylorModelType operator/(C1TaylorModelType f, NumericType c) { _imul(f,rec(c)); return std::move(f); }
    friend C1TaylorModelType operator+(C1TaylorModelType f1, C1TaylorModelType f2) { return _add(f1,f2); }
    friend C1TaylorModelType operator*(C1TaylorModelType f1, C1TaylorModelType f2) { return _mul(f1,f2); }

    friend C1Errors<FE> c1_seminorms(C1TaylorModel<F,FE> const& f) { return _c1_seminorms(f); }
    friend NumericType evaluate(C1TaylorModelType f, Vector<NumericType> x) { return _evaluate(f,x); }
    friend C1TaylorModelType compose(UnivariateC1TaylorModelType const& f, C1TaylorModelType g) { return _compose(f,g); }
    friend C1TaylorModelType compose(C1TaylorModelType f, Vector<C1TaylorModelType> g) { return _compose(f,g); }
    friend OutputStream& operator<< (OutputStream& os, const C1TaylorModelType& f) { return _write(os,f); }
  public:
    static C1Errors<FE> _c1_seminorms(C1TaylorModel<F,FE> const&);
    static C1TaylorModel<F,FE>& _iadd(C1TaylorModel<F,FE>&, Bounds<F> const&);
    static C1TaylorModel<F,FE>& _imul(C1TaylorModel<F,FE>&, Bounds<F> const&);
    static C1TaylorModel<F,FE> _add(C1TaylorModel<F,FE> const&, C1TaylorModel<F,FE> const&);
    static C1TaylorModel<F,FE> _mul(C1TaylorModel<F,FE> const&, C1TaylorModel<F,FE> const&);
    static Bounds<F> _evaluate(C1TaylorModel<F,FE> const&, Vector<Bounds<F>> const&);
    static C1TaylorModel<F,FE> _compose(UnivariateC1TaylorModel<F,FE> const&, C1TaylorModel<F,FE> const&);
    static C1TaylorModel<F,FE> _compose(C1TaylorModel<F,FE> const&, Vector<C1TaylorModel<F,FE>> const&);
    static OutputStream& _write(OutputStream& os, const C1TaylorModel<F,FE>& f);
};

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_HPP
