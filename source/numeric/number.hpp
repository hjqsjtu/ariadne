/***************************************************************************
 *            numeric/number.hpp
 *
 *  Copyright 2013-17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file numeric/number.hpp
 *  \brief Generic numbers
 */



#ifndef ARIADNE_NUMBER_HPP
#define ARIADNE_NUMBER_HPP

#include "../utility/handle.hpp"
#include "../numeric/paradigm.hpp"
#include "../utility/prototype.hpp"

#include "logical.decl.hpp"
#include "number.decl.hpp"
#include "float.decl.hpp"

#include "number_interface.hpp"

#include "arithmetic.hpp"
#include "integer.hpp"
#include "rational.hpp"
#include "real.hpp"

#include "number_interface.hpp"

namespace Ariadne {

/************ Number *********************************************************/

template<class X> struct IsNumericType;

class NumberInterface;

template<class P> class Number;
template<class P> struct IsNumericType<Number<P>> : True { };
template<class P> struct IsNumericType<LowerNumber<P>> : True { };
template<class P> struct IsNumericType<UpperNumber<P>> : True { };

template<class R> struct IsConcreteNumericType : IsConvertible<R,Real> { };

struct DispatchException : public std::runtime_error {
    using std::runtime_error::runtime_error;
};


template<class P> Positive<Number<P>> cast_positive(Number<P> y);
template<class P> Positive<UpperNumber<P>> cast_positive(UpperNumber<P> y);


class DeclareNumberOperators {
    friend ApproximateNumber operator+(ApproximateNumber const&, ApproximateNumber const&);
    friend ApproximateNumber operator-(ApproximateNumber const&, ApproximateNumber const&);
    friend ApproximateNumber operator*(ApproximateNumber const&, ApproximateNumber const&);
    friend ApproximateNumber operator/(ApproximateNumber const&, ApproximateNumber const&);

    friend ValidatedLowerNumber operator+(ValidatedLowerNumber const& y1, ValidatedLowerNumber const& y2);
    friend ValidatedLowerNumber operator-(ValidatedLowerNumber const& y1, ValidatedUpperNumber const& y2);
    friend ValidatedUpperNumber operator+(ValidatedUpperNumber const& y1, ValidatedUpperNumber const& y2);
    friend ValidatedUpperNumber operator-(ValidatedUpperNumber const& y1, ValidatedLowerNumber const& y2);

    template<class N, class D, EnableIf<IsGenericNumericType<N>> =dummy, EnableIf<IsSame<D,Dbl>> =dummy> friend auto
    operator+(N const& y1, D const& d2) -> decltype(y1+Number<ApproximateTag>(d2)) { return y1+Number<ApproximateTag>(d2); }
    template<class N, class D, EnableIf<IsGenericNumericType<N>> =dummy, EnableIf<IsSame<D,Dbl>> =dummy> friend auto
    operator-(N const& y1, D const& d2) -> decltype(y1-Number<ApproximateTag>(d2)) { return y1-Number<ApproximateTag>(d2); }
    template<class N, class D, EnableIf<IsGenericNumericType<N>> =dummy, EnableIf<IsSame<D,Dbl>> =dummy> friend auto
    operator*(N const& y1, D const& d2) -> decltype(y1*Number<ApproximateTag>(d2)) { return y1*Number<ApproximateTag>(d2); }
    template<class N, class D, EnableIf<IsGenericNumericType<N>> =dummy, EnableIf<IsSame<D,Dbl>> =dummy> friend auto
    operator/(N const& y1, D const& d2) -> decltype(y1/Number<ApproximateTag>(d2)) { return y1/Number<ApproximateTag>(d2); }
    template<class N, class D, EnableIf<IsGenericNumericType<N>> =dummy, EnableIf<IsSame<D,Dbl>> =dummy> friend auto
    operator+(D const& d1, N const& y2) -> decltype(Number<ApproximateTag>(d1)+y2) { return Number<ApproximateTag>(d1)+y2; }
    template<class N, class D, EnableIf<IsGenericNumericType<N>> =dummy, EnableIf<IsSame<D,Dbl>> =dummy> friend auto
    operator-(D const& d1, N const& y2) -> decltype(Number<ApproximateTag>(d1)-y2) { return Number<ApproximateTag>(d1)-y2; }
    template<class N, class D, EnableIf<IsGenericNumericType<N>> =dummy, EnableIf<IsSame<D,Dbl>> =dummy> friend auto
    operator*(D const& d1, N const& y2) -> decltype(Number<ApproximateTag>(d1)*y2) { return Number<ApproximateTag>(d1)*y2; }
    template<class N, class D, EnableIf<IsGenericNumericType<N>> =dummy, EnableIf<IsSame<D,Dbl>> =dummy> friend auto
    operator/(D const& d1, N const& y2) -> decltype(Number<ApproximateTag>(d1)/y2) { return Number<ApproximateTag>(d1)/y2; }

    template<class R, class P, EnableIf<IsConcreteNumericType<R>> =dummy> friend auto
    operator+(R const& r1, Number<P> const& y2) -> decltype(Number<Paradigm<R>>(r1)+y2) { return Number<Paradigm<R>>(r1)+y2; }
    template<class R, class P, EnableIf<IsConcreteNumericType<R>> =dummy> friend auto
    operator+(Number<P> const& y1, R const& r2) -> decltype(y1+Number<Paradigm<R>>(r2)) { return y1+Number<Paradigm<R>>(r2); }

};

template<class X, class P=Void> struct HasOperatorNumber {
    template<class XX, class PP, class=decltype(declval<XX>().operator Number<PP>())> static True test(int);
    template<class XX, class PP> static False test(...);
    static const bool value = decltype(test<X,P>(1))::value;
};

template<class X> struct HasOperatorNumber<X,Void> {
    template<class XX, class=decltype(declval<XX>().operator Number<Paradigm<XX>>())> static True test(int);
    template<class XX> static False test(...);
    static const bool value = decltype(test<X>(1))::value;
};

//! \ingroup NumericModule
//! \brief Generic numbers with computational paradigm \a P, which may be %EffectiveTag, %ValidatedTag, %UpperTag, %LowerTag or %ApproximateTag.
template<class P> class Number
    : public DeclareNumberOperators
{
    static_assert(IsParadigm<P>::value,"P must be a paradigm");
    static_assert(IsSame<P,ExactTag>::value or IsSame<P,EffectiveTag>::value or IsSame<P,ValidatedTag>::value or IsSame<P,ApproximateTag>::value);
    
    template<class PP> friend class Number;
    template<class PP> friend class LowerNumber;
    template<class PP> friend class UpperNumber;
    
    template<class X> using IsGettableAs = And<IsNumericType<X>,IsWeaker<typename X::Paradigm,P>,Not<IsSame<typename X::Paradigm,ExactTag>>>;
  private: public:
    Handle<NumberInterface> _handle;
  public:
    explicit Number(NumberInterface* p) : _handle(p) { }
  private: public:
    explicit Number(Handle<NumberInterface> h) : _handle(h) { }
    Handle<NumberInterface> handle() const { return this->_handle; }
    NumberInterface const& _ref() const { return this->_handle.reference(); }
  public:
    typedef P Paradigm;
    typedef Number<P> NumericType;

    Number() : Number(Integer(0)) { }

    // Construct from a Number of a stronger paradigm
    template<class SP, EnableIf<IsStronger<SP,P>> = dummy> Number(const Number<SP>& y) : Number<P>(y.handle()) { }

    // Construct from a builtin integer
    template<class N, EnableIf<IsBuiltinIntegral<N>> =dummy> Number(const N& n) : Number<P>(Integer(n)) { }
    // Construct from a builtin floating-point number
    template<class X, EnableIf<And<IsSame<P,ApproximateTag>,IsBuiltinFloatingPoint<X>>> =dummy> 
        Number(const X& x) : Number<P>(FloatType<P,DoublePrecision>(x)) { }

    // Construct from a type which is convertible to Real.
    template<class X, EnableIf<IsWeaker<P,ParadigmTag<X>>> =dummy,
                               EnableIf<IsConvertible<X,Real>> =dummy>
        Number<P>(X const & x) : Number<P>(x.operator Number<ParadigmTag<X>>()) { }

    // Construct from a type which is convertible to another Number type.
    // TODO: Decide conversion properties from concrete type to Number<P>
    template<class X, EnableIf<IsWeaker<P,ParadigmTag<X>>> =dummy,
                      DisableIf<IsConvertible<X,Real>> =dummy,
                      EnableIf<IsConvertible<X,Number<ParadigmTag<X>>>> =dummy>
        explicit Number<P>(X const & x) : Number<P>(x.operator Number<ParadigmTag<X>>()) { }

    //! \brief Convert to an UpperNumber, losing lower-bound information.
    UpperNumber<P> upper() const;

    //! \brief Get the value of the number as a double-precision floating-point type
    template<class WP, EnableIf<IsWeaker<WP,P>> =dummy>
    FloatType<WP,DoublePrecision> get(WP par, DoublePrecision const& prec) const { return this->_ref()._get(WP(),prec); }
    //! \brief Get the value of the number as a multiple-precision floating-point type
    template<class WP, EnableIf<IsWeaker<WP,P>> =dummy>
    FloatType<WP,MultiplePrecision> get(WP par, MultiplePrecision const& prec) const { return this->_ref()._get(WP(),prec); }
    //! \brief Get the value of the number as a double-precision floating-point type
    template<class PR, EnableIf<IsSame<PR,DoublePrecision>> =dummy>
    FloatType<P,PR> get(PR pr) const { return this->_ref()._get(P(),pr); }
    //! \brief Get the value of the number as a multiple-precision floating-point type
    template<class PR, EnableIf<IsSame<PR,MultiplePrecision>> =dummy>
    FloatType<P,PR> get(PR pr) const { return this->_ref()._get(P(),pr); }

    //! \brief Get the value of the number as a multiple-prision floating-point ball.
    template<class WP, EnableIf<IsWeaker<WP,P>> =dummy, EnableIf<IsSame<WP,MetricTag>> =dummy>
    FloatType<WP,DoublePrecision> get(WP par, DoublePrecision const& pr, DoublePrecision const& pre) const { return this->_ref()._get(WP(),pr,pre); }
    //! \brief Get the value of the number as a multiple-prision floating-point ball, with double-prision error.
    template<class WP, EnableIf<IsWeaker<WP,P>> =dummy, EnableIf<IsSame<WP,MetricTag>> =dummy>
    FloatType<WP,MultiplePrecision,DoublePrecision> get(WP par, MultiplePrecision const& pr, DoublePrecision const& pre) const { return this->_ref()._get(WP(),pr,pre); }
    //! \brief Get the value of the number as a multiple-prision floating-point ball.
    template<class WP, EnableIf<IsWeaker<WP,P>> =dummy, EnableIf<IsSame<WP,MetricTag>> =dummy>
    FloatType<WP,MultiplePrecision,MultiplePrecision> get(WP par, MultiplePrecision const& pr, MultiplePrecision const& pre) const { return this->_ref()._get(WP(),pr,pre); }

    template<class PR> FloatBounds<PR> get(PR const& pr) const { return this->_ref()._get(P(),pr); }
    template<class PR, class PRE> FloatBall<PR> get(PR const& pr, PRE const& pre) const { return this->_ref()._get(P(),pr,pre); }

/*
    //! \brief Get the value of the number as a double-prision floating-point type  DEPRECATED
    FloatBounds<DoublePrecision> get() const { return this->_ref()._get(P()); }
    //! \brief Get the value of the number as a double-precision floating-point type
    FloatBounds<DoublePrecision> get(DoublePrecision const& prec) const { return this->_ref()._get(P(),prec); }
    //! \brief Get the value of the number as a multiple-precision floating-point type
    FloatBounds<MultiplePrecision> get(MultiplePrecision const& prec) const { return this->_ref()._get(P(),prec); }

    FloatBall<DoublePrecision,DoublePrecision> get(DoublePrecision const& pr, DoublePrecision const& pre) const { 
        return this->_ref()._get(MetricTag(),pr,pre); }
    FloatBall<DoublePrecision,DoublePrecision> get(MultiplePrecision const& pr, DoublePrecision const& pre) const { 
        return this->_ref()._get(MetricTag(),pr,pre); }
    FloatBall<MultiplePrecision,MultiplePrecision> get(MultiplePrecision const& pr, MultiplePrecision const& pre) const { 
        return this->_ref()._get(MetricTag(),pr,pre); }
*/

    template<class X> X extract() const;

    friend Number<P> operator+(Number<P> const& y) { return pos(y); }
    friend Number<P> operator-(Number<P> const& y) { return neg(y); }
    friend Number<P> operator+(Number<P> const& y1, Number<P> const& y2) { return add(y1,y2); }
    friend Number<P> operator-(Number<P> const& y1, Number<P> const& y2) { return sub(y1,y2); }
    friend Number<P> operator*(Number<P> const& y1, Number<P> const& y2) { return mul(y1,y2); }
    friend Number<P> operator/(Number<P> const& y1, Number<P> const& y2) { return div(y1,y2); }
    friend Number<P>& operator+=(Number<P>& y1, Number<P> const& y2) { return y1=y1+y2; }
    friend Number<P>& operator-=(Number<P>& y1, Number<P> const& y2) { return y1=y1-y2; }
    friend Number<P>& operator*=(Number<P>& y1, Number<P> const& y2) { return y1=y1*y2; }
    friend Number<P>& operator/=(Number<P>& y1, Number<P> const& y2) { return y1=y1/y2; }

    friend Number<P> pos(Number<P> const& y) { return Number<P>(y._ref()._apply(Pos())); }
    friend Number<P> neg(Number<P> const& y) { return Number<P>(y._ref()._apply(Neg())); }
    friend Positive<Number<P>> sqr(Number<P> const& y) { return cast_positive(Number<P>(y._ref()._apply(Sqr()))); }
    friend Number<P> rec(Number<P> const& y) { return Number<P>(y._ref()._apply(Rec())); }
    friend Number<P> add(Number<P> const& y1, Number<P> const& y2) { return Number<P>(y1._ref()._apply(Add(),&y2._ref())); }
    friend Number<P> sub(Number<P> const& y1, Number<P> const& y2) { return Number<P>(y1._ref()._apply(Sub(),&y2._ref())); }
    friend Number<P> mul(Number<P> const& y1, Number<P> const& y2) { return Number<P>(y1._ref()._apply(Mul(),&y2._ref())); }
    friend Number<P> div(Number<P> const& y1, Number<P> const& y2) { return Number<P>(y1._ref()._apply(Div(),&y2._ref())); }

    friend Number<P> sqrt(Number<P> const& y) { return Number<P>(y._ref()._apply(Sqrt())); }
    friend Number<P> exp(Number<P> const& y) { return Number<P>(y._ref()._apply(Exp())); }
    friend Number<P> log(Number<P> const& y) { return Number<P>(y._ref()._apply(Log())); }
    friend Number<P> sin(Number<P> const& y) { return Number<P>(y._ref()._apply(Sin())); }
    friend Number<P> cos(Number<P> const& y) { return Number<P>(y._ref()._apply(Cos())); }
    friend Number<P> tan(Number<P> const& y) { return Number<P>(y._ref()._apply(Tan())); }
    friend Number<P> atan(Number<P> const& y) { return Number<P>(y._ref()._apply(Atan())); }

    friend Number<P> pow(Number<P> const& y, Nat m) { return Number<P>(y._ref()._apply(PowOf<Int>(m))); }
    friend Number<P> pow(Number<P> const& y, Int n) { return Number<P>(y._ref()._apply(PowOf<Int>(n))); }

    friend Positive<Number<P>> abs(Number<P> const& y) { return cast_positive(Number<P>(y._ref()._apply(Abs()))); }
    friend Number<P> max(Number<P> const& y1, Number<P> const& y2) { return Number<P>(y1._ref()._apply(Max(),&y2._ref())); }
    friend Number<P> min(Number<P> const& y1, Number<P> const& y2) { return Number<P>(y1._ref()._apply(Min(),&y2._ref())); }

    friend EqualityLogicalType<P> operator==(Number<P> const& y1, Number<P> const& y2) {
        return EqualityLogicalType<P>(y1._ref()._apply(Equal(),&y2._ref())); }
    friend InequalityLogicalType<P> operator!=(Number<P> const& y1, Number<P> const& y2) { return !(y1==y2); }
    friend LogicalType<P> operator< (Number<P> const& y1, Number<P> const& y2) {
        return LogicalType<P>(y1._ref()._apply(Less(),&y2._ref())); }
    friend LogicalType<P> operator> (Number<P> const& y1, Number<P> const& y2) { return (y2<y1); }
    friend LogicalType<P> operator<=(Number<P> const& y1, Number<P> const& y2) { return !(y1>y2); }
    friend LogicalType<P> operator>=(Number<P> const& y1, Number<P> const& y2) { return !(y1<y2); }

    String class_name() const { return this->_ref()._class_name(); }

    friend OutputStream& operator<<(OutputStream& os, Number<P> const& y) { return y._ref()._write(os); }
};


//! \ingroup NumericModule
//! \brief Generic lower (real) numbers with computational paradigm \a P, which may be %EffectiveTag or %ValidatedTag.
template<class P> class LowerNumber
{
    static_assert(IsSame<P,EffectiveTag>::value or IsSame<P,ValidatedTag>::value,"P must be a paradigm");
    friend class UpperNumber<P>;
  private: public:
    Handle<NumberInterface> _handle;
  private:
    explicit LowerNumber(Handle<NumberInterface> h) : _handle(h) { }
    Handle<NumberInterface> handle() const { return this->_handle; }
    NumberInterface const& _ref() const { return this->_handle.reference(); }
  public:
    typedef P Paradigm;
    typedef LowerNumber<P> NumericType;

    LowerNumber() : LowerNumber(Integer(0)) { }

    //! \brief Construct from a LowerNumber of a stronger paradigm
    template<class SP, EnableIf<IsStronger<SP,P>> = dummy> LowerNumber(const LowerNumber<SP>& y) : LowerNumber<P>(y.handle()) { }
    //! \brief Construct from a type convertible to a Number.
    template<class X, EnableIf<IsConvertible<X,Number<P>>> = dummy> LowerNumber(const X& x) : LowerNumber<P>(Number<P>(x).handle()) { }
    
    template<class PR> FloatLowerBound<PR> get(LowerTag, PR pr) const { return this->_ref()._get(LowerTag(),pr); }
    template<class PR> FloatLowerBound<PR> get(PR pr) const { return this->_ref()._get(LowerTag(),pr); }

    template<class X> X extract() const;

    friend LowerNumber<P> operator+(LowerNumber<P> const& y) { return pos(y); }
    friend UpperNumber<P> operator-(LowerNumber<P> const& y) { return neg(y); }
    friend LowerNumber<P> operator+(LowerNumber<P> const& y1, LowerNumber<P> const& y2) { return add(y1,y2); }
    friend LowerNumber<P> operator-(LowerNumber<P> const& y1, UpperNumber<P> const& y2) { return sub(y1,y2); }
    friend LowerNumber<P>& operator+=(LowerNumber<P>& y1, LowerNumber<P> const& y2) { return y1=y1+y2; }
    friend LowerNumber<P>& operator-=(LowerNumber<P>& y1, UpperNumber<P> const& y2) { return y1=y1-y2; }
    
    friend LowerNumber<P> pos(LowerNumber<P> const& y) { return LowerNumber<P>(y._ref()._apply(Pos())); }
    friend UpperNumber<P> neg(LowerNumber<P> const& y) { assert(false); } //return UpperNumber<P>(y._ref()._apply(Neg())); }
    friend LowerNumber<P> add(LowerNumber<P> const& y1, LowerNumber<P> const& y2) { return LowerNumber<P>(y1._ref()._apply(Add(),&y2._ref())); }
    friend LowerNumber<P> sub(LowerNumber<P> const& y1, UpperNumber<P> const& y2) { assert(false); } //return LowerNumber<P>(y1._ref()._apply(Sub(),&y2._ref())); }
    
    friend LowerNumber<P> sqrt(LowerNumber<P> const& y) { return LowerNumber<P>(y._ref()._apply(Sqrt())); }
    friend LowerNumber<P> exp(LowerNumber<P> const& y) { return LowerNumber<P>(y._ref()._apply(Exp())); }
    friend LowerNumber<P> log(LowerNumber<P> const& y) { return LowerNumber<P>(y._ref()._apply(Log())); }
    friend LowerNumber<P> atan(LowerNumber<P> const& y) { return LowerNumber<P>(y._ref()._apply(Atan())); }

    friend LowerNumber<P> max(LowerNumber<P> const& y1, LowerNumber<P> const& y2) { return LowerNumber<P>(y1._ref()._apply(Max(),&y2._ref())); }
    friend LowerNumber<P> min(LowerNumber<P> const& y1, LowerNumber<P> const& y2) { return LowerNumber<P>(y1._ref()._apply(Min(),&y2._ref())); }

    friend UpperLogicalType<P> operator==(LowerNumber<P> const& y1, UpperNumber<P> const& y2) {
        return UpperLogicalType<P>(y1._ref()._equals(y2._ref())); }
    friend LowerLogicalType<P> operator!=(LowerNumber<P> const& y1, UpperNumber<P> const& y2) {
        return not (y1 == y2); }
    friend UpperLogicalType<P> operator< (LowerNumber<P> const& y1, UpperNumber<P> const& y2) {
        return UpperLogicalType<P>(y1._ref()._apply(Less(),y2._ref())); }
    friend LowerLogicalType<P> operator> (LowerNumber<P> const& y1, UpperNumber<P> const& y2) { 
        return y2 <  y1; }
 
    String class_name() const { return this->_ref()._class_name(); }

    friend OutputStream& operator<<(OutputStream& os, LowerNumber<P> const& y) { return y._ref()._write(os); }
};


//! \ingroup NumericModule
//! \brief Generic upper (real) numbers with computational paradigm \a P, which may be %EffectiveTag or %ValidatedTag.
template<class P> class UpperNumber
{
    static_assert(IsSame<P,EffectiveTag>::value or IsSame<P,ValidatedTag>::value,"P must be a paradigm");
    friend class LowerNumber<P>;
  private: public:
    Handle<UpperNumberInterface> _handle;
  private: public:
    explicit UpperNumber(Handle<UpperNumberInterface> h) : _handle(h) { }
    Handle<UpperNumberInterface> handle() const { return this->_handle; }
    UpperNumberInterface const& _ref() const { return this->_handle.reference(); }
  public:
    typedef P Paradigm;
    typedef UpperNumber<P> NumericType;

    UpperNumber() : UpperNumber(Integer(0)) { }

    //! \brief Construct from a UpperNumber of a stronger paradigm
    template<class SP, EnableIf<IsStronger<SP,P>> = dummy> 
    UpperNumber(const UpperNumber<SP>& y) : UpperNumber<P>(y.handle()) { }
    //! \brief Construct from a Number of a stronger paradigm
    template<class SP, EnableIf<IsStronger<SP,P>> = dummy> 
    UpperNumber(const Number<SP>& y) : UpperNumber(nullptr) { assert(false); }
    //! \brief Construct from a type convertible to a Number.
    template<class X, EnableIf<IsConvertible<X,Number<P>>> = dummy, DisableIf<IsConvertible<X,UpperNumber<P>>> =dummy> 
    UpperNumber(const X& x) : UpperNumber<P>(Number<P>(x)) { }

    // Construct from a type which is convertible to another Number type.
    // TODO: Decide conversion properties from concrete type to Number<P>
    template<class X, EnableIf<IsConvertible<X,UpperNumber<P>>> =dummy>
        explicit UpperNumber<P>(X const & x) : UpperNumber<P>(x.operator UpperNumber<P>()) { }

    template<class PR> FloatUpperBound<PR> get(UpperTag, PR pr) const { return this->_ref()._get(pr); }
    template<class PR> FloatUpperBound<PR> get(PR pr) const { return this->_ref()._get(pr); }

    template<class X> X extract() const { return dynamic_cast<Wrapper<X,NumberInterface> const&>(this->_ref); }

    friend UpperNumber<P> operator+(UpperNumber<P> const& y) { return pos(y); }
    friend LowerNumber<P> operator-(UpperNumber<P> const& y) { return neg(y); }
    friend UpperNumber<P> operator+(UpperNumber<P> const& y1, UpperNumber<P> const& y2) { return add(y1,y2); }
    friend UpperNumber<P> operator-(UpperNumber<P> const& y1, LowerNumber<P> const& y2) { return sub(y1,y2); }
    friend UpperNumber<P>& operator+=(UpperNumber<P>& y1, UpperNumber<P> const& y2) { return y1=y1+y2; }
    friend UpperNumber<P>& operator-=(UpperNumber<P>& y1, LowerNumber<P> const& y2) { return y1=y1-y2; }
    
    friend UpperNumber<P> pos(UpperNumber<P> const& y) { return UpperNumber<P>(y._ref()._apply(Pos())); }
    friend LowerNumber<P> neg(UpperNumber<P> const& y) { return LowerNumber<P>(y._ref()._apply(Neg())); }
    friend UpperNumber<P> add(UpperNumber<P> const& y1, UpperNumber<P> const& y2) { return UpperNumber<P>(y1._ref()._apply(Add(),&y2._ref())); }
//    friend UpperNumber<P> sub(UpperNumber<P> const& y1, LowerNumber<P> const& y2) { return UpperNumber<P>(y1._ref()._apply(Sub(),&y2._ref())); }
    friend UpperNumber<P> sub(UpperNumber<P> const& y1, LowerNumber<P> const& y2) { return add(y1,neg(y2)); }
    
    friend UpperNumber<P> sqrt(UpperNumber<P> const& y) { return UpperNumber<P>(y._ref()._apply(Sqrt())); }
    friend UpperNumber<P> exp(UpperNumber<P> const& y) { return UpperNumber<P>(y._ref()._apply(Exp())); }
    friend UpperNumber<P> log(UpperNumber<P> const& y) { return UpperNumber<P>(y._ref()._apply(Log())); }
    friend UpperNumber<P> atan(UpperNumber<P> const& y) { return UpperNumber<P>(y._ref()._apply(Atan())); }

    friend UpperNumber<P> max(UpperNumber<P> const& y1, UpperNumber<P> const& y2) { return UpperNumber<P>(y1._ref()._apply(Max(),&y2._ref())); }
    friend UpperNumber<P> min(UpperNumber<P> const& y1, UpperNumber<P> const& y2) { return UpperNumber<P>(y1._ref()._apply(Min(),&y2._ref())); }

    friend UpperLogicalType<P> operator==(UpperNumber<P> const& y1, LowerNumber<P> const& y2) {
        return UpperLogicalType<P>(y1._ref()._apply(Equal(),y2._ref())); }
    friend LowerLogicalType<P> operator!=(UpperNumber<P> const& y1, LowerNumber<P> const& y2) {
        return not (y1 == y2); }
    friend LowerLogicalType<P> operator< (UpperNumber<P> const& y1, LowerNumber<P> const& y2) {
        assert(false); } //return LowerLogicalType<P>(y1._ref()._apply(Less(),y2._ref())); }
    friend UpperLogicalType<P> operator> (UpperNumber<P> const& y1, LowerNumber<P> const& y2) { 
        return y2 <  y1; }
        
    String class_name() const { return this->_ref()._class_name(); }

    friend OutputStream& operator<<(OutputStream& os, UpperNumber<P> const& y) { return y._ref()._write(os); }
};



template<class P> class Positive<Number<P>> : public Number<P> {
    friend Number<P> const& unsign(Positive<Number<P>> const& y) { return y; }
  public:
    Positive<Number<P>>() : Number<P>() { }
    explicit Positive<Number<P>>(Number<P> const& y)
        : Number<P>(y) { }
    template<class N, EnableIf<IsConvertible<N,ExactNumber>> =dummy> 
        Positive<Number<P>>(const Positive<N>& n) : Number<P>(ExactNumber(static_cast<N const&>(n))) { }
    template<class N, EnableIf<IsConstructible<ExactNumber,N>> =dummy> 
        explicit Positive<Number<P>>(const N& n) : Number<P>(ExactNumber(n)) { }
    explicit operator Number<P> () const { return *this; }

    friend Positive<Number<P>> operator+(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) { 
        return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> operator*(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) { 
        return cast_positive(mul(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> add(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) { 
        return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> mul(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) { 
        return cast_positive(mul(unsign(y1),unsign(y2))); }
    friend Positive<Number<P>> max(Positive<Number<P>> const& y1, Positive<Number<P>> const& y2) { 
        return cast_positive(max(unsign(y1),unsign(y2))); }
};

template<class P> class Positive<UpperNumber<P>> : public UpperNumber<P> {
    friend UpperNumber<P> const& unsign(Positive<UpperNumber<P>> const& y) { return y; }
  public:
    Positive<UpperNumber<P>>() : UpperNumber<P>() { }
    explicit Positive<UpperNumber<P>>(UpperNumber<P> const& y) : UpperNumber<P>(y) { }
    template<class N, EnableIf<IsConstructible<ExactNumber,N>> =dummy> 
        Positive<UpperNumber<P>>(const Positive<N>& n) : UpperNumber<P>(ExactNumber(static_cast<N const&>(n))) { }
    template<class N, EnableIf<IsConstructible<ExactNumber,N>> =dummy> 
        Positive<UpperNumber<P>>(const N& n) : UpperNumber<P>(ExactNumber(n)) { }
    explicit operator UpperNumber<P> () const { return *this; }

    friend UpperNumber<P> mul(UpperNumber<P> const& y1, UpperNumber<P> const& y2);
        
    friend Positive<UpperNumber<P>> operator+(Positive<UpperNumber<P>> const& y1, Positive<UpperNumber<P>> const& y2) { 
        return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<UpperNumber<P>> operator*(Positive<UpperNumber<P>> const& y1, Positive<UpperNumber<P>> const& y2) {
        return cast_positive(mul(unsign(y1),unsign(y2))); }
    friend Positive<UpperNumber<P>> add(Positive<UpperNumber<P>> const& y1, Positive<UpperNumber<P>> const& y2) {
        return cast_positive(add(unsign(y1),unsign(y2))); }
    friend Positive<UpperNumber<P>> mul(Positive<UpperNumber<P>> const& y1, Positive<UpperNumber<P>> const& y2) { 
        return cast_positive(Number<P>(y1._ref()._apply(Mul(),y2._ref()))); }
    friend Positive<UpperNumber<P>> max(Positive<UpperNumber<P>> const& y1, Positive<UpperNumber<P>> const& y2) { 
        return cast_positive(max(unsign(y1),unsign(y2))); }
};

template<class X> decltype(auto) cast_generic(X const& x) { return x.generic(); }

template<class P> Positive<Number<P>> cast_positive(Number<P> y) { return Positive<Number<P>>(y); }
template<class P> Positive<UpperNumber<P>> cast_positive(UpperNumber<P> y) { return Positive<UpperNumber<P>>(y); }

template<class P> Positive<ExactNumber> cast_exact(Positive<Number<P>> const& y) { return cast_positive(cast_exact(unsign(y))); }
template<class P> Positive<ExactNumber> cast_exact(Positive<UpperNumber<P>> const& y) { return cast_positive(cast_exact(unsign(y))); }


} // namespace Ariadne

#endif
