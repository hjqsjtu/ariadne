/***************************************************************************
 *            numeric/floatdp.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file numeric/floatdp.hpp
 *  \brief RawTag floating-point number class based on double-precision floats.
 */

#ifndef ARIADNE_FLOAT64_HPP
#define ARIADNE_FLOAT64_HPP

#include <iosfwd> // For std::floor std::ceil etc
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "../utility/declarations.hpp"
#include "../numeric/operators.hpp"
#include "../numeric/rounding.hpp"
#include "../numeric/sign.hpp"
#include "../numeric/number.decl.hpp"
#include "../numeric/float.decl.hpp"

namespace Ariadne {

class FloatDP;
typedef FloatDP RawFloatDP;

class Rational;
enum class Comparison : char;

//! \ingroup NumericModule
//! \brief The precision of a FloatDP object. Since this is fixed, the class is only a tag; all objects are equal.
//! \relates FloatDP
class DoublePrecision;
using DP = DoublePrecision;

// Correctly rounded functions
double nul_rnd(double x);
double pos_rnd(double x);
double neg_rnd(double x);
double hlf_rnd(double x);
double sqr_rnd(double x);
double rec_rnd(double x);
double add_rnd(double x1, double x2);
double sub_rnd(double x1, double x2);
double mul_rnd(double x1, double x2);
double div_rnd(double x1, double x2);
double fma_rnd(double x1, double x2, double x3);
double pow_rnd(double x, Nat n);
double pow_rnd(double x, int n);
double sqrt_rnd(double x);
double exp_rnd(double x);
double log_rnd(double x);
double sin_rnd(double x);
double cos_rnd(double x);
double tan_rnd(double x);
double asin_rnd(double x);
double acos_rnd(double x);
double atan_rnd(double x);
double neg_rec_rnd(double x);
double atan_rnd_series(double x);
double pi_rnd();

double texp(double x);

double pi_opp();
double add_opp(double x, double y);
double sub_opp(double x, double y);
double mul_opp(double x, double y);
double div_opp(double x, double y);
double neg_rec_opp(double x);

template<class FLT> class Rounded;

using RoundedFloatDP = Rounded<FloatDP>;

template<> class Rounded<FloatDP>
{
    volatile double dbl;
  public:
    typedef RawTag Paradigm;
    typedef FloatDP FloatType;
    typedef Rounded<FloatDP> NumericType;
    typedef DoublePrecision PrecisionType;
    typedef BuiltinRoundingModeType RoundingModeType;
  public:
    static RoundingModeType get_rounding_mode() { return FloatDP::get_rounding_mode(); }
    static Void set_rounding_mode(RoundingModeType rnd) { FloatDP::set_rounding_mode(rnd); }
    static Void set_rounding_downward() { FloatDP::set_rounding_downward(); }
    static Void set_rounding_upward() { FloatDP::set_rounding_upward(); }
    static Void set_rounding_to_nearest() { FloatDP::set_rounding_to_nearest(); }
    static Void set_rounding_toward_zero() { FloatDP::set_rounding_toward_zero(); }
  public:
    Rounded() : dbl(0.0) { }
    Rounded(double d) : dbl(d) { }
    double data() const { return dbl; }
    PrecisionType precision() const { return DoublePrecision(); }

    Rounded(Dyadic const& w, PrecisionType pr) : Rounded(FloatDP(w,pr)) { }
    template<class Y, EnableIf<IsConstructible<FloatType,Y,RoundingModeType,PrecisionType>> =dummy>
        Rounded(Y const& y, PrecisionType pr) : Rounded(FloatDP(y,FloatDP::get_rounding_mode(),pr)) { }

    Rounded<FloatDP>& operator=(double d) { this->dbl=d; return *this; }

    inline explicit Rounded(FloatType x, PrecisionType) : dbl(x.dbl) { }
    inline explicit Rounded(Rounded<FloatType> x, PrecisionType) : dbl(x.dbl) { }
    inline explicit Rounded(PrecisionType pr) : Rounded(FloatType(pr)) { }
    inline explicit Rounded(FloatType x) : dbl(x.dbl) { }
    inline explicit operator FloatType() const { return FloatType(this->dbl); }
    inline FloatType raw() const { return FloatDP(this->dbl); }

    explicit Rounded(double d, PrecisionType pr) : dbl(d) { }
    Rounded(Value<FloatDP> const& x);
    Rounded(Approximation<FloatDP> const& x);
    operator Approximation<FloatDP> () const;

    // Integer conversions
    friend Rounded<FloatDP> floor(Rounded<FloatDP> x) { return Rounded<FloatDP>(std::floor(x.dbl)); }
    friend Rounded<FloatDP> ceil(Rounded<FloatDP> x) { return Rounded<FloatDP>(std::ceil(x.dbl)); }
    friend Rounded<FloatDP> round(Rounded<FloatDP> x) { return Rounded<FloatDP>(std::round(x.dbl)); }

    // Correctly rounded arithmetic
    friend Rounded<FloatDP> nul(Rounded<FloatDP> x) { return Rounded<FloatDP>(nul_rnd(x.dbl)); }
    friend Rounded<FloatDP> pos(Rounded<FloatDP> x) { return Rounded<FloatDP>(pos_rnd(x.dbl)); }
    friend Rounded<FloatDP> neg(Rounded<FloatDP> x) { return Rounded<FloatDP>(neg_rnd(x.dbl)); }
    friend Rounded<FloatDP> sqr(Rounded<FloatDP> x) { return Rounded<FloatDP>(sqr_rnd(x.dbl)); }
    friend Rounded<FloatDP> hlf(Rounded<FloatDP> x) { return Rounded<FloatDP>(hlf_rnd(x.dbl)); }
    friend Rounded<FloatDP> rec(Rounded<FloatDP> x) { return Rounded<FloatDP>(rec_rnd(x.dbl)); }
    friend Rounded<FloatDP> add(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(add_rnd(x1.dbl,x2.dbl)); }
    friend Rounded<FloatDP> sub(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(sub_rnd(x1.dbl,x2.dbl)); }
    friend Rounded<FloatDP> mul(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(mul_rnd(x1.dbl,x2.dbl)); }
    friend Rounded<FloatDP> div(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(div_rnd(x1.dbl,x2.dbl)); }
    friend Rounded<FloatDP> fma(Rounded<FloatDP> x1, Rounded<FloatDP> x2, Rounded<FloatDP> x3) { return Rounded<FloatDP>(fma_rnd(x1.dbl,x2.dbl,x3.dbl)); }
    friend Rounded<FloatDP> pow(Rounded<FloatDP> x, Int n) { return Rounded<FloatDP>(pow_rnd(x.dbl,n)); }
    friend Rounded<FloatDP> sqrt(Rounded<FloatDP> x) { return Rounded<FloatDP>(sqrt_rnd(x.dbl)); }
    friend Rounded<FloatDP> exp(Rounded<FloatDP> x) { return Rounded<FloatDP>(exp_rnd(x.dbl)); }
    friend Rounded<FloatDP> log(Rounded<FloatDP> x) { return Rounded<FloatDP>(log_rnd(x.dbl)); }
    friend Rounded<FloatDP> sin(Rounded<FloatDP> x) { return Rounded<FloatDP>(sin_rnd(x.dbl)); }
    friend Rounded<FloatDP> cos(Rounded<FloatDP> x) { return Rounded<FloatDP>(cos_rnd(x.dbl)); }
    friend Rounded<FloatDP> tan(Rounded<FloatDP> x) { return Rounded<FloatDP>(tan_rnd(x.dbl)); }
    friend Rounded<FloatDP> asin(Rounded<FloatDP> x) { return Rounded<FloatDP>(asin_rnd(x.dbl)); }
    friend Rounded<FloatDP> acos(Rounded<FloatDP> x) { return Rounded<FloatDP>(acos_rnd(x.dbl)); }
    friend Rounded<FloatDP> atan(Rounded<FloatDP> x) { return Rounded<FloatDP>(atan_rnd(x.dbl)); }
    static Rounded<FloatDP> pi(PrecisionType pr) { return Rounded<FloatDP>(pi_rnd()); }

    friend Rounded<FloatDP> abs(Rounded<FloatDP> x) { return Rounded<FloatDP>(std::abs(x.dbl)); }
    friend Rounded<FloatDP> max(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(std::max(x1.dbl,x2.dbl)); }
    friend Rounded<FloatDP> min(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(std::min(x1.dbl,x2.dbl)); }
    friend Rounded<FloatDP> mag(Rounded<FloatDP> x) { return Rounded<FloatDP>(std::abs(x.dbl)); }
    friend Rounded<FloatDP> mig(Rounded<FloatDP> x) { return Rounded<FloatDP>(std::abs(x.dbl)); }

    friend Rounded<FloatDP> med(Rounded<FloatDP> x1, Rounded<FloatDP> x2) {
        return Rounded<FloatDP>(hlf_rnd(add_rnd(x1.dbl,x2.dbl))); }
    friend Rounded<FloatDP> rad(Rounded<FloatDP> x1, Rounded<FloatDP> x2) {
        return Rounded<FloatDP>(hlf_rnd(sub_rnd(x2.dbl,x1.dbl))); }

    friend Rounded<FloatDP> operator+(Rounded<FloatDP> x) { return Rounded<FloatDP>(pos_rnd(x.dbl)); }
    friend Rounded<FloatDP> operator-(Rounded<FloatDP> x) { return Rounded<FloatDP>(neg_rnd(x.dbl)); }
    friend Rounded<FloatDP> operator+(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(x1.dbl+x2.dbl); }
    friend Rounded<FloatDP> operator-(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(x1.dbl-x2.dbl); }
    friend Rounded<FloatDP> operator*(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(x1.dbl*x2.dbl); }
    friend Rounded<FloatDP> operator/(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return Rounded<FloatDP>(x1.dbl/x2.dbl); }
    friend Rounded<FloatDP>& operator+=(Rounded<FloatDP>& x1, Rounded<FloatDP> x2) { x1.dbl=x1.dbl+x2.dbl; return x1; }
    friend Rounded<FloatDP>& operator-=(Rounded<FloatDP>& x1, Rounded<FloatDP> x2) { x1.dbl=x1.dbl-x2.dbl; return x1; }
    friend Rounded<FloatDP>& operator*=(Rounded<FloatDP>& x1, Rounded<FloatDP> x2) { x1.dbl=x1.dbl*x2.dbl; return x1; }
    friend Rounded<FloatDP>& operator/=(Rounded<FloatDP>& x1, Rounded<FloatDP> x2) { x1.dbl=x1.dbl/x2.dbl; return x1; }

    friend Comparison cmp(Rounded<FloatDP> x1, Rounded<FloatDP> x2);
    friend Bool operator==(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1.dbl == x2.dbl; }
    friend Bool operator!=(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1.dbl != x2.dbl; }
    friend Bool operator<=(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1.dbl <= x2.dbl; }
    friend Bool operator>=(Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1.dbl >= x2.dbl; }
    friend Bool operator< (Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1.dbl <  x2.dbl; }
    friend Bool operator> (Rounded<FloatDP> x1, Rounded<FloatDP> x2) { return x1.dbl >  x2.dbl; }

    friend OutputStream& operator<<(OutputStream& os, Rounded<FloatDP> const& x) { return os << x.data(); }
};

template<class FLT> class Rounded
{
    FLT _flt;
  public:
    typedef RawTag Paradigm;
    typedef FLT FloatType;
    typedef Rounded<FloatType> NumericType;
    typedef typename FloatType::PrecisionType PrecisionType;
    typedef typename FloatType::RoundingModeType RoundingModeType;
  public:
    static RoundingModeType get_rounding_mode() { return FloatType::get_rounding_mode(); }
    static Void set_rounding_mode(RoundingModeType rnd) { FloatType::set_rounding_mode(rnd); }
    static Void set_rounding_downward() { FloatType::set_rounding_downward(); }
    static Void set_rounding_upward() { FloatType::set_rounding_upward(); }
    static Void set_rounding_to_nearest() { FloatType::set_rounding_to_nearest(); }
    static Void set_rounding_toward_zero() { FloatType::set_rounding_toward_zero(); }
  public:
    Rounded() : _flt(0.0) { }
    Rounded(FloatType flt) : _flt(flt) { }
    double data() const { return _flt; }
    PrecisionType precision() const { return this->_flt.precision(); }

    Rounded(Dyadic const& w, PrecisionType pr) : Rounded(FloatType(w,pr)) { }
    template<class Y, EnableIf<IsConstructible<FloatType,Y,RoundingModeType,PrecisionType>> =dummy>
        Rounded(Y const& y, PrecisionType pr) : Rounded(FloatType(y,FloatType::get_rounding_mode(),pr)) { }

    Rounded<FloatType>& operator=(FLT flt) { this->_flt=d; return *this; }

    inline explicit Rounded(FloatType x, PrecisionType) : _flt(x._flt) { }
    inline explicit Rounded(Rounded<FloatType> x, PrecisionType) : _flt(x._flt) { }
    inline explicit Rounded(PrecisionType pr) : Rounded(FloatType(pr)) { }
    inline explicit Rounded(FloatType x) : _flt(x._flt) { }
    inline explicit operator FloatType() const { return FloatType(this->_flt); }
    inline FloatType raw() const { return FloatType(this->_flt); }

    explicit Rounded(double d, PrecisionType pr) : _flt(d) { }
    Rounded(Value<FloatType> const& x);
    Rounded(Approximation<FloatType> const& x);
    operator Approximation<FloatType> () const;

    // Integer conversions
    friend Rounded<FloatType> floor(Rounded<FloatType> x) { return Rounded<FloatType>(std::floor(x._flt)); }
    friend Rounded<FloatType> ceil(Rounded<FloatType> x) { return Rounded<FloatType>(std::ceil(x._flt)); }
    friend Rounded<FloatType> round(Rounded<FloatType> x) { return Rounded<FloatType>(std::round(x._flt)); }

    // Correctly rounded arithmetic
    friend Rounded<FloatType> nul(Rounded<FloatType> x) { return Rounded<FloatType>(nul_rnd(x._flt)); }
    friend Rounded<FloatType> pos(Rounded<FloatType> x) { return Rounded<FloatType>(pos_rnd(x._flt)); }
    friend Rounded<FloatType> neg(Rounded<FloatType> x) { return Rounded<FloatType>(neg_rnd(x._flt)); }
    friend Rounded<FloatType> sqr(Rounded<FloatType> x) { return Rounded<FloatType>(sqr_rnd(x._flt)); }
    friend Rounded<FloatType> hlf(Rounded<FloatType> x) { return Rounded<FloatType>(hlf_rnd(x._flt)); }
    friend Rounded<FloatType> rec(Rounded<FloatType> x) { return Rounded<FloatType>(rec_rnd(x._flt)); }
    friend Rounded<FloatType> add(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(add_rnd(x1._flt,x2._flt)); }
    friend Rounded<FloatType> sub(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(sub_rnd(x1._flt,x2._flt)); }
    friend Rounded<FloatType> mul(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(mul_rnd(x1._flt,x2._flt)); }
    friend Rounded<FloatType> div(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(div_rnd(x1._flt,x2._flt)); }
    friend Rounded<FloatType> fma(Rounded<FloatType> x1, Rounded<FloatType> x2, Rounded<FloatType> x3) { return Rounded<FloatType>(fma_rnd(x1._flt,x2._flt,x3._flt)); }
    friend Rounded<FloatType> pow(Rounded<FloatType> x, Int n) { return Rounded<FloatType>(pow_rnd(x._flt,n)); }
    friend Rounded<FloatType> sqrt(Rounded<FloatType> x) { return Rounded<FloatType>(sqrt_rnd(x._flt)); }
    friend Rounded<FloatType> exp(Rounded<FloatType> x) { return Rounded<FloatType>(exp_rnd(x._flt)); }
    friend Rounded<FloatType> log(Rounded<FloatType> x) { return Rounded<FloatType>(log_rnd(x._flt)); }
    friend Rounded<FloatType> sin(Rounded<FloatType> x) { return Rounded<FloatType>(sin_rnd(x._flt)); }
    friend Rounded<FloatType> cos(Rounded<FloatType> x) { return Rounded<FloatType>(cos_rnd(x._flt)); }
    friend Rounded<FloatType> tan(Rounded<FloatType> x) { return Rounded<FloatType>(tan_rnd(x._flt)); }
    friend Rounded<FloatType> asin(Rounded<FloatType> x) { return Rounded<FloatType>(asin_rnd(x._flt)); }
    friend Rounded<FloatType> acos(Rounded<FloatType> x) { return Rounded<FloatType>(acos_rnd(x._flt)); }
    friend Rounded<FloatType> atan(Rounded<FloatType> x) { return Rounded<FloatType>(atan_rnd(x._flt)); }
    static Rounded<FloatType> pi(PrecisionType pr) { return Rounded<FloatType>(pi_rnd()); }

    friend Rounded<FloatType> abs(Rounded<FloatType> x) { return Rounded<FloatType>(std::abs(x._flt)); }
    friend Rounded<FloatType> max(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(std::max(x1._flt,x2._flt)); }
    friend Rounded<FloatType> min(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(std::min(x1._flt,x2._flt)); }
    friend Rounded<FloatType> mag(Rounded<FloatType> x) { return Rounded<FloatType>(std::abs(x._flt)); }
    friend Rounded<FloatType> mig(Rounded<FloatType> x) { return Rounded<FloatType>(std::abs(x._flt)); }

    friend Rounded<FloatType> med(Rounded<FloatType> x1, Rounded<FloatType> x2) {
        return Rounded<FloatType>(hlf_rnd(add_rnd(x1._flt,x2._flt))); }
    friend Rounded<FloatType> rad(Rounded<FloatType> x1, Rounded<FloatType> x2) {
        return Rounded<FloatType>(hlf_rnd(sub_rnd(x2._flt,x1._flt))); }

    friend Rounded<FloatType> operator+(Rounded<FloatType> x) { return Rounded<FloatType>(pos_rnd(x._flt)); }
    friend Rounded<FloatType> operator-(Rounded<FloatType> x) { return Rounded<FloatType>(neg_rnd(x._flt)); }
    friend Rounded<FloatType> operator+(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(x1._flt+x2._flt); }
    friend Rounded<FloatType> operator-(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(x1._flt-x2._flt); }
    friend Rounded<FloatType> operator*(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(x1._flt*x2._flt); }
    friend Rounded<FloatType> operator/(Rounded<FloatType> x1, Rounded<FloatType> x2) { return Rounded<FloatType>(x1._flt/x2._flt); }
    friend Rounded<FloatType>& operator+=(Rounded<FloatType>& x1, Rounded<FloatType> x2) { x1._flt=x1._flt+x2._flt; return x1; }
    friend Rounded<FloatType>& operator-=(Rounded<FloatType>& x1, Rounded<FloatType> x2) { x1._flt=x1._flt-x2._flt; return x1; }
    friend Rounded<FloatType>& operator*=(Rounded<FloatType>& x1, Rounded<FloatType> x2) { x1._flt=x1._flt*x2._flt; return x1; }
    friend Rounded<FloatType>& operator/=(Rounded<FloatType>& x1, Rounded<FloatType> x2) { x1._flt=x1._flt/x2._flt; return x1; }

    friend Comparison cmp(Rounded<FloatType> x1, Rounded<FloatType> x2);
    friend Bool operator==(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt == x2._flt; }
    friend Bool operator!=(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt != x2._flt; }
    friend Bool operator<=(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt <= x2._flt; }
    friend Bool operator>=(Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt >= x2._flt; }
    friend Bool operator< (Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt <  x2._flt; }
    friend Bool operator> (Rounded<FloatType> x1, Rounded<FloatType> x2) { return x1._flt >  x2._flt; }

    friend OutputStream& operator<<(OutputStream& os, Rounded<FloatType> const& x) { return os << x.data(); }
};




} // namespace Ariadne

#endif
