/***************************************************************************
 *            test_rounding_float.cpp
 *
 *  Copyright  2006-20  Alberto Casagrande, Pieter Collins
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

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>

#include "config.hpp"
#include "numeric/rounding.hpp"
#include "numeric/float.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/numeric.hpp"
#include "numeric/dyadic.hpp"

#include "../test.hpp"

namespace Ariadne {
/*
template<> String class_name<DoublePrecision>() { return "DoublePrecision"; }
template<> String class_name<MultiplePrecision>() { return "MultiplePrecision"; }
FloatDP operator+(FloatDP x1, FloatDP x2);
FloatDP operator-(FloatDP x1, FloatDP x2);
FloatDP operator*(FloatDP x1, FloatDP x2);
FloatDP operator/(FloatDP x1, FloatDP x2);
FloatDP& operator+=(FloatDP& x1, FloatDP x2);
FloatDP& operator-=(FloatDP& x1, FloatDP x2);
FloatDP& operator*=(FloatDP& x1, FloatDP x2);
FloatDP& operator/=(FloatDP& x1, FloatDP x2);
FloatDP sqr(FloatDP x) { return sqr_rnd(x.dbl); }
FloatDP rec(FloatDP x) { return rec_rnd(x.dbl); }
FloatDP add(FloatDP x1, FloatDP x2) { return add_rnd(x1.dbl,x2.dbl); }
FloatDP sub(FloatDP x1, FloatDP x2) { return sub_rnd(x1.dbl,x2.dbl); }
FloatDP mul(FloatDP x1, FloatDP x2) { return mul_rnd(x1.dbl,x2.dbl); }
FloatDP div(FloatDP x1, FloatDP x2) { return div_rnd(x1.dbl,x2.dbl); }
FloatDP fma(FloatDP x1, FloatDP x2, FloatDP x3);
FloatDP pow(FloatDP x, Int n) { return pow_rnd(x.dbl,n); }
FloatDP sqrt(FloatDP x) { return sqrt_rnd(x.dbl); }
FloatDP exp(FloatDP x) { return exp_rnd(x.dbl); }
FloatDP log(FloatDP x) { return log_rnd(x.dbl); }
FloatDP sin(FloatDP x) { return sin_rnd(x.dbl); }
FloatDP cos(FloatDP x) { return cos_rnd(x.dbl); }
FloatDP tan(FloatDP x) { return tan_rnd(x.dbl); }
//FloatDP asin(FloatDP x) { return asin_rnd(x.dbl); }
//FloatDP acos(FloatDP x) { return acos_rnd(x.dbl); }
FloatDP atan(FloatDP x) { return atan_rnd(x.dbl); }
FloatMP operator+(FloatMP const& x1, FloatMP const& x2);
FloatMP operator-(FloatMP const& x1, FloatMP const& x2);
FloatMP operator*(FloatMP const& x1, FloatMP const& x2);
FloatMP operator/(FloatMP const& x1, FloatMP const& x2);
FloatMP& operator+=(FloatMP& x1, FloatMP const& x2);
FloatMP& operator-=(FloatMP& x1, FloatMP const& x2);
FloatMP& operator*=(FloatMP& x1, FloatMP const& x2);
FloatMP& operator/=(FloatMP& x1, FloatMP const& x2);
FloatMP sqr(FloatMP const& x);
FloatMP rec(FloatMP const& x);
FloatMP add(FloatMP const& x1, FloatMP const& x2);
FloatMP sub(FloatMP const& x1, FloatMP const& x2);
FloatMP mul(FloatMP const& x1, FloatMP const& x2);
FloatMP div(FloatMP const& x1, FloatMP const& x2);
FloatMP fma(FloatMP const& x1, FloatMP const& x2, FloatMP const& x3);
FloatMP pow(FloatMP const& x, Int n);
FloatMP sqrt(FloatMP const& x);
FloatMP exp(FloatMP const& x);
FloatMP log(FloatMP const& x);
FloatMP sin(FloatMP const& x);
FloatMP cos(FloatMP const& x);
FloatMP tan(FloatMP const& x);
FloatMP asin(FloatMP const& x);
FloatMP acos(FloatMP const& x);
FloatMP atan(FloatMP const& x);
*/
} // namespace Ariadne

using namespace std;
using namespace Ariadne;

template<class FLT>
class TestRounded
{
    using FloatType = FLT;
    using RoundedFloatType = Rounded<FLT>;
    using PR=typename FLT::PrecisionType;
    PR precision;

  public:
    TestRounded(PR prec);
    Void test();
  private:
    Void test_concept();
    Void test_class();
    Void test_conversion_from_to();
    Void test_stream();
    Void test_comparison();
    Void test_rounding();
    Void test_arithmetic();
    Void test_cosine();
    Void test_arctan();
    Void test_function();
};


Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

#warning
    ARIADNE_TEST_CLASS(Rounded<FloatDP>,TestRounded<FloatDP>(dp));
//    ARIADNE_TEST_CLASS(Rounded<FloatMP>,TestRounded<FloatMP>(MultiplePrecision(64_bits)));
//    ARIADNE_TEST_CLASS(Rounded<FloatMP>,TestRounded<FloatMP>(MultiplePrecision(192_bits)));


    return ARIADNE_TEST_FAILURES;
}


template<class FLT>
TestRounded<FLT>::TestRounded(PR pr)
    : precision(pr)
{
}

template<class FLT> Void
TestRounded<FLT>::test()
{
    ARIADNE_TEST_CALL(test_class());
    ARIADNE_TEST_CALL(test_conversion_from_to());
    ARIADNE_TEST_CALL(test_stream());
    ARIADNE_TEST_CALL(test_comparison());
    ARIADNE_TEST_CALL(test_rounding());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_function());
    ARIADNE_TEST_CALL(test_cosine());
    ARIADNE_TEST_CALL(test_arctan());
}


// Test that the type implements all operations of
// the RoundedFloatType concept without testing correctness
template<class FLT> Void
TestRounded<FLT>::test_concept()
{
    Bool b=true;
    Int n=1;
    Nat m=1;
    double d=1;
    RoundedFloatType x=1, x2=1;

    // Constructors
    x=RoundedFloatType(); x=RoundedFloatType(n); x=RoundedFloatType(m); x=RoundedFloatType(d); x=RoundedFloatType(x);

    // Assignment
    x=n; x=m; x=d; x=x2;

    // Conversion
    d=x.get_d();

    // Maximum and minimum and absolute value
    x=max(x,x); x=min(x,x); x=abs(x);


    // ExactTag operations
    x=nul(x); x=pos(x); x=neg(x); x=hlf(x);

    // Rounded arithmetic operations
    x=add(near,x,x); x=add(approx,x,x); x=add(down,x,x); x=add(up,x,x); // x=add_chop(x,x);
    x=sub(near,x,x); x=add(approx,x,x); x=sub(down,x,x); x=sub(up,x,x); // x=sub_chop(x,x);
    x=mul(near,x,x); x=add(approx,x,x); x=mul(down,x,x); x=mul(up,x,x); // x=mul_chop(x,x);
    x=div(near,x,x); x=add(approx,x,x); x=div(down,x,x); x=div(up,x,x); // x=div_chop(x,x);
    x=fma(approx,x,x,x); x=fma(down,x,x,x); x=fma(up,x,x,x);
    x=pow(approx,x,n); x=pow(down,x,n); x=pow(up,x,n); // x=pow_chop(x,n);
    x=pow(approx,x,m); x=pow(down,x,m); x=pow(up,x,m); // x=pow_chop(x,m);

    // Non-exact operations
    x=sqr(approx,x); x=sqr(down,x); x=sqr(up,x); // x=sqr_chop(x);
    x=rec(approx,x); x=rec(down,x); x=rec(up,x); // x=rec_chop(x);
    x=sqrt(approx,x); x=sqrt(down,x); x=sqrt(up,x); // x=sqrt_chop(x);
    x=exp(approx,x); x=exp(down,x); x=exp(up,x); // x=exp_chop(x);
    x=log(approx,x); x=log(down,x); x=log(up,x); // x=log_chop(x);
    x=sin(approx,x); x=sin(down,x); x=sin(up,x); // x=sin_chop(x);
    x=cos(approx,x); x=cos(down,x); x=cos(up,x); // x=cos_chop(x);
    x=tan(approx,x); x=tan(down,x); x=tan(up,x); // x=tan_chop(x);
    x=asin(approx,x); x=asin(down,x); x=asin(up,x); // x=asin_chop(x);
    x=acos(approx,x); x=acos(down,x); x=acos(up,x); // x=acos_chop(x);
    x=atan(approx,x); x=atan(down,x); x=atan(up,x); // x=atan_chop(x);

    x=med(near,x,x); x=rad(up,x,x);

    // Mixed RoundedFloatType/Int arithmetic
    x=mul(approx,n,x); x=mul(down,n,x); x=mul(up,n,x); // x=mul_chop(n,x);
    x=mul(approx,m,x); x=mul(down,m,x); x=mul(up,m,x); // x=mul_chop(m,x);
    x=mul(approx,x,n); x=mul(down,x,n); x=mul(up,x,n); // x=mul_chop(x,n);
    x=mul(approx,x,m); x=mul(down,x,m); x=mul(up,x,m); // x=mul_chop(x,m);
    x=div(approx,x,n); x=div(down,x,n); x=div(up,x,n); // x=div_chop(x,n);
    x=div(approx,x,m); x=div(down,x,m); x=div(up,x,m); // x=div_chop(x,m);

    // Mixed RoundedFloatType/double arithmetic
    x=mul(approx,d,x); x=mul(approx,x,d); x=div(approx,x,d);

    // Reset x to zero
    x=0; x=0.0;

    // Reset x to 1
    x=1; x=1.0;

    // Operators in rounding mode
    x=+x; x=-x;
    x=x+x; x=x-x; x=x*x; x=x/x;
    x+x; x-=x2; x*=x2; x/=x2;

    // Comparisons
    b=(x==n); b=(x!=n); b=(x<=n); b=(x>=n); b=(x<n); b=(x>n);
    b=(n==x); b=(n!=x); b=(n<=x); b=(n>=x); b=(n<x); b=(n>x);
    b=(x==m); b=(x!=m); b=(x<=m); b=(x>=m); b=(x<m); b=(x>m);
    b=(m==x); b=(m!=x); b=(m<=x); b=(m>=x); b=(m<x); b=(m>x);
    b=(x==d); b=(x!=d); b=(x<=d); b=(x>=d); b=(x<d); b=(x>d);
    b=(d==x); b=(d!=x); b=(d<=x); b=(d>=x); b=(d<x); b=(d>x);
    b=(x==x); b=(x!=x); b=(x<=x); b=(x>=x); b=(x<x); b=(x>x);

    // Rounded mode
    RoundedFloatType::set_rounding_to_nearest();
    RoundedFloatType::set_rounding_downward();
    RoundedFloatType::set_rounding_upward();
    RoundedFloatType::set_rounding_toward_zero();

    RoundedFloatType::set_rounding_mode(to_nearest);
    RoundedFloatType::set_rounding_mode(downward);
    RoundedFloatType::set_rounding_mode(upward);
    RoundedFloatType::set_rounding_mode(toward_zero);

    RoundedFloatType::set_rounding_mode(near);
    RoundedFloatType::set_rounding_mode(down);
    RoundedFloatType::set_rounding_mode(up);

    typename RoundedFloatType::RoundedModeType rnd=RoundedFloatType::get_rounding_mode();
    RoundedFloatType::set_rounding_mode(rnd);

    // DoublePrecision
    typename RoundedFloatType::PrecisionType pr=RoundedFloatType::get_default_precision();
    pr=x.precision();
    x.set_precision(pr);

}


template<class FLT> Void
TestRounded<FLT>::test_class()
{
    cout << __PRETTY_FUNCTION__ << endl;
    // Construct from an Int
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f1,(2));
    ARIADNE_TEST_EQUALS(f1,2);
    // Construct from a double
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f2,(1.25));
    ARIADNE_TEST_EQUALS(f2,1.25);
    // Copy constructor
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f3,(f2));
    ARIADNE_TEST_EQUAL(f3,f2);


    // Assign from an Int
    ARIADNE_TEST_EXECUTE(f1=3);
    ARIADNE_TEST_EQUALS(f1,3);
    // Assign from a double
    ARIADNE_TEST_EXECUTE(f2=2.25);
    ARIADNE_TEST_EQUALS(f2,2.25);
    // Copy assignment
    ARIADNE_TEST_EXECUTE(f3=f2);
    ARIADNE_TEST_EQUAL(f3,f2);
    // Self-assignment
    ARIADNE_TEST_EXECUTE(f3=f3);
    ARIADNE_TEST_EQUAL(f3,f2);

    if(not (precision==FloatType::get_default_precision())) {
        PR low_precision=min(precision,FloatType::get_default_precision());
        PR high_precision=max(precision,FloatType::get_default_precision());
        ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f4,(2.25,high_precision));
        ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f5,(3.75,low_precision));
        ARIADNE_TEST_EXECUTE(f5=f4);
        ARIADNE_TEST_EQUAL(f5.precision(),high_precision);
        ARIADNE_TEST_EQUAL(f5,f4);
        ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f6,(low_precision));
        ARIADNE_TEST_EXECUTE(f6=std::move(f4));
        ARIADNE_TEST_EQUAL(f6.precision(),high_precision);
        ARIADNE_TEST_EQUAL(f6,f5);
    }


}



template<class FLT> Void
TestRounded<FLT>::test_conversion_from_to()
{
    // Convert from dyadic
    RoundedFloatType f1(Dyadic(5,2u),precision);
    ARIADNE_TEST_EQUAL(f1,1.25);

    // Convert from self
    RoundedFloatType f2(f1,precision);
    ARIADNE_TEST_EQUAL(f1,f2);
    RoundedFloatType f3(f1,precision);
    ARIADNE_TEST_EQUAL(f1,f3);

    // Convert from integers
    Int n;
    n=std::numeric_limits<Int>::min();
    ARIADNE_TEST_ASSERT(RoundedFloatType(n)==n);
    n=std::numeric_limits<Int>::max();
    ARIADNE_TEST_ASSERT(RoundedFloatType(n)==n);
    Nat nn = std::numeric_limits<Nat>::max();
    ARIADNE_TEST_ASSERT(RoundedFloatType(nn)==nn);

    // Convert to a rational
    ARIADNE_TEST_EQUAL(Rational(RoundedFloatType(1.0).raw()),Rational(1));
    ARIADNE_TEST_EQUAL(Rational(RoundedFloatType(2.25).raw()),Rational(9,4));

    // Convert from a rational
    Int num=1; Int den=3;
    Rational q(1,3);
    RoundedFloatType::set_rounding_downward();
    ARIADNE_TEST_COMPARE(Rational(RoundedFloatType(q,precision).raw()),<,q);
    RoundedFloatType::set_rounding_upward();
    ARIADNE_TEST_COMPARE(Rational(RoundedFloatType(Rational(num,den),precision).raw()),>,q);

    // Convert from a negative rational
    num=-2; den=5;
    RoundedFloatType::set_rounding_downward();
    ARIADNE_TEST_COMPARE(Rational(RoundedFloatType(q,precision).raw()),<,q);
    RoundedFloatType::set_rounding_upward();
    ARIADNE_TEST_COMPARE(Rational(RoundedFloatType(q,precision).raw()),>,q);
}

template<class FLT> Void
TestRounded<FLT>::test_stream()
{
#warning
/*
    cout << __PRETTY_FUNCTION__ << endl;

    stringstream ss("1.25 -2.25 42 2.375e1 2.35e1");
    RoundedFloatType f1,f2,f3,f4,f5;
    ss >> f1;
    cout << f1 << endl;
    ARIADNE_TEST_EQUALS(f1,1.25);
    ss >> f2;
    cout << f2 << endl;
    ARIADNE_TEST_EQUALS(f2,-2.25);
    ss >> f3;
    cout << f3 << endl;
    ARIADNE_TEST_EQUALS(f3,42.0);
    try {
        ss >> f4;
        cout << f4 << endl;
        if(f4!=23.75) {
            ARIADNE_TEST_WARN("Cannot create RoundedFloatType<"<<class_name<PR>()<<"> from string literal in exponential form 2.375e1");
        }
    }
    catch(const std::exception& e) {
        cerr << e.what() << endl;
    }

    try {
        ss >> f5;
        cout << f5 << endl;
        if(f4!=23.5) {
            ARIADNE_TEST_WARN("Cannot create RoundedFloatType<"<<class_name<PR>()<<"> from string literal in exponential form 2.35e1");
        }
    }
    catch(const std::exception& e) {
        cerr << e.what() << endl;
    }
*/
}


template<class FLT> Void
TestRounded<FLT>::test_comparison()
{
#warning
/*
    cout << __PRETTY_FUNCTION__ << endl;

    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f0,(3.00));
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f1,(1.25));
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f2,(-1.25));
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f3,(-2.25));
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,f4,(1.25));

    // Test comparison of two equal numbers
    ARIADNE_TEST_ASSERT(f1==f4); ARIADNE_TEST_ASSERT(!(f1!=f4));
    ARIADNE_TEST_ASSERT(f1<=f4); ARIADNE_TEST_ASSERT(!(f1> f4));
    ARIADNE_TEST_ASSERT(f1>=f4); ARIADNE_TEST_ASSERT(!(f1< f4));


    // Test comparison of two different numbers
    ARIADNE_TEST_ASSERT(!(f1==f2)); ARIADNE_TEST_ASSERT(f1!=f2);
    ARIADNE_TEST_ASSERT(!(f1<=f2)); ARIADNE_TEST_ASSERT(f1> f2);
    ARIADNE_TEST_ASSERT(f1>=f2); ARIADNE_TEST_ASSERT(!(f1< f2));

    // Test comparison of two negative numbers
    ARIADNE_TEST_ASSERT(!(f2==f3)); ARIADNE_TEST_ASSERT(f2!=f3);
    ARIADNE_TEST_ASSERT(!(f2<=f3)); ARIADNE_TEST_ASSERT(f2> f3);
    ARIADNE_TEST_ASSERT(f2>=f3); ARIADNE_TEST_ASSERT(!(f2< f3));

    // Test comparison with in int
    Int i2=1;
    ARIADNE_TEST_ASSERT(!(f1==i2)); ARIADNE_TEST_ASSERT(f1!=i2);
    ARIADNE_TEST_ASSERT(!(f1<=i2)); ARIADNE_TEST_ASSERT(f1> i2);
    ARIADNE_TEST_ASSERT(f1>=i2); ARIADNE_TEST_ASSERT(!(f1< i2));

    Int i1=1;
    ARIADNE_TEST_ASSERT(!(i1==f2)); ARIADNE_TEST_ASSERT(i1!=f2);
    ARIADNE_TEST_ASSERT(!(i1<=f2)); ARIADNE_TEST_ASSERT(i1> f2);
    ARIADNE_TEST_ASSERT(i1>=f2); ARIADNE_TEST_ASSERT(!(i1< f2));

    // Test comparison with a double
    double x2=1.0;
    ARIADNE_TEST_ASSERT(!(f1==x2)); ARIADNE_TEST_ASSERT(f1!=x2);
    ARIADNE_TEST_ASSERT(!(f1<=x2)); ARIADNE_TEST_ASSERT(f1> x2);
    ARIADNE_TEST_ASSERT(f1>=x2); ARIADNE_TEST_ASSERT(!(f1< x2));

    double x1=1.0;
    ARIADNE_TEST_ASSERT(!(x1==f2)); ARIADNE_TEST_ASSERT(x1!=f2);
    ARIADNE_TEST_ASSERT(!(x1<=f2)); ARIADNE_TEST_ASSERT(x1> f2);
    ARIADNE_TEST_ASSERT(x1>=f2); ARIADNE_TEST_ASSERT(!(x1< f2));

    // Test comparison with an integer
    ARIADNE_TEST_CONSTRUCT(Integer,z0,(-3));
    ARIADNE_TEST_ASSERT(!(f1==z0)); ARIADNE_TEST_ASSERT(f1!=z0);
    ARIADNE_TEST_ASSERT(!(f1<=z0)); ARIADNE_TEST_ASSERT(f1> z0);
    ARIADNE_TEST_ASSERT(f1>=z0); ARIADNE_TEST_ASSERT(!(f1< z0));
    ARIADNE_TEST_ASSERT(!(z0==f1)); ARIADNE_TEST_ASSERT(z0<f1);

    // Test comparison with a dyadic
    ARIADNE_TEST_CONSTRUCT(Dyadic,w2,(-5,2u));
    ARIADNE_TEST_ASSERT(!(f1==w2)); ARIADNE_TEST_ASSERT(f1!=w2);
    ARIADNE_TEST_ASSERT(!(w2==f1)); ARIADNE_TEST_ASSERT(w2<f1);
    ARIADNE_TEST_ASSERT(!(f1<=w2)); ARIADNE_TEST_ASSERT(f1> w2);
    ARIADNE_TEST_ASSERT(f1>=w2); ARIADNE_TEST_ASSERT(!(f1< w2));
    ARIADNE_TEST_ASSERT(f2==w2); ARIADNE_TEST_ASSERT(!(f2!=w2));

    // Test comparison with a rational
    ARIADNE_TEST_CONSTRUCT(Rational,q2,(-5,4));
    ARIADNE_TEST_ASSERT(!(f1==q2)); ARIADNE_TEST_ASSERT(f1!=q2);
    ARIADNE_TEST_ASSERT(!(f1<=q2)); ARIADNE_TEST_ASSERT(f1> q2);
    ARIADNE_TEST_ASSERT(f1>=q2); ARIADNE_TEST_ASSERT(!(f1< q2));
    ARIADNE_TEST_ASSERT(f2==q2); ARIADNE_TEST_ASSERT(!(f2!=q2));
*/
}

template<> Void
TestRounded<FloatDP>::test_rounding()
{
    volatile double one   = 1;
    volatile double two_  = 2;
    volatile double three = 3;
    volatile double five  = 5;
    const double onethirddown    = 0.33333333333333331483;
    const double onethirdup      = 0.33333333333333337034;
    const double onethirdchop    = 0.33333333333333331483;
    const double onethirdnearest = 0.33333333333333331483;
    const double twofifthsdown   = 0.39999999999999996669;
    const double twofifthsup     = 0.40000000000000002220;
    const double twofifthschop   = 0.39999999999999996669;
    const double twofifthsnearest= 0.40000000000000002220;

    Ariadne::set_rounding_mode(downward);
    double onethirdrounddown=one/three;
    ARIADNE_TEST_EQUAL(onethirdrounddown, onethirddown);
    Ariadne::set_rounding_mode(upward);
    double onethirdroundup=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundup, onethirdup);
    Ariadne::set_rounding_mode(toward_zero);
    double onethirdroundchop=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundchop, onethirdchop);
    Ariadne::set_rounding_mode(to_nearest);
    double onethirdroundnearest=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundnearest, onethirdnearest);

    Ariadne::set_rounding_downward();
    double twofifthsrounddown=two_/five;
    ARIADNE_TEST_EQUAL(twofifthsrounddown, twofifthsdown);
    Ariadne::set_rounding_upward();
    double twofifthsroundup=two_/five;
    ARIADNE_TEST_EQUAL(twofifthsroundup, twofifthsup);
    Ariadne::set_rounding_toward_zero();
    double twofifthsroundchop=two_/five;
    ARIADNE_TEST_EQUAL(twofifthsroundchop, twofifthschop);
    Ariadne::set_rounding_to_nearest();
    double twofifthsroundnearest=two_/five;
    ARIADNE_TEST_EQUAL(twofifthsroundnearest, twofifthsnearest);
}

template<class FLT> Void
TestRounded<FLT>::test_rounding()
{
    volatile double one   = 1;
    volatile double two_  = 2;
    volatile double three = 3;
    volatile double five  = 5;
    const double onethirddown    = 0.33333333333333331483;
    const double onethirdup      = 0.33333333333333337034;
    const double onethirdchop    = 0.33333333333333331483;
    const double onethirdnearest = 0.33333333333333331483;
    const double twofifthsdown   = 0.39999999999999996669;
    const double twofifthsup     = 0.40000000000000002220;
    const double twofifthschop   = 0.39999999999999996669;
    const double twofifthsnearest= 0.40000000000000002220;

    Ariadne::set_rounding_mode(Ariadne::downward);
    double onethirdrounddown=one/three;
    ARIADNE_TEST_EQUAL(onethirdrounddown, onethirddown);
    Ariadne::set_rounding_mode(Ariadne::upward);
    double onethirdroundup=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundup, onethirdup);
    Ariadne::set_rounding_mode(Ariadne::toward_zero);
    double onethirdroundchop=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundchop, onethirdchop);
    Ariadne::set_rounding_mode(Ariadne::to_nearest);
    double onethirdroundnearest=one/three;
    ARIADNE_TEST_EQUAL(onethirdroundnearest, onethirdnearest);

    Ariadne::set_rounding_downward();
    double twofifthsrounddown=two_/five;
    ARIADNE_TEST_EQUAL(twofifthsrounddown, twofifthsdown);
    Ariadne::set_rounding_upward();
    double twofifthsroundup=two_/five;
    ARIADNE_TEST_EQUAL(twofifthsroundup, twofifthsup);
    Ariadne::set_rounding_toward_zero();
    double twofifthsroundchop=two_/five;
    ARIADNE_TEST_EQUAL(twofifthsroundchop, twofifthschop);
    Ariadne::set_rounding_to_nearest();
    double twofifthsroundnearest=two_/five;
    ARIADNE_TEST_EQUAL(twofifthsroundnearest, twofifthsnearest);
}

template<class FLT> Void
TestRounded<FLT>::test_arithmetic()
{
    cout << __PRETTY_FUNCTION__ << endl;

    static bool full_precision_warning = false;
    RoundedFloatType::set_rounding_to_nearest();
    const RoundedFloatType pi_near=RoundedFloatType::pi(precision);
    const RoundedFloatType four   =RoundedFloatType(4,precision);
    RoundedFloatType::set_rounding_upward();
    if( mul(pi_near,4.0) != mul(pi_near,four) ) {
        if(not full_precision_warning) {
            full_precision_warning=true;
            ARIADNE_TEST_WARN("Mixed RoundedFloatType/BuiltinTag operations may not use full precision.");
        }
    }

    static const double eps = std::numeric_limits<double>::epsilon();

    // Set next_up some variables
    RoundedFloatType f1(1.25); RoundedFloatType f2(2.25); RoundedFloatType f3(-3.25); RoundedFloatType f4; RoundedFloatType f5;

    // Minimum (this should always remain exact)
    f4=min(f1,f2);
    cout << "min(" << f1 << "," << f2 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f1);
    f4=min(f1,f3);
    cout << "min(" << f1 << "," << f3 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f3);

    // Maximum (this should always remain exact)
    f4=max(f1,f2);
    cout << "max(" << f1 << "," << f2 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f2);
    f4=max(f1,f3);
    cout << "max(" << f1 << "," << f3 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==f1);

    // Absolute value (this should always remain exact)
    f4=abs(f1);
    cout << "abs(" << f1 << ") = " << f4 << endl;
    ARIADNE_TEST_ASSERT(f4==1.25);
    f5=abs(f3);
    cout << "abs(" << f3 << ") = " << f5 << endl;
    ARIADNE_TEST_ASSERT(f5==3.25);

    // Median (this should remain exact here)
    f3=med(f1,f2);
    cout << f1 << " <= med(" << f1 << "," << f2 << ")=" << f3 << " <= " << f2 << endl;
    ARIADNE_TEST_ASSERT(f1<=f3); ARIADNE_TEST_ASSERT(f3<=f2);
    ARIADNE_TEST_ASSERT(f3==1.75);

    // Negation (this should always remain exact)
    f3=neg(f1);
    cout << "neg(" << f1 << ") = " << f3 << endl;
    f3=-f1;
    cout << "- " << f1 << " = " << f3 << endl;
    ARIADNE_TEST_ASSERT(f3==-1.25);

    // Addition (this should remain exact here)
    f3=add(f1,f2);
    f4=add(f1,f2);
    cout << f3 << " <= " << f1 << " + " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=3.5); ARIADNE_TEST_ASSERT(f4>=3.5);
    // Addition
    ARIADNE_TEST_EXECUTE(RoundedFloatType::set_rounding_upward());
    ARIADNE_TEST_COMPARE(add(RoundedFloatType(1.0),RoundedFloatType(eps/2)),>,1.0);
    ARIADNE_TEST_EXECUTE(RoundedFloatType::set_rounding_downward());
    ARIADNE_TEST_COMPARE(add(RoundedFloatType(1.0),RoundedFloatType(eps)),>,1.0);
    if(add(RoundedFloatType(1.0),RoundedFloatType(eps/2)) != RoundedFloatType(1.0)) {
        ARIADNE_TEST_WARN("Results of floating-point operations stored to higher-precision in registers than memory.");
        ARIADNE_TEST_COMPARE(add(RoundedFloatType(1.0),RoundedFloatType(eps/(1<<11))),>,RoundedFloatType(1.0));
    }
    ARIADNE_TEST_COMPARE(add(RoundedFloatType(1.0),RoundedFloatType(eps/(1<<12))),==,RoundedFloatType(1.0));
    ARIADNE_TEST_CONSTRUCT(RoundedFloatType,one_add_down_half_epsilon,(add(RoundedFloatType(1.0),RoundedFloatType(eps/2))));
    ARIADNE_TEST_COMPARE(one_add_down_half_epsilon,==,RoundedFloatType(1.0));
    ARIADNE_TEST_EQUAL(RoundedFloatType(1.0),1.0);
    ARIADNE_TEST_COMPARE(RoundedFloatType(1.0),==,1.0);


    // Subtraction (this should remain exact here)
    RoundedFloatType::set_rounding_downward();
    f3=sub(f1,f2);
    RoundedFloatType::set_rounding_upward();
    f4=sub(f1,f2);
    cout << f3 << " <= " << f1 << " - " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=-1); ARIADNE_TEST_ASSERT(f4>=-1);

    // Multiplication (this should remain exact here)
    RoundedFloatType::set_rounding_downward();
    f3=mul(f1,f2);
    RoundedFloatType::set_rounding_upward();
    f4=mul(f1,f2);
    cout << f3 << " <= " << f1 << " * " << f2 << " <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=2.8125); ARIADNE_TEST_ASSERT(f4>=2.8125);

    // Division (not exact; should catch errors here)
    RoundedFloatType::set_rounding_downward();
    f3=div(f1,f2);
    RoundedFloatType::set_rounding_upward();
    f4=div(f1,f2);
    RoundedFloatType::set_rounding_to_nearest();
    f5=div(f1,f2);
    cout << f3 << " <= " << f1 << " / " << f2 << " <= " << f4 << endl;
    RoundedFloatType five = 5;
    RoundedFloatType nine = 9;

    RoundedFloatType expected_five_ninths_up=0.55555555555555558023;
    RoundedFloatType::set_rounding_downward();
    ARIADNE_TEST_COMPARE(five/nine,<,expected_five_ninths_up);
    ARIADNE_TEST_COMPARE(div(five,nine),<,expected_five_ninths_up);
    RoundedFloatType five_divby_nine_down=five;
    five_divby_nine_down/=nine;
    ARIADNE_TEST_COMPARE(five_divby_nine_down,<,expected_five_ninths_up);

    RoundedFloatType::set_rounding_downward();
    RoundedFloatType five_ninths_down = div(five,nine);
    RoundedFloatType five_down = mul(five_ninths_down,nine);
    RoundedFloatType::set_rounding_upward();
    RoundedFloatType five_ninths_up = div(five,nine);
    RoundedFloatType five_up = mul(five_ninths_up,nine);
    ARIADNE_TEST_COMPARE(five_ninths_down, < , five_ninths_up);
    ARIADNE_TEST_COMPARE(five_down, < , five);
    ARIADNE_TEST_COMPARE(five_up, > , five);

    RoundedFloatType::set_rounding_downward();
    RoundedFloatType pi_down=RoundedFloatType::pi(precision);
    ARIADNE_TEST_EQUAL(pi_down/2,div(pi_down,2));
    ARIADNE_TEST_EQUAL(pi_down/2,hlf(pi_down));


    // Power (not exact; should catch errors here)
    RoundedFloatType::set_rounding_downward();
    f3=pow(f1,3);
    RoundedFloatType::set_rounding_upward();
    f4=pow(f1,3);
    cout << f3 << " <= pow(" << f1 << ",3) <= " << f4 << endl;
    ARIADNE_TEST_ASSERT(f3<=1.953125); ARIADNE_TEST_ASSERT(f4>=1.953125);

    RoundedFloatType::set_rounding_downward();
    f3=pow(f1,-2);
    RoundedFloatType::set_rounding_upward();
    f4=pow(f1,-2);
    cout << f3 << " <= pow(" << f1 << ",-2) <= " << f4 << endl;
    //ARIADNE_TEST_ASSERT(Rational(f3)<Rational(16,25));
    //ARIADNE_TEST_ASSERT(Rational(f4)>Rational(16,25));

    // Floor and ceiling
    f2=RoundedFloatType(-3.25); f3=RoundedFloatType(-2);

    ARIADNE_TEST_ASSERT(floor(f1)==1); ARIADNE_TEST_ASSERT(ceil(f1)==2);
    ARIADNE_TEST_ASSERT(floor(f2)==-4); ARIADNE_TEST_ASSERT(ceil(f2)==-3);
    ARIADNE_TEST_ASSERT(floor(f3)==-2); ARIADNE_TEST_ASSERT(ceil(f3)==-2);

    // Conversion to integer types
    Int i3,i4;
    i3=integer_cast<Int>(floor(f1));
    i4=integer_cast<Int>(ceil(f1));
    cout << i3 << " < " << f1 << " < " << i4 << endl;
    ARIADNE_TEST_ASSERT(i3==1); ARIADNE_TEST_ASSERT(i4==2);
    i3=integer_cast<Int>(floor(f2));
    i4=integer_cast<Int>(ceil(f2));
    cout << i3 << " < " << f2 << " < " << i4 << endl;
    ARIADNE_TEST_ASSERT(i3==-4); ARIADNE_TEST_ASSERT(i4==-3);

    // Check interval conversions
    RoundedFloatType z(0); RoundedFloatType o(1); RoundedFloatType t(3);

    RoundedFloatType odtd=div(down,o,t);
    RoundedFloatType odtu=div(up,o,t);
    cout << odtd << " <= 1/3 <= " << odtu << endl;
    ARIADNE_TEST_COMPARE(odtd,<,odtu);
    // Regression test to catch errors when RoundedFloatType result is not assigned to a variable
    cout << div(down,o,t) << " <= 1/3 <= " << div(up,o,t) << endl;
    ARIADNE_TEST_COMPARE(div(down,o,t),<,div(up,o,t));

    ARIADNE_TEST_COMPARE(med(near,RoundedFloatType(2),RoundedFloatType(3)),==,2.5);
    ARIADNE_TEST_COMPARE(rad(up,RoundedFloatType(2),RoundedFloatType(3)),>=,0.5);
    ARIADNE_TEST_COMPARE(rad(up,RoundedFloatType(2),RoundedFloatType(3)),<=,0.5000000000000002);

    // The following line should not compile
    // f5=f1+f2;

}


template<class FLT> Void
TestRounded<FLT>::test_function()
{
    cout << __PRETTY_FUNCTION__ << endl;

    cout << setprecision(20);

    // Set up some variables
    RoundedFloatType x; RoundedFloatType ra; RoundedFloatType rl; RoundedFloatType ru;

    x=1;
    ra=exp(x);
    rl=next(down,ra);
    ru=next(up,ra);
    ARIADNE_TEST_PRINT(rl);
    ARIADNE_TEST_PRINT(ru);
    ARIADNE_TEST_ASSERT(rl<ru);
    ARIADNE_TEST_ASSERT(2.71<rl);
    ARIADNE_TEST_ASSERT(ru<2.72);


    //The following don't work as rounded operators not exported.
    //test_inverse_pair("sin",&sin_down,&sin_up,&asin_down,&asin_up);
    //The following don't work as acos is decreasing
    //test_inverse_pair("cos",&cos_down,&cos_up,&acos_down,&acos_up);
    //test_inverse_pair("cos",&cos_up,&cos_down,&acos_down,&acos_up);
    //test_inverse_pair("tan",&tan_down,&tan_up,&atan_down,&atan_up);
    //test_inverse_pair("sinh",&sinh_down,&sinh_up,&asinh_down,&asinh_up);
    //test_inverse_pair("cosh",&cosh_down,&cosh_up,&acosh_down,&acosh_up);
    //test_inverse_pair("tanh",&tanh_down,&tanh_up,&atanh_down,&atanh_up);
}

template<class FLT> Void
TestRounded<FLT>::test_cosine()
{
    //3.14159265358979323846264338327950288419716939937510
    const RoundedFloatType pi_down=FloatType::pi(down,precision);
    const RoundedFloatType pi_up  =FloatType::pi(up,precision);
    const RoundedFloatType three  =FloatType(3,precision);

    static const RoundedFloatType third_pi_down=div(down,pi_down,three);
    static const RoundedFloatType third_pi_up  =div(up  ,pi_up  ,three);

    ARIADNE_TEST_COMPARE(cos(up,div(down,pi_down,2)),>,0.0);

    CurrentRoundingMode rndm;

    RoundedFloatType::set_rounding_mode(upward);
    ARIADNE_TEST_COMPARE(cos(hlf(pi_down)),>,0.0);
    ARIADNE_TEST_EQUAL(cos(RoundedFloatType(0.0)),1.0);
    ARIADNE_TEST_COMPARE(cos(third_pi_down),>,0.5);
    ARIADNE_TEST_COMPARE(sqr(cos(pi_down/4)),>,0.5);
    ARIADNE_TEST_COMPARE(cos(pi_down/2),>,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_up/2),<=,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_down),>,-1.0);
    ARIADNE_TEST_COMPARE(cos(pi_up),>,-1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_down),==,1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_up),==,1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_down),>,-1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_up),>,-1.0);

    RoundedFloatType::set_rounding_mode(downward);
    ARIADNE_TEST_EQUAL(cos(RoundedFloatType(0.0)),1.0);
    ARIADNE_TEST_COMPARE(cos(third_pi_up),<,0.5);
    ARIADNE_TEST_COMPARE(sqr(cos(pi_up/4)),<,0.5);
    ARIADNE_TEST_COMPARE(cos(pi_down/2),>=,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_up/2),<,0.0);
    ARIADNE_TEST_COMPARE(cos(pi_down),==,-1.0);
    ARIADNE_TEST_COMPARE(cos(pi_up),==,-1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_down),<,1.0);
    ARIADNE_TEST_COMPARE(cos(2*pi_up),<,1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_down),==,-1.0);
    ARIADNE_TEST_COMPARE(cos(3*pi_up),==,-1.0);

    static const RoundedFloatType one=1.0;
    RoundedFloatType::set_rounding_mode(upward);
    RoundedFloatType sin_rnd_up_one=sin(one);
    RoundedFloatType::set_rounding_mode(downward);
    RoundedFloatType sin_rnd_down_one=sin(one);
    RoundedFloatType::set_rounding_mode(to_nearest);
    RoundedFloatType sin_rnd_approx_one=sin(one);
    ARIADNE_TEST_COMPARE(sin_rnd_down_one,<=,sin_rnd_approx_one);
    ARIADNE_TEST_COMPARE(sin_rnd_approx_one,<=,sin_rnd_up_one);
    ARIADNE_TEST_COMPARE(sin_rnd_down_one,< ,sin_rnd_up_one);
}

template<> Void
TestRounded<FloatDP>::test_arctan()
{
    PR pr;

    static const RoundedFloatType pi_down=3.1415926535897931;
    //static const FloatType pi_near=3.1415926535897931;
    static const RoundedFloatType pi_up  =3.1415926535897936;

    static const RoundedFloatType atan_quarter_down=0.244978663126864143;
    //static const FloatType atan_quarter_near=0.244978663126864154;
    static const RoundedFloatType atan_quarter_up  =0.244978663126864171;

    static const RoundedFloatType sqrt_three_down=1.732050807568877193;
    //static const FloatType sqrt_three_near=1.732050807568877294;
    static const RoundedFloatType sqrt_three_up  =1.732050807568877415;

    static const RoundedFloatType zero(0.0,pr);
    static const RoundedFloatType one(1.0,pr);
    static const RoundedFloatType eps(FloatType::eps(pr));

    RoundedFloatType::set_rounding_mode(upward);
    ARIADNE_TEST_EQUAL(atan(zero),zero);
    ARIADNE_TEST_COMPARE(atan(one)*4,>=,pi_up);
    ARIADNE_TEST_COMPARE(atan(-one)*4,>=,-pi_down);
    ARIADNE_TEST_COMPARE(atan(sqrt_three_up)*3,>=,pi_up);
    ARIADNE_TEST_COMPARE(atan(-sqrt_three_down)*3,>=,-pi_down);
    ARIADNE_TEST_COMPARE(atan(1/sqrt_three_down)*6,>=,pi_up);
    ARIADNE_TEST_COMPARE(atan(-1/sqrt_three_up)*6,>=,-pi_down);
    ARIADNE_TEST_COMPARE(atan(one/4),>=,atan_quarter_up);
    ARIADNE_TEST_COMPARE(atan(-one/4),>=,-atan_quarter_down);

    RoundedFloatType::set_rounding_mode(downward);
    ARIADNE_TEST_EQUAL(atan(zero),zero);
    ARIADNE_TEST_COMPARE(atan(+one)*4,<=,pi_down);
    ARIADNE_TEST_COMPARE(atan(-one)*4,<=,-pi_up);
    ARIADNE_TEST_COMPARE(atan(+sqrt_three_down)*3,<=,pi_down);
    ARIADNE_TEST_COMPARE(atan(-sqrt_three_up)*3,<=,-pi_up);
    ARIADNE_TEST_COMPARE(atan(1/sqrt_three_up)*6,<=,pi_down);
    ARIADNE_TEST_COMPARE(atan(-1/sqrt_three_down)*6,<=,-pi_up);
    ARIADNE_TEST_COMPARE(atan(+one/4),<=,atan_quarter_down);
    ARIADNE_TEST_COMPARE(atan(-one/4),<=,-atan_quarter_up);

    ARIADNE_TEST_COMPARE(atan(one)*4-pi_up,<=,4*eps);
    ARIADNE_TEST_COMPARE(atan(sqrt_three_up)*3-pi_up,<=,4*eps);
    ARIADNE_TEST_COMPARE(atan(1/sqrt_three_down)*6-pi_up,<=,4*eps);
    ARIADNE_TEST_COMPARE(atan(one/4)-atan_quarter_up,<=,eps/2);

}


template<class FLT> Void
TestRounded<FLT>::test_arctan()
{
    //pi~=3.14159265358979323846264338327950288419716939937510
    const RoundedFloatType pi_down=RoundedFloatType::pi(down,precision);
    const RoundedFloatType pi_near=RoundedFloatType::pi(near,precision);
    const RoundedFloatType pi_up  =RoundedFloatType::pi(up,precision);

    if constexpr(std::is_same<FloatType,FloatDP>::value) {
        ARIADNE_TEST_EQUAL(pi_down.raw(),3.1415926535897931);
        ARIADNE_TEST_EQUAL(pi_near.raw(),3.1415926535897931);
        ARIADNE_TEST_EQUAL(pi_up.raw(),  3.1415926535897936);
    }

    RoundedFloatType three=RoundedFloatType(3,precision);
    const RoundedFloatType sqrt_three_down=sqrt(down,three);
    const RoundedFloatType sqrt_three_up  =sqrt(up,three);

    const RoundedFloatType zero(0,precision);
    const RoundedFloatType one(1,precision);
    const RoundedFloatType four(4,precision);
    const RoundedFloatType six(6,precision);
    const RoundedFloatType eps=RoundedFloatType::eps(precision);

    RoundedFloatType::set_rounding_mode(upward);
    ARIADNE_TEST_EQUAL(atan(zero),zero);
    ARIADNE_TEST_COMPARE(mul(atan(+one),four),>=,pi_up);
    ARIADNE_TEST_COMPARE(mul(atan(-one),four),>=,-pi_down);
    ARIADNE_TEST_COMPARE(mul(atan(+sqrt_three_up),three),>=,pi_up);
    ARIADNE_TEST_COMPARE(mul(atan(-sqrt_three_down),three),>=,-pi_down);
    ARIADNE_TEST_COMPARE(mul(atan(+1/sqrt_three_down),six),>=,pi_up);
    ARIADNE_TEST_COMPARE(mul(atan(-1/sqrt_three_up),six),>=,-pi_down);

    RoundedFloatType::set_rounding_mode(downward);
    ARIADNE_TEST_EQUAL(atan(zero),zero);
    ARIADNE_TEST_COMPARE(mul(atan(+one),four),<=,pi_down);
    ARIADNE_TEST_COMPARE(mul(atan(-one),four),<=,-pi_up);
    ARIADNE_TEST_COMPARE(mul(atan(+sqrt_three_down),three),<=,pi_down);
    ARIADNE_TEST_COMPARE(mul(atan(-sqrt_three_up),three),<=,-pi_up);
    ARIADNE_TEST_COMPARE(mul(atan(+1/sqrt_three_up),six),<=,pi_down);
    ARIADNE_TEST_COMPARE(mul(atan(-1/sqrt_three_down),six),<=,-pi_up);

    ARIADNE_TEST_COMPARE(sub(mul(atan(one),four),pi_up),<=,4*eps);
    ARIADNE_TEST_COMPARE(sub(mul(atan(sqrt_three_up),three),pi_up),<=,4*eps);
    ARIADNE_TEST_COMPARE(sub(mul(atan(1/sqrt_three_down),six),pi_up),<=,4*eps);

    // Regression test
    ARIADNE_TEST_COMPARE(pi_down/4,<,pi_up/4);

}
