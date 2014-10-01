/***************************************************************************
 *            constraint.h
 *
 *  Copyright 2009  Pieter Collins
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

/*! \file constraint.h
 *  \brief Linear and nonlinear constraints.
 */

#ifndef ARIADNE_CONSTRAINT_H
#define ARIADNE_CONSTRAINT_H

#include "numeric.h"
#include "function.h"
#include "metaprogramming.h"

namespace Ariadne {


template<class F, class R>
class Constraint {
  public:
    typedef F FunctionType;
    typedef R BoundType;
    typedef R UpperBoundType;
    typedef decltype(-std::declval<R>()) LowerBoundType;
//    typedef typename IntervalOfType<Real>::Type IntervalBoundsType;
  public:
    Constraint(LowerBoundType const& l, FunctionType const& f, UpperBoundType const& u)
        : _function(f), _lower_bound(l), _upper_bound(u) { ARIADNE_ASSERT(l<=u); }

    Constraint(FunctionType const& f, BoundType const& x)
        : _function(f), _lower_bound(x), _upper_bound(x) { }

    template< class FF, class RR, class=typename EnableIf<IsConvertible<FF,F>>::Type, class=typename EnableIf<IsConvertible<RR,R>>::Type >
    Constraint(const Constraint<FF,RR>& c)
        : _function(static_cast<F>(c.function())), _lower_bound(c.lower_bound()), _upper_bound(c.upper_bound()) { }

    void set_function(const FunctionType& f) { this->_function = f; }
    FunctionType& function() { return this->_function; }
    FunctionType const& function() const { return this->_function; }
    Nat argument_size() const { return this->_function.argument_size(); }
    LowerBoundType const& lower_bound() const { return this->_lower_bound; }
    UpperBoundType const& upper_bound() const { return this->_upper_bound; }
    const Interval bounds() const { return Interval(this->_lower_bound,this->_upper_bound); }
  private:
    F _function;
    R _lower_bound;
    R _upper_bound;
};

typedef Constraint<RealScalarFunction,Real> RealConstraint;
typedef Constraint<EffectiveScalarFunction,EffectiveNumberType> EffectiveConstraint;
typedef Constraint<ValidatedScalarFunction,ValidatedNumberType> ValidatedConstraint;
typedef Constraint<ValidatedScalarFunction,ExactNumberType> ValidatedExactConstraint;

inline EffectiveConstraint operator<=(const EffectiveNumberType& c, const EffectiveScalarFunction& f) {
    return EffectiveConstraint(c,f,infinity);
}

inline EffectiveConstraint operator>=(const EffectiveNumberType& c, const EffectiveScalarFunction& f) {
    return EffectiveConstraint(-infinity,f,c);
}

inline EffectiveConstraint operator<=(const EffectiveScalarFunction& f, const EffectiveNumberType& c) {
    return EffectiveConstraint(-infinity,f,c);
}

inline EffectiveConstraint operator>=(const EffectiveScalarFunction& f, const EffectiveNumberType& c) {
    return EffectiveConstraint(c,f,infinity);
}

inline EffectiveConstraint operator==(const EffectiveScalarFunction& f, const EffectiveNumberType& c) {
    return EffectiveConstraint(f,c);
}

inline EffectiveConstraint operator<=(const EffectiveScalarFunction& f, double c) {
    return EffectiveConstraint(-infinity,f,Real(c));
}

inline EffectiveConstraint operator>=(const EffectiveScalarFunction& f, double c) {
    return EffectiveConstraint(Real(c),f,infinity);
}

inline EffectiveConstraint operator==(const EffectiveScalarFunction& f, double c) {
    return EffectiveConstraint(f,Real(c));
}

inline EffectiveConstraint operator<=(const EffectiveConstraint& nc, const EffectiveNumberType& c) {
    ARIADNE_ASSERT(static_cast<ApproximateFloatType>(nc.upper_bound()).get_d()==inf);
    return EffectiveConstraint(nc.lower_bound(),nc.function(),c);
}


inline ValidatedExactConstraint operator<=(const ExactNumberType& c, const ValidatedScalarFunction& f) {
    return ValidatedExactConstraint(c,f,+infty);
}

inline ValidatedExactConstraint operator<=(const ValidatedScalarFunction& f, const ExactNumberType& c) {
    return ValidatedExactConstraint(-infty,f,c);
}

inline ValidatedExactConstraint operator>=(const ValidatedScalarFunction& f, const ExactNumberType& c) {
    return ValidatedExactConstraint(c,f,+infty);
}

inline ValidatedExactConstraint operator==(const ValidatedScalarFunction& f, const ExactNumberType& c) {
    return ValidatedExactConstraint(c,f,c);
}

inline ValidatedExactConstraint operator<=(const ValidatedExactConstraint& nc, const ExactNumberType& c) {
    ARIADNE_ASSERT(nc.upper_bound()==infty);
    return ValidatedExactConstraint(nc.lower_bound(),nc.function(),c);
}

/*
inline ValidatedConstraint operator<=(const ValidatedScalarFunction& f, double c) {
    return f<=ExactNumberType(c);
}

inline ValidatedConstraint operator>=(const ValidatedScalarFunction& f, double c) {
    return f>=ExactNumberType(c);
}

inline ValidatedConstraint operator==(const ValidatedScalarFunction& f, double c) {
    return f==ExactNumberType(c);
}
*/

inline ValidatedExactConstraint operator<=(const ValidatedScalarFunction& f1, const ValidatedScalarFunction& f2) {
    return (f1-f2) <= ExactNumberType(0);
}

inline ValidatedExactConstraint operator>=(const ValidatedScalarFunction& f1, const ValidatedScalarFunction& f2) {
    return (f1-f2) >= ExactNumberType(0);
}


inline ValidatedConstraint operator<=(const ValidatedNumberType& c, const ValidatedScalarFunction& f) {
    return ValidatedConstraint(c,f,+infty);
}

inline ValidatedConstraint operator<=(const ValidatedConstraint& nc, const ValidatedNumberType& c) {
    ARIADNE_ASSERT(nc.upper_bound()==infty);
    return ValidatedConstraint(nc.lower_bound(),nc.function(),c);
}


template<class X, class R> std::ostream& operator<<(std::ostream& os, const Constraint<X,R>& c) {
    return os << c.lower_bound() << "<=" << c.function() << "<=" << c.upper_bound();
}


} //namespace Ariadne

#endif
