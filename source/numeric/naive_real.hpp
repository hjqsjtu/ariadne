/***************************************************************************
 *            numeric/naive_real.hpp
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

/*! \file numeric/naive_real.hpp
 *  \brief "Naive" real numbers
 */



#ifndef ARIADNE_NAIVE_REAL_HPP
#define ARIADNE_NAIVE_REAL_HPP

#include <functional>

#include "../utility/typedefs.hpp"
#include "../utility/pointer.hpp"

#include "../numeric/logical.decl.hpp"
#include "../numeric/number.decl.hpp"
#include "../numeric/float.decl.hpp"

#include "paradigm.hpp"
#include "../numeric/arithmetic.hpp"
#include "../numeric/sequence.hpp"

namespace Ariadne {

class Real;
class LowerReal;
class UpperReal;
class NaiveReal;
class PositiveReal;
class PositiveLowerReal;
class PositiveUpperReal;
class PositiveNaiveReal;

class ValidatedReal;
class ValidatedLowerReal;
class ValidatedUpperReal;
class ApproximateReal;

class RealInterface;

//! \ingroup NumericModule
//! \brief %Real number type defined as limits of convergent sequences of (dyadic) rationals, without bounds on the convergence rate.
//! It is possible to compute arbitrarily accurate (dyadic) rational approximations, but the convergence rate is not known,
//! so at no point in the the computation can anything concrete be deduced about the value of the number.
//! \details In principle useless for rigorous computation, but quickly-computed approximations to naive_real numbers may be useful for preconditioning rigorous algorithms.
//! \sa Real, LowerReal, UpperReal
class NaiveReal
    : public DeclareRealOperations<NaiveReal,PositiveNaiveReal>
    , public DeclareAnalyticFieldOperations<NaiveReal>
    , public DeclareLatticeOperations<NaiveReal,PositiveNaiveReal>
    , public DeclareComparisonOperations<NaiveReal,ApproximateKleenean>
    , public DefineFieldOperators<NaiveReal>
{
  private: public:
    using Interface = RealInterface;
    using InterfaceType = RealInterface;
  private: public:
    SharedPointer<Interface> _ptr;
  public:
    typedef ApproximateTag Paradigm;
    typedef NaiveReal NumericType;
  public:
    //@{
    //! \name Constructors
    NaiveReal(); //!< Default constructor yields the integer \c 0 as a naive_real number.
    explicit NaiveReal(SharedPointer<NaiveReal::InterfaceType>); //!< Construct from any class implementing the naive_real number interface.
    //explicit NaiveReal(ConvergentSequence<DyadicBounds> const&); //!< Construct from a sequence of dyadic bounds converging to a singleton intersection.
    //@}

    //@{
    //! \name Conversion operators
    NaiveReal(Dbl d); //!< Construction from a builtin double-precision approximation.
    NaiveReal(Integer const& n); //!< Construct from an integer.
    NaiveReal(Dyadic const& d); //!< Construct from a dyadic number \a d=p/2<sup>q</sup> for integers \a p, \a q.
    NaiveReal(Decimal const& d); //!< Construct from a decimal number, given by its decimal expansion.
    NaiveReal(Rational const& q); //!< Construct from a rational number.

    NaiveReal(Real const& r); //!< Construct from an (effective) naive_real number.
    NaiveReal(LowerReal const& lr); //!< Construct from a lower naive_real number.
    NaiveReal(UpperReal const& ur); //!< Construct from an upper naive_real number.

//    operator Number<Tag>() const; //!< Convert to an effective number for computations.
    //@}

    //@{
    //! \name Computation of approximations to the number.
/*
    //! Compute an approximation with no guarantees on the error. The error converges to \a 0 as \a eff approaches infinity.
    //! The time taken should be roughly polynomial in \a eff.
    ApproximateReal compute(Effort eff) const;
*/
    //! Compute a concrete approximation using double-precision.
    FloatDPApproximation get(DoublePrecision pr) const;
    //! Compute a concrete approximation using the given precision.
    FloatMPApproximation get(MultiplePrecision pr) const;
    //@}

    //@{
    //! \name Standard arithmetic operators
    friend NaiveReal operator+(NaiveReal const& r); //!< Unary plus.
    friend NaiveReal operator-(NaiveReal const& r); //!< Unary minus.
    friend NaiveReal operator+(NaiveReal const& r1, NaiveReal const& r2); //!< Plus.
    friend NaiveReal operator-(NaiveReal const& r1, NaiveReal const& r2); //!< Minus.
    friend NaiveReal operator*(NaiveReal const& r1, NaiveReal const& r2); //!< Times.
    friend NaiveReal operator/(NaiveReal const& r1, NaiveReal const& r2); //!< Divides.
    friend NaiveReal& operator+=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace plus.
    friend NaiveReal& operator-=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace minus.
    friend NaiveReal& operator*=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace times.
    friend NaiveReal& operator/=(NaiveReal& r1, NaiveReal const& r2); //!< Inplace divides.
    //@}

    //@{
    //! \name Named arithmetical functions
    friend NaiveReal pos(NaiveReal const& r); //!< Identity \a +r.
    friend NaiveReal neg(NaiveReal const& r); //!< Negative \a -r.
    friend NaiveReal hlf(NaiveReal const& r); //!< Half \a r÷2.
    friend NaiveReal sqr(NaiveReal const& r); //!< Square \a r<sup>2</sup>.
    friend NaiveReal rec(NaiveReal const& r); //!< Reciprocal \a 1/r.
    friend NaiveReal add(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Sum \a r1+r2.
    friend NaiveReal sub(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Difference \a r1-r2.
    friend NaiveReal mul(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Product \a r1×r2.
    friend NaiveReal div(NaiveReal const& r1, NaiveReal const& r2); //!< \brief Quotient \a r1÷r2.
    friend NaiveReal fma(NaiveReal const& r1, NaiveReal const& r2, NaiveReal const& r3); //!< \brief Fused multiply-and-add \a r1×r2+r3.
    friend NaiveReal pow(NaiveReal const& r, Int n); //!< \brief Power \a r<sup>n</sup>.
    //@}

    //@{
    //! \name Algebraic and transcendental functions
    friend NaiveReal sqrt(NaiveReal const& r); //!< The square root of \a r, √\a r. Requires \a r ≥ 0.
    friend NaiveReal exp(NaiveReal const& r); //!< The natural exponent of \a r, \em e<sup>r</sup>.
    friend NaiveReal log(NaiveReal const& r); //!< The natural logarithm of \a r. Requires \a r ≥ 0.
    friend NaiveReal sin(NaiveReal const& r); //!< The sine of \a r.
    friend NaiveReal cos(NaiveReal const& r); //!< The cosine of \a r.
    friend NaiveReal tan(NaiveReal const& r); //!< The tangent of \a r, sin(\a r)/cos(\a r) \f$.
    friend NaiveReal atan(NaiveReal const& r); //!< The arc-tangent of \a r.
    //@}

    //@{
    //! \name Lattice operations
    friend PositiveNaiveReal abs(NaiveReal const&); //!< Absolute value \a |r|.
    friend NaiveReal min(NaiveReal const& r1, NaiveReal const& r2); //!< The mimimum of \a r1 and \a r2.
    friend NaiveReal max(NaiveReal const& r1, NaiveReal const& r2); //!< The maximum of \a r1 and \a r2.
    //@}

    //@{
    //! Operations based on the metric structure.
    friend PositiveNaiveReal dist(NaiveReal const& r1, NaiveReal const& r2);
        //< The distance |\a r <sub>1</sub>-\a r <sub>2</sub>| between \a r<sub>1</sub> and \a r<sub>2</sub>.
    //@}

    //@{
    //! \name Comparison operations and operators.
    friend ApproximateKleenean leq(NaiveReal const& r1, NaiveReal const& r2);
        //!< Returns \c likely if \a r1<r2, \c unlikely if \a r1\>r2, and \c indetermiate if \a r1==r2.
    //@}

    //@{
    //! \name Limit operations.
    friend NaiveReal limit(ConvergentSequence<Dyadic> const& qas);
        //!< Create a naive_real number from an converging sequence of dyadic approximants.
    //@}

    //@{
    //! \name Input/output operations
    friend OutputStream& operator<<(OutputStream& os, NaiveReal const& r); //< Write to an output stream.
    //@}
};

} // namespace Ariadne

#endif
