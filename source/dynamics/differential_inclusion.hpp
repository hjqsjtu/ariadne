/***************************************************************************
 *            differential_inclusion.hpp
 *
 *  Copyright  2008-17  Pieter Collins, Sanja Zivanovic
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

/*! \file differential_inclusion.hpp
 *  \brief Methods for computing solutions of differential inclusions.
 */

#ifndef ARIADNE_DIFFERENTIAL_INCLUSION_HPP
#define ARIADNE_DIFFERENTIAL_INCLUSION_HPP

#include "utility/typedefs.hpp"
#include "utility/attribute.hpp"
#include "algebra/sweeper.hpp"
#include "algebra/algebra.hpp"
#include "function/domain.hpp"
#include "function/function_model.hpp"
#include "function/formula.hpp"

namespace Ariadne {

class Real;

struct StepSize : public Attribute<Float64> { };
struct NumberOfStepsBetweenSimplifications : public Attribute<Nat> { };
struct NumberOfVariablesToKeep : public Attribute<Nat> { };

Generator<StepSize> step_size;
Generator<NumberOfStepsBetweenSimplifications> number_of_steps_between_simplifications;
Generator<NumberOfVariablesToKeep> number_of_variables_to_keep;

using ValidatedScalarFunctionModelType = ValidatedScalarFunctionModel64;
using ValidatedVectorFunctionModelType = ValidatedVectorFunctionModel64;

using ThresholdSweeper64 = ThresholdSweeper<Float64>;
using Sweeper64 = Sweeper<Float64>;

using ApproximateTimeStepType = PositiveFloat64Approximation;
using ExactTimeStepType = PositiveFloat64Value;

Vector<Float64Value> const& cast_exact(Vector<Float64Error> const& v) {
    return reinterpret_cast<Vector<Float64Value>const&>(v); }

Bool refines(Vector<UpperIntervalType> const& v1, UpperBoxType const& bx2) {
    return refines(v1,static_cast<Vector<UpperIntervalType>const&>(bx2)); }

Box<Interval<Float64Value>> over_approximation(Box<Interval<Real>> const&);

template<class F1, class F2, class F3, class... FS> decltype(auto) combine(F1 const& f1, F2 const& f2, F3 const& f3, FS const& ... fs) {
    return combine(combine(f1,f2),f3,fs...); }
template<class F1, class F2, class F3, class... FS> decltype(auto) join(F1 const& f1, F2 const& f2, F3 const& f3, FS const& ... fs) {
    return join(join(f1,f2),f3,fs...); }

BoxDomainType error_domain(SizeType n, Float64Error e) {
    Float64Value const& v=cast_exact(e);
    return BoxDomainType(n,IntervalDomainType(-v,+v));
}

//! \brief The function dexp(x)=(exp(x)-1)/x.
//! Note that the function is positive and monotone increasing.
template<class PR> PositiveFloatUpperBound<PR> dexp(FloatUpperBound<PR> const& x) {
    if(x.raw()>=0) { return cast_positive(exp(x)-1)/cast_positive(cast_exact(x)); }
    else { return cast_positive(cast_exact(1-exp(x)))/cast_positive(-x); }
}


class Reconditioner {
  public:
    virtual Void simplify(ValidatedVectorFunctionModelType& phi) const = 0;
    virtual ValidatedVectorFunctionModelType expand_errors(ValidatedVectorFunctionModelType Phi) const = 0;
};

class InclusionIntegratorInterface {
  public:
    virtual List<ValidatedVectorFunctionModelType> flow(EffectiveVectorFunction f, BoxDomainType V, BoxDomainType X0, Real T) const = 0;

    virtual Pair<ExactTimeStepType,UpperBoxType> flow_bounds(ValidatedVectorFunction f, UpperBoxType V, ExactBoxType D, ApproximateTimeStepType hsug) const = 0;
    virtual ValidatedVectorFunctionModelType compute_step(EffectiveVectorFunction f, BoxDomainType V, BoxDomainType D, ExactTimeStepType h, UpperBoxType B) const = 0;
};

class InclusionIntegratorBase : public virtual InclusionIntegratorInterface, public Loggable {
  protected:
    SharedPointer<Reconditioner> _reconditioner;
    Sweeper64 _sweeper;
    Float64 _step_size;
    Nat _number_of_steps_between_simplifications;
    Nat _number_of_variables_to_keep;
    InclusionIntegratorBase(Sweeper64 sweeper, StepSize step_size);
    template<class... AS> InclusionIntegratorBase(Sweeper64 sweeper, StepSize step_size, AS... attributes)
        : InclusionIntegratorBase(sweeper,step_size) { this->set(attributes...); }
  public:
    InclusionIntegratorBase& set(NumberOfStepsBetweenSimplifications n) { _number_of_steps_between_simplifications=n; return *this; }
    InclusionIntegratorBase& set(NumberOfVariablesToKeep n) { _number_of_variables_to_keep=n; return *this; }
    template<class A, class... AS> InclusionIntegratorBase& set(A a, AS... as) { this->set(a); this->set(as...); return *this; }

    ValidatedVectorFunctionModelType expand_errors(ValidatedVectorFunctionModelType Phi) const;
    ValidatedVectorFunctionModelType simplify(ValidatedVectorFunctionModelType Phi) const;
    List<ValidatedVectorFunctionModelType> flow(EffectiveVectorFunction f, BoxDomainType V, BoxDomainType X0, Real T) const;

    Pair<ExactTimeStepType,UpperBoxType> flow_bounds(ValidatedVectorFunction f, UpperBoxType V, ExactBoxType D, ApproximateTimeStepType hsug) const;

    virtual ValidatedVectorFunctionModelType compute_step(EffectiveVectorFunction f, BoxDomainType V, BoxDomainType D, ExactTimeStepType h, UpperBoxType B) const = 0;
};

class InclusionIntegrator3rdOrder : public InclusionIntegratorBase {
  public:
    template<class... AS> InclusionIntegrator3rdOrder(Sweeper64 sweeper, StepSize step_size, AS... attributes)
        : InclusionIntegratorBase(sweeper,step_size) { this->set(attributes...); }
    Tuple<Float64Error,Float64Error,Float64Error,Float64UpperBound> compute_norms(EffectiveVectorFunction const& f, UpperBoxType const& B) const;
    virtual ValidatedVectorFunctionModelType compute_step(EffectiveVectorFunction f, BoxDomainType V, BoxDomainType D, ExactTimeStepType h, UpperBoxType B) const override;
};



class InclusionIntegrator2ndOrder : public InclusionIntegratorBase {
  public:
    template<class... AS> InclusionIntegrator2ndOrder(Sweeper64 sweeper, StepSize step_size, AS... attributes)
        : InclusionIntegratorBase(sweeper,step_size) { this->set(attributes...); }
    Tuple<Float64Error,Float64Error,Float64UpperBound> compute_norms(EffectiveVectorFunction const& f, UpperBoxType const& B) const;
    virtual ValidatedVectorFunctionModelType compute_step(EffectiveVectorFunction f, BoxDomainType V, BoxDomainType D, PositiveFloat64Value h, UpperBoxType B) const override;
};



class LohnerReconditioner : public Reconditioner, public Loggable {
    Sweeper64 _sweeper;
    SizeType _number_of_variables_to_keep;
  public:
    LohnerReconditioner(Sweeper64 sweeper, NumberOfVariablesToKeep number_of_variables_to_keep);
    virtual ValidatedVectorFunctionModelType expand_errors(ValidatedVectorFunctionModelType Phi) const override;
    virtual Void simplify(ValidatedVectorFunctionModelType& phi) const override;
};

} // namespace Ariadne;

#endif // ARIADNE_DIFFERENTIAL_INCLUSION_HPP