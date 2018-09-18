/***************************************************************************
 *            measurable_function.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file measurable_function.hpp
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_MEASURABLE_FUNCTION_HPP
#define ARIADNE_MEASURABLE_FUNCTION_HPP

#include "../numeric/number.hpp"

#include "domain.hpp"
#include "function_model.hpp"

namespace Ariadne {

template<class P> class PositiveUpperNumber;
template<> class PositiveUpperNumber<ValidatedTag> : public Positive<Number<ValidatedUpperTag>> { };
using ValidatedPositiveUpperNumber = Positive<Number<ValidatedUpperTag>>;
using EffectivePositiveUpperNumber = Positive<Number<EffectiveUpperTag>>;

template<class P> class MeasurableFunction;
using EffectiveMeasurableFunction = MeasurableFunction<EffectiveTag>;
using ValidatedMeasurableFunction = MeasurableFunction<ValidatedTag>;

template<class P> using ContinuousFunction = ScalarUnivariateFunction<P>;
using EffectiveContinuousFunction = ContinuousFunction<EffectiveTag>;
using ValidatedContinuousFunction = ContinuousFunction<ValidatedTag>;

//using ValidatedScalarUnivariateFunctionModel

template<Nat p> class Norm {
    template<class P> PositiveUpperNumber<P> operator() (ContinuousFunction<P> const&, IntervalDomainType dom);
};

template<> class Norm<2u> {
    ValidatedPositiveUpperNumber operator() (ValidatedContinuousFunction const& f, IntervalDomainType dom) {
        ValidatedScalarUnivariateFunctionModel mf;
        auto imf=antiderivative(mf*mf);
        sqr_nrm = imf(dom.upper_bound())-imf(dom.lower_bound());
        return sqrt(sqr_nrm);
    }
    
    EffectivePositiveUpperNumber operator() (EffectiveContinuousFunction const& f, IntervalDomainType dom);
};

class LowerMeasurableSet {
    
    EffectivePositiveLowerNumber measure() const;
    
    friend LowerMeasurableSet intersection(LowerMeasurableSet, OpenSet)
    
};



template<> class MeasurableFunction<EffectiveTag> {
    LowerMeasurableSet preimage(OpenSet);
};


} // namespace Ariadne

#endif
