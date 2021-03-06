/***************************************************************************
 *            differential.decl.hpp
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

/*! \file differential.decl.hpp
 *  \brief Class declarations and type synonyms for differential classes
 */

#ifndef ARIADNE_DIFFERENTIAL_DECL_HPP
#define ARIADNE_DIFFERENTIAL_DECL_HPP

#include "../numeric/number.decl.hpp"
#include "../numeric/float.decl.hpp"
#include "../geometry/interval.decl.hpp"

namespace Ariadne {

template<class X> class Vector;

template<class X> class Differential;
template<class X> class UnivariateDifferential;
template<class X> class NonAssignableDifferential;

template<class X> class Series;


//@{
//! \relates Differential
//! \name Type synonyms
template<class X> using DifferentialVector = Vector<Differential<X>>; //!< .

using FloatDPDifferential = Differential<FloatDP>; //!< DEPRECATED.
using FloatDPApproximationDifferential = Differential<FloatDPApproximation>; //!< .
using FloatDPBoundsDifferential = Differential<FloatDPBounds>; //!< .

using FloatMPDifferential = Differential<FloatMP>; //!< DEPRECATED.
using FloatMPApproximationDifferential = Differential<FloatMPApproximation>; //!< .
using FloatMPBoundsDifferential = Differential<FloatMPBounds>; //!< .

using FloatDPApproximationDifferentialVector = DifferentialVector<FloatDPApproximation>; //!< .
using FloatDPBoundsDifferentialVector = DifferentialVector<FloatDPBounds>; //!< .
using FloatMPApproximationDifferentialVector = DifferentialVector<FloatMPApproximation>; //!< .
using FloatMPBoundsDifferentialVector = DifferentialVector<FloatMPBounds>; //!< .

using FloatDPUpperIntervalDifferential = Differential<FloatDPUpperInterval>; //!< .
using FloatMPUpperIntervalDifferential = Differential<FloatMPUpperInterval>; //!< .
//@}

} //namespace Ariadne

#endif // ARIADNE_DIFFERENTIAL_DECL_HPP
