/***************************************************************************
 *            differential.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "numeric/numeric.hpp"
#include "geometry/interval.hpp"

#include "algebra/differential.hpp"

#include "operations.hpp"

#include "differential.tpl.hpp"

namespace Ariadne {

template class UnivariateDifferential<Float64>;
template class UnivariateDifferential<Float64Approximation>;
template class UnivariateDifferential<Float64Bounds>;
template class UnivariateDifferential<UpperIntervalType>;

template class UnivariateDifferential<FloatMPApproximation>;
template class UnivariateDifferential<FloatMPBounds>;


template class Differential<Float64>;
template class Differential<Float64Bounds>;
template class Differential<Float64Approximation>;
template class Differential<UpperIntervalType>;

template class AlgebraOperations<Differential<Float64>>;
template class AlgebraOperations<Differential<Float64Approximation>>;
template class AlgebraOperations<Differential<Float64Bounds>>;
template class GradedAlgebraOperations<Differential<Float64>>;
template class GradedAlgebraOperations<Differential<Float64Approximation>>;
template class GradedAlgebraOperations<Differential<Float64Bounds>>;

template class Vector<Differential<Float64>>;
template class Vector<Differential<Float64Bounds>>;
template class Vector<Differential<Float64Approximation>>;
//template class Vector<Differential<UpperIntervalType>>;

template class Differential<FloatMPBounds>;
template class Differential<FloatMPApproximation>;
template class AlgebraOperations<Differential<FloatMPApproximation>>;
template class AlgebraOperations<Differential<FloatMPBounds>>;
template class GradedAlgebraOperations<Differential<FloatMPApproximation>>;
template class GradedAlgebraOperations<Differential<FloatMPBounds>>;

template class Vector<Differential<FloatMPBounds>>;
template class Vector<Differential<FloatMPApproximation>>;

}