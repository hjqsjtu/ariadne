/***************************************************************************
 *            function/function_patch.cpp
 *
 *  Copyright  2020  Pieter Collins
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

#include "../function/function_patch.hpp"

#include "../algebra/algebra.hpp"

#include "../function/formula.hpp"

namespace Ariadne {

template class FunctionPatchFactory<ValidatedTag>;
template class FunctionPatchFactory<ApproximateTag>;

template class FunctionPatchCreator<FunctionPatchFactory<ValidatedTag>,IntervalDomainType>;
template class FunctionPatchCreator<FunctionPatchFactory<ValidatedTag>,BoxDomainType>;

template class FunctionPatch<ValidatedTag,IntervalDomainType,IntervalDomainType>;
template class FunctionPatch<ValidatedTag,IntervalDomainType,BoxDomainType>;
template class FunctionPatch<ValidatedTag,BoxDomainType,IntervalDomainType>;
template class FunctionPatch<ValidatedTag,BoxDomainType,BoxDomainType>;

template class FunctionPatch<ApproximateTag,IntervalDomainType,IntervalDomainType>;
template class FunctionPatch<ApproximateTag,IntervalDomainType,BoxDomainType>;
template class FunctionPatch<ApproximateTag,BoxDomainType,IntervalDomainType>;
template class FunctionPatch<ApproximateTag,BoxDomainType,BoxDomainType>;

} // namespace Ariadne
