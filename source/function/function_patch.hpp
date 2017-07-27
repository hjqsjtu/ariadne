/***************************************************************************
 *            function_patch.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file function_patch.hpp
 *  \brief Over-approximations of functions on box domains.
 */

#ifndef ARIADNE_FUNCTION_PATCH_HPP
#define ARIADNE_FUNCTION_PATCH_HPP

#include "function/function_interface.hpp"
#include "function/function.hpp"

namespace Ariadne {

template<class Y> class Upper;

template<> class Upper<ValidatedNumber> : ValidatedUpperNumber { };

/*! \ingroup FunctionModelSubModule
 *  \brief A FunctionPatch is a function defined on an interval, box, or other compact domain.
 *   It supports the supremum \a norm() method.
 */
template<class P, class D, class C> class FunctionPatch
    : public Function<P,D,C>
{
    //! \brief Default constructor.
    FunctionPatch(Function<P,D,C>);

    Positive<Upper<Number<P>>> norm() const;
};

} // namespace Ariadne

#endif // ARIADNE_FUNCTION_PATCH_HPP
