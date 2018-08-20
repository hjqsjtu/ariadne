/***************************************************************************
 *            map.hpp
 *
 *  Copyright  2004-8  Alberto Casagrande, Pieter Collins
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

/*! \file map.hpp
 *  \brief Main continuous dynamics system class.
 */

#ifndef ARIADNE_MAP_HPP
#define ARIADNE_MAP_HPP

#include <memory>

#include "../function/function.hpp"
#include "../geometry/set_interface.hpp"
#include "../geometry/grid.hpp"

namespace Ariadne {

/*! \brief An iterated function system in Euclidean space.
 */
class IteratedMap
{
  public:
    //! \brief The type used to represent time.
    typedef Integer TimeType;
    //! \brief The type used to represent real numbers.
    typedef Real RealType ;
    //! \brief The type used to describe the state space.
    typedef EuclideanSpace StateSpaceType;
  public:
    IteratedMap(const EffectiveVectorFunction& f) : _function(f) { }
    virtual IteratedMap* clone() const { return new IteratedMap(*this); }
    virtual ~IteratedMap() = default;
    const EffectiveVectorFunction& function() const { return _function; }
    Grid grid() const { return Grid(_function.argument_size()); }
  private:
    EffectiveVectorFunction _function;
};

inline OutputStream& operator<<(OutputStream& os, const IteratedMap& vf) {
    return os << "IteratedMap( " << vf.function() << " )";
}


} // namespace Ariadne

#endif // ARIADNE_MAP_HPP
