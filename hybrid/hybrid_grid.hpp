/***************************************************************************
 *            hybrid_grid.h
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

/*! \file hybrid_grid.h
 *  \brief Grids in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_GRID_H
#define ARIADNE_HYBRID_GRID_H

#include <iostream>

#include "expression/space.h"
#include "hybrid/hybrid_space.h"
#include "hybrid/hybrid_scaling.h"

#include "geometry/grid.h"

namespace Ariadne {

typedef OutputStream OutputStream;
typedef Bool Bool;
typedef Vector<Float64Value> ExactFloatVector;
class HybridGridTreeSet;

//! \ingroup HybridModule
//! \brief A grid in a hybrid space
class HybridGrid
{
    // NOTE: The use of the system is to allow the "Grid" of a compositional
    // hybrid automaton to be computed on-the-fly since it might not be
    // feasible to compute the reachable states.
    HybridSpace _space;
    HybridScaling _scaling;
  public:
    //! Test whether the grid has location \a q.
    Bool has_location(const DiscreteLocation& q) const {
        return this->_space.has_location(q); }
    //! The grid in location \a loc.
    Grid operator[](const DiscreteLocation& loc) const;
    //! The underlying hybrid space.
    const HybridSpace& space() const { return this->_space; }
    //! The underlying real space in location \a loc.
    RealSpace space(const DiscreteLocation& loc) const { return this->_space[loc]; }
    //! The variable scalings used.
    HybridScalingInterface& scalings() { return this->_scaling; }
    //! The grid in location \a loc.
    Grid grid(const DiscreteLocation& loc) const { return this->operator[](loc); }
  public:
    //! The grid corresponding to unit resolution on each dimension.
    HybridGrid(const HybridSpace& hsp) : _space(hsp), _scaling(SimpleHybridScaling()) { }
    //! The grid corresponding to scaling each variable in each location according to the policy \a hsc.
    HybridGrid(const HybridSpace& hsp, const HybridScaling& hsc)
        : _space(hsp), _scaling(hsc) { }
    HybridGrid* clone() const { return new HybridGrid(*this); }
    //!
    friend OutputStream& operator<<(OutputStream& os, const HybridGrid& hgrid);
};

inline Grid HybridGrid::operator[](const DiscreteLocation& loc) const
{
    RealSpace continuous_space = this->_space[loc];
    Vector<RawFloat64> lengths(continuous_space.size());
    for(Nat i=0; i!=continuous_space.size(); ++i) {
        lengths[i] = (this->_scaling.scaling(loc,continuous_space.variable(i))).raw();
    }
    return Grid(lengths);
}

inline OutputStream& operator<<(OutputStream& os, const HybridGrid& hgrid) {
    return os << "HybridGrid( space="<<hgrid._space << ", scalings=" << hgrid._scaling << " )";
}

}

#endif // ARIADNE_HYBRID_GRID_H
