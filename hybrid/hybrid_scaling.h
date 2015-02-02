/***************************************************************************
 *            hybrid_scaling.h
 *
 *  Copyright 2008-12  Pieter Collins
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

/*! \file hybrid_scaling.h
 *  \brief Scalings for variables in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SCALING_H
#define ARIADNE_HYBRID_SCALING_H

#include <iostream>
#include <map>

#include "utility/container.h"
#include "utility/stlio.h"

#include "expression/variables.h"

typedef Void Void;
typedef OutputStream OutputStream;

namespace Ariadne {


//! \ingroup HybridModule
//! \brief A class which defines the state space grid to use in location \a loc given the continuous state variables \a spc.
class HybridScalingInterface
{
  public:
    //!
    virtual HybridScalingInterface* clone() const = 0;
    virtual ExactFloat scaling(const DiscreteLocation& loc, const RealVariable& var) const = 0;
    virtual Void write(OutputStream& os) const = 0;
};
inline OutputStream& operator<<(OutputStream& os, const HybridScalingInterface& hsc) { hsc.write(os); return os; }


//! \ingroup HybridModule
//! \brief A class which defines the state space grid to use in location \a loc given the continuous state variables \a spc.
class HybridScaling
{
    std::shared_ptr<HybridScalingInterface> _ptr;
  public:
    HybridScaling(const HybridScalingInterface& ref) : _ptr(ref.clone()) { }
    operator const HybridScalingInterface& () const { return *this->_ptr; }
    operator HybridScalingInterface& () { return *this->_ptr; }
    ExactFloat scaling(const DiscreteLocation& loc, const RealVariable& var) const { return this->_ptr->scaling(loc,var); }
  public:
    HybridScaling(const Map<Identifier,ExactFloat>& scalings);
    HybridScaling(const InitializerList<Pair<RealVariable,Float>>& scalings);
};

class SimpleHybridScaling
    : public HybridScalingInterface
{
    Map<Identifier,ExactFloat> _scalings;
  public:
    SimpleHybridScaling() : _scalings() { }
    SimpleHybridScaling(const Map<Identifier,ExactFloat>& scalings) : _scalings(scalings) { }
    SimpleHybridScaling(const InitializerList<Pair<RealVariable,Float>>& scalings);
    Void set_scaling(const RealVariable& var, ExactFloat res) { ARIADNE_ASSERT(res>0.0); _scalings[var.name()]=res; }
    virtual SimpleHybridScaling* clone() const { return new SimpleHybridScaling(*this); }
    virtual ExactFloat scaling(const DiscreteLocation& loc, const RealVariable& var) const {
        return (this->_scalings.has_key(var.name())) ? this->_scalings[var.name()] : ExactFloat(1.0); }
    virtual Void write(OutputStream& os) const { os << "HybridScaling( " << this->_scalings << " )"; }
};

inline Pair<RealVariable,Float> operator|(const RealVariable& var, Float scal) {
    return Pair<RealVariable,Float>(var,scal); }

inline HybridScaling::HybridScaling(const Map<Identifier,ExactFloat>& scalings)
    : _ptr(new SimpleHybridScaling(scalings)) { }

inline HybridScaling::HybridScaling(const InitializerList<Pair<RealVariable,Float>>& scalings)
    : _ptr(new SimpleHybridScaling(scalings)) { }

inline SimpleHybridScaling::SimpleHybridScaling(const InitializerList<Pair<RealVariable,Float>>& scalings) {
    for(auto iter=scalings.begin(); iter!=scalings.end(); ++iter) {
        this->_scalings.insert(iter->first.name(),ExactFloat(iter->second));
    }
}

}

#endif // ARIADNE_HYBRID_SCALING_H

