/***************************************************************************
 *            hybrid_enclosure.h
 *
 *  Copyright  2009-10  Pieter Collins
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

/*! \file hybrid_enclosure.h
 *  \brief Enclosure sets for hybrid systems
 */

#ifndef ARIADNE_HYBRID_ENCLOSURE_H
#define ARIADNE_HYBRID_ENCLOSURE_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>
#include "discrete_location.h"
#include "discrete_event.h"
#include "taylor_set.h"
#include "graphics_interface.h"
#include "container.h"
#include "logging.h"

namespace Ariadne {

class Interval;
template<class X> class Vector;
template<class X> class LinearProgram;
class ScalarFunction;
class VectorFunction;
class ScalarTaylorFunction;
class VectorTaylorFunction;
class TaylorConstrainedImageSet;
class HybridEnclosure;
class Box;
class Grid;
class GridTreeSet;
class AffineSet;
class DiscreteEvent;
class Figure;
class CanvasInterface;

template<class ES> class ListSet;
template<class ES> class HybridListSet;
template<> class ListSet<HybridEnclosure>;

template<class BS> class HybridBasicSet;
typedef HybridBasicSet<Box> HybridBox;
class HybridGridTreeSet;

class HybridEnclosure
    : public DrawableInterface
{
    friend class ConstraintHybridEvolver;
  public:
    typedef TaylorConstrainedImageSet ContinuousStateSetType;
    private:
    DiscreteLocation _location;
    List<DiscreteEvent> _events;
    List< List<DiscreteEvent> > _constraint_events;
    ContinuousStateSetType _set;
    ScalarTaylorFunction _time;
  public:
    HybridEnclosure(const DiscreteLocation&, const Box&);
    HybridEnclosure(const std::pair<DiscreteLocation,ContinuousStateSetType>&);
    HybridEnclosure(const DiscreteLocation&, const ContinuousStateSetType&);
    ~HybridEnclosure();
    HybridEnclosure* clone() const;

    const DiscreteLocation& location() const;
    const ContinuousStateSetType& continuous_state_set() const;

    void apply_reset(DiscreteEvent, DiscreteLocation, VectorFunction);
    void apply_flow(VectorFunction, Interval);
    void apply_flow(VectorTaylorFunction, Float);
    void apply_flow(VectorTaylorFunction, Interval);

    void set_maximum_time(DiscreteEvent,Float);
    void set_step_time(Float);
    void set_time(DiscreteEvent,Float);
    void set_dwell_time(DiscreteEvent,Float);
    void set_dwell_time(DiscreteEvent,ScalarFunction);
    void set_dwell_time(DiscreteEvent,ScalarTaylorFunction);

    void new_invariant(DiscreteEvent,ScalarFunction,ScalarFunction);
    void new_activation(DiscreteEvent,ScalarFunction,ScalarFunction);
    void new_guard(DiscreteEvent,ScalarFunction,ScalarFunction);
    void new_time_bound(DiscreteEvent,ScalarFunction);

    uint dimension() const;
    tribool empty() const;

    HybridBox bounding_box() const;
    tribool disjoint(const HybridBox& bx) const;
    void adjoin_outer_approximation_to(HybridGridTreeSet& paving, int depth) const;
    Pair<HybridEnclosure,HybridEnclosure> split(uint dim) const;

    void draw(CanvasInterface&) const;
    std::ostream& write(std::ostream&) const;
};

inline std::ostream& operator<<(std::ostream& os, const HybridEnclosure& s) { return s.write(os); }

}

#include "hybrid_set.h"

namespace Ariadne {

template<>
class ListSet<HybridEnclosure>
    : public HybridListSet<HybridEnclosure::ContinuousStateSetType>
{
  public:
    ListSet() { }
    ListSet(const HybridEnclosure& hes) { this->adjoin(hes); }
    using HybridListSet<HybridEnclosure::ContinuousStateSetType>::adjoin;
    void adjoin(const HybridEnclosure& hes) { this->adjoin(hes.location(),hes.continuous_state_set()); }
};

inline std::ostream& operator<<(std::ostream& os, const ListSet<HybridEnclosure>& hls) {
    return os << static_cast<const HybridListSet<HybridEnclosure::ContinuousStateSetType>&>(hls);
}

} // namespace Ariadne

#endif // ARIADNE_HYBRID_ENCLOSURE_H
