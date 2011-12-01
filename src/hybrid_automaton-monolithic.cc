/***************************************************************************
 *            hybrid_automaton.cc
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

#include <map>

#include "config.h"

#include "macros.h"
#include "stlio.h"
#include "function.h"
#include "hybrid_time.h"
#include "hybrid_space.h"
#include "hybrid_automaton-monolithic.h"

namespace Ariadne {

typedef uint DimensionType;



MonolithicHybridAutomaton::Mode::Mode(DiscreteLocation q, RealSpace s, RealVectorFunction f)
    : _location(q), _variable_names(s.variable_names()), _dynamic(f)
{
}


MonolithicHybridAutomaton::~MonolithicHybridAutomaton()
{
}

MonolithicHybridAutomaton::MonolithicHybridAutomaton()
{
}

MonolithicHybridAutomaton::MonolithicHybridAutomaton(const String& name)
{
    ARIADNE_DEPRECATED("MonolithicHybridAutomaton(String)","MonolithicHybridAutomaton class instances have no name.");
}





void
MonolithicHybridAutomaton::new_mode(DiscreteLocation location,
                                    RealVectorFunction dynamic)
{
    List<Identifier> names;
    for(uint i=0; i!=dynamic.result_size(); ++i) {
        std::stringstream ss;
        ss << "x" << i;
        names.append(ss.str());
    }
    this->new_mode(location,RealSpace(names),dynamic);
}


void
MonolithicHybridAutomaton::new_mode(DiscreteLocation location,
                                    RealSpace space,
                                    RealVectorFunction dynamic)
{
    if(this->has_mode(location)) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_mode(location,dynamic)",
                      "The hybrid automaton already has a mode with location label "<<location<<"\n");
    }
    if(dynamic.result_size()!=dynamic.argument_size()) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_mode(location,dynamic)",
            "The dynamic has argument size " << dynamic.argument_size()
                << " and result size "<<dynamic.result_size()<<", so does not define a vector field.");
    }
    if(space.size()!=dynamic.result_size()) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_mode(location,space,dynamic)",
                      "The number of variables in the state space "<<space<<" does not match the number of variables defined by the dynamic "<<dynamic<<"\n");
    }
    this->_modes.insert(location,Mode(location,space,dynamic));
}


void
MonolithicHybridAutomaton::new_invariant(DiscreteLocation location,
                                         DiscreteEvent event,
                                         RealScalarFunction invariant,
                                         EventKind kind)
{
    if(!this->has_mode(location)) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_invariant(location,label,invariant,kind)",
                      "The location "<<location<<" of the invariant must be in the automaton.");
    }
    if(this->has_guard(location,event)) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_invariant(location,label,invariant,kind)",
                      "The automaton already has a guard or invariant in location "<<location<<" with event label "<<event<<"\n");
    }
    if(invariant.argument_size()!=this->dimension(location)) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_invariant(location,label,invariant,kind)",
            "The invariant has argument size " << invariant.argument_size()
                << " but the mode has state-space dimension " << this->dimension(location));
    }

    this->_modes.value(location)._invariants.insert(event,Invariant(location,event,invariant,kind));
}



void
MonolithicHybridAutomaton::new_transition(DiscreteLocation source,
                                          DiscreteEvent event,
                                          DiscreteLocation target,
                                          RealVectorFunction reset,
                                          RealScalarFunction guard,
                                          EventKind kind)
{
    if(!this->has_mode(source)) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_transition(...)",
                      "The source location "<<source<<" of transition event "<<event<<" must be in the automaton\n");
    }
    if(!this->has_mode(target)) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_transition(...)",
                      "The target location "<<target<<" of transition event "<<event<<" from "<<source<<" must be in the automaton\n");
    }
    if(this->has_guard(source,event)) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_transition(...)",
                      "The automaton already has an invariant or transition in location "<<source<<" with event label "<<event<<"\n");
    }

    this->_modes.value(source)._transitions.insert(event,Transition(source,event,target,reset,guard,kind));
}







bool
MonolithicHybridAutomaton::has_mode(DiscreteLocation location) const
{
    return this->_modes.has_key(location);
}


bool
MonolithicHybridAutomaton::has_guard(DiscreteLocation source, DiscreteEvent event) const
{
    if(!this->has_mode(source)) { return false; }
    const Mode& mode = this->_modes[source];
    return mode._transitions.has_key(event);
}

bool
MonolithicHybridAutomaton::has_invariant(DiscreteLocation source, DiscreteEvent event) const
{
    if(!this->has_mode(source)) { return false; }
    const Mode& mode = this->_modes[source];
    return mode._invariants.has_key(event);
}

bool
MonolithicHybridAutomaton::has_transition(DiscreteLocation source, DiscreteEvent event) const
{
    if(!this->has_mode(source)) { return false; }
    const Mode& mode = this->_modes[source];
    return mode._transitions.has_key(event);
}

const MonolithicHybridAutomaton::Mode&
MonolithicHybridAutomaton::mode(DiscreteLocation location) const
{
    ARIADNE_ASSERT(this->has_mode(location));
    return this->_modes[location];
}


const MonolithicHybridAutomaton::Invariant&
MonolithicHybridAutomaton::invariant(DiscreteLocation source, DiscreteEvent event) const
{
    ARIADNE_ASSERT(this->has_invariant(source,event));
    return this->_modes[source]._invariants[event];
}

const MonolithicHybridAutomaton::Transition&
MonolithicHybridAutomaton::transition(DiscreteLocation source, DiscreteEvent event) const
{
    ARIADNE_ASSERT(this->has_transition(source,event));
    return this->_modes[source]._transitions[event];
}


uint
MonolithicHybridAutomaton::dimension(DiscreteLocation location) const
{
    return this->mode(location).dimension();
}

RealVectorFunction
MonolithicHybridAutomaton::dynamic_function(DiscreteLocation location) const
{
    return this->mode(location)._dynamic;
}

Set<DiscreteLocation>
MonolithicHybridAutomaton::locations() const
{
    return this->_modes.keys();
}

Set<DiscreteEvent>
MonolithicHybridAutomaton::events(DiscreteLocation location) const
{
    return join(this->mode(location)._invariants.keys(), this->mode(location)._transitions.keys());
    //return this->mode(location)._transitions.keys();
}

EventKind
MonolithicHybridAutomaton::event_kind(DiscreteLocation location, DiscreteEvent event) const
{
    ARIADNE_ASSERT_MSG(this->has_guard(location,event) || this->has_invariant(location,event),"No event "<<event<<" in location "<<location);
    const Mode& mode = this->_modes[location];
    if(mode._invariants.has_key(event)) { return mode._invariants[event]._kind; }
    else { return mode._transitions[event]._kind; }
}

RealScalarFunction
MonolithicHybridAutomaton::invariant_function(DiscreteLocation location, DiscreteEvent event) const
{
    ARIADNE_ASSERT_MSG(this->has_invariant(location,event),"No invariant "<<event<<" in location "<<location);
    const Mode& mode = this->_modes[location];
    return mode._invariants[event]._guard;
}


RealScalarFunction
MonolithicHybridAutomaton::guard_function(DiscreteLocation location, DiscreteEvent event) const
{
    ARIADNE_ASSERT_MSG(this->has_guard(location,event),"No guard "<<event<<" in location "<<location);
    const Mode& mode = this->_modes[location];
    return mode._transitions[event]._guard;
}


DiscreteLocation
MonolithicHybridAutomaton::target(DiscreteLocation source, DiscreteEvent event) const
{
    return this->transition(source,event)._target;
}

RealVectorFunction
MonolithicHybridAutomaton::reset_function(DiscreteLocation source, DiscreteEvent event) const
{
    return this->transition(source,event)._reset;
}





RealSpace
MonolithicHybridAutomaton::continuous_state_space(DiscreteLocation location) const
{
    return this->mode(location)._variable_names;
}


HybridSpace
MonolithicHybridAutomaton::state_space() const
{
    MonolithicHybridSpace result;
    for(Map<DiscreteLocation,Mode>::const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        result.new_location(mode_iter->first,continuous_state_space(mode_iter->first));
    }
    return result;
}



std::ostream&
MonolithicHybridAutomaton::write(std::ostream& os) const
{
    MonolithicHybridAutomaton const& automaton = *this;
    typedef MonolithicHybridAutomaton::Mode Mode;
    typedef MonolithicHybridAutomaton::Transition Transition;
    typedef MonolithicHybridAutomaton::Invariant Invariant;

    os << "MonolithicHybridAutomaton( \n";
    for(Map<DiscreteLocation,Mode>::const_iterator mode_iter=automaton._modes.begin();
        mode_iter!=automaton._modes.end(); ++mode_iter)
    {
        const Mode& mode = mode_iter->second;
        os << "  " << mode._location << ": " << mode._dynamic << ";\n";
        for(Map<DiscreteEvent,Invariant>::const_iterator invariant_iter=mode._invariants.begin();
            invariant_iter!=mode._invariants.end(); ++invariant_iter)
        {
            const Invariant& invariant = invariant_iter->second;
            os << "    " << invariant._event << ": "  << invariant._kind << ", " << invariant._guard << "<=0,\n";
        }
        for(Map<DiscreteEvent,Transition>::const_iterator transition_iter=mode._transitions.begin();
            transition_iter!=mode._transitions.end(); ++transition_iter)
        {
            const Transition& transition = transition_iter->second;
            os << "    " << transition._event << ": " << transition._kind << ", " << transition._guard << ">=0, "
               << transition._target << ", " << transition._reset << ";\n";
        }
    }
    return os << ")\n";
}




}
