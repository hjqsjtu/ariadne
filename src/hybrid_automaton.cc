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

#include "macros.h"
#include "stlio.h"
#include "function.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "grid_set.h"
#include "discrete_event.h"
#include "space.h"


namespace Ariadne {

typedef uint DimensionType;

class HybridSet {};


uint
DiscreteMode::
dimension() const
{
    return this->_dynamic.argument_size();
}


DiscreteMode::
DiscreteMode(DiscreteState location,
             const VectorFunction& dynamic)
    :  _location(location), _dynamic(dynamic), _invariants(), _grid(new Grid(dynamic.argument_size()))
{
}


/*
DiscreteMode::
DiscreteMode(DiscreteState location,
             const boost::shared_ptr< const VectorFunction > dynamic,
             const std::map< DiscreteEvent, boost::shared_ptr< const VectorFunction > >& invariants)
    :  _location(location), _dynamic(dynamic), _invariants(invariants), _grid(new Grid(dynamic->argument_size()))
{
    ARIADNE_ASSERT(dynamic->result_size()==dynamic->argument_size());
    for(uint i=0; i!=invariants.size(); ++i) {
        ARIADNE_ASSERT(invariants[i]->argument_size()==dynamic->argument_size());
        ARIADNE_ASSERT(invariants[i]->result_size()==1u);
    }
}
*/


std::ostream&
operator<<(std::ostream& os, const DiscreteMode& mode)
{
    return os << "DiscreteMode( "
              << "location=" << mode.location() << ", "
              << "dynamic=" << mode.dynamic() << ", "
              << "invariants=" << mode.invariants() << ", "
              << "grid=" << mode.grid() << " )";
}





DiscreteTransition::
DiscreteTransition(DiscreteEvent event,
                   const DiscreteMode& source,
                   const DiscreteMode& target,
                   const VectorFunction& reset,
                   const VectorFunction& activation,
                   bool forced)
    : _event(event), _source(source.location()), _target(target.location()),
      _activation(activation), _reset(reset), _forced(forced)
{
    ARIADNE_ASSERT(activation.result_size()==1);
    ARIADNE_ASSERT(activation.argument_size()==source.dimension());
    ARIADNE_ASSERT(reset.argument_size()==source.dimension());
    ARIADNE_ASSERT(reset.result_size()==target.dimension());
}




std::ostream&
operator<<(std::ostream& os, const DiscreteTransition& transition)
{
    return os << "DiscreteTransition( "
              << "event=" << transition.event() << ", "
              << "source=" << transition.source() << ", "
              << "target=" << transition.target() << ", "
              << "reset=" << transition.reset() << ", "
              << "activation=" << transition.activation() << ", "
              << (transition.forced() ? "forced" : "unforced") << " )";
}




HybridAutomaton::~HybridAutomaton()
{
}

HybridAutomaton::HybridAutomaton()
{
}

HybridAutomaton::HybridAutomaton(const std::string& name)
    : _name(name)
{
}





const DiscreteMode&
HybridAutomaton::new_mode(DiscreteState location,
                          const VectorFunction& dynamic)
{
    if(this->has_mode(location)) {
        throw std::runtime_error("The hybrid automaton already has a mode with the given id");
    }
    if(dynamic.result_size()!=dynamic.argument_size()) {
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_mode(location,dynamic)",
            "The dynamic has argument size " << dynamic.argument_size()
                << " and result size " << dynamic.result_size() << ", so does not define a vector field.");
    }
    this->_modes.push_back(DiscreteMode(location,dynamic));
    return this->mode(location);
}


const DiscreteMode&
HybridAutomaton::new_invariant(DiscreteState location,
                               const ScalarFunction& invariant)
{
    return this->new_invariant(location, VectorFunction(1u,invariant));
}


const DiscreteMode&
HybridAutomaton::new_invariant(DiscreteState location,
                               const VectorFunction& invariant)
{
    if(!this->has_mode(location)) {
        throw std::runtime_error("The location of the invariant must be in the automaton.");
    }
    DiscreteMode& mode=const_cast<DiscreteMode&>(this->mode(location));
    if(invariant.argument_size()!=mode.dimension()) {
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_invariant(location,invariant)",
            "The invariant has argument size " << invariant.argument_size()
                << " but the mode has state-space dimension " << mode.dimension());
    }
    if(invariant.result_size()!=1u) {
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_invariant(location,invariant)",
            "The invariant has result size " << invariant.result_size()
                << " but only scalar invariants are currently supported.");
    }
    DiscreteEvent invariant_event(std::string("invariant_"+to_str(mode._invariants.size())));
    mode._invariants[invariant_event]=invariant;
    return mode;
}

const DiscreteTransition&
HybridAutomaton::new_transition(DiscreteEvent event,
                                DiscreteState source,
                                DiscreteState target,
                                const VectorFunction &reset,
                                const VectorFunction &activation,
                                bool forced)
{
    if(!event.is_transition()) {
        throw std::runtime_error("Transition event names cannot start with \"invariant\".");    
    }
    if(this->has_transition(event,source)) {
        throw std::runtime_error("The automaton already has a transition with the given id and source id.");
    }
    if(!this->has_mode(source)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given source id");
    }
    if(!this->has_mode(target)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given desitination id");
    }

    const DiscreteMode& source_mode=this->mode(source);
    const DiscreteMode& target_mode=this->mode(target);
    this->_transitions.push_back(DiscreteTransition(event,source_mode,target_mode,reset,activation,forced));
    return this->transition(event,source);
}

const DiscreteTransition&
HybridAutomaton::
new_transition(DiscreteEvent event,
               DiscreteState source,
               DiscreteState target,
               const VectorFunction &reset,
               const ScalarFunction &activation,
               bool forced)
{
    return this->new_transition(event,source,target,reset,VectorFunction(1u,activation),forced);
}

const DiscreteTransition&
HybridAutomaton::
new_transition(DiscreteEvent event,
               const DiscreteMode &source,
               const DiscreteMode &target,
               const VectorFunction &reset,
               const VectorFunction &activation,
               bool forced)
{
    DiscreteState source_id=source.location();
    DiscreteState target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,activation,forced);
}

const DiscreteTransition&
HybridAutomaton::
new_transition(DiscreteEvent event,
               const DiscreteMode &source,
               const DiscreteMode &target,
               const VectorFunction &reset,
               const ScalarFunction &activation,
               bool forced)
{
    DiscreteState source_id=source.location();
    DiscreteState target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,VectorFunction(1u,activation),forced);
}

const DiscreteTransition&
HybridAutomaton::
new_forced_transition(DiscreteEvent event,
                      DiscreteState source,
                      DiscreteState target,
                      const VectorFunction &reset,
                      const ScalarFunction &activation)
{
    return this->new_transition(event,source,target,reset,VectorFunction(1u,activation),true);
}

const DiscreteTransition&
HybridAutomaton::
new_forced_transition(DiscreteEvent event,
                      DiscreteState source,
                      DiscreteState target,
                      const VectorFunction &reset,
                      const VectorFunction &activation)
{
    return this->new_transition(event,source,target,reset,activation,true);
}

const DiscreteTransition&
HybridAutomaton::
new_forced_transition(DiscreteEvent event,
                      const DiscreteMode &source,
                      const DiscreteMode &target,
                      const VectorFunction &reset,
                      const VectorFunction &activation)
{
    DiscreteState source_id=source.location();
    DiscreteState target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,activation,true);
}

const DiscreteTransition&
HybridAutomaton::
new_forced_transition(DiscreteEvent event,
                      const DiscreteMode &source,
                      const DiscreteMode &target,
                      const VectorFunction &reset,
                      const ScalarFunction &activation)
{
    DiscreteState source_id=source.location();
    DiscreteState target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,VectorFunction(1u,activation),true);
}


const DiscreteTransition&
HybridAutomaton::
new_unforced_transition(DiscreteEvent event,
                        DiscreteState source,
                        DiscreteState target,
                        const VectorFunction &reset,
                        const ScalarFunction &activation)
{
    return this->new_transition(event,source,target,reset,VectorFunction(1u,activation),false);
}

const DiscreteTransition&
HybridAutomaton::
new_unforced_transition(DiscreteEvent event,
                        DiscreteState source,
                        DiscreteState target,
                        const VectorFunction &reset,
                        const VectorFunction &activation)
{
    return this->new_transition(event,source,target,reset,activation,false);

}

const DiscreteTransition&
HybridAutomaton::
new_unforced_transition(DiscreteEvent event,
                        const DiscreteMode &source,
                        const DiscreteMode &target,
                        const VectorFunction &reset,
                        const VectorFunction &activation)
{
    DiscreteState source_id=source.location();
    DiscreteState target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,activation,false);
}

const DiscreteTransition&
HybridAutomaton::
new_unforced_transition(DiscreteEvent event,
                        const DiscreteMode &source,
                        const DiscreteMode &target,
                        const VectorFunction &reset,
                        const ScalarFunction &activation)
{
    DiscreteState source_id=source.location();
    DiscreteState target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,VectorFunction(1u,activation),false);
}



bool
HybridAutomaton::has_mode(DiscreteState state) const
{
    // FIXME: This is a hack since we use std::list which cannot be searched by id.
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            if(mode_iter->location()==state) {
                return true;
            }
        }
    return false;
}


bool
HybridAutomaton::has_transition(DiscreteEvent event, DiscreteState source) const
{
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->event()==event && transition_iter->source()==source) {
                return true;
            }
        }
    return false;
}


bool
HybridAutomaton::has_accessible_constant(const RealConstant& con) const
{
    return _accessible_constants.find(con) != _accessible_constants.end();
}


HybridSet
HybridAutomaton::invariant() const
{
    ARIADNE_NOT_IMPLEMENTED;
}



HybridSpace
HybridAutomaton::state_space() const
{
    HybridSpace result;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            result[mode_iter->location()]=mode_iter->dimension();
        }
    return result;
}


const std::list< DiscreteMode >&
HybridAutomaton::modes() const
{
    return this->_modes;
}



const DiscreteMode&
HybridAutomaton::mode(DiscreteState state) const
{
    // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given discrete state.
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            if(mode_iter->location()==state) {
                return *mode_iter;
            }
        }
    throw std::runtime_error("The hybrid automaton does not have a mode with the given id.");
}


const std::list< DiscreteTransition >&
HybridAutomaton::transitions() const
{
    return this->_transitions;
}


RealConstantSet
HybridAutomaton::accessible_constants() const
{
	RealConstantSet result;

	for (RealConstantSet::const_iterator constant_it = _accessible_constants.begin();
										 constant_it != _accessible_constants.end();
									   ++constant_it)
			result.insert(*constant_it);

	return result;
}


std::list< DiscreteTransition >
HybridAutomaton::transitions(DiscreteState source) const
{
    std::list< DiscreteTransition > result;
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->source()==source) {
                result.push_back(*transition_iter);
            }
        }
    return result;
}


std::map<DiscreteEvent,VectorFunction>
HybridAutomaton::blocking_guards(DiscreteState source) const
{
    std::map<DiscreteEvent,VectorFunction> result;
    const DiscreteMode& mode=this->mode(source);
    for(invariant_const_iterator invariant_iter=mode._invariants.begin();
        invariant_iter!=mode._invariants.end(); ++invariant_iter)
    {
        const DiscreteEvent event=invariant_iter->first;
        const VectorFunction invariant=invariant_iter->second;
        result[event]=invariant;
    }

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source()==source && transition_iter->forced()) {
            const DiscreteEvent event=transition_iter->event();
            const VectorFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}


std::map<DiscreteEvent,VectorFunction>
HybridAutomaton::permissive_guards(DiscreteState source) const
{
    std::map<DiscreteEvent,VectorFunction> result;

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source()==source && !transition_iter->forced()) {
            const DiscreteEvent event=transition_iter->event();
            const VectorFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}



const DiscreteTransition&
HybridAutomaton::transition(DiscreteEvent event, DiscreteState source) const
{
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->event()==event && transition_iter->source()==source) {
                return *transition_iter;
            }
        }
    ARIADNE_THROW(std::runtime_error, "HybridAutomaton::transition(event, source)",
        "The hybrid automaton does not have a transition with event \"" << event << "\" and source \"" << source << "\".");
}


const std::string&
HybridAutomaton::name() const
{
    return this->_name;
}

HybridAutomaton& 
HybridAutomaton::operator=(const std::pair< HybridAutomaton, RealSpace >& pair) 
{
    //std::cout << pair.first << std::endl;
    // Handle self-assignment
    if(this == &(pair.first)) return *this;
    this->_name = pair.first.name();
    this->_modes = pair.first.modes();
    this->_transitions = pair.first.transitions();
    return *this;
}


std::ostream&
operator<<(std::ostream& os, const HybridAutomaton& ha)
{
    return os << "HybridAutomaton( name=" << ha.name() <<
    							", accessible constants=" << ha.accessible_constants() <<
    							", modes=" << ha.modes() <<
    							", transitions=" << ha.transitions() <<
    							")";
}

void 
HybridAutomaton::substitute(Constant<Real> con, const Real& c)
{
	RealConstantSet::iterator constant_it = _accessible_constants.find(con);
	if (constant_it != _accessible_constants.end()) {
		_accessible_constants.erase(*constant_it);
		_accessible_constants.insert(RealConstant(con.name(),c));
	}
	else {
		ARIADNE_FAIL_MSG("The constant to substitute is not registered as accessible.");
	}

	// Substitutes on the modes and transitions
	for (std::list<DiscreteMode>::iterator modes_it=this->_modes.begin();modes_it!=this->_modes.end();modes_it++)
		modes_it->substitute(con,c);
	for (std::list<DiscreteTransition>::iterator trans_it=this->_transitions.begin();trans_it!=this->_transitions.end();trans_it++)
		trans_it->substitute(con,c);
}

void
HybridAutomaton::substitute(const RealConstantSet& cons)
{
	// Restores the system
	for (RealConstantSet::const_iterator const_it = cons.begin(); const_it != cons.end(); ++const_it)
		substitute(*const_it);
}


void
HybridAutomaton::substitute(const RealConstantSet& cons, bool use_midpoint)
{
	for (RealConstantSet::const_iterator cons_it = cons.begin(); cons_it != cons.end(); ++cons_it) {
		if (use_midpoint)
			substitute(*cons_it,cons_it->value().midpoint());
		else
			substitute(*cons_it);
	}
}

void
HybridAutomaton::register_accessible_constant(RealConstant con)
{
	// Silently accepts duplicated constant registrations
	if (!this->has_accessible_constant(con))
		this->_accessible_constants.insert(con);
}

const Real&
HybridAutomaton::accessible_constant_value(const String& name) const
{
	for (RealConstantSet::const_iterator constants_it=this->_accessible_constants.begin(); constants_it!=this->_accessible_constants.end(); constants_it++)
		if (constants_it->name() == name)
			return constants_it->value();

	ARIADNE_FAIL_MSG("The constant is not registered as accessible.");
}

RealConstantSet
HybridAutomaton::nonsingleton_accessible_constants() const
{
	RealConstantSet result;

	const RealConstantSet& constants = _accessible_constants;

	for (RealConstantSet::const_iterator constant_it = constants.begin();
												 constant_it != constants.end();
												 ++constant_it) {
		if (!constant_it->value().singleton())
			result.insert(*constant_it);
	}

	return result;
}


}
