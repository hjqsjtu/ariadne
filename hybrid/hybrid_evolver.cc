/***************************************************************************
 *            hybrid_evolver.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include <typeinfo>
#include "config.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "function/polynomial.h"
#include "function/function.h"
#include "function/function_model.h"
#include "geometry/grid_set.h"
#include "hybrid/hybrid_time.h"
#include "hybrid/hybrid_automaton_interface.h"
#include "hybrid/hybrid_evolver.h"
#include "dynamics/orbit.h"

#include "solvers/integrator.h"
#include "solvers/solver.h"

#include "hybrid/hybrid_evolver.h"

namespace {

} // namespace

namespace Ariadne {

static const ExactFloat zero=0;

ScalarTaylorFunction unchecked_compose(const ValidatedScalarFunctionModel& f, const VectorTaylorFunction& g);

inline auto operator+(int n, ValidatedFloat x) -> decltype(ExactFloat(n)+x) { return ExactFloat(n)+x; }
inline auto operator-(int n, ValidatedFloat x) -> decltype(ExactFloat(n)-x) { return ExactFloat(n)-x; }
inline auto operator/(ExactFloat x, uint n) -> decltype(x/ExactFloat(n)) { return x/ExactFloat(n); }

typedef Integer Natural;

static const DiscreteEvent final_event("_tmax_");
static const DiscreteEvent step_event("_h_");

typedef Vector<ExactFloat> ExactFloatVector;
typedef Vector<ExactInterval> ExactIntervalVector;

std::ostream& operator<<(std::ostream& os, const HybridTerminationCriterion& termination) {
    return os << "HybridTerminationCriterion( maximum_time=" << termination.maximum_time()
              << ", maximum_steps="<<termination.maximum_steps()
              << ", terminating_events="<<termination.terminating_events() << " )";
}

std::ostream& operator<<(std::ostream& os, const CrossingKind& crk) {
    switch(crk) {
        case CrossingKind::DEGENERATE: os<<"DEGENERATE"; break;
        case CrossingKind::POSITIVE: os<<"POSITIVE"; break;
        case CrossingKind::NEGATIVE: os<<"NEGATIVE"; break;
        case CrossingKind::INCREASING: os<<"INCREASING"; break;
        case CrossingKind::DECREASING: os<<"DECREASING"; break;
        case CrossingKind::CONVEX: os<<"CONVEX"; break;
        case CrossingKind::CONCAVE: os<<"CONCAVE"; break;
        case CrossingKind::TRANSVERSE: os<<"TRANSVERSE"; break;
        case CrossingKind::GRAZING: os<<"GRAZING"; break;
        default: os << "unknown"; break;
    } return os;
}

std::ostream& operator<<(std::ostream& os, const DirectionKind& dir) {
    switch(dir) {
        case DirectionKind::POSITIVE: os<<"POSITIVE"; break;
        case DirectionKind::NEGATIVE: os<<"NEGATIVE"; break;
        case DirectionKind::INCREASING: os<<"INCREASING"; break;
        case DirectionKind::DECREASING: os<<"DECREASING"; break;
        case DirectionKind::CONVEX: os<<"CONVEX"; break;
        case DirectionKind::CONCAVE: os<<"CONCAVE"; break;
        case DirectionKind::INDETERMINATE: os<<"INDETERMINATE"; break;
        default: os << "UNKNOWN"; break;
    } return os;
}

std::ostream& operator<<(std::ostream& os, const StepKind& stpk) {
    switch(stpk) {
        case StepKind::CONSTANT_EVOLUTION_TIME: os<<"CONSTANT_EVOLUTION_TIME"; break;
        case StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME: os<<"SPACETIME_DEPENDENT_EVOLUTION_TIME"; break;
        case StepKind::PARAMETER_DEPENDENT_EVOLUTION_TIME: os<<"PARAMETER_DEPENDENT_EVOLUTION_TIME"; break;
        case StepKind::PARAMETER_DEPENDENT_FINISHING_TIME: os<<"PARAMETER_DEPENDENT_FINISHING_TIME"; break;
        case StepKind::SPACETIME_DEPENDENT_FINISHING_TIME: os<<"SPACETIME_DEPENDENT_FINISHING_TIME"; break;
        case StepKind::CONSTANT_FINISHING_TIME: os<<"CONSTANT_FINISHING_TIME"; break;
        default: os << "UNKNOWN"; break;
    } return os;
}

std::ostream& operator<<(std::ostream& os, const FinishingKind& fnshk) {
    switch(fnshk) {
        case FinishingKind::BEFORE_FINAL_TIME: os<<"BEFORE_FINAL_TIME"; break;
        case FinishingKind::AT_FINAL_TIME: os<<"AT_FINAL_TIME"; break;
        case FinishingKind::AFTER_FINAL_TIME: os<<"AFTER_FINAL_TIME"; break;
        case FinishingKind::STRADDLE_FINAL_TIME: os<<"STRADDLE_FINAL_TIME"; break;
        default: os << "UNKNOWN"; break;
    } return os;
}

std::ostream& operator<<(std::ostream& os, const TransitionData& transition) {
    return os << "kind="<<transition.event_kind<<", guard="<<transition.guard_function<<", "
                 "target="<<transition.target<<", reset="<<transition.reset_function;
}

std::ostream& operator<<(std::ostream& os, const TimingData& timing) {
    os << "step_kind="<<timing.step_kind<<", finishing_kind="<<timing.finishing_kind<<", step_size="<<timing.step_size<<", "
       << "final_time="<<timing.final_time;
    if(timing.step_kind==StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME) {
        os <<", spacial_evolution_time="<<timing.spacetime_dependent_evolution_time;
    } else if(timing.step_kind==StepKind::PARAMETER_DEPENDENT_FINISHING_TIME) {
        os << ", parameter_dependent_finishing_time="<<timing.parameter_dependent_finishing_time;
    }
    os <<", parameter_dependent_evolution_time="<<timing.parameter_dependent_evolution_time;
    return os;
}

std::ostream& operator<<(std::ostream& os, const CrossingData& crossing_data) {
    os << "{kind="<<crossing_data.crossing_kind;
    if(crossing_data.crossing_kind==CrossingKind::TRANSVERSE) {
        os << ", crossing_time="<<crossing_data.crossing_time;
    }
    if(crossing_data.crossing_kind==CrossingKind::GRAZING) {
        os << ", critical_time="<<crossing_data.critical_time;
    }
    os << "}";
    return os;
}

// Test if an event 'blocks' continuous evolution.
bool is_blocking(EventKind evk) {
    switch(evk) {
        case INVARIANT: case PROGRESS: case URGENT: case IMPACT:
            return true;
        case PERMISSIVE:
            return false;
        default:
            ARIADNE_FAIL_MSG("EventKind "<<evk<<" not recognised by is_blocking(...) predicate.");
    }
}

// Test if an event 'activates' a discrete transition.
bool is_activating(EventKind evk) {
    switch(evk) {
        case PERMISSIVE: case URGENT: case IMPACT:
            return true;
        case INVARIANT: case PROGRESS:
            return false;
        default:
            ARIADNE_FAIL_MSG("EventKind "<<evk<<" not recognised by is_activating(...) predicate.");
    }
}

// Extract the blocking events.
Set<DiscreteEvent> blocking_events(const Map<DiscreteEvent,TransitionData>& transitions) {
    Set<DiscreteEvent> events;
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        if(is_blocking(transition_iter->second.event_kind)) {
            events.insert(transition_iter->first);
        }
    }
    return events;
}

// Extract the activating events.
Set<DiscreteEvent> activating_events(const Map<DiscreteEvent,TransitionData>& transitions) {
    Set<DiscreteEvent> events;
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        if(is_activating(transition_iter->second.event_kind)) {
            events.insert(transition_iter->first);
        }
    }
    return events;
}





Orbit<HybridEnclosure>
HybridEvolverBase::
orbit(const HybridBox& initial_box,
      const HybridTerminationCriterion& termination,
      Semantics semantics) const
{
    ARIADNE_LOG(2,"HybridEvolverBase::orbit(HybridAutomaton, HybridBox, HybridTime, Semantics)\n");
    ARIADNE_LOG(3,"factory="<<this->function_factory()<<"\n");
    ARIADNE_LOG(3,"initial_box="<<initial_box<<"\n");
    HybridEnclosure initial_enclosure(initial_box,this->function_factory());
    ARIADNE_LOG(3,"initial_enclosure="<<initial_enclosure<<"\n");
    return this->orbit(initial_enclosure,termination,semantics);
}

Orbit<HybridEnclosure>
HybridEvolverBase::
orbit(const HybridSet& initial_set,
      const HybridTerminationCriterion& termination,
      Semantics semantics) const
{
    ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
    HybridEnclosure initial_enclosure(initial_set,_sys_ptr->continuous_state_space(initial_set.location()),this->function_factory());
    ARIADNE_LOG(3,"initial_enclosure="<<initial_enclosure<<"\n");
    return this->orbit(initial_enclosure,termination,semantics);
}


Orbit<HybridEnclosure>
HybridEvolverBase::
orbit(const HybridEnclosure& initial,
      const HybridTerminationCriterion& termination,
      Semantics semantics) const
{
    ARIADNE_LOG(2,"\nHybridEvolverBase::orbit(...): verbosity="<<verbosity<<"\n");

    EvolutionData evolution_data;
    evolution_data.semantics=semantics;
    evolution_data.initial_sets.push_back(HybridEnclosure(initial));
    while(!evolution_data.initial_sets.empty()) {
        this->_evolution_in_mode(evolution_data,termination);
    }
    ARIADNE_ASSERT(evolution_data.initial_sets.empty());
    ARIADNE_ASSERT(evolution_data.working_sets.empty());

    Orbit<HybridEnclosure> orbit(initial);
    orbit.adjoin_intermediate(ListSet<HybridEnclosure>(evolution_data.intermediate_sets));
    orbit.adjoin_reach(evolution_data.reach_sets);
    orbit.adjoin_final(evolution_data.final_sets);
    return orbit;
}

HybridEvolverBase::FunctionFactoryType* make_taylor_function_factory();

HybridEvolverBase::HybridEvolverBase(const SystemType& system)
{
    this->_create(system,make_taylor_function_factory());
}

HybridEvolverBase::HybridEvolverBase(const SystemType& system,
                                     const FunctionFactoryType& factory)
{
    this->_create(system,factory.clone());
}

void
HybridEvolverBase::_create(
        const SystemType& system,
        FunctionFactoryType* factory)
{
    this->_sys_ptr=std::shared_ptr<SystemType>(system.clone());
    this->_function_factory_ptr=std::shared_ptr<FunctionFactoryType>(factory);
    this->_solver_ptr=std::shared_ptr<SolverInterface>(new IntervalNewtonSolver(1e-8,12));
    this->ALLOW_CREEP=true;
    this->ALLOW_UNWIND=false;
    //this->_configuration_ptr=std::shared_ptr<ConfigurationType>(new ConfigurationType());
}


const HybridEvolverBase::SystemType&
HybridEvolverBase::system() const
{
    return *this->_sys_ptr;
}

HybridEvolverBase::ConfigurationType&
HybridEvolverBase::configuration()
{
    return *this->_configuration_ptr;
}

const HybridEvolverBase::ConfigurationType&
HybridEvolverBase::configuration() const
{
    return *this->_configuration_ptr;
}

void
HybridEvolverBase::set_function_factory(const FunctionFactoryType& factory)
{
    this->_function_factory_ptr=std::shared_ptr<FunctionFactoryType>(factory.clone());
}

const HybridEvolverBase::FunctionFactoryType&
HybridEvolverBase::function_factory() const
{
    return*this->_function_factory_ptr;
}

void
HybridEvolverBase::set_integrator(const IntegratorInterface& integrator)
{
    this->_integrator_ptr=std::shared_ptr<IntegratorInterface>(integrator.clone());
}

void
HybridEvolverBase::set_solver(const SolverInterface& solver)
{
    this->_solver_ptr=std::shared_ptr<SolverInterface>(solver.clone());
}


HybridEvolverBase::EnclosureType
HybridEvolverBase::enclosure(const HybridBox& initial_box) const
{
    return HybridEnclosure(initial_box,this->function_factory());
}

HybridEvolverBase::EnclosureType
HybridEvolverBase::enclosure(const HybridSet& initial_set) const
{
    return HybridEnclosure(initial_set,_sys_ptr->continuous_state_space(initial_set.location()),this->function_factory());
}


void
HybridEvolverBase::
_evolution(ListSet<HybridEnclosure>& final,
           ListSet<HybridEnclosure>& reachable,
           ListSet<HybridEnclosure>& intermediate,
           HybridEnclosure const& initial_set,
           HybridTerminationCriterion const& termination_criterion,
           Semantics semantics,
           bool reach) const
{
    EvolutionData evolution_data;
    evolution_data.semantics=semantics;
    evolution_data.initial_sets.push_back(HybridEnclosure(initial_set));

    while(!evolution_data.initial_sets.empty()) {
        this->_evolution_in_mode(evolution_data,termination_criterion);
    }

    final=evolution_data.final_sets;
    reachable=evolution_data.reach_sets;
    intermediate=evolution_data.intermediate_sets;
}

struct EvolutionStepData {
    EvolutionStepData() : progress(false), finishing(false), events() { }
    bool progress;
    bool finishing;
    Set<DiscreteEvent> events;
};

void
HybridEvolverBase::
_log_summary(const EvolutionData& evolution_data, HybridEnclosure const& starting_set) const
{
    UpperBox starting_bounding_box=starting_set.space_bounding_box();
    UpperInterval starting_time_range=starting_set.time_range();
    UpperInterval starting_dwell_time_range=starting_set.dwell_time_range();
    int old_precision = std::clog.precision();
    if(verbosity>=1) { std::clog<<"\n"; }
    ARIADNE_LOG(1,(verbosity==1?"\r":"")
            <<"#w="<<std::setw(4)<<std::left<<evolution_data.initial_sets.size()+1u
            <<"#r="<<std::setw(5)<<std::left<<evolution_data.reach_sets.size()
            <<"#f="<<std::setw(4)<<std::left<<evolution_data.final_sets.size()
            <<"#e="<<std::setw(3)<<std::left<<starting_set.previous_events().size()
            <<" #p="<<std::setw(2)<<std::left<<starting_set.number_of_parameters()
            <<" #c="<<std::setw(1)<<std::left<<starting_set.number_of_constraints()
            <<" t=["<<std::setw(6)<<std::setprecision(3)<<std::left<<std::fixed<<starting_time_range.lower()
            <<","<<std::setw(6)<<std::left<<std::fixed<<starting_time_range.upper()<<"]"<<std::flush
            <<" dwt=["<<std::setw(6)<<std::setprecision(3)<<std::left<<std::fixed<<starting_dwell_time_range.lower()
            <<","<<std::setw(6)<<std::left<<std::fixed<<starting_dwell_time_range.upper()<<"]"<<std::flush
            <<" c="<<starting_bounding_box.centre()
            <<" r="<<std::setw(4)<<std::fixed<<std::setprecision(3)<<starting_bounding_box.radius()
            <<" te="<<std::setw(7)<<std::scientific<<std::setprecision(1)<<starting_set.time_function().error()<<std::flush
            <<" se="<<std::setw(7)<<std::scientific<<std::setprecision(1)<<sup_norm(starting_set.space_function().errors())<<std::fixed<<std::flush
            <<" l="<<std::left<<starting_set.location()
            <<" e="<<starting_set.previous_events()
            <<"                      \n"<<std::setprecision(old_precision));
    ARIADNE_LOG(4,"\r    \r"<<starting_set<<"\n");

}

Map<DiscreteEvent,TransitionData>
HybridEvolverBase::
_extract_transitions(DiscreteLocation const& location) const
{
    ARIADNE_LOG(3,"HybridEvolverBase::_extract_transitions(...)\n");
    const EffectiveVectorFunction& dynamic=_sys_ptr->dynamic_function(location);

    Map<DiscreteEvent,TransitionData> transitions;
    Set<DiscreteEvent> events = _sys_ptr->events(location);
    for(Set<DiscreteEvent>::const_iterator event_iter=events.begin();
        event_iter!=events.end(); ++event_iter)
    {
        DiscreteEvent event=*event_iter;
        EventKind event_kind=_sys_ptr->event_kind(location,event);
        ARIADNE_LOG(5,"event="<<event<<", kind="<<event_kind<<"\n");
        EffectiveScalarFunction constraint_function;
        if(is_activating(event_kind)) {
            constraint_function=_sys_ptr->guard_function(location,event);
        } else {
            constraint_function=_sys_ptr->invariant_function(location,event);
        }
        ARIADNE_LOG(5,"constraint_function="<<constraint_function<<"\n");
        EffectiveScalarFunction constraint_flow_derivative_function=lie_derivative(constraint_function,dynamic);
        EffectiveVectorFunction reset_function; DiscreteLocation target; RealSpace target_space;
        if(is_activating(event_kind)) {
            reset_function=_sys_ptr->reset_function(location,event);
            target=_sys_ptr->target(location,event);
            target_space=_sys_ptr->continuous_state_space(target);
        }
        TransitionData transition_data={event,event_kind,constraint_function,constraint_flow_derivative_function,target,reset_function,target_space};
        transitions.insert(event,transition_data);
    }
    ARIADNE_LOG(3,"transitions="<<transitions<<"\n");
    return transitions;
}

void
HybridEvolverBase::
_apply_invariants(HybridEnclosure& initial_set,
                  Map<DiscreteEvent,TransitionData> const& transitions) const
{
    ARIADNE_LOG(3,"HybridEvolverBase::_apply_invariants(...)\n");
    HybridEnclosure& invariant_set=initial_set;

    // Apply restrictions due to invariants
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        DiscreteEvent event=transition_iter->first;
        TransitionData const & transition=transition_iter->second;
        if(transition.event_kind==INVARIANT) {
            if (possibly(initial_set.satisfies(transition.guard_function>=0))) {
                invariant_set.new_invariant(event,transition.guard_function);
            }
        }
    }
}

void
HybridEvolverBase::
_process_starting_events(EvolutionData& evolution_data,
                        HybridEnclosure const& initial_set,
                        Map<DiscreteEvent,TransitionData> const& transitions) const
{
    ARIADNE_LOG(3,"HybridEvolverBase::_process_starting_events(...)\n");
    ARIADNE_ASSERT(evolution_data.working_sets.empty());
    HybridEnclosure invariant_set=initial_set;

    // Apply restrictions due to invariants
    this->_apply_invariants(invariant_set,transitions);

    // Set the flowable set, storing the invariant set as a base for jumps
    HybridEnclosure flowable_set = invariant_set;

    // Compute possibly initially active events
    Set<DiscreteEvent> events=transitions.keys();
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        DiscreteEvent event=transition_iter->first;
        TransitionData const & transition=transition_iter->second;
        // Assume that an impact cannot occur immediately after any other event.
        // This is true after the impact, since $L_{f}g$ changes sign.
        // After a different event, this is not so clear, though it should be
        // a modelling error for the set to intersect the interior of the guard
        // FIXME: Need to consider impact which really can occur immediately due to mapping to boundary.
        // TODO: The condition for an impact to occur is $L_{f}g>0$.
        //       Since this is now available, we should implement this!
        if(transition.event_kind!=INVARIANT && transition.event_kind!=IMPACT) {
            if(possibly(initial_set.satisfies(transition.guard_function>=0))) {
                if(transition.event_kind!=PROGRESS) {
                    HybridEnclosure immediate_jump_set=invariant_set;
                    immediate_jump_set.new_activation(event,transition.guard_function);
                    if(!definitely(immediate_jump_set.empty())) {
                        // Put the immediate jump set in the reached sets, since it does not occur in the flowable set
                        ARIADNE_LOG(2,event<<": "<<transition.event_kind<<", immediate\n");
                        evolution_data.reach_sets.append(immediate_jump_set);
                        immediate_jump_set.apply_reset(event,transition.target,transition.target_space,transition.reset_function);
                        ARIADNE_LOG(4,"immediate_jump_set="<<immediate_jump_set<<"\n");
                        evolution_data.intermediate_sets.append(immediate_jump_set);
                        evolution_data.initial_sets.append(immediate_jump_set);
                    }
                }
                if(transition.event_kind!=PERMISSIVE) {
                    flowable_set.new_invariant(event,transition.guard_function);
                }
            }
        }
    }

    // Put the flowable set in the starting sets for ordinary evolution
    if(!definitely(flowable_set.empty())) {
        evolution_data.working_sets.append(flowable_set);
    }
}

ValidatedVectorFunctionModel
HybridEvolverBase::
_compute_flow(EffectiveVectorFunction dynamic,
              ExactBox const& initial_box,
              const ExactFloat& maximum_step_size) const
{
    ARIADNE_LOG(7,"HybridEvolverBase::_compute_flow(...)\n");

    IntegratorInterface& integrator=*this->_integrator_ptr;

    // Compute flow and actual time step size used
    //
    // The Integrator classes compute the flow as a function on a symmetrical
    // time domain [-h,+h], since this means the time is centred at 0.
    // We then restrict to the time domain [0,h] since this can make evaluation
    // more accurate, and the time domain might be used explicitly for the domain
    // of the resulting set.
    ExactFloat step_size=maximum_step_size;
    ValidatedVectorFunctionModel flow_model=integrator.flow_step(dynamic,initial_box,step_size.raw());

    ARIADNE_LOG(6,"twosided_flow_model="<<flow_model<<"\n");
    ExactIntervalVector flow_domain=flow_model.domain();
    ARIADNE_ASSERT(step_size==flow_domain[flow_domain.size()-1u].upper());
    flow_domain[flow_domain.size()-1u]=ExactInterval(zero,step_size);
    flow_model=restrict(flow_model,flow_domain);
    ARIADNE_LOG(6,"flow_model="<<flow_model<<"\n");
    ARIADNE_LOG(2,"flow_model: step_size="<<step_size<<", errors="<<std::scientific<<flow_model.errors()<<", range="<<std::fixed<<flow_model.range()<<"\n");
    return flow_model;
}

Set<DiscreteEvent>
HybridEvolverBase::
_compute_active_events(EffectiveVectorFunction const& dynamic,
                       Map<DiscreteEvent,EffectiveScalarFunction> const& guards,
                       ValidatedVectorFunctionModel const& flow,
                       HybridEnclosure const& starting_set) const
{
    ARIADNE_LOG(7,"HybridEvolverBase::_compute_active_events(...)\n");
    Set<DiscreteEvent> events=guards.keys();
    Set<DiscreteEvent> active_events;
    HybridEnclosure reach_set=starting_set;
    UpperBox flow_bounds=flow.range();
    ARIADNE_LOG(8,"flow_bounds="<<flow_bounds<<"\n");
    reach_set.apply_full_reach_step(flow);
    ARIADNE_LOG(8,"reach_set="<<reach_set<<"\n");
    for(Set<DiscreteEvent>::iterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
        const DiscreteEvent event=*event_iter;
        const EffectiveScalarFunction& guard_function=guards[event];
        ARIADNE_LOG(8,"event="<<event<<", guard="<<guard_function<<", flow_derivative="<<lie_derivative(guard_function,dynamic)<<"\n");
        // First try a simple test based on the bounding box
        UpperInterval guard_range=apply(guard_function,reach_set.space_bounding_box());
        if(guard_range.upper()>=zero) {
            // Now make a set containing the complement of the constraint,
            // and test for emptiness. If the set is empty, then the guard is
            // not satisfied anywhere.
            HybridEnclosure test_set=reach_set;
            test_set.new_activation(event,guard_function);
            if(!definitely(test_set.empty())) {
                active_events.insert(*event_iter);
/*
                // FIXME: Need to allow permissive events with strictly decreasing guard.
                // Test direction of guard increase
                EffectiveScalarFunction flow_derivative = lie_derivative(guard_function,dynamic);
                ExactInterval flow_derivative_range = flow_derivative(flow_bounds);
                if(flow_derivative_range.upper()>zero) {
                    active_events.insert(*event_iter);
                }
*/
            }
        }
    }
    return active_events;
}


Map<DiscreteEvent,CrossingData>
HybridEvolverBase::
_compute_crossings(Set<DiscreteEvent> const& active_events,
                   EffectiveVectorFunction const& dynamic,
                   Map<DiscreteEvent,EffectiveScalarFunction> const& guards,
                   FlowFunctionModel const& flow,
                   HybridEnclosure const& initial_set) const
{
    ARIADNE_LOG(7,"HybridEvolverBase::_compute_crossings(...)\n");

    const SolverInterface& solver=*this->_solver_ptr;

    Map<DiscreteEvent,CrossingData> crossings;
    crossings.clear();

    ExactIntervalVector flow_spacial_domain=project(flow.domain(),range(0,flow.argument_size()-1u));
    ExactInterval flow_time_domain=flow.domain()[flow.argument_size()-1u];
    UpperBox flow_bounds=flow.range();
    for(Set<DiscreteEvent>::const_iterator event_iter=active_events.begin();
        event_iter!=active_events.end(); ++event_iter)
    {
        const DiscreteEvent event=*event_iter;
        EffectiveScalarFunction const& guard=guards[event];

        // Compute the derivative of the guard function g along flow lines of $\dot(x)=f(x)$
        // This is given by the Lie derivative at a point x, defined as $L_{f}g(x) = (\nabla g\cdot f)(x)$
        EffectiveScalarFunction derivative=lie_derivative(guard,dynamic);
        UpperInterval derivative_range=apply(derivative,flow_bounds);
        if(derivative_range.lower()>zero) {
            // If the derivative $L_{f}g$is strictly positive over the bounding box for the flow,
            // then the guard function is strictly increasing.
            // There is at most one crossing with the guard, and the time of this
            // crossing must be the time of the event along the trajectory.
            // The crossing time $\gamma(x_0)$ given the initial state can usually be computed
            // by solving the equation $g(\phi(x_0,\gamma(x_0))) = 0$
            ValidatedScalarFunctionModel crossing_time;
            try {
                crossing_time=solver.implicit(compose(guard,flow),flow_spacial_domain,flow_time_domain);
                if(crossing_time.error()>1e-8) { ARIADNE_LOG(2,event<<": crossing_time: error="<<crossing_time.error()<<", range="<<crossing_time.range()<<"\n"); }
                crossings[event]=CrossingData(CrossingKind::TRANSVERSE,crossing_time);
                ARIADNE_LOG(8,"crossing_time="<<crossing_time<<"\n");
            }
            catch(const UnknownSolutionException& e) {
                // If the crossing time cannot be computed, then we can still
                // use the fact that the crossing occurs as soon as $g(x(t))=0$.
                ARIADNE_LOG(2,"Error in computing crossing time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
                crossings[event]=CrossingData(CrossingKind::INCREASING);
            }
            catch(const SolverException& e) {
                // If the crossing time cannot be computed, then we can still
                // use the fact that the crossing occurs as soon as $g(x(t))=0$.
                ARIADNE_LOG(2,"Error in computing crossing time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
                crossings[event]=CrossingData(CrossingKind::INCREASING);
            }
            catch(const std::runtime_error& e) {
                ARIADNE_LOG(0,"Unexpected error in computing crossing time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
                ARIADNE_FAIL_MSG("ERROR!!");
                crossings[event]=CrossingData(CrossingKind::INCREASING);
            }
        } else if(derivative_range.upper()<zero) {
            // If the derivative is strictly negative over the bounding box for the flow,
            // then the guard function is strictly decreasing.
            // This means that the event is either initially active, or does not occur.
            // There is no need to compute a crossing time.
            crossings[event]=CrossingData(CrossingKind::DECREASING);
        } else {
            // If the derivative of the guard function along flow lines cannot be shown
            // to have a definite sign over the entire flow box, then try to compute
            // the sign of the second derivative $L_{f}^{2}g(x)=L_{f}L_{f}g(x)$.
            ValidatedScalarFunction second_derivative=lie_derivative(derivative,dynamic);
            UpperInterval second_derivative_bounds_range=apply(second_derivative,flow_bounds);
            UpperInterval second_derivative_flow_range=compose(second_derivative,flow).range();
            UpperInterval second_derivative_range=intersection(second_derivative_bounds_range,second_derivative_flow_range);
            if(second_derivative_range.lower()>zero) {
                // If the second derivative is positive, then either
                //    (i) the event is immediately active
                //   (ii) the event is never active, or
                //  (iii) the event is initially inactive, but becomes active
                //        due to a transverse crossing.
                //   (iv) the initial state is on the boundary of the guard
                //        set, possibly with the flow tangent to this set
                // We cannot compute the crossing time, even in case (iii),
                // due to the singularity due to the tangency in (iv). However,
                // we do know that in (iii), the event occurs when $t>0$ and
                // $g(\phi(x_0,t))=0$. The crossing time is not computed.
                crossings[event]=CrossingData(CrossingKind::CONVEX);
            } else if(second_derivative_range.upper()<zero) {
                // If the second derivative is negative, then the guard
                // values $g(x(t))$ are concave along flow lines. There are
                // four main cases:
                //   (i) The event is initially active.
                //  (ii) The event is not initially active, but later becomes active.
                // (iii) The event is never active, but would become active if
                //       flowing backwards in time.
                //  (iv) The event is never active, and the maximum value of
                //       the guard along the flow lines is zero.
                // Additionally, there is the degenerate case
                //   (v) At some point in the (forward) flow, the state touches
                //       the guard set at a point of tangency.
                // Due to the presence of the tangency, the event time is
                // not a smooth function of the initial state. Further, since
                // some trajectories cross the boundary of the guard set twice,
                // the condition $g(x(t))=0$ is not sufficient for determining
                // the jump time. A necessary and sufficient condition,
                // assuming the event is not initially active, is that
                // $g(x)=0$ and $L_{f}g(x)\geq 0$. Alternatively, we can compute
                // the <em>critical time</em> $|mu(x_0) at which the guard value
                // reaches a maximum. A necessary and sufficient condition
                // is then $g(\phi(x_0,t))=0$ and $t<= \mu(x_0)$.
                //   Note that while $g(\phi(x_0,t))=0$ and $L_{f}f(\phi(x_0,t))\geq0$ is a
                // necessary and sufficient condition for a crossing, there
                // is no necessary and sufficient condition for no crossing
                // which does not involve the critical time. A necessary and
                // sufficient condition for no crossing involving the critical
                // time is $(g(\phi(x_0,t))<=0 /\ t<=\mu(x_0)) \/ g(\phi(x_0,\mu(x_0)))<=0$
                try {
                    ValidatedScalarFunctionModel critical_time=solver.implicit(compose(derivative,flow),flow_spacial_domain,flow_time_domain);
                    if(critical_time.error()>1e-8) { ARIADNE_LOG(2,event<<": critical_time: error="<<critical_time.error()<<", range="<<critical_time.range()<<"\n"); }
                    crossings[event]=CrossingData(CrossingKind::GRAZING);
                    crossings[event].critical_time=critical_time;
                }
                catch(const SolverException& e) {
                    ARIADNE_LOG(2,"Error in computing critical time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
                    crossings[event]=CrossingData(CrossingKind::CONCAVE);
                }
                catch(const std::runtime_error& e) {
                    ARIADNE_LOG(0,"Unexpected error in computing critical time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
                    crossings[event]=CrossingData(CrossingKind::CONCAVE);
                }
            } else {
                // The crossing cannot be shown to be one of the kinds mentioned
                // above. A theoretically exact expression for the crossing
                // set is generally not available.
                crossings[event]=CrossingData(CrossingKind::DEGENERATE);
            }
        }
    }
    if(!crossings.empty()) {
        Map<DiscreteEvent,Tuple<CrossingKind,UpperInterval,ErrorFloat> > crossing_log_data;
        for (Map<DiscreteEvent,CrossingData>::const_iterator crossing_iter=crossings.begin(); crossing_iter!=crossings.end(); ++crossing_iter) {
            DiscreteEvent event=crossing_iter->first;
            CrossingKind crossing_kind=crossing_iter->second.crossing_kind;
            UpperInterval crossing_time_range(-infty,+infty);
            ErrorFloat crossing_time_error=0u;
            if (crossing_kind==CrossingKind::TRANSVERSE) {
                crossing_time_range=crossing_iter->second.crossing_time.range();
                crossing_time_error=crossing_iter->second.crossing_time.error();
            }
            crossing_log_data.insert(event,make_tuple(crossing_kind,crossing_time_range,crossing_time_error));
        }
        ARIADNE_LOG(2,"crossings: "<<crossing_log_data<<"\n");
    }
    return crossings;
}






void
HybridEvolverBase::
_recondition(HybridEnclosure& set) const
{
    set.recondition();
}



void
HybridEvolverBase::
_apply_reach_step(HybridEnclosure& set,
                  ValidatedVectorFunctionModel const& flow,
                  TimingData const& timing_data) const
{
    set.apply_reach_step(flow,timing_data.parameter_dependent_evolution_time);
}

void
HybridEvolverBase::
_apply_evolve_step(HybridEnclosure& set,
                  ValidatedVectorFunctionModel const& flow,
                  TimingData const& timing_data) const
{

    switch(timing_data.step_kind) {
        case StepKind::CONSTANT_EVOLUTION_TIME:
        case StepKind::SPACE_DEPENDENT_EVOLUTION_TIME:
        case StepKind::TIME_DEPENDENT_EVOLUTION_TIME:
        case StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME:
        case StepKind::PARAMETER_DEPENDENT_EVOLUTION_TIME:
            set.apply_evolve_step(flow,timing_data.parameter_dependent_evolution_time);
            break;
        case StepKind::PARAMETER_DEPENDENT_FINISHING_TIME:
        case StepKind::SPACETIME_DEPENDENT_FINISHING_TIME:
        case StepKind::CONSTANT_FINISHING_TIME:
            set.apply_finishing_evolve_step(flow,timing_data.parameter_dependent_finishing_time);
            break;
        default:
            ARIADNE_FAIL_MSG("Unhandled step kind "<<timing_data.step_kind);
    }
}

void
HybridEvolverBase::
_apply_guard_step(HybridEnclosure& set,
                  EffectiveVectorFunction const& dynamic,
                  ValidatedVectorFunctionModel const& flow,
                  TimingData const& timing_data,
                  TransitionData const& transition_data,
                  CrossingData const& crossing_data,
                  const Semantics semantics) const
{
    ARIADNE_LOG(4,"HybridEvolverBase::_apply_guard_step(...)\n");
    // Compute flow to guard set up to evolution time.
    HybridEnclosure& jump_set=set;
    const DiscreteEvent event=transition_data.event;
    ValidatedVectorFunctionModel starting_state=set.space_function();
    ValidatedVectorFunctionModel reach_starting_state=embed(starting_state,timing_data.evolution_time_domain);
    ValidatedScalarFunctionModel reach_step_time=embed(starting_state.domain(),timing_data.evolution_time_coordinate);
    ValidatedScalarFunctionModel step_time;

    switch(transition_data.event_kind) {
        case PERMISSIVE:
            // The continuous evolution is just the same as a reachability step,
            // so we need to embed the starting state and the step time function into one higher dimension.
            jump_set.apply_reach_step(flow,timing_data.parameter_dependent_evolution_time);
            jump_set.new_activation(event,transition_data.guard_function);
            break;
        case URGENT: case IMPACT:
            switch(crossing_data.crossing_kind) {
                case CrossingKind::TRANSVERSE:
                    step_time=unchecked_compose(crossing_data.crossing_time,starting_state);
                    // If the jump step might occur after the final evolution time, then introduce constraint that this does not happen
                    if(timing_data.step_kind!=StepKind::CONSTANT_EVOLUTION_TIME && (step_time-timing_data.parameter_dependent_evolution_time).range().upper()>zero) {
                        jump_set.new_parameter_constraint(step_event,step_time<=timing_data.parameter_dependent_evolution_time);
                    }
                    jump_set.apply_evolve_step(flow,unchecked_compose(crossing_data.crossing_time,starting_state));
                    break;
                case CrossingKind::INCREASING: case CrossingKind::CONVEX:
                    jump_set.apply_reach_step(flow,timing_data.parameter_dependent_evolution_time);
                    jump_set.new_guard(event,transition_data.guard_function);
                    break;
                case CrossingKind::CONCAVE: case CrossingKind::GRAZING:
                    jump_set.apply_reach_step(flow,timing_data.parameter_dependent_evolution_time);
                    jump_set.new_guard(event,transition_data.guard_function);
                    jump_set.new_invariant(event,-lie_derivative(transition_data.guard_function,dynamic));
                    break;
                case CrossingKind::DEGENERATE: // Just check positive derivative in this case; NOT EXACT
                    if(semantics==UPPER_SEMANTICS) {
                        jump_set.apply_reach_step(flow,timing_data.parameter_dependent_evolution_time);
                        jump_set.new_guard(event,transition_data.guard_function);
                        jump_set.new_invariant(event,-lie_derivative(transition_data.guard_function,dynamic));
                    } else {
                        // Make the empty set
                        jump_set.new_invariant(event,transition_data.guard_function*0+1);
                    }
                    break;
                case CrossingKind::NEGATIVE:
                case CrossingKind::POSITIVE:
                case CrossingKind::DECREASING:
                    // Since guard must be satisfied on impact, the jump set is empty
                    if(transition_data.event_kind==IMPACT) {
                        jump_set.new_invariant(event,transition_data.guard_function*0+1);
                    } else {
                        ARIADNE_WARN("Crossing "<<crossing_data<<" does not introduce additional restrictions on flow\n");
                    }
                    break;
                default:
                    ARIADNE_FAIL_MSG("Unhandled crossing kind in "<<crossing_data<<"\n");
            }
            break;
        default:
            ARIADNE_FAIL_MSG("Invalid event kind "<<transition_data.event_kind<<" for transition.");
    }

}


// Apply guard to a single set.
// In the case of concave crossings, splits the set into two, one part
// corresponding to points which actually hit the set (and stop on first crossing)
// the other part corresponding to points which miss the set.
void HybridEvolverBase::
_apply_guard(List<HybridEnclosure>& sets,
             const HybridEnclosure& starting_set,
             const ValidatedVectorFunctionModel& flow,
             const ValidatedScalarFunctionModel& elapsed_time,
             const TransitionData& transition_data,
             const CrossingData crossing_data,
             const Semantics semantics) const
{
    ARIADNE_LOG(3,"HybridEvolverBase::_apply_guard(...)\n");
    static const uint SUBDIVISIONS_FOR_DEGENERATE_CROSSING = 2;
    const DiscreteEvent event=transition_data.event;
    const ValidatedScalarFunction& guard_function=transition_data.guard_function;
    //const ValidatedScalarFunction& guard_flow_derivative_function=transition_data.guard_flow_derivative_function;
    ValidatedVectorFunctionModel starting_state=starting_set.space_function();
    if(elapsed_time.argument_size()>starting_state.argument_size()) {
        starting_state=embed(starting_state,elapsed_time.domain()[elapsed_time.domain().size()-1]);
    }
    ARIADNE_ASSERT(starting_state.domain()==elapsed_time.domain());

    List<HybridEnclosure>::iterator end=sets.end();
    for(List<HybridEnclosure>::iterator iter=sets.begin(); iter!=end; ++iter) {
        HybridEnclosure& set=*iter;

        switch(crossing_data.crossing_kind) {
            case CrossingKind::TRANSVERSE:
                //set.new_state_constraint(event, guard_function <= zero);
                set.new_invariant(event, guard_function);
                // Alternatively:
                // set.new_parameter_constraint(event, elapsed_time <= compose(crossing_data.crossing_time,starting_state) );
                break;
            case CrossingKind::CONVEX: case CrossingKind::INCREASING:
                //set.new_state_constraint(event, guard_function <= zero);
                set.new_invariant(event, guard_function);
                break;
            case CrossingKind::GRAZING: {
                ValidatedScalarFunctionModel critical_time = unchecked_compose(crossing_data.critical_time,starting_state);
                ValidatedScalarFunctionModel final_guard
                    = compose( guard_function, unchecked_compose( flow, join(starting_state, elapsed_time) ) );
                ValidatedScalarFunctionModel maximal_guard
                    = compose( guard_function, unchecked_compose( flow, join(starting_state, critical_time) ) );
                ARIADNE_ASSERT_MSG(starting_state.argument_size()==set.parameter_domain().size(),
                                   starting_state<<" "<<set);
                ARIADNE_ASSERT_MSG(critical_time.argument_size()==set.parameter_domain().size(),
                                   critical_time<<" "<<set);
                ARIADNE_ASSERT_MSG(maximal_guard.argument_size()==set.parameter_domain().size(),
                                   maximal_guard<<" "<<set);
                ARIADNE_ASSERT_MSG(final_guard.argument_size()==set.parameter_domain().size(),
                                   final_guard<<" "<<set);
                // If no points in the set arise from trajectories which will later leave the progress set,
                // then we only need to look at the maximum value of the guard.
                HybridEnclosure eventually_hitting_set=set;
                eventually_hitting_set.new_parameter_constraint( event, elapsed_time <= critical_time );
                eventually_hitting_set.new_parameter_constraint( event, maximal_guard >= zero);
                if(definitely(eventually_hitting_set.empty())) {
                    set.new_parameter_constraint(event, maximal_guard <= zero);
                    break;
                }
                // If no points in the set arise from trajectories which leave the progress set and
                // later return, then we only need to look at the guard at the final value
                HybridEnclosure returning_set=set;
                returning_set.new_parameter_constraint( event, elapsed_time >= critical_time );
                returning_set.new_parameter_constraint( event, final_guard <= zero );
                if(definitely(returning_set.empty())) {
                    set.new_parameter_constraint(event, final_guard <= zero);
                    break;
                }
                // Split the set into two components, one corresponding to
                // points which miss the guard completely, the other to points which
                // eventually hit the guard, ensuring that the set stops at the first crossing
                HybridEnclosure extra_set=set;
                set.new_parameter_constraint(event,maximal_guard<=zero);
                extra_set.new_parameter_constraint( event, elapsed_time <= critical_time );
                extra_set.new_state_constraint(event,guard_function<=zero);
                sets.append(extra_set);
                break;
                // Code below is always exact, but uses two sets
                // set1.new_parameter_constraint(event,final_guard <= zero);
                // set2.new_parameter_constraint(event, elapsed_time <= critical_time);
                // set1.new_parameter_constraint(event, maximal_guard <= zero);
                // set2.new_parameter_constraint(event, elapsed_time >= critical_time);
            }
            case CrossingKind::DEGENERATE: case CrossingKind::CONCAVE: {
                // The crossing with the guard set is not one of the kinds handled above.
                // We obtain an over-approximation by testing at finitely many time points
                const uint n=SUBDIVISIONS_FOR_DEGENERATE_CROSSING;
                switch(semantics) {
                    case UPPER_SEMANTICS:
                        for(uint i=0; i!=n; ++i) {
                            ValidatedFloat alpha=ExactFloat(i+1)/n;
                            ValidatedScalarFunctionModel intermediate_guard
                                = compose( guard_function, unchecked_compose( flow, join(starting_state, alpha*elapsed_time) ) );
                            set.new_parameter_constraint(event, intermediate_guard <= zero);
                        }
                        break;
                    case LOWER_SEMANTICS:
                        // Can't continue the evolution, so set a trivially-falsified constraint
                        set.new_parameter_constraint(event, this->function_factory().create_constant(set.parameter_domain(),1) <= zero);
                        break;
                }
                break;
            }
            case CrossingKind::POSITIVE:
                // No need to do anything since all points are initially
                // active and should have been handled already
            case CrossingKind::NEGATIVE:
                // No points are active
            case CrossingKind::DECREASING:
                // No need to do anything, since only initially active points
                // become active during the evolution, and these have been
                // handled already.
                break;
            default:
                ARIADNE_FAIL_MSG("Unhandled crossing "<<crossing_data<<"\n");
        }
    }
}



void
HybridEvolverBase::
_evolution_in_mode(EvolutionData& evolution_data,
                   HybridTerminationCriterion const& termination_criterion) const
{

    //  Select a working set and evolve this in the current location until either
    // all initial points have undergone a discrete transition (possibly to
    // the same location) or the final time is reached.
    //   Evolving within one location avoids having to re-extract event sets,
    // and means that initially active events are tested for only once.
    ARIADNE_LOG(3,"HybridEvolverBase::_evolution_in_mode\n");

    typedef Map<DiscreteEvent,ValidatedScalarFunction>::const_iterator constraint_iterator;
    typedef Set<DiscreteEvent>::const_iterator event_iterator;

    const Real final_time=termination_criterion.maximum_time();
    const Natural maximum_steps=termination_criterion.maximum_steps();
    const Set<DiscreteEvent>& terminating_events=termination_criterion.terminating_events();

    // Routine check for emptiness
    if(evolution_data.initial_sets.empty()) { return; }

    // Get the initial set for this round of evolution
    HybridEnclosure initial_set=evolution_data.initial_sets.back(); evolution_data.initial_sets.pop_back();
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n\n");

    if(initial_set.time_range().lower()>=final_time) {
        ARIADNE_WARN("initial_set.time_range()="<<initial_set.time_range()<<" which exceeds final time="<<final_time<<"\n");
        return;
    }

    // Extract starting location
    const DiscreteLocation location=initial_set.location();

    // Cache dynamic and constraint functions
    EffectiveVectorFunction dynamic=_sys_ptr->dynamic_function(location);
    Map<DiscreteEvent,TransitionData> transitions = this->_extract_transitions(location);
    Set<DiscreteEvent> events = transitions.keys();

    ARIADNE_LOG(4,"dynamic="<<dynamic<<"\n");
    ARIADNE_LOG(4,"transitions="<<transitions<<"\n\n");

    // Test if maximum number of steps has been exceeded; if so, the set should be discarded.
    // NOTE: We could also place a test for the maximum number of steps being reaches which computing jump sets
    // This is not done since the maximum_steps information is not passed to the _apply_evolution_step(...) method.
    if(initial_set.previous_events().size()>=maximum_steps) {
        ARIADNE_LOG(4,"initial_set "<<initial_set<<" has undergone more than maximum number of events "<<maximum_steps<<"\n");
        this->_apply_invariants(initial_set,transitions);
        evolution_data.final_sets.append(initial_set);
        return;
    }

    // Test if a terminating event has been reached.
    if(initial_set.previous_events().size()>=1 && terminating_events.contains(initial_set.previous_events().back())) {
        ARIADNE_LOG(4,"initial_set "<<initial_set<<" has undergone event "<<initial_set.previous_events().back()<<"\n");
        this->_apply_invariants(initial_set,transitions);
        evolution_data.final_sets.append(initial_set);
        return;
    }

    // Process the initially active events; cut out active points to leave initial flowable set.
    this->_process_starting_events(evolution_data, initial_set,transitions);
    ARIADNE_ASSERT(evolution_data.working_sets.size()<=1);

    // NOTE: Uncomment the lines below to stop evolution immediately after the maximum event, without further flow
    //if(initial_set.previous_events().size()>maximum_steps) {
    //    evolution_data.final_sets.append(initial_set);
    //    return;
    //}

    while(!evolution_data.working_sets.empty()) {
        this->_evolution_step(evolution_data,dynamic,transitions,final_time);
    }
}

void
HybridEvolverBase::
_evolution_step(EvolutionData& evolution_data,
                EffectiveVectorFunction const& dynamic,
                Map<DiscreteEvent,TransitionData> const& transitions,
                Real const& final_time) const
{
    ARIADNE_LOG(3,"HybridEvolverBase::_evolution_step\n");
    HybridEnclosure starting_set=evolution_data.working_sets.back(); evolution_data.working_sets.pop_back();

    ARIADNE_LOG(4,"starting_set="<<starting_set<<"\n");
    ARIADNE_LOG(4,"starting_time="<<starting_set.time_function()<<"\n");
    if(definitely(starting_set.empty())) {
        ARIADNE_LOG(4,"Empty starting_set "<<starting_set<<"\n");
        return;
    }

    if(starting_set.time_range().lower() >= final_time) {
        ARIADNE_WARN("starting_set.time_range()="<<starting_set.time_range()<<" which exceeds final time="<<final_time<<"\n");
        return;
    }

    if(verbosity==1 || verbosity==2) { _log_summary(evolution_data,starting_set); }
    ARIADNE_LOG(2,"starting_set: bounding_box="<<starting_set.bounding_box()<<", time_range="<<starting_set.time_function().range()<<"    \n");


    // Compute the bounding box of the enclosure
    const ExactBox starting_bounding_box=make_exact_box(starting_set.space_bounding_box());
    ARIADNE_LOG(4,"starting_bounding_box="<<starting_bounding_box<<"\n");

    // Test to see if set requires reconditioning
    if(this->_configuration_ptr->enable_reconditioning() &&
            norm(starting_set.space_function().errors()) > this->_configuration_ptr->maximum_spacial_error()) {
        HybridEnclosure reconditioned_set=starting_set;
        reconditioned_set.recondition();
        evolution_data.working_sets.append(reconditioned_set);
        return;
    }

    // Handle a set that is too large, based on semantics
    if (starting_bounding_box.radius() > this->_configuration_ptr->maximum_enclosure_radius()) {
        if (evolution_data.semantics == LOWER_SEMANTICS) {
            ARIADNE_LOG(1,"\r  too large, discarding\n");
            return;
        } else if (this->_configuration_ptr->enable_subdivisions()) {
            ARIADNE_LOG(1,"\r  too large, splitting\n");
            List<HybridEnclosure> split_sets = starting_set.split();
            for(uint i=0; i!=split_sets.size(); ++i) {
                if(!definitely(split_sets[i].empty())) { evolution_data.working_sets.append(split_sets[i]); }
            }
            return;
        }
    }

    Map<DiscreteEvent,EffectiveScalarFunction> guard_functions;
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        guard_functions.insert(transition_iter->first,transition_iter->second.guard_function);
    }
    ARIADNE_LOG(4,"guards="<<guard_functions<<"\n");

    // Compute flow and actual time step size used
    const FlowFunctionModel flow_model=this->_compute_flow(dynamic,starting_bounding_box,ExactFloat(this->configuration().maximum_step_size()));
    ARIADNE_LOG(4,"flow_model.domain()="<<flow_model.domain()<<" flow_model.range()="<<flow_model.range()<<"\n");

    // Compute possibly active urgent events with increasing guards, and crossing times
    Set<DiscreteEvent> active_events =
        this->_compute_active_events(dynamic,guard_functions,flow_model,starting_set);
    if(!active_events.empty()) {
        Map<DiscreteEvent,EventKind> active_events_log_data;
        for(Set<DiscreteEvent>::const_iterator event_iter=active_events.begin(); event_iter!=active_events.end(); ++event_iter) {
            active_events_log_data.insert(*event_iter,transitions[*event_iter].event_kind);
        }
        ARIADNE_LOG(2,"active_events: "<<active_events_log_data<<"\n");
    }

    // Compute the kind of crossing (increasing, convex, etc);
    Map<DiscreteEvent,CrossingData> crossings =
        this->_compute_crossings(active_events,dynamic,guard_functions,flow_model,starting_set);
    ARIADNE_LOG(4,"crossings="<<crossings<<"\n");

    // Compute end conditions for flow
    TimingData timing_data = this->_estimate_timing(active_events,Real(final_time),flow_model,crossings,transitions,starting_set);
    ARIADNE_LOG(2,"timing_data: "<<timing_data<<"\n");

    // Apply the time step
    HybridEnclosure reach_set, evolve_set;
    this->_apply_evolution_step(evolution_data,starting_set,flow_model,timing_data,crossings,dynamic,transitions);
}


void HybridEvolverBase::
_apply_evolution_step(EvolutionData& evolution_data,
                      HybridEnclosure const& starting_set,
                      ValidatedVectorFunctionModel const& flow,
                      TimingData const& timing_data,
                      Map<DiscreteEvent,CrossingData> const& crossings,
                      EffectiveVectorFunction const& dynamic,
                      Map<DiscreteEvent,TransitionData> const& transitions) const
{
    ARIADNE_LOG(3,"GeneralHybridEvolver::_apply_evolution_step(...)\n");
    //ARIADNE_LOG(4,"parameter_dependent_evolution_time="<<parameter_dependent_evolution_time.range()<<" final_time="<<final_time<<"\n");

    EvolutionStepData _step_data;
    HybridEnclosure starting_set_copy=starting_set;
    Tribool starting_set_empty=starting_set_copy.empty();

    // Counters for number of sucessor sets
    Nat jump_sets = 0;
    Bool progress = false;

    if(definitely(starting_set_empty)) {
        ExactIntervalVector reduced_domain=starting_set.continuous_set().reduced_domain();
        ARIADNE_WARN("empty starting_set "<<representation(starting_set)<<"\n");
        return;
    }
    //ARIADNE_ASSERT_MSG(!definitely(starting_set.empty()),"starting_set="<<repr(starting_set)<<"\n");

    Semantics semantics = evolution_data.semantics;

    // Compute events enabling transitions and events blocking continuous evolution
    Set<DiscreteEvent> active_events=crossings.keys();
    Set<DiscreteEvent> activating_events=intersection(Ariadne::activating_events(transitions),active_events);
    Set<DiscreteEvent> blocking_events=intersection(Ariadne::blocking_events(transitions),active_events);

    ARIADNE_LOG(6,"activating_events="<<activating_events<<"\n");
    ARIADNE_LOG(6,"blocking_events="<<blocking_events<<"\n");

    ValidatedScalarFunctionModel const& evolve_step_time=timing_data.parameter_dependent_evolution_time;
    ValidatedScalarFunctionModel reach_step_time=embed(starting_set.parameter_domain(),timing_data.evolution_time_coordinate);

    ARIADNE_LOG(8,"evolve_step_time="<<evolve_step_time<<"\n")
    ARIADNE_LOG(8,"reach_step_time="<<reach_step_time<<"\n")

    // Compute the reach and evolve sets, without introducing bounds due to the final time.
    List<HybridEnclosure> reach_sets={starting_set};
    List<HybridEnclosure> evolve_sets={starting_set};

    _apply_reach_step(reach_sets.front(),flow,timing_data);
    _apply_evolve_step(evolve_sets.front(),flow,timing_data);

    ARIADNE_LOG(8,"flow_reach_set="<<reach_sets.front()<<"\n")
    ARIADNE_LOG(8,"flow_evolve_set="<<evolve_sets.front()<<"\n")

    // Apply constraints on reach and evolve sets due to invariants and urgent guards
    for(Set<DiscreteEvent>::const_iterator event_iter=blocking_events.begin();
        event_iter!=blocking_events.end(); ++event_iter)
    {
        const TransitionData& transition_data=transitions[*event_iter];
        const CrossingData& crossing_data=crossings[*event_iter];

        this->_apply_guard(reach_sets,starting_set,flow,reach_step_time,
                           transition_data,crossing_data,semantics);

        switch(crossing_data.crossing_kind) {
            case CrossingKind::INCREASING: case CrossingKind::TRANSVERSE:
                // Delay applying guard until all splittings due to non-increasing crossings have been processed.
                // We can then test emptiness on sets with crossings, but not actually apply crossing immediately,
                // since it will be introduced at the next time step anyway; this avoids duplication
                break;
            default:
                this->_apply_guard(evolve_sets,starting_set,flow,evolve_step_time,
                                   transition_data,crossing_data,semantics);
        }
    }

    // Make copy of evolve sets without transverse and increasing crossings which can be used for future evolution
    List<HybridEnclosure> next_working_sets = evolve_sets;

    for(Set<DiscreteEvent>::const_iterator event_iter=blocking_events.begin();
        event_iter!=blocking_events.end(); ++event_iter)
    {
        const TransitionData& transition_data=transitions[*event_iter];
        const CrossingData& crossing_data=crossings[*event_iter];
        switch (crossing_data.crossing_kind) {
            case CrossingKind::INCREASING: case CrossingKind::TRANSVERSE:
                this->_apply_guard(evolve_sets,starting_set,flow,timing_data.parameter_dependent_evolution_time,
                                   transition_data,crossing_data,semantics);
                break;
            default:
                break; // Constraint has already been handled
        }
    }

    // Compute final set depending on whether the finishing kind is exactly AT_FINAL_TIME.
    // Insert sets into evolution_data as appropriate
    if(timing_data.finishing_kind==FinishingKind::AT_FINAL_TIME) {
        for(List<HybridEnclosure>::const_iterator evolve_set_iter=evolve_sets.begin();
            evolve_set_iter!=evolve_sets.end(); ++evolve_set_iter)
        {
            HybridEnclosure const& evolve_set=*evolve_set_iter;
            if(!definitely(evolve_set.empty())) {
                ARIADNE_LOG(4,"final_set="<<evolve_set<<"\n");
                evolution_data.final_sets.append(evolve_set);
                _step_data.finishing = true;
            }
        }
        for(List<HybridEnclosure>::const_iterator reach_set_iter=reach_sets.begin();
            reach_set_iter!=reach_sets.end(); ++reach_set_iter)
        {
            HybridEnclosure const& reach_set=*reach_set_iter;
            evolution_data.reach_sets.append(reach_set);
        }
    } else { // (timing_data.finishing_kind!=AT_FINAL_TIME)
        List<HybridEnclosure>::iterator next_working_set_iter=next_working_sets.begin();
        for(List<HybridEnclosure>::iterator evolve_set_iter=evolve_sets.begin();
            evolve_set_iter!=evolve_sets.end(); ++evolve_set_iter, ++next_working_set_iter )
        {
            HybridEnclosure& evolve_set=*evolve_set_iter;
            HybridEnclosure& next_working_set=*next_working_set_iter;
            UpperInterval evolve_set_time_range=evolve_set.time_range();
            if(evolve_set_time_range.lower()>timing_data.final_time) {
                // Do nothing, since evolve set is past final time is definitely empty
            } else if(evolve_set_time_range.upper()<=timing_data.final_time) {
                // No need to introduce timing constraints
                if(!definitely(evolve_set.empty())) {
                    ARIADNE_LOG(4,"evolve_set="<<evolve_set<<"\n");
                    ARIADNE_LOG(4,"next_working_set="<<next_working_set<<"\n");
                    evolution_data.working_sets.append(next_working_set);
                    evolution_data.intermediate_sets.append(evolve_set);
                    _step_data.progress = true;
                }
            } else {
                // Bound time if necessary
                if(evolve_set_time_range.upper()>timing_data.final_time) {
                    evolve_set.bound_time(timing_data.final_time);
                }
                // Only continue evolution if the time-bounded evolve set is nonempty;
                // However, continue evolution without adding time constraint,
                // since this will be introduced in the next step
                if(!definitely(evolve_set.empty())) {
                    ARIADNE_LOG(4,"evolve_set="<<evolve_set<<"\n");
                    ARIADNE_LOG(4,"next_working_set="<<next_working_set<<"\n");
                    evolution_data.working_sets.append(next_working_set);
                    evolution_data.intermediate_sets.append(evolve_set);
                    _step_data.progress = true;
                }
            }
        }
        for(List<HybridEnclosure>::iterator reach_set_iter=reach_sets.begin();
            reach_set_iter!=reach_sets.end(); ++reach_set_iter)
        {
            HybridEnclosure& reach_set=*reach_set_iter;
            //reach_set.continuous_set().reduce();
            if(reach_set.time_range().upper()>timing_data.final_time) {
                HybridEnclosure final_set=reach_set;
                final_set.set_time(timing_data.final_time);
                if(timing_data.finishing_kind!=FinishingKind::BEFORE_FINAL_TIME && !definitely(final_set.empty())) {
                //if(!definitely(final_set.empty())) {
                    ARIADNE_LOG(4,"final_set="<<final_set<<"\n");
                    evolution_data.final_sets.append(final_set);
                    _step_data.finishing = true;
                }
            }
            reach_set.bound_time(timing_data.final_time);
            evolution_data.reach_sets.append(reach_set);
        }
    }




    // Compute jump sets
    for(Set<DiscreteEvent>::const_iterator event_iter=activating_events.begin();
        event_iter!=activating_events.end(); ++event_iter)
    {
        DiscreteEvent event=*event_iter;


        // Compute active set
        List<HybridEnclosure> jump_sets={starting_set};
        HybridEnclosure& jump_set=jump_sets.front();
        ARIADNE_LOG(3,"  "<<event<<": "<<transitions[event].event_kind<<", "<<crossings[event].crossing_kind<<"\n");
        _apply_guard_step(jump_set,dynamic,flow,timing_data,transitions[event],crossings[event],semantics);
        ValidatedScalarFunctionModel jump_step_time;
        switch(crossings[event].crossing_kind) {
            case CrossingKind::INCREASING: case CrossingKind::CONVEX: case CrossingKind::CONCAVE: case CrossingKind::DEGENERATE:
                jump_step_time=reach_step_time;
                break;
            case CrossingKind::GRAZING:
                jump_step_time=unchecked_compose(crossings[event].critical_time,starting_set.space_function());
                break;
            case CrossingKind::TRANSVERSE:
                jump_step_time=unchecked_compose(crossings[event].crossing_time,starting_set.space_function());
                break;
            case CrossingKind::DECREASING: case CrossingKind::POSITIVE: case CrossingKind::NEGATIVE:
                // No need to set jump set time;
                break;
            default:
                ARIADNE_FAIL_MSG("Unknown crossing kind in "<<crossings[event]);
        }

        // Apply maximum time bound, as this will be applied after the next flow step
        for(List<HybridEnclosure>::iterator jump_set_iter=jump_sets.begin(); jump_set_iter!=jump_sets.end(); ++jump_set_iter) {
            if(jump_set.time_range().upper()>timing_data.final_time) {
                jump_set.bound_time(timing_data.final_time);
                if(definitely(not jump_set.empty())) { ARIADNE_WARN("Explicitly bounding time in jump set\n"); }
            }
        }

        // Apply blocking conditions for other active events
        for(Set<DiscreteEvent>::const_iterator other_event_iter=blocking_events.begin();
            other_event_iter!=blocking_events.end(); ++other_event_iter)
        {
            DiscreteEvent other_event=*other_event_iter;
            if(other_event!=event) {
                const TransitionData& other_transition_data=transitions[other_event];
                const CrossingData& other_crossing_data=crossings[other_event];
                _apply_guard(jump_sets,starting_set,flow,jump_step_time,
                             other_transition_data,other_crossing_data,semantics);
            }
        }

        ARIADNE_LOG(3, "  "<<event<<": "<<transitions[event].event_kind<<", "<<crossings[event].crossing_kind<<"\n");
        // Apply reset
        for(List<HybridEnclosure>::iterator jump_set_iter=jump_sets.begin(); jump_set_iter!=jump_sets.end(); ++jump_set_iter) {
            HybridEnclosure& jump_set=*jump_set_iter;
            if(!definitely(jump_set.empty())) {
                jump_set.apply_reset(event,transitions[event].target,transitions[event].target_space,transitions[event].reset_function);
                evolution_data.initial_sets.append(jump_set);
                _step_data.events.insert(event);
                ARIADNE_LOG(6, "jump_set="<<jump_set<<"\n");
            }
        }
    }

    if(verbosity>=1 && (_step_data.finishing || !_step_data.progress || !_step_data.events.empty()) ) {
        ARIADNE_LOG(3,timing_data.finishing_kind<<" "<<_step_data.events<<(_step_data.progress?" progress":" no progress")<<"\n");
    }

}





TimingData
HybridEvolverBase::
_estimate_timing(Set<DiscreteEvent>& active_events,
                Real final_time,
                FlowFunctionModel const& flow,
                Map<DiscreteEvent,CrossingData>& crossings,
                Map<DiscreteEvent,TransitionData> const& transitions,
                HybridEnclosure const& initial_set) const
{
    // Compute the evolution time for the given step.
    ARIADNE_LOG(7,"HybridEvolverBase::_estimate_timing(...)\n");
    const ExactFloat step_size=static_cast<ExactFloat>(flow.domain()[flow.domain().size()-1].upper());
    TimingData result;
    result.step_kind=StepKind::CONSTANT_EVOLUTION_TIME;
    result.finishing_kind=FinishingKind::STRADDLE_FINAL_TIME;
    result.step_size=step_size;
    result.final_time=final_time;
    result.evolution_time_domain=ExactInterval(zero,step_size);
    result.evolution_time_coordinate=this->function_factory().create_identity(result.evolution_time_domain);
    result.parameter_dependent_evolution_time=this->function_factory().create_constant(initial_set.parameter_domain(),ExactFloat(result.step_size));
    ARIADNE_LOG(8,"  timing_data="<<result<<"\n");
    return result;
}


HybridEvolverBaseConfiguration::HybridEvolverBaseConfiguration(HybridEvolverBase& evolver)
    : _evolver(evolver)
{
    set_flow_accuracy(1e-5);
    set_maximum_step_size(1.0);
    set_maximum_enclosure_radius(100.0);
    set_maximum_spacial_error(1e-2);
    set_enable_reconditioning(true);
    set_enable_subdivisions(true);
}

void
HybridEvolverBaseConfiguration::set_flow_accuracy(const RealType value)
{
    _evolver._integrator_ptr=std::shared_ptr<TaylorSeriesIntegrator>(new TaylorSeriesIntegrator(MaximumError(value)));
    _flow_accuracy = value;
}


std::ostream&
HybridEvolverBaseConfiguration::write(std::ostream& os) const
{
    os << "HybridEvolverBaseConfiguration"
       << ",\n  flow_accuracy=" << flow_accuracy()
       << ",\n  maximum_step_size=" << maximum_step_size()
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n  maximum_spacial_error=" << maximum_spacial_error()
       << ",\n  enable_reconditioning=" << enable_reconditioning()
       << ",\n  enable_subdivisions=" << enable_subdivisions()
       << "\n)\n";
    return os;
}

GeneralHybridEvolver::GeneralHybridEvolver(const SystemType& system)
    : HybridEvolverBase(system)
{
    this->_configuration_ptr.reset(new GeneralHybridEvolverConfiguration(*this));
}


GeneralHybridEvolver::GeneralHybridEvolver(
        const SystemType& system,
        const ValidatedFunctionModelFactoryInterface& factory)
    : HybridEvolverBase(system,factory)
{
    this->_configuration_ptr.reset(new GeneralHybridEvolverConfiguration(*this));
}


TimingData
GeneralHybridEvolver::
_estimate_timing(Set<DiscreteEvent>& active_events,
                Real final_time,
                FlowFunctionModel const& flow,
                Map<DiscreteEvent,CrossingData>& crossings,
                Map<DiscreteEvent,TransitionData> const& transitions,
                HybridEnclosure const& initial_set) const
{
    // Compute the evolution time for the given step.
    ARIADNE_LOG(7,"GeneralHybridEvolver::_estimate_timing(...)\n");

    const uint n = flow.result_size();
    const ExactFloat step_size=flow.domain()[flow.domain().size()-1].upper();

    TimingData result;
    result.step_size=flow.step_size();
    result.final_time=final_time;

    ExactBox space_domain = make_exact_box(initial_set.space_bounding_box());
    ExactInterval time_domain = make_exact_interval(initial_set.time_range()+ExactInterval(zero,step_size));
    ExactBox spacetime_domain = join(space_domain,time_domain);

    //ValidatedVectorFunctionModel space_coordinates=this->function_factory().create_identity(space_domain);
    ValidatedScalarFunctionModel time_coordinate=this->function_factory().create_coordinate(spacetime_domain,n);
    ValidatedScalarFunctionModel time_identity=this->function_factory().create_identity(time_domain);

    result.evolution_time_domain=ExactInterval(zero,step_size);
    result.evolution_time_coordinate=this->function_factory().create_identity(result.evolution_time_domain);

    ExactBox flow_space_domain = ExactBox(project(flow.domain(),range(0,n)));
    if(!subset(space_domain,flow_space_domain)) {
        ARIADNE_WARN(std::setprecision(17)<<"Bounding box "<<space_domain<<" is not subset of the flow spacial domain "<<flow_space_domain<<"\n");
        space_domain=hull(space_domain,flow_space_domain);
    }

    // NOTE: The starting time function may be negative or greater than the final time
    // over part of the parameter domain.
    ValidatedVectorFunctionModel const& starting_space_function=initial_set.space_function();
    ValidatedScalarFunctionModel const& starting_time_function=initial_set.time_function();
    UpperInterval starting_time_range=initial_set.time_range();
    UpperInterval remaining_time_range=static_cast<UpperInterval>(final_time)-starting_time_range;

    ARIADNE_LOG(7,std::fixed<<"starting_time_range="<<starting_time_range<<" step_size="<<step_size<<" final_time="<<final_time<<"\n");


    // The time-dependent part of the evolution time
    ValidatedScalarFunctionModel temporal_evolution_time=this->function_factory().create_zero(ExactIntervalVector(1u,time_domain));

    if(remaining_time_range.lower()<zero) {
        // Some of the points may already have reached the final time.
        // Don't try anything fancy, just do a simple constant time step.
        if(remaining_time_range.upper()<=step_size) {
            if(false) {
                // Within one time step we can go beyond final time
                result.step_kind=StepKind::CONSTANT_EVOLUTION_TIME;
                result.finishing_kind=FinishingKind::AFTER_FINAL_TIME;
                temporal_evolution_time=ValidatedFloat(step_size); //   remaining_time_range.upper();
            } else {
                result.step_kind=StepKind::CONSTANT_FINISHING_TIME;
                result.finishing_kind=FinishingKind::AT_FINAL_TIME;
                temporal_evolution_time=ValidatedFloat(final_time)-time_identity;
            }
        } else {
            result.step_kind=StepKind::CONSTANT_EVOLUTION_TIME;
            result.finishing_kind=FinishingKind::STRADDLE_FINAL_TIME;
            temporal_evolution_time=ValidatedFloat(step_size);
        }
    } else if(remaining_time_range.upper()<=result.step_size) {
        // The rest of the evolution can be computed within a single time step.
        // The finishing kind is given as AT_FINAL_TIME so that the evolution algorithm
        // knows that the evolved set does not need to be evolved further.
        // This knowledge is required to be given combinarially, since
        // specifying the final time as a constant Function is not
        // exact if the final_time parameter is not exactly representable as
        // a ExactFloat
        result.step_kind=StepKind::CONSTANT_FINISHING_TIME;
        result.finishing_kind=FinishingKind::AT_FINAL_TIME;
        temporal_evolution_time=ValidatedFloat(final_time)-time_identity;
    } else if(remaining_time_range.lower()<=step_size && ALLOW_CREEP) {
        // Some of the evolved points can be evolved to the final time in a single step
        // The evolution is performed over a step size which moves points closer to the final time, but does not cross.

        // Using the final_time as a guide, set the finishing time to closer to the final time.
        // This method ensures that points do not pass the final time after the transition.
        result.step_kind=StepKind::SPACETIME_DEPENDENT_FINISHING_TIME;
        result.finishing_kind=FinishingKind::BEFORE_FINAL_TIME;
        ExactFloat sf=1;
        while(remaining_time_range.upper()*sf>step_size) { sf = half(sf); }
        temporal_evolution_time= ValidatedFloat(sf)*(ValidatedFloat(final_time)-time_identity);
    } else { // remaining_time_range.lower()>step_size)
        // As far as timing goes, perform the evolution over a full time step
        result.step_kind=StepKind::CONSTANT_EVOLUTION_TIME;
        result.finishing_kind=FinishingKind::BEFORE_FINAL_TIME;
        temporal_evolution_time=ExactFloat(result.step_size);
    }

    ARIADNE_LOG(7,"finishing_kind="<<result.finishing_kind<<"\n");
    ARIADNE_LOG(7,"temporal_evolution_time="<<temporal_evolution_time<<"\n");


    ValidatedScalarFunctionModel spacial_evolution_time=this->function_factory().create_constant(space_domain,ExactFloat(step_size));

    // Select one of GUARD_CREEP or TIME_CREEP
    static const bool GUARD_CREEP=true;
    static const bool TIME_CREEP=false;

    bool creep=false;

    // Test for creep step
    if(ALLOW_CREEP && !crossings.empty() && (result.finishing_kind == FinishingKind::BEFORE_FINAL_TIME || result.finishing_kind == FinishingKind::STRADDLE_FINAL_TIME) ) {
        // If an event is possible, but only some points reach the guard set
        // after a full time step, then terminating the evolution here will
        // cause a splitting of the enclosure set. To prevent this, attempt
        // to "creep" up to the event guard boundary, so that in the next step,
        // all points can be made to cross.
        // Test to see if a creep step is needed
        ARIADNE_LOG(6,"Possible creep step; h="<<step_size<<"\n");

        HybridEnclosure evolve_set=initial_set;
        evolve_set.apply_fixed_evolve_step(flow,flow.step_size());

        for(Map<DiscreteEvent,CrossingData>::iterator crossing_iter=crossings.begin();
            crossing_iter!=crossings.end(); ++crossing_iter)
        {
            DiscreteEvent event=crossing_iter->first;
            EventKind event_kind=transitions[event].event_kind;
            CrossingKind crossing_kind=crossing_iter->second.crossing_kind;
            ARIADNE_LOG(6,"  Event "<<event<<": "<<event_kind<<": "<<crossing_kind<<"\n");
            if(event_kind!=PERMISSIVE) {
                evolve_set.new_invariant(event,transitions[event].guard_function);
            }
            // FIXME: When using permissive crossings, jumps in the step after the crossing time are lost.
            // A hack to fix this is to only creep on non-permissive events. Check that evolution is correct in this case.
            // FIXME: What should we do on increasing but non-transverse crossings?
            if((crossing_kind==CrossingKind::TRANSVERSE ) // || crossing_kind==CrossingKind::INCREASING)
                    && event_kind!=PERMISSIVE)
            {
                ARIADNE_LOG(6,"  crossing_time_range="<<crossing_iter->second.crossing_time.range()<<"\n");
                const ValidatedScalarFunctionModel& crossing_time=crossing_iter->second.crossing_time;
                UpperInterval crossing_time_range=crossing_time.range();
                if(Ariadne::is_blocking(event_kind) && crossing_time_range.upper()<step_size) {
                    // NOTE: Use strict comparison here so that guard is fully crossed
                    // In principle this is not necessary, but this would involve testing
                    // to ensure that the evolved set is not propagated

                    // This event ensures that the evolve set is empty after a full step, so use this.
                    ARIADNE_LOG(6,std::setprecision(18)<<"crossing_time_range="<<crossing_time_range<<", crossing_time_range .upper()="<<crossing_time_range.upper()<<", step_size="<<step_size<<"\n");
                    EffectiveScalarFunction guard=transitions[event].guard_function;
                    ValidatedVectorFunctionModel identity=this->function_factory().create_identity(crossing_time.domain());
                    ValidatedScalarFunctionModel step_time=crossing_time*zero+ExactFloat(step_size);
                    ARIADNE_LOG(6,"full flow="<<compose(flow,join(identity,step_time))<<"\n");
                    ARIADNE_LOG(6,"guard range at crossing time="<<compose(guard,compose(flow,join(initial_set.space_function(),compose(crossing_time,initial_set.space_function())))).range()<<"\n");
                    ARIADNE_LOG(6,"guard range at crossing time="<<compose(guard,compose(flow,join(identity,crossing_time))).range());
                    ARIADNE_LOG(6,"No creep; event "<<crossing_iter->first<<" completely taken\n");
                    creep=false;
                    break;
                } else if(crossing_time_range.lower()<=zero) {
                    // This event is already partially active, so carry on with a full step
                    ARIADNE_LOG(6,"No creep; event "<<crossing_iter->first<<" already partially active\n");
                    creep=false;
                    break;
                } else if(crossing_time_range.lower()>=step_size) {
                    ARIADNE_LOG(6,"Event "<<crossing_iter->first<<" is not actually reached\n");
                    //crossings.erase(pending_erase_crossing_iter);
                } else {
                    ARIADNE_LOG(6,"Event "<<crossing_iter->first<<" can be crept up to.\n");
                    creep=true;
                }
            }
        }
        // If the evolved set is definitely empty, no creeping occurs
        if(definitely(evolve_set.empty())) {
            ARIADNE_LOG(6,"No creep; evolve set is empty");
            creep=false;
        }
    }

    if(ALLOW_CREEP && creep==false) { ARIADNE_LOG(6,"No creep\n"); }

    if(TIME_CREEP && creep==true) {
        // Compute reduced evolution time; note that for every remaining increasing crossing, we have
        // crossing_time_range.lower()>zero and crossing_time_range.upper()>step_size.

        for(Map<DiscreteEvent,CrossingData>::const_iterator crossing_iter=crossings.begin();
            crossing_iter!=crossings.end(); ++crossing_iter)
        {
            if(crossing_iter->second.crossing_kind==CrossingKind::TRANSVERSE) {
                // Modify the crossing time function to be the smallest possible; this ensures that the evaluation time is
                // essentially exact
                ValidatedScalarFunctionModel lower_crossing_time=crossing_iter->second.crossing_time;
                ErrorFloat crossing_time_error=lower_crossing_time.error();
                lower_crossing_time.set_error(0u);
                lower_crossing_time-=ExactFloat(crossing_time_error.raw());

                // One possibility is to use quadratic restrictions
                //   If 0<=x<=2h, then x(1-x/4h)<=min(x,h)
                //   If 0<=x<=4h, then x(1-x/8h)<=min(x,2h)
                // Formula below works if x<=2h
                //   evolution_time=evolution_time*(crossing_time/result.step_size)*(1.0-crossing_time/(4*result.step_size));

                // Prefer simpler linear restrictions.
                // Multiply evolution time by crossing_time/max_crossing_time
                spacial_evolution_time=spacial_evolution_time*lower_crossing_time/make_exact(lower_crossing_time.range().upper());
            }
        }
        // Erase increasing transverse crossings since these cannot occur
        for(Map<DiscreteEvent,CrossingData>::iterator crossing_iter=crossings.begin();
            crossing_iter!=crossings.end(); )
        {
            if(crossing_iter->second.crossing_kind==CrossingKind::TRANSVERSE) {
                crossings.erase(crossing_iter++);
            } else {
                ++crossing_iter;
            }
        }
        ARIADNE_LOG(6,"Creep step: spacial_evolution_time="<<spacial_evolution_time<<"\n");
        result.step_kind=StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME;
        result.finishing_kind=FinishingKind::BEFORE_FINAL_TIME;
    }

    if(GUARD_CREEP && creep==true) {
        // If an event is possible, but only some points reach the guard set
        // after a full time step, then terminating the evolution here will
        // cause a splitting of the enclosure set. To prevent this, attempt
        // to "creep" up to the event guard boundary, so that in the next step,
        // all points can be made to cross.
        // Test to see if a creep step is needed
        ARIADNE_LOG(6,"Possible creep step; h="<<step_size<<"\n");
        ARIADNE_LOG(6,"crossings="<<crossings<<"\n");

        const SolverInterface& solver=*this->_solver_ptr;

        HybridEnclosure evolve_set=initial_set;
        evolve_set.apply_fixed_evolve_step(flow,flow.step_size());
        EffectiveVectorFunction dynamic=this->_sys_ptr->dynamic_function(initial_set.location());
        ExactIntervalVector flow_spacial_domain=project(flow.domain(),range(0,flow.argument_size()-1u));
        ExactInterval flow_time_domain=flow.domain()[flow.argument_size()-1u];
        ValidatedScalarFunctionModel zero_function=flow[0].create_zero();
        ValidatedVectorFunctionModel identity_function=flow.create_identity();
        ValidatedVectorFunctionModel space_projection=flow*zero;
        for(uint i=0; i!=n; ++i) { space_projection[i]=space_projection[i]+identity_function[i]; }

        //static const ExactFloat CREEP_MAXIMUM=ExactFloat(1.0);
        static const ExactFloat CREEP_MAXIMUM=ExactFloat(15.0/16);
        spacial_evolution_time=this->function_factory().create_constant(flow.space_domain(),flow.step_size()*CREEP_MAXIMUM);

        for(Map<DiscreteEvent,CrossingData>::iterator crossing_iter=crossings.begin();
            crossing_iter!=crossings.end(); )
        {
            DiscreteEvent event=crossing_iter->first;
            EventKind event_kind=transitions[event].event_kind;
            CrossingKind crossing_kind=crossing_iter->second.crossing_kind;
            EffectiveScalarFunction guard_function=transitions[event].guard_function;
            ARIADNE_LOG(6,"  Event "<<event<<": "<<event_kind<<": "<<crossing_kind<<"\n");
            if(event_kind!=PERMISSIVE) {
                UpperInterval guard_range = compose(guard_function,flow).range();
                ARIADNE_ASSERT(guard_range.lower()<zero);
                UpperInterval guard_derivative_range = compose(lie_derivative(guard_function,dynamic),flow).range();

                //ExactFloat alpha=numeric_cast<ExactFloat>(1+flow.step_size()*guard_derivative_range.lower()/guard_range.lower());
                ValidatedFloat alpha_val=(1+flow.step_size()*make_exact(guard_derivative_range.lower())/make_exact(guard_range.lower()));
                ExactFloat alpha=make_exact(alpha_val);
                assert(alpha_val.value()==alpha);
                ARIADNE_LOG(6,"  step_size: "<<flow.step_size()<<", guard_range: "<<guard_range<<", guard_derivative_range: "<<guard_derivative_range<<", alpha: "<<alpha<<"\n");
                if(alpha>0 && alpha<=1) {
                    ValidatedScalarFunctionModel guard_creep_time;
                    bool sucessfully_computed_guard_creep_time=false;
                    try {
                        guard_creep_time=solver.implicit(compose(guard_function,flow)-alpha*compose(guard_function,space_projection),
                                                        flow_spacial_domain,flow_time_domain);
                        ARIADNE_LOG(6,"  guard_creep_time= "<<guard_creep_time<<"\n");
                        ARIADNE_LOG(6,"  guard_creep_time.range()="<<guard_creep_time.range()<<"\n");
                        sucessfully_computed_guard_creep_time=true;
                        ARIADNE_LOG(9,"  sucessfully_computed_guard_creep_time="<<sucessfully_computed_guard_creep_time<<"\n");
                    }
                    catch(...) {
                        ARIADNE_LOG(6,"  Error in computing guard creep time\n");
                    }
                    if(sucessfully_computed_guard_creep_time) {
                        spacial_evolution_time = spacial_evolution_time * (guard_creep_time/flow.step_size());
                        ARIADNE_LOG(9,"  spacial_evolution_time="<<spacial_evolution_time<<"\n");
                        ARIADNE_LOG(9,"  crossings before erasing="<<crossings<<"\n");
                        crossings.erase(crossing_iter++);
                        ARIADNE_LOG(9,"  crossings after erasing="<<crossings<<"\n");
                    } else {
                      ++crossing_iter;
                    }
                } else {
                  ++crossing_iter;
                }
            } else {
              ++crossing_iter;
            }
        }
        spacial_evolution_time.set_error(0u);


        ARIADNE_LOG(6,"Creep step: spacial_evolution_time="<<spacial_evolution_time<<"\n");
        ARIADNE_LOG(6,"  spacial_evolution_time.range()="<<spacial_evolution_time.range()<<"\n");
        ARIADNE_LOG(6,"  remaining crossings="<<crossings<<"\n");
        result.step_kind=StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME;
        result.finishing_kind=FinishingKind::BEFORE_FINAL_TIME;
    }


    ValidatedScalarFunctionModel evolution_time = embed(spacial_evolution_time,time_domain) * embed(space_domain,temporal_evolution_time/ExactFloat(step_size));
    ValidatedScalarFunctionModel finishing_time=evolution_time+time_coordinate;

    ARIADNE_LOG(7,"evolution_time="<<(evolution_time)<<"\n");
    ARIADNE_LOG(7,"finishing_time="<<(finishing_time)<<"\n");
    result.spacetime_dependent_evolution_time=evolution_time;
    result.spacetime_dependent_finishing_time=finishing_time;
    result.parameter_dependent_evolution_time=unchecked_compose(evolution_time,join(starting_space_function,starting_time_function));
    result.parameter_dependent_finishing_time=unchecked_compose(finishing_time,join(starting_space_function,starting_time_function));

    // Test to see if it is possible to unwind crossings
    if(this->ALLOW_UNWIND && crossings.empty() && (starting_time_range.lower()<starting_time_range.upper())
            && result.step_kind==StepKind::CONSTANT_EVOLUTION_TIME && result.finishing_kind==FinishingKind::BEFORE_FINAL_TIME) {
        ARIADNE_LOG(6,"Possible unwind step; starting_time_range="<<starting_time_range<<", step_size="<<step_size<<"\n");
        // Try to unwind the evolution time to a constant
        result.step_kind=StepKind::PARAMETER_DEPENDENT_FINISHING_TIME;
        result.finishing_kind=FinishingKind::BEFORE_FINAL_TIME;
        if(starting_time_range.width()*2<step_size) {
            result.parameter_dependent_finishing_time=this->function_factory().create_constant(initial_set.parameter_domain(),make_exact(starting_time_range.lower())+step_size);
        } else {
            // Try to reduce the time interval by half the step size
            // Corresponds to setting omega(smin)=tau(smin)+h, omega(smax)=tau(smax)+h/2
            // Taking omega(s)=a tau(s) + b, we obtain
            //   a=1-h/2(tmax-tmin);  b=h(tmax-tmin/2)/(tmax-tmin) = (2tmax-tmin)a
            ExactFloat h=result.step_size;
            ExactFloat tmin=make_exact(starting_time_range.lower());
            ExactFloat tmax=make_exact(starting_time_range.upper());
            ValidatedFloat a=1-(half(h)/(tmax-tmin));
            ValidatedFloat b=h*(tmax-half(tmin))/(tmax-tmin);
            result.parameter_dependent_finishing_time=a*starting_time_function+b;
        }
        ARIADNE_LOG(7,"Unwinding to time "<<result.parameter_dependent_finishing_time<<"\n");
        result.parameter_dependent_evolution_time=result.parameter_dependent_finishing_time-starting_time_function;
    }

    ARIADNE_LOG(7,"step_kind="<<result.step_kind<<", finishing_kind="<<result.finishing_kind<<"\n");
    ARIADNE_LOG(7,"parameter_dependent_evolution_time="<<result.parameter_dependent_evolution_time<<"\n");
    ARIADNE_LOG(7,"parameter_dependent_finishing_time="<<result.parameter_dependent_finishing_time<<"\n\n");
    return result;
}


GeneralHybridEvolverConfiguration::GeneralHybridEvolverConfiguration(GeneralHybridEvolver& evolver)
    : HybridEvolverBaseConfiguration(evolver)
{
}

GeneralHybridEvolverFactory::GeneralHybridEvolverFactory()
    : _function_factory(make_taylor_function_factory())
{
}

GeneralHybridEvolverFactory::GeneralHybridEvolverFactory(const ValidatedFunctionModelFactoryInterface& factory)
    : _function_factory(factory.clone())
{
}


GeneralHybridEvolver*
GeneralHybridEvolverFactory::create(const HybridAutomatonInterface& system) const
{
    return new GeneralHybridEvolver(system,*_function_factory);
}


} // namespace Ariadne