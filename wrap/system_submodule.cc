/***************************************************************************
 *            system_submodule.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include <iostream>
#include <iomanip>

#include "function.h"
#include "real.h"
#include "formula.h"
#include "hybrid_automaton.h"
#include "hybrid_system.h"
#include "hybrid_time.h"
#include "hybrid_set.h"

#include "utilities.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;

namespace Ariadne {



ScalarFunction discrete_transition_guard(const DiscreteTransition& transition) {
    ARIADNE_ASSERT(transition.activation().result_size()==1);
    return transition.activation()[0];
}

std::map<DiscreteEvent,ScalarFunction> discrete_mode_invariants(const DiscreteMode& mode) {
    std::map<DiscreteEvent,ScalarFunction> scalar_invariants;
    std::map<DiscreteEvent,VectorFunction>const& vector_invariants=mode.invariants();
    
    for(std::map<DiscreteEvent,VectorFunction>::const_iterator iter=vector_invariants.begin();
        iter!=vector_invariants.end(); ++iter)
    {
        ARIADNE_ASSERT(iter->second.result_size()==1);
        scalar_invariants.insert(std::make_pair(iter->first,iter->second[0]));
    }

    return scalar_invariants;
}



template<class T>
struct from_python< Space<T> > {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Space<T> >()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        boost::python::list elements=extract<boost::python::list>(obj_ptr);
        Space<T>* spc_ptr = new (storage) Space<T>();
        for(int i=0; i!=len(elements); ++i) {
            extract<String> xs(elements[i]);
            if(xs.check()) { spc_ptr->append(Variable<T>(xs())); }
            else { Variable<T> v=extract< Variable<T> >(elements[i]); spc_ptr->append(v); }
        }
        data->convertible = storage;
    }
};

template<>
struct from_python< EventSet > {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< EventSet >()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((converter::rvalue_from_python_storage<EventSet>*)data)->storage.bytes;
        boost::python::list elements=extract<boost::python::list>(obj_ptr);
        EventSet* evnts_ptr = new (storage) EventSet();
        for(int i=0; i!=len(elements); ++i) {
            extract<String> xs(elements[i]);
            if(xs.check()) { evnts_ptr->insert(Event(xs())); }
            else { Event e=extract< Event >(elements[i]); evnts_ptr->insert(e); }
        }
        data->convertible = storage;
    }
};

template<class T> uint __hash__(const T&);
template<> uint __hash__<DiscreteEvent>(const DiscreteEvent& e) { return reinterpret_cast<const uint&>(e); }
template<> uint __hash__<DiscreteState>(const DiscreteState& q) { return reinterpret_cast<const uint&>(q); }

RealExpression var(const std::string& s) { return RealExpression(RealVariable(s)); }
RealExpression operator+(const RealVariable& v) { return +RealExpression(v); }
RealExpression operator-(const RealVariable& v) { return -RealExpression(v); }
RealExpression operator+(const RealVariable& v, const RealExpression& e) { return RealExpression(v)+e; }
RealExpression operator-(const RealVariable& v, const RealExpression& e) { return RealExpression(v)-e; }
RealExpression operator*(const RealVariable& v, const RealExpression& e) { return RealExpression(v)*e; }
RealExpression operator/(const RealVariable& v, const RealExpression& e) { return RealExpression(v)/e; }
RealExpression operator+(const RealExpression& e, const RealVariable& v) { return e+RealExpression(v); }
RealExpression operator-(const RealExpression& e, const RealVariable& v) { return e-RealExpression(v); }
RealExpression operator*(const RealExpression& e, const RealVariable& v) { return e*RealExpression(v); }
RealExpression operator/(const RealExpression& e, const RealVariable& v) { return e/RealExpression(v); }

} // namespace Ariadne


namespace Ariadne { int length(const array<std::string>& a) { return a.size(); } }

void export_formula()
{
    implicitly_convertible<Event,EventSet>();

    implicitly_convertible<String,StringExpression>();

    implicitly_convertible<int,IntegerExpression>();
    implicitly_convertible<Integer,IntegerExpression>();
    implicitly_convertible<IntegerVariable,IntegerExpression>();

    implicitly_convertible<double,RealExpression>();
    implicitly_convertible<RealVariable,RealExpression>();

    to_python< List<RealExpression> >();

    class_<Event> event_class("Event", init<std::string>());
    event_class.def(self_ns::str(self));

    class_<EventSet> event_set_class("EventSet", init<EventSet>());
    event_set_class.def(init<>());
    event_set_class.def("__invert__", &__not__<EventSet,EventSet>);
    event_set_class.def(self_ns::str(self));

    from_python<EventSet>();


    // TODO: These interval conversions are dangerous since they are applied when they sometimes should not be.
    //implicitly_convertible<double,RealExpression>();
    //implicitly_convertible<Interval,RealExpression>();

    class_<StringVariable> string_variable_class("StringVariable", init<std::string>());
    string_variable_class.def("__eq__", &__eq__<Expression<bool>,StringVariable,std::string>);
    string_variable_class.def("__ne__", &__ne__<Expression<bool>,StringVariable,std::string>);
    string_variable_class.def(self_ns::str(self));

    class_<StringNextVariable> string_next_variable_class("StringNextVariable", no_init);
    string_next_variable_class.def("__lshift__", (StringUpdate(StringNextVariable::*)(const String&)const) &StringNextVariable::operator=);
    string_next_variable_class.def(self_ns::str(self));

    def("next", (StringNextVariable(*)(const StringVariable&)) &next);


    class_<IntegerVariable> integer_variable_class("IntegerVariable", init<std::string>());
    integer_variable_class.def("__lshift__", (IntegerAssignment(IntegerVariable::*)(const IntegerExpression&)const) &IntegerVariable::operator=);
    integer_variable_class.def("__pos__", &__pos__<IntegerExpression,IntegerVariable>);
    integer_variable_class.def("__neg__", &__neg__<IntegerExpression,IntegerVariable>);
    integer_variable_class.def("__add__", &__add__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__sub__", &__sub__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__mul__", &__mul__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__radd__", &__radd__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__rsub__", &__rsub__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__rmul__", &__rmul__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__eq__", &__eq__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__ne__", &__ne__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__le__", &__le__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__ge__", &__ge__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__lt__", &__lt__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__gt__", &__gt__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def(self_ns::str(self));

    class_<IntegerNextVariable> integer_next_variable_class("IntegerNextVariable", no_init);
    integer_next_variable_class.def("__lshift__", (IntegerUpdate(IntegerNextVariable::*)(const IntegerExpression&)const) &IntegerNextVariable::operator=);
    integer_next_variable_class.def(self_ns::str(self));

    def("next", (IntegerNextVariable(*)(const IntegerVariable&)) &next);

    class_<IntegerExpression> integer_expression_class("IntegerExpression", init<IntegerExpression>());
    integer_expression_class.def("__pos__", &__pos__<IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__neg__", &__neg__<IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__add__", &__add__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__sub__", &__sub__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__mul__", &__mul__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__radd__", &__radd__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__rsub__", &__rsub__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__rmul__", &__rmul__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__eq__", &__eq__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__ne__", &__ne__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__le__", &__le__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__ge__", &__ge__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__lt__", &__lt__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__gt__", &__gt__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def(self_ns::str(self));


    class_<RealVariable> real_variable_class("RealVariable", init<std::string>());
    real_variable_class.def("__pos__", &__pos__<RealExpression,RealVariable>);
    real_variable_class.def("__neg__", &__neg__<RealExpression,RealVariable>);
    real_variable_class.def("__add__", &__add__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__sub__", &__sub__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__mul__", &__mul__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__div__", &__div__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__radd__", &__radd__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__rsub__", &__rsub__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__rmul__", &__rmul__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__rdiv__", &__rdiv__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__le__", &__le__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__ge__", &__ge__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__lt__", &__lt__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__gt__", &__gt__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__lshift__", (RealAssignment(RealVariable::*)(const RealExpression&)const) &RealVariable::operator=);
    real_variable_class.def("eq", (RealAssignment(RealVariable::*)(const RealExpression&)const) &RealVariable::operator=);
    real_variable_class.def(self_ns::str(self));

    class_<RealDottedVariable> real_dotted_variable_class("RealDottedVariable", no_init);
    real_dotted_variable_class.def("__lshift__", (RealDynamic(RealDottedVariable::*)(const RealExpression&)const) &RealDottedVariable::operator=);
    real_dotted_variable_class.def(self_ns::str(self));

    def("dot", (RealDottedVariable(*)(const RealVariable&)) &dot);

    class_<RealNextVariable> real_next_variable_class("RealNextVariable", no_init);
    real_next_variable_class.def("__lshift__", (RealUpdate(RealNextVariable::*)(const RealExpression&)const) &RealNextVariable::operator=);
    real_next_variable_class.def(self_ns::str(self));

    def("next", (RealNextVariable(*)(const RealVariable&)) &next);

    class_<RealSpace> real_space_class("RealSpace", init<RealSpace>());
    real_space_class.def("dimension", &RealSpace::dimension);
    real_space_class.def("variable", &RealSpace::variable, return_value_policy<reference_existing_object>());
    real_space_class.def("index", &RealSpace::index);
    real_space_class.def(self_ns::str(self));

    def("variable",(RealVariable(*)(const String& s)) &variable<Real>);
    def("variables",(RealSpace(*)(const List<String>& s)) &variables<Real>);

    from_python<RealSpace>();

    class_<RealExpression> real_expression_class("RealExpression", init<RealExpression>());
    real_expression_class.def("name", &RealExpression::operator_name);
    real_expression_class.def("subexpressions", &RealExpression::subexpressions);
    real_expression_class.def("substitute", &RealExpression::substitute<Real>);
    real_expression_class.def("simplify", &RealExpression::simplify);
    real_expression_class.def("__pos__", &__pos__<RealExpression,RealExpression>);
    real_expression_class.def("__neg__", &__neg__<RealExpression,RealExpression>);
    real_expression_class.def("__add__", &__add__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__sub__", &__sub__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__mul__", &__mul__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__div__", &__div__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__radd__", &__radd__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__rsub__", &__rsub__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__rmul__", &__rmul__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__rdiv__", &__rdiv__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__le__", &__le__<ContinuousPredicate,RealExpression,RealExpression>);
    real_expression_class.def("__ge__", &__ge__<ContinuousPredicate,RealExpression,RealExpression>);
    real_expression_class.def("__lt__", &__lt__<ContinuousPredicate,RealExpression,RealExpression>);
    real_expression_class.def("__gt__", &__gt__<ContinuousPredicate,RealExpression,RealExpression>);
    //real_expression_class.def("__cmp__", &__cmp__<ContinuousPredicate,RealExpression,RealExpression>);
    real_expression_class.def(self_ns::str(self));

    def("neg", (RealExpression(*)(RealExpression)) &neg);
    def("rec", (RealExpression(*)(RealExpression)) &rec);
    def("sqr", (RealExpression(*)(RealExpression)) &sqr);
    //def("pow", (RealExpression(*)(RealExpression,int)) &pow);
    def("sqrt", (RealExpression(*)(RealExpression)) &sqrt);
    def("exp", (RealExpression(*)(RealExpression)) &exp);
    def("log", (RealExpression(*)(RealExpression)) &log);
    def("sin", (RealExpression(*)(RealExpression)) &sin);
    def("cos", (RealExpression(*)(RealExpression)) &cos);
    def("tan", (RealExpression(*)(RealExpression)) &tan);

    class_<RealAssignment> real_assignment_class("RealAssignment",no_init);
    real_assignment_class.def(self_ns::str(self));
    class_<RealDynamic> real_dynamic_class("RealDynamic",no_init);
    real_dynamic_class.def(self_ns::str(self));
    class_<RealUpdate> real_update_class("RealUpdate",no_init);
    real_update_class.def(self_ns::str(self));
    class_<StringUpdate> string_update_class("StringUpdate",no_init);
    string_update_class.def(self_ns::str(self));
    class_<IntegerAssignment> integer_assignment_class("IntegerAssignment",no_init);
    integer_assignment_class.def(self_ns::str(self));
    class_<IntegerUpdate> integer_update_class("IntegerUpdate",no_init);
    integer_update_class.def(self_ns::str(self));

    typedef Variable<Tribool> TriboolVariable;
    typedef Expression<Tribool> TriboolExpression;

    to_python< List<TriboolExpression> >();

    class_<DiscretePredicate> discrete_predicate_class("DiscretePredicate", init<DiscretePredicate>());
    discrete_predicate_class.def(init<bool>());
    discrete_predicate_class.def("__and__", &__and__<DiscretePredicate,DiscretePredicate,DiscretePredicate>);
    discrete_predicate_class.def("__or__", &__or__<DiscretePredicate,DiscretePredicate,DiscretePredicate>);
    discrete_predicate_class.def("__invert__", &__not__<DiscretePredicate,DiscretePredicate>);
    discrete_predicate_class.def(self_ns::str(self));

    class_<TriboolVariable> tribool_variable_class("TriboolVariable", init<std::string>());
    tribool_variable_class.def(self_ns::str(self));

    class_<TriboolExpression> tribool_expression_class("TriboolExpression",init<TriboolExpression>());
    tribool_expression_class.def("name", &TriboolExpression::operator_name);
    tribool_expression_class.def("subexpressions", &TriboolExpression::subexpressions);
    tribool_expression_class.def("substitute", &TriboolExpression::substitute<Real>);
    tribool_expression_class.def("substitute", &TriboolExpression::substitute<Tribool>);
    tribool_expression_class.def("simplify", &TriboolExpression::simplify);
    tribool_expression_class.def("__and__", &__and__<TriboolExpression,TriboolExpression,TriboolExpression>);
    tribool_expression_class.def("__or__", &__or__<TriboolExpression,TriboolExpression,TriboolExpression>);
    tribool_expression_class.def("__neg__", &__not__<TriboolExpression,TriboolExpression>);
    tribool_expression_class.def(self_ns::str(self));

    implicitly_convertible<TriboolVariable,TriboolExpression>();

    def("sgn", (TriboolExpression(*)(RealExpression)) &sgn);



    //class_<RealVariable> real_variable_class("RealVariable", init<std::string>());

}


void export_hybrid_automaton()
{
    // Don't use return_value_policy<copy_const_reference> since reference lifetime should not exceed automaton lifetime

    to_python< std::map<DiscreteEvent,VectorFunction> >();
    to_python< std::map<DiscreteEvent,ScalarFunction> >();
    to_python< std::set<DiscreteMode> >();
    to_python< std::set<DiscreteTransition> >();

    class_<DiscreteState> discrete_state_class("DiscreteState",init<DiscreteState>());
    discrete_state_class.def("__eq__", &__eq__<bool,DiscreteState,DiscreteState>);
    discrete_state_class.def("__ne__", &__ne__<bool,DiscreteState,DiscreteState>);
    discrete_state_class.def("__hash__", &__hash__<DiscreteState>);
    discrete_state_class.def(self_ns::str(self));
    implicitly_convertible<int,DiscreteState>();

    class_<DiscreteEvent> discrete_event_class("DiscreteEvent",init<DiscreteEvent>());
    discrete_event_class.def("__eq__", &__eq__<bool,DiscreteEvent,DiscreteEvent>);
    discrete_event_class.def("__ne__", &__ne__<bool,DiscreteEvent,DiscreteEvent>);
    discrete_event_class.def("__hash__", &__hash__<DiscreteEvent>);
    discrete_event_class.def(self_ns::str(self));
    implicitly_convertible<int,DiscreteEvent>();

    class_<DiscreteMode, shared_ptr<DiscreteMode> > discrete_mode_class("DiscreteMode",no_init);
    discrete_mode_class.def("location",&DiscreteMode::location);
    discrete_mode_class.def("dynamic",&DiscreteMode::dynamic,return_value_policy<reference_existing_object>());
    //discrete_mode_class.def("invariants",&DiscreteMode::invariants,return_value_policy<reference_existing_object>());
    discrete_mode_class.def("invariants",&DiscreteMode::invariants,return_value_policy<copy_const_reference>());
    discrete_mode_class.def("invariants",&discrete_mode_invariants);
    discrete_mode_class.def(self_ns::str(self));

    class_<DiscreteTransition, shared_ptr<DiscreteTransition> > discrete_transition_class("DiscreteTransition",no_init);
    discrete_transition_class.def("event",&DiscreteTransition::event);
    discrete_transition_class.def("source",&DiscreteTransition::source,return_value_policy<reference_existing_object>());
    discrete_transition_class.def("target",&DiscreteTransition::target,return_value_policy<reference_existing_object>());
    discrete_transition_class.def("reset",&DiscreteTransition::reset,return_value_policy<reference_existing_object>());
    discrete_transition_class.def("activation",&DiscreteTransition::activation,return_value_policy<reference_existing_object>());
    discrete_transition_class.def("guard",&discrete_transition_guard);
    discrete_transition_class.def("urgency",&DiscreteTransition::forced);
    discrete_transition_class.def(self_ns::str(self));

    class_<HybridTime> hybrid_time_class("HybridTime",init<double,int>());
    //hybrid_time_class.def("continuous_time",&HybridTime::continuous_time);
    //hybrid_time_class.def("discrete_time",&HybridTime::discrete_time);

    class_<HybridAutomaton> hybrid_automaton_class("HybridAutomaton",init<>());
    hybrid_automaton_class.def("mode",&HybridAutomaton::mode,return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def("transition",&HybridAutomaton::transition,return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def("modes",&HybridAutomaton::modes,return_value_policy<copy_const_reference>());
    hybrid_automaton_class.def("transitions",(const std::set<DiscreteTransition>&(HybridAutomaton::*)()const) &HybridAutomaton::transitions,return_value_policy<copy_const_reference>());
    hybrid_automaton_class.def("transitions",(std::set<DiscreteTransition>(HybridAutomaton::*)(DiscreteState)const) &HybridAutomaton::transitions);
    hybrid_automaton_class.def("blocking_guards",(std::map<DiscreteEvent,VectorFunction>(HybridAutomaton::*)(DiscreteState)const) &HybridAutomaton::blocking_guards);
    hybrid_automaton_class.def("new_mode",(const DiscreteMode&(HybridAutomaton::*)(DiscreteState,const VectorFunction&)) &HybridAutomaton::new_mode, return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def("new_invariant",(const DiscreteMode&(HybridAutomaton::*)(DiscreteState,const ScalarFunction&)) &HybridAutomaton::new_invariant, return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def("new_transition",(const DiscreteTransition&(HybridAutomaton::*)(DiscreteEvent,DiscreteState,DiscreteState,const VectorFunction&,const ScalarFunction&,bool)) &HybridAutomaton::new_transition, return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def(self_ns::str(self));

}

void export_hybrid_system()
{
    class_<DiscreteValuation> discrete_valuation_class("DiscreteValuation",init<DiscreteValuation>());
    discrete_valuation_class.def(init<>());
    discrete_valuation_class.def("set", (void(DiscreteValuation::*)(const StringVariable&, const String&)) &DiscreteValuation::set);
    discrete_valuation_class.def("set", (void(DiscreteValuation::*)(const IntegerVariable&, const int&)) &DiscreteValuation::set);
    discrete_valuation_class.def(self_ns::str(self));


    class_<HybridSystem> hybrid_system_class("HybridSystem",init<>());
    hybrid_system_class.def("new_equation",(void(HybridSystem::*)(DiscretePredicate,RealAssignment)) &HybridSystem::new_equation);
    hybrid_system_class.def("new_dynamic",(void(HybridSystem::*)(DiscretePredicate,RealDynamic)) &HybridSystem::new_dynamic);
    hybrid_system_class.def("new_transition",(void(HybridSystem::*)(EventSet,DiscretePredicate,StringUpdate)) &HybridSystem::new_transition);
    hybrid_system_class.def("new_reset",(void(HybridSystem::*)(EventSet,DiscretePredicate,RealUpdate)) &HybridSystem::new_reset);
    hybrid_system_class.def("new_guard",(void(HybridSystem::*)(EventSet,DiscretePredicate,ContinuousPredicate)) &HybridSystem::new_guard);
    hybrid_system_class.def("new_invariant",(void(HybridSystem::*)(DiscretePredicate,ContinuousPredicate)) &HybridSystem::new_invariant);
    hybrid_system_class.def("new_disabled_events",(void(HybridSystem::*)(EventSet,DiscretePredicate))  &HybridSystem::new_disabled_events);
    hybrid_system_class.def("events",(EventSet(HybridSystem::*)(const DiscreteValuation&)const)  &HybridSystem::events);
    hybrid_system_class.def("dynamic",(VectorFunction(HybridSystem::*)(const DiscreteValuation&)const)  &HybridSystem::dynamic);
    hybrid_system_class.def("target",(DiscreteValuation(HybridSystem::*)(const Event&,const DiscreteValuation&)const)  &HybridSystem::target);
    hybrid_system_class.def("reset",(VectorFunction(HybridSystem::*)(const Event&,const DiscreteValuation&)const)  &HybridSystem::reset);
    hybrid_system_class.def("guard",(ScalarFunction(HybridSystem::*)(const Event&,const DiscreteValuation&)const)  &HybridSystem::guard);
    hybrid_system_class.def(self_ns::str(self));
}

void system_submodule() {
    export_formula();
    export_hybrid_automaton();
    export_hybrid_system();
}
