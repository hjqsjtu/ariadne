/***************************************************************************
 *            inclusion_evolver.cpp
 *
 *  Copyright  2008-20  Luca Geretti, Pieter Collins, Sanja Zivanovic
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

#include "../function/taylor_function.hpp"
#include "../solvers/integrator.hpp"
#include "../solvers/bounder.hpp"
#include "../algebra/expansion.inl.hpp"
#include "inclusion_evolver.hpp"

namespace Ariadne {
/*
FloatDP volume(Vector<IntervalValidatedRangeType> const& box) {
    FloatDP result = 1.0;
    for (auto i: range(box.size())) {
        result *= box[i].width().raw();
    }
    return result;
}
*/
BoxDomainType initial_ranges_to_box(RealVariablesBox const& var_ranges) {
    auto vars = var_ranges.variables();
    List<IntervalDomainType> result;
    for (auto v : vars) {
        result.push_back(cast_exact(widen(IntervalDomainType(var_ranges[v].lower().get_d(),var_ranges[v].upper().get_d()))));
    }
    return Vector<IntervalDomainType>(result);
}

struct ScheduledApproximator
{
    Nat step;
    InclusionIntegratorHandle approximator;

    ScheduledApproximator(Nat s, InclusionIntegratorHandle a) : step(s), approximator(a) {}
};

inline OutputStream& operator<<(OutputStream& os, ScheduledApproximator const& sa) {
    return os << "(" << sa.step << ":" << sa.approximator << ")"; }

struct ScheduledApproximatorComparator {
    inline bool operator() (ScheduledApproximator const& sa1, ScheduledApproximator const& sa2) {
        return (sa1.step > sa2.step); }
};

class InclusionEvolverState {
  private:
    Nat _step;
    List<ScheduledApproximator> _schedule;
    Map<InclusionIntegratorHandle,Nat> _approximator_check_delay;
    Map<InclusionIntegratorHandle,Nat> _approximator_global_optima_count;
    Map<InclusionIntegratorHandle,Nat> _approximator_local_optima_count;
  public:
    InclusionEvolverState(InclusionVectorField const& ivf, List<InputApproximation> const& approximations, IntegratorInterface const& integrator)
        : _step(0u)
    {
        InclusionIntegratorFactory factory(integrator);
        for (auto appro : approximations) {
            InclusionIntegratorHandle approximator = factory.create(ivf.function(),ivf.inputs(),appro);
            _schedule.push_back(ScheduledApproximator(0u,approximator));
            _approximator_global_optima_count[approximator] = 0;
            _approximator_local_optima_count[approximator] = 0;
            _approximator_check_delay[approximator] = 0;
        }
    }

    Nat step() const { return _step; }
    List<ScheduledApproximator> const& schedule() const { return _schedule; }
    Map<InclusionIntegratorHandle,Nat> const& check_delay() const { return _approximator_check_delay; }
    Map<InclusionIntegratorHandle,Nat> const& global_optima_count() const { return _approximator_global_optima_count; }
    Map<InclusionIntegratorHandle,Nat> const& local_optima_count() const { return _approximator_local_optima_count; }

    Void reset_local_optima_count() {
        for (auto entry : _approximator_local_optima_count)
            _approximator_local_optima_count[entry.first] = 0;
    }

    Void update_with_best(InclusionIntegratorHandle const& best) {
        for (auto appro : approximators_to_use()) {
            if (best == appro) _reset_check_delay(appro);
            else _increase_check_delay(appro);

            Nat offset = 1u<<_approximator_check_delay[appro];
            _schedule.push_back(ScheduledApproximator(_step+offset,appro));
        }
        std::sort(_schedule.begin(),_schedule.end(),ScheduledApproximatorComparator());

        _increase_optima_count(best);
    }

    Void append_to_schedule(Nat const& step, InclusionIntegratorHandle const& approximator) {
        _schedule.push_back(ScheduledApproximator(step,approximator));
    }

    List<InclusionIntegratorHandle> approximators_to_use() const {
        List<InclusionIntegratorHandle> result;
        for (auto i = _schedule.size(); i > 0; --i) {
            auto entry = _schedule[i-1];
            if (entry.step == _step) {
                result.push_back(entry.approximator);
            } else if (entry.step > _step) {
                break;
            }
        }
        return result;
    }

    Void next_step() {
        while (!_schedule.empty()) {
            auto entry = _schedule.back();
            if (entry.step == _step) {
                _schedule.pop_back();
            } else if (entry.step > _step) {
                break;
            }
        }
        _step++;
    }

  private:
    Void _reset_check_delay(InclusionIntegratorHandle const& approximator) { _approximator_check_delay[approximator] = 0; }
    Void _increase_check_delay(InclusionIntegratorHandle const& approximator) { _approximator_check_delay[approximator]++; }
    Void _increase_optima_count(InclusionIntegratorHandle const& approximator) {
        _approximator_global_optima_count[approximator]++;
        _approximator_local_optima_count[approximator]++;
    }
};


inline char activity_symbol(SizeType step) {
    switch (step % 4) {
    case 0: return '\\';
    case 1: return '|';
    case 2: return '/';
    default: return '-';
    }
}


inline Map<InclusionIntegratorHandle,FloatDP> convert_to_percentages(Map<InclusionIntegratorHandle,Nat> const& approximation_global_frequencies) {

    Nat total_steps(0);
    for (auto entry: approximation_global_frequencies) {
        total_steps += entry.second;
    }

    Map<InclusionIntegratorHandle,FloatDP> result;
    for (auto entry: approximation_global_frequencies) {
        result[entry.first] = 1.0/total_steps*entry.second;
    }

    return result;
}

InclusionEvolver::InclusionEvolver(SystemType const& system, SweeperDP const& sweeper, IntegratorInterface const& integrator, ReconditionerHandle const& reconditioner)
    : _system(system)
    , _sweeper(sweeper)
    , _integrator(integrator.clone())
    , _reconditioner(reconditioner)
    , _configuration(new ConfigurationType())
{
    assert(system.inputs().size() > 0);
}

Void InclusionEvolver::_recondition_and_update(ValidatedVectorMultivariateFunctionModelType& function, InclusionEvolverState& state) {
    if (_reconditioner.must_reduce_parameters(state)) {
        _reconditioner.update_from(state);
        _reconditioner.reduce_parameters(function);
        state.reset_local_optima_count();
    }

    if (_reconditioner.must_incorporate_errors(state)) {
        function = _reconditioner.incorporate_errors(function);
    }
}

List<ValidatedVectorMultivariateFunctionModelDP> InclusionEvolver::reach(BoxDomainType const& initial, Real const& tmax) {

    ARIADNE_LOG(2,"System: "<<_system<<"\n");
    ARIADNE_LOG(2,"Initial: "<<initial<<"\n");

    StepSizeType hsug(_configuration->maximum_step_size());

    ValidatedVectorMultivariateFunctionModelDP evolve_function = ValidatedVectorMultivariateTaylorFunctionModelDP::identity(initial,this->_sweeper);

    TimeStepType t;

    InclusionEvolverState state(_system,_configuration->approximations(),*_integrator);

    List<ValidatedVectorMultivariateFunctionModelDP> result;

    while (possibly(t<lower_bound(tmax))) {

        if (verbosity == 1)
            std::cout << "\r[" << activity_symbol(state.step()) << "] " << static_cast<int>(std::round(100*t.get_d()/tmax.get_d())) << "% " << std::flush;

        ARIADNE_LOG(3,"step#:"<<state.step()<<", t:"<<t<<", hsug:"<<hsug << "\n");

        ARIADNE_LOG(4,"n. of parameters="<<evolve_function.argument_size()<<"\n");

        auto approximators_to_use = state.approximators_to_use();

        ARIADNE_LOG(4,"approximators to use="<<approximators_to_use<<"\n");

        auto domx = cast_exact_box(evolve_function.range());

        UpperBoxType B;
        StepSizeType h;
        std::tie(h,B)=approximators_to_use.at(0).flow_bounds(domx,_system.inputs(),hsug);
        ARIADNE_LOG(3,"flow bounds = "<<B<<" (using h = " << h << ")\n");

        TimeStepType new_t = lower_bound(t+h);

        List<ValidatedVectorMultivariateFunctionModelDP> reach_functions;
        List<ValidatedVectorMultivariateFunctionModelDP> best_reach_functions;
        ValidatedVectorMultivariateFunctionModelDP best_evolve_function;
        InclusionIntegratorHandle best = approximators_to_use.at(0);
        FloatDPApproximation best_volume(std::numeric_limits<double>::infinity());

        for (auto approximator : approximators_to_use) {
            ARIADNE_LOG(5,"checking "<<approximator<<" approximator\n");

            auto current_reach=approximator.reach(domx,evolve_function,B,t,h);
            auto current_evolve=approximator.evolve(current_reach.at(current_reach.size()-1u),new_t);

            FloatDPApproximation current_volume = volume(current_evolve.range());
            if (possibly(current_volume < best_volume)) {
                best = approximator;
                ARIADNE_LOG(6,"best approximator: " << best << "\n");
                best_reach_functions = current_reach;
                best_evolve_function = current_evolve;
                best_volume = current_volume;
            }
        }

        if (approximators_to_use.size() > 1)
            ARIADNE_LOG(4,"chosen approximator: " << best << "\n");

        state.update_with_best(best);
        ARIADNE_LOG(4,"updated schedule: " << state.schedule() << "\n");

        reach_functions = best_reach_functions;
        evolve_function = best_evolve_function;

        ARIADNE_LOG(3,"evolve bounds="<<evolve_function.range()<<"\n");

        result.concatenate(reach_functions);

        this->_recondition_and_update(evolve_function,state);

        state.next_step();
        t=new_t;
    }

    ARIADNE_LOG(2,"approximation % ="<<convert_to_percentages(state.global_optima_count())<<"\n");

    return result;
}

struct IndexedFloatDPError
{
    SizeType index;
    FloatDPError value;

    IndexedFloatDPError() : index(0), value(FloatDPError()) {}
};

inline OutputStream& operator<<(OutputStream& os, IndexedFloatDPError const& ifl) {
    return os << "(" << ifl.index << ":" << std::scientific << ifl.value.raw() << std::fixed << ")"; }

struct IndexedFloatDPErrorComparator
{
    inline bool operator() (const IndexedFloatDPError& ifl1, const IndexedFloatDPError& ifl2)
    {
        return (ifl1.value.raw() < ifl2.value.raw());
    }
};

Bool LohnerReconditioner::must_reduce_parameters(InclusionEvolverState const& state) const {
    if (_number_of_steps_between_simplifications == 0) return false;
    return (state.step()%_number_of_steps_between_simplifications == _number_of_steps_between_simplifications-1);
}

Bool LohnerReconditioner::must_incorporate_errors(InclusionEvolverState const& state) const {
    return true;
}

Void LohnerReconditioner::update_from(InclusionEvolverState const& state) {
    auto freq = _number_of_steps_between_simplifications;
    auto n = _number_of_variables;
    auto m = _number_of_inputs;

    FloatDP npk = 0;
    FloatDP rho = _ratio_of_parameters_to_keep;
    for (auto entry: state.local_optima_count()) {
        SizeType ppi = entry.first.num_params_per_input();
        FloatDP partial = n + rho*(n+2*m) + (freq-1)*m*(2 - ppi);
        npk += partial*entry.second/freq;
    }

    _number_of_parameters_to_keep = static_cast<Nat>(round(npk).get_d());
}

ValidatedVectorMultivariateFunctionModelDP LohnerReconditioner::incorporate_errors(ValidatedVectorMultivariateFunctionModelDP const& f) const {

    ValidatedVectorMultivariateTaylorFunctionModelDP const& tf = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(f.reference());

    BoxDomainType domain=f.domain();
    BoxDomainType errors=cast_exact(cast_exact(f.errors())*FloatDPUpperInterval(-1,+1)); // FIXME: Avoid cast;

    ARIADNE_LOG(6,"Uniform errors:"<<errors<<"\n");

    ValidatedVectorMultivariateFunctionModelDP error_function=ValidatedVectorMultivariateTaylorFunctionModelDP::identity(errors,tf.properties());
    ValidatedVectorMultivariateFunctionModelDP result = embed(f,errors)+embed(domain,error_function);
    for(SizeType i=0; i!=result.result_size(); ++i) { result[i].clobber(); }
    return result;
}

Void LohnerReconditioner::reduce_parameters(ValidatedVectorMultivariateFunctionModelDP& f) const {
    ARIADNE_LOG(6,"simplifying\n");
    ARIADNE_LOG(6,"f="<<f<<"\n");

    auto m=f.argument_size();
    auto n=f.result_size();

    ARIADNE_LOG(6,"num.parameters="<<m<<", to keep="<< this->_number_of_parameters_to_keep <<"\n");

    ValidatedVectorMultivariateTaylorFunctionModelDP& tf = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP&>(f.reference());

    auto sweeper = tf.properties();

    // Compute effect of error terms, but not of original variables;
    Matrix<FloatDPError> C(m,n);
    for (auto i : range(n)) {
        auto p=tf[i].model().expansion();

        for (auto ac : p) {
            UniformConstReference<MultiIndex> a=ac.index();
            UniformReference<FloatDPValue> c=ac.coefficient();
            for (auto j : range(m)) {
                if (a[j]!=0) {
                    C[j][i] += mag(c);
                }
            }
        }
    }

    ARIADNE_LOG(6,"C"<<C<<"\n");

    Array<IndexedFloatDPError> Ce(m);
    for (auto j : range(m)) {
        Ce[j].index = j;
        for (auto i : range(n)) {
            Ce[j].value += C[j][i];
        }
    }
    ARIADNE_LOG(6,"Ce:"<<Ce<<"\n");
    auto SCe=Ce;
    std::sort(SCe.begin(),SCe.end(),IndexedFloatDPErrorComparator());
    ARIADNE_LOG(6,"SortedCe:"<<SCe<<"\n");
    List<SizeType> keep_indices;
    List<SizeType> remove_indices;

    if (m <= this->_number_of_parameters_to_keep) {
        ARIADNE_LOG(6, "Insufficient number of variables, not simplifying\n");
        return;
    }

    Nat number_of_variables_to_remove = m - this->_number_of_parameters_to_keep;
    ARIADNE_LOG(6, "Number of parameters to remove:" << _number_of_parameters_to_keep<<"\n");

    for (auto j : range(number_of_variables_to_remove)) {
        remove_indices.append(SCe[j].index);
    }

    for (auto j : range(number_of_variables_to_remove,m)) {
        keep_indices.append(SCe[j].index);
    }

    ARIADNE_LOG(2,"number of kept parameters: " << keep_indices.size() << "/" << m << "\n");

    ARIADNE_LOG(6,"keep_indices:"<<keep_indices<<"\n");
    ARIADNE_LOG(6,"remove_indices:"<<remove_indices<<"\n");

    for (auto i : range(n)) {
        FloatDPError error = tf[i].error();
        for(SizeType k=0; k!=remove_indices.size(); ++k) {
            error += mag(C[remove_indices[k]][i]);
        }
        tf[i].set_error(error);
    }

    auto old_domain=f.domain();
    auto new_domain=BoxDomainType(Vector<IntervalDomainType>(keep_indices.size(),[&old_domain,&keep_indices](SizeType j){return old_domain[keep_indices[j]];}));
    auto projection=ValidatedVectorMultivariateTaylorFunctionModelDP(m,new_domain,sweeper);
    for (auto i : range(new_domain.size())) { projection[keep_indices[i]]=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(new_domain,i,sweeper); }
    for (auto i : range(remove_indices.size())) {
        auto j=remove_indices[i]; auto cj=old_domain[j].midpoint();
        projection[j]=ValidatedScalarMultivariateTaylorFunctionModelDP::constant(new_domain,cj,sweeper); }
    f=compose(f,projection);
}


InclusionEvolverConfiguration::InclusionEvolverConfiguration()
{
    maximum_step_size(1);
    maximum_enclosure_radius(100.0);
    enable_parameter_reduction(true);
    approximations({ZeroApproximation(),ConstantApproximation(),AffineApproximation(),SinusoidalApproximation(),PiecewiseApproximation()});
}


OutputStream&
InclusionEvolverConfiguration::_write(OutputStream& os) const
{
    os << "InclusionEvolverConfiguration"
       << ",\n  maximum_step_size=" << maximum_step_size()
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n  enable_parameter_reduction=" << enable_parameter_reduction()
       << ",\n  approximations=" << approximations()
       << "\n)\n";
    return os;
}

} // namespace Ariadne;
