/***************************************************************************
 *            function_model.hpp
 *
 *  Copyright 2011-17  Pieter Collins
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

/*! \file function_model.hpp
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_MODEL_HPP
#define ARIADNE_FUNCTION_MODEL_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../function/function.decl.hpp"
#include "../function/function_model_interface.hpp"

#include "../numeric/operators.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/operations.hpp"
#include "../function/domain.hpp"

#include "../function/function_interface.hpp"
#include "../function/function_mixin.hpp"
#include "../function/function.hpp"

namespace Ariadne {

template<class P, class D, class TR> struct AlgebraOperations<ScalarFunctionModel<P,D,TR>>;

// FIXME: Extend with univariate case
template<class P, class TR> class FunctionModelFactory {
    SharedPointer<const FunctionModelFactoryInterface<P,TR>> _ptr;
    typedef IntervalDomainType SD;
    typedef BoxDomainType VD;
  public:
    typedef P Paradigm;
    typedef TR Traits;
    
    typedef SD ScalarDomainType;
    typedef VD VectorDomainType;

    typedef Number<P> NumberType;
    typedef ScalarFunction<P,VD> ScalarMultivariateFunctionType;
    typedef VectorFunction<P,VD> VectorMultivariateFunctionType;

    typedef typename TR::NumericType NumberModelType;
    typedef ScalarFunctionModel<P,VD,TR> ScalarMultivariateFunctionModelType;
    typedef VectorFunctionModel<P,VD,TR> VectorMultivariateFunctionModelType;

    typedef FunctionModelFactoryInterface<P,TR> FunctionModelFactoryInterfaceType;
    typedef FunctionModelFactory<P,TR> FunctionModelFactoryType;
    
  public:    
    operator const FunctionModelFactoryInterfaceType& () const { return *_ptr; }

    explicit FunctionModelFactory(const FunctionModelFactoryInterfaceType* p) : _ptr(p) { }
    explicit FunctionModelFactory(SharedPointer<const FunctionModelFactoryInterfaceType> p) : _ptr(p) { }

    NumberModelType create(NumberType const& c) const {
        return NumberModelType(this->_ptr->_create(c)); }
    ScalarMultivariateFunctionModelType create(VectorDomainType const& dom, ScalarMultivariateFunctionType const& f) const {
        return ScalarMultivariateFunctionModelType(this->_ptr->_create(dom,f)); }
    VectorMultivariateFunctionModelType create(VectorDomainType const& dom, VectorMultivariateFunctionType const& f) const {
        return VectorMultivariateFunctionModelType(this->_ptr->_create(dom,f)); }

    NumberModelType create_number(NumberType const& c) const {
        return NumberModelType(this->_ptr->_create(c)); }
    ScalarMultivariateFunctionModelType create_zero(VectorDomainType const& dom) const {
        return ScalarMultivariateFunctionModelType(this->_ptr->_create_zero(dom)); }
    ScalarMultivariateFunctionModelType create_constant(VectorDomainType const& dom, NumberType const& c) const {
        return ScalarMultivariateFunctionModelType(this->_ptr->_create_constant(dom,c)); }
    ScalarMultivariateFunctionModelType create_constant(VectorDomainType const& dom, ExactNumber const& c) const {
        return this->create_constant(dom,NumberType(c)); }
#warning 
//    ScalarMultivariateFunctionModelType create_constant(VectorDomainType const& dom, NumberModelType const& c) const {
//        return ScalarMultivariateFunctionModelType(this->_ptr->_create_constant(dom,c)); }
    ScalarMultivariateFunctionModelType create_coordinate(VectorDomainType const& dom, SizeType j) const {
        return ScalarMultivariateFunctionModelType(this->_ptr->_create_coordinate(dom,j)); }
    VectorMultivariateFunctionModelType create_zeros(SizeType n, VectorDomainType const& dom) const {
        return VectorMultivariateFunctionModelType(this->_ptr->_create_zeros(n,dom)); }
    VectorMultivariateFunctionModelType create_identity(VectorDomainType const& dom) const {
        return VectorMultivariateFunctionModelType(this->_ptr->_create_identity(dom)); }

    // FIXME: Should return a univariate model
    ScalarMultivariateFunctionModelType create_identity(ScalarDomainType const& dom) const {
        return ScalarMultivariateFunctionModelType(this->_ptr->_create_coordinate(VectorDomainType(1u,dom),0u)); }
    VectorMultivariateFunctionModelType create_zeros(SizeType n, ScalarDomainType const& dom) const {
        return VectorMultivariateFunctionModelType(this->_ptr->_create_zeros(n,VectorDomainType(1u,dom))); }

    friend OutputStream& operator<<(OutputStream& os, FunctionModelFactoryType const& factory) { return factory._ptr->_write(os); }
};

template<class P, class TR> auto inline
FunctionModelFactoryInterface<P,TR>::create_identity(IntervalDomainType const& dom) const -> ScalarMultivariateFunctionModelType {
    return ScalarMultivariateFunctionModelType(this->create_coordinate(BoxDomainType(1u,dom),0u)); }

template<class FCTRY, class D> class FunctionModelCreator {
    using P = typename FCTRY::Paradigm;
    using TR = typename FCTRY::Traits;
  public:
    typedef FCTRY FactoryType;
    typedef D DomainType;
    typedef P Paradigm;

    typedef Number<P> NumberType;
    typedef ScalarFunction<P,D> ScalarFunctionType;
    typedef VectorFunction<P,D> VectorFunctionType;
    
    typedef typename FactoryType::NumberModelType NumberModelType;
    typedef ScalarFunctionModel<P,D,TR> ScalarFunctionModelType;
    typedef VectorFunctionModel<P,D,TR> VectorFunctionModelType;
    
    explicit FunctionModelCreator(DomainType domain, FactoryType factory) : _factory(factory), _domain(domain) { }

    NumberModelType create(NumberType const& c) const { return this->_factory.create(c); }
    // FIXME: Declaring this to return ScalarFunctionModel causes ambiguity in add(ScalarTaylorFunctionModel,ScalarFunction)
    decltype(auto) create(ScalarFunctionType const& f) { return this->_factory.create(this->_domain,f); }
    VectorFunctionModelType create(VectorFunctionType const& f) { return this->_factory.create(this->_domain,f); }
    ScalarFunctionModelType create_zero() { return this->_factory.create_zero(this->_domain); }
    VectorFunctionModelType create_zeros(SizeType n) { return this->_factory.create_zeros(n,this->_domain); }
    VectorFunctionModelType create_identity() { return this->_factory.create_identity(this->_domain); }

    template<class X, EnableIf<IsConvertible<X,NumberModelType>> =dummy, DisableIf<IsConstructible<X,NumberType>> =dummy>
        NumberModelType const& create(X const& c) const { return c; }
    ScalarFunctionModelType const& create(ScalarFunctionModelType const& f) const { return f; }
    VectorFunctionModelType const& create(VectorFunctionModelType const& f) const { return f; }
  protected:
    FactoryType _factory;
    DomainType _domain;
};

// FIXME: Merge with multivariate case
template<class FCTRY> class FunctionModelCreator<FCTRY,IntervalDomainType> {
    using P = typename FCTRY::Paradigm;
    using TR = typename FCTRY::Traits;
    using D = IntervalDomainType;
  public:
    typedef FCTRY FactoryType;
    typedef D DomainType;
    typedef P Paradigm;

    typedef Number<P> NumberType;
    typedef ScalarFunction<P,D> ScalarFunctionType;
    typedef VectorFunction<P,D> VectorFunctionType;
    
    typedef typename FactoryType::NumberModelType NumberModelType;
    typedef ScalarFunctionModel<P,D,TR> ScalarFunctionModelType;
    typedef VectorFunctionModel<P,D,TR> VectorFunctionModelType;
    
    explicit FunctionModelCreator(DomainType domain, FactoryType factory) : _factory(factory), _domain(domain) { }

    NumberModelType create(NumberType const& c) const { return this->_factory.create(c); }
    ScalarFunctionModelType create(ScalarFunctionType const& f) { ARIADNE_NOT_IMPLEMENTED; }
    VectorFunctionModelType create(VectorFunctionType const& f) { ARIADNE_NOT_IMPLEMENTED; }
    ScalarFunctionModelType create_zero() { ARIADNE_NOT_IMPLEMENTED; }
    VectorFunctionModelType create_zeros(SizeType n) { ARIADNE_NOT_IMPLEMENTED; }
    ScalarFunctionModelType create_identity() { ARIADNE_NOT_IMPLEMENTED; }

#warning
//    NumberModelType const& create(NumberModelType const& c) const { return c; }
  protected:
    FactoryType _factory;
    DomainType _domain;
};

#warning Should have concrete-generic operators if NumericType is concrete
//! \ingroup FunctionModelSubModule
//! \brief Generic scalar functions on singleton domains.
template<class P, class D, class TR> class FunctionModel<P,D,IntervalDomainType,TR>
//    : public DispatchTranscendentalAlgebraOperations<ScalarFunctionModel<P,D>, Number<P>>
    : public DispatchTranscendentalAlgebraOperations<ScalarFunctionModel<P,D,TR>, Number<P>>
    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P,D,TR>,ScalarMultivariateFunction<P>>
//    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P,D,TR>,Number<P>>
{
    static_assert(IsSame<D,IntervalDomainType>::value or IsSame<D,BoxDomainType>::value,"");
    typedef IntervalDomainType C;
  public:
    typedef ScalarFunction<P,D> GenericType;
    typedef D DomainType;
    typedef C CodomainType;

    typedef TR Traits;
    typedef typename Traits::CoefficientType CoefficientType;
    typedef typename Traits::ErrorType ErrorType;
    typedef typename Traits::NumericType NumericType;
    typedef typename Traits::NormType NormType;
    typedef typename Traits::RangeType RangeType;

    typedef Number<P> NumberType;
    typedef typename TR::NumericType NumberModelType;
    typedef Function<P,D,C> FunctionType;
    typedef FunctionModelInterface<P,D,C,TR> FunctionModelInterfaceType;
    typedef FunctionModel<P,D,C,TR> FunctionModelType;

    typedef FunctionModel<ValidatedTag,D,C,TR> ValidatedFunctionModelType;
    
    typedef FunctionModelFactory<P,TR> FunctionModelFactoryType;
    typedef FunctionModelCreator<FunctionModelFactoryType,D> FunctionModelCreatorType;

    typedef FunctionModel<P,D,IntervalDomainType,TR> ScalarFunctionModelType;
    typedef FunctionModel<P,D,BoxDomainType,TR> VectorFunctionModelType;
    
  public:
    clone_on_copy_ptr< FunctionModelInterfaceType > _ptr;
  public:
    FunctionModel() : _ptr() { }
    explicit FunctionModel(FunctionModelInterfaceType* p) : _ptr(p) { }
    FunctionModel(const SharedPointer<const FunctionModelInterfaceType> p) : _ptr(p->_clone()) { }
    FunctionModel(const FunctionModelType& f) : _ptr(f._ptr) { }
    FunctionModel(const FunctionModelInterfaceType& f) : _ptr(f._clone()) { }
    FunctionModel(const FunctionType& f) : _ptr(dynamic_cast<FunctionModelInterfaceType*>(f.raw_pointer()->_clone())) { }
    operator FunctionType() const { return FunctionType(this->_ptr->_clone()); }
    operator FunctionModelInterfaceType& () { return *_ptr; }
    operator const FunctionModelInterfaceType& () const { return *_ptr; }
    const FunctionModelInterfaceType* raw_pointer() const { return _ptr.operator->(); }
    FunctionModelInterfaceType& reference() { return *_ptr; }
    const FunctionModelInterfaceType& reference() const { return *_ptr; }

    FunctionModelType& operator=(const NumberType& c);
    // FIXME: Useful if NumberModelType and NumberType are different
//    FunctionModelType& operator=(const NumberModelType& c);
    FunctionModelType& operator=(const FunctionType& f);
    FunctionModelType& operator=(const FunctionModelInterfaceType& f);
//    FunctionModelType& operator=(const ValidatedScalarMultivariateTaylorFunctionModelDP& f);

    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    template<class X> X operator() (const Vector<X>& x) const {
        return this->_ptr->_evaluate(x); }
    template<class X> X evaluate(const Vector<X>& x) const {
        return this->_ptr->_evaluate(x); }
    inline DomainType const domain() const { return this->_ptr->domain(); }
    inline CodomainType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const { return this->_ptr->_range(); }

    inline CoefficientType value() const { return this->_ptr->_value(); }
    inline CoefficientType gradient_value(SizeType j) const { return this->_ptr->_gradient_value(j); }
    inline ErrorType error() const { return this->_ptr->_error(); }

    inline Void set_error(const ErrorType& e) { return this->_ptr->_set_error(e); }
#warning    inline Void set_error(Nat e) { return this->_ptr->_set_error(ErrorType(e,this->error().precision())); }
    inline Void clobber() { return this->_ptr->_clobber(); }

    inline FunctionModelType apply(Operator op) const { return FunctionModelType(this->_ptr->_apply(op)); }
    inline Void restrict(const DomainType& d) { *this=restriction(*this,d); }
  public:
    friend FunctionModelCreatorType factory(FunctionModelType const& f) {
        FunctionModelFactoryType factory(f._ptr->_factory()); return FunctionModelCreatorType(f.domain(),factory); }
  public:
  public:
    friend FunctionModelType partial_evaluate(const FunctionModelType& f, SizeType j, const NumberModelType& c) {
        return f._ptr->_partial_evaluate(j,c); }

    friend NormType norm(const FunctionModelType& f) {
        return f._ptr->_norm(); }
    friend FunctionModelType derivative(const FunctionModelType& f, SizeType j) {
        return FunctionModelType(f._ptr->_derivative(j)); }
    friend FunctionModelType antiderivative(const FunctionModelType& f, SizeType j) {
        return FunctionModelType(f._ptr->_antiderivative(j)); }
    friend FunctionModelType antiderivative(const FunctionModelType& f, SizeType j, NumberModelType c) {
        return FunctionModelType(f._ptr->_antiderivative(j,c)); }
    friend FunctionModelType antiderivative(const FunctionModelType& f, SizeType j, const NumberType& c) {
        return antiderivative(f,j,NumberModelType(c,f.value().precision())); }

    friend FunctionModelType embed(const DomainType& d1, const FunctionModelType& f, const DomainType& d2) {
        return FunctionModelType(f._ptr->_embed(d1,d2)); }
    friend FunctionModelType embed(const DomainType& d, const FunctionModelType& f) {
        return embed(d,f,DomainType()); }
    friend FunctionModelType embed(const FunctionModelType& f, const BoxDomainType& d) {
        return embed(DomainType(),f,d); }
    friend FunctionModelType embed(const FunctionModelType& f, const IntervalDomainType& d) {
        return embed(f,DomainType(1,d)); }
    friend FunctionModelType restrict(const FunctionModelType& f, const DomainType& d) {
        return FunctionModelType(f._ptr->_restriction(d)); }
    friend FunctionModelType restriction(const FunctionModelType& f, const DomainType& d) {
        return FunctionModelType(f._ptr->_restriction(d)); }

    friend VectorFunctionModelType join(const ScalarFunctionModelType& f1, const ScalarFunctionModelType& f2) {
        return join(VectorFunctionModelType(1,f1),f2); }
    friend VectorFunctionModelType combine(const ScalarFunctionModelType& f1, const ScalarFunctionModelType& f2);
  public:
    friend ValidatedFunctionModelType refinement(const ValidatedFunctionModelType& f1, const ValidatedFunctionModelType& f2) {
        return ValidatedFunctionModelType(f1._ptr->_refinement(f2)); }
    friend Boolean inconsistent(const ValidatedFunctionModelType& f1, const ValidatedFunctionModelType& f2) {
        return f1._ptr->_inconsistent(f2); }
    friend Boolean refines(const ValidatedFunctionModelType& f1, const ValidatedFunctionModelType& f2) {
        return f1._ptr->_refines(f2); }
  public:
    friend OutputStream& operator<<(OutputStream& os, const FunctionModelType& f) {
        return os <<  f.operator ScalarMultivariateFunction<P>(); }
};

template<class P, class D, class TR> struct AlgebraOperations<ScalarFunctionModel<P,D,TR>> {
    typedef Number<P> NumberType;
    typedef typename TR::NumericType NumberModelType;
    typedef ScalarFunctionModel<P,D,TR> FunctionModelType;
    
    static FunctionModelType apply(Nul, FunctionModelType f) {
        f._ptr->_imul(NumberModelType(0)); return std::move(f); }
    static FunctionModelType apply(Neg, FunctionModelType f) {
        f._ptr->_imul(NumberModelType(-1)); return std::move(f); }
    static FunctionModelType apply(Add, FunctionModelType f1, const FunctionModelType& f2) {
        f1._ptr->_isma(NumberModelType(+1),f2); return std::move(f1); }
    static FunctionModelType apply(Sub, FunctionModelType f1, const FunctionModelType& f2) {
        f1._ptr->_isma(NumberModelType(-1),f2); return std::move(f1); }
    static FunctionModelType apply(Mul, const FunctionModelType& f1, const FunctionModelType& f2) {
        FunctionModelType r=factory(f1).create_zero(); r._ptr->_ifma(f1,f2); return r; }
    static FunctionModelType apply(Add, FunctionModelType f1, const NumberModelType& c2) {
        f1._ptr->_iadd(c2); return std::move(f1); }
    static FunctionModelType apply(Mul, FunctionModelType f1, const NumberModelType& c2) {
        f1._ptr->_imul(c2); return std::move(f1); }
    static FunctionModelType apply(Add, const NumberModelType& c1, FunctionModelType f2) {
        f2._ptr->_iadd(c1); return std::move(f2); }
    static FunctionModelType apply(Mul, const NumberModelType& c1, FunctionModelType f2) {
        f2._ptr->_imul(c1); return std::move(f2); }

    static FunctionModelType apply(Sub, FunctionModelType f1, const NumberModelType& c2) {
        return add(std::move(f1),neg(c2)); }
    static FunctionModelType apply(Div, FunctionModelType f1, const NumberModelType& c2) {
        return mul(std::move(f1),rec(c2)); }
    static FunctionModelType apply(Sub, const NumberModelType& c1, FunctionModelType f2) {
        return add(neg(std::move(f2)),c1); }
    static FunctionModelType apply(Div, const NumberModelType& c1, FunctionModelType f2) {
        return mul(rec(std::move(f2)),c1); }

/*
    static FunctionModelType apply(Add, FunctionModelType f1, const NumberType& c2) {
        NumberModelType s2=factory(f1).create(c2); return add(f1,s2); }
    static FunctionModelType apply(Sub, FunctionModelType f1, const NumberType& c2) {
        return add(f1,neg(c2)); }
    static FunctionModelType apply(Mul, FunctionModelType f1, const NumberType& c2) {
        NumberModelType s2=factory(f1).create(c2); return mul(f1,s2); }
    static FunctionModelType apply(Div, FunctionModelType f1, const NumberType& c2) {
        return mul(f1,rec(c2)); }
    static FunctionModelType apply(Add, const NumberType& c1, FunctionModelType f2) {
        return add(f2,c1); }
    static FunctionModelType apply(Sub, const NumberType& c1, FunctionModelType f2) {
        return add(neg(f2),c1); }
    static FunctionModelType apply(Mul, const NumberType& c1, FunctionModelType f2) {
        return mul(f2,c1); }
    static FunctionModelType apply(Div, const NumberType& c1, FunctionModelType f2) {
        return mul(rec(f2),c1); }
*/

    static FunctionModelType apply(Rec, const FunctionModelType& f) {
        return f.apply(Rec()); }
    static FunctionModelType apply(Div, const FunctionModelType& f1, const FunctionModelType& f2) {
        return mul(f1,rec(f2)); }
    static FunctionModelType apply(Pow, const FunctionModelType& f1, Int n2) {
        return generic_pow(f1,n2); }

    template<class OP> static FunctionModelType apply(OP op, const FunctionModelType& f) {
        ARIADNE_NOT_IMPLEMENTED; }

};


    // FIXME: Should not be needed since ScalarMultivariateFunctionModel has a representation
template<class P> inline ScalarMultivariateFunctionModel<P> embed(const ScalarMultivariateFunction<P>& f, const IntervalDomainType& d) {
    return embed(ScalarMultivariateFunctionModel<P>(f),d); }

//template<class P, class D, class TR> inline auto
//ScalarFunctionModel<P,D,TR>::operator=(const NumberModelType& c) -> FunctionModelType& {
//    (*this)*=NumberType(0); (*this)+=c; return *this; }
template<class P, class D, class TR> inline auto
ScalarFunctionModel<P,D,TR>::operator=(const NumberType& c) -> FunctionModelType& {
    return (*this)=factory(*this).create(c); }
template<class P, class D, class TR> inline auto
ScalarFunctionModel<P,D,TR>::operator=(const FunctionType& f) -> FunctionModelType& {
    return (*this)=factory(*this).create(f); }
template<class P, class D, class TR> inline auto
ScalarFunctionModel<P,D,TR>::operator=(const FunctionModelInterfaceType& f) -> FunctionModelType& {
    return (*this)=ScalarFunctionModelType(f._clone()); }



template<class V> struct Element;

template<class M> class ScaledFunctionPatch;
template<class M> class VectorScaledFunctionPatch;
template<class M> struct Element<VectorScaledFunctionPatch<M>> { typedef ScaledFunctionPatch<M> Type; };

typedef ScaledFunctionPatch<ValidatedTaylorModelDP> ValidatedScalarMultivariateTaylorFunctionModelDP;

template<class P, class D, class TR> class VectorFunctionModelElement
    : public DispatchTranscendentalAlgebraOperations<ScalarFunctionModel<P,D,TR>, Number<P>>
    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P,D,TR>>
{
    typedef ScalarFunctionModel<P,D,TR> ScalarFunctionModelType;
    typedef VectorFunctionModel<P,D,TR> VectorFunctionModelType;
    typedef VectorFunctionModelElement<P,D,TR> VectorFunctionModelElementType;

    typedef ScalarFunctionModel<ValidatedTag,D,TR> ValidatedScalarFunctionModelType;
  private:
     VectorFunctionModelType* _p; SizeType _i;
  public:
    typedef typename ScalarFunctionModelType::GenericType GenericType;
    typedef typename ScalarFunctionModelType::ErrorType ErrorType;
    
    operator const ScalarFunctionModelType () const;
    VectorFunctionModelElement(VectorFunctionModelType* p, SizeType i) : _p(p), _i(i) { }
    VectorFunctionModelElementType& operator=(const ScalarFunctionModelType& sf) {
        _p->set(_i,sf); return *this; }
    VectorFunctionModelElementType& operator=(const VectorFunctionModelElementType& sf) {
        return this->operator=(static_cast<ScalarFunctionModelType const>(sf)); }
    Void clobber() { ScalarFunctionModelType sf=_p->get(_i); sf.clobber(); _p->set(_i,sf); }
    decltype(auto) model() const { ScalarFunctionModelType sf=_p->get(_i); return sf.model(); }
    const ErrorType error() const { ScalarFunctionModelType sf=_p->get(_i); return sf.error(); }
    Void set_error(ErrorType e) const { ScalarFunctionModelType sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    Void set_error(Nat e) const { ScalarFunctionModelType sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    friend Boolean refines(const ValidatedScalarFunctionModelType& f1, const ValidatedScalarFunctionModelType& f2);
    friend Boolean inconsistent(const ValidatedScalarFunctionModelType& f1, const ValidatedScalarFunctionModelType& f2);
    friend ScalarFunctionModelType antiderivative(VectorFunctionModelElementType const& f, SizeType k) {
        return antiderivative(ScalarFunctionModelType(f),k); }
    friend inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModelElementType& function) {
        return os << static_cast< const ScalarFunctionModelType >(function); }
};

//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on singleton domains.
template<class P, class D, class TR> class FunctionModel<P,D,BoxDomainType,TR>
{
    static_assert(IsSame<D,IntervalDomainType>::value or IsSame<D,BoxDomainType>::value,"");
    typedef BoxDomainType C;
  public:
    clone_on_copy_ptr< VectorFunctionModelInterface<P,D> > _ptr;
  public:
    typedef D DomainType;
    typedef C CodomainType;

    typedef TR Traits;
    typedef typename Traits::CoefficientType CoefficientType;
    typedef typename Traits::ErrorType ErrorType;
    typedef typename Traits::NumericType NumericType;
    typedef typename Traits::NormType NormType;
    typedef Box<typename Traits::RangeType> RangeType;

    typedef Number<P> NumberType;
    typedef typename Traits::NumericType NumberModelType;
    typedef Function<P,D,C> FunctionType;
    typedef FunctionModelInterface<P,D,C,TR> FunctionModelInterfaceType;
    typedef FunctionModel<P,D,C,TR> FunctionModelType;

    typedef FunctionModel<ValidatedTag,D,C,TR> ValidatedFunctionModelType;
    
    typedef FunctionModelFactory<P,TR> FunctionModelFactoryType;
    typedef FunctionModelCreator<FunctionModelFactoryType,D> FunctionModelCreatorType;

    typedef VectorFunctionModelElement<P,D,TR> VectorFunctionModelElementType;
    typedef FunctionModel<P,D,IntervalDomainType,TR> ScalarFunctionModelType;
    typedef FunctionModel<P,D,BoxDomainType,TR> VectorFunctionModelType;

    
    typedef ScalarMultivariateFunction<P> ScalarMultivariateFunctionType;
    typedef VectorMultivariateFunction<P> VectorMultivariateFunctionType;
    typedef FunctionModelInterface<P,BoxDomainType,IntervalDomainType,TR> ScalarMultivariateFunctionModelInterfaceType;
    typedef FunctionModelInterface<P,BoxDomainType,BoxDomainType,TR> VectorMultivariateFunctionModelInterfaceType;
    typedef FunctionModel<P,BoxDomainType,IntervalDomainType,TR> ScalarMultivariateFunctionModelType;
    typedef FunctionModel<P,BoxDomainType,BoxDomainType,TR> VectorMultivariateFunctionModelType;

  public:
    inline FunctionModel() : _ptr() { }
    inline FunctionModel(SharedPointer<const FunctionModelInterfaceType> vfp)
        : _ptr(vfp->_clone()) { }
    inline FunctionModel(SizeType n, const ScalarFunctionModelType& sf) {
        FunctionModelCreatorType factry(factory(sf)); *this=factry.create_zeros(n);
        for(SizeType i=0; i!=n; ++i) { (*this)[i]=sf; } }
    inline FunctionModel(Array<ScalarFunctionModelType> const& asf)
        : FunctionModel(asf.size(),asf[0]) { for(SizeType i=0; i!=asf.size(); ++i) { (*this)[i]=asf[i]; } }
    inline FunctionModel(List<ScalarFunctionModelType> const& lsf)
        : FunctionModel(lsf.size(),lsf[0]) { for(SizeType i=0; i!=lsf.size(); ++i) { (*this)[i]=lsf[i]; } }
    inline explicit FunctionModel(FunctionModelInterfaceType* p) : _ptr(p) { }
    inline FunctionModel(const FunctionModelInterfaceType& f) : _ptr(f._clone()) { }
    inline FunctionModel(const FunctionModelType& f) : _ptr(f._ptr) { }
    inline operator const FunctionModelInterfaceType& () const { return *_ptr; }
    inline operator FunctionType () const { return FunctionType(*_ptr); }
    inline const FunctionModelInterfaceType* raw_pointer() const { return _ptr.operator->(); }
    inline const FunctionModelInterfaceType& reference() const { return *_ptr; }
    inline FunctionModelInterfaceType& reference() { return *_ptr; }

    inline SizeType result_size() const { return this->_ptr->result_size(); }
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    inline SizeType size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }

    inline ScalarFunctionModelType const get(SizeType i) const { return ScalarFunctionModelType(this->_ptr->_get(i)); }
    inline Void set(SizeType i, ScalarFunctionModelType const& sf) { this->_ptr->_set(i,sf); }
    inline ScalarFunctionModelType const operator[](SizeType i) const { return this->get(i); }
    inline VectorFunctionModelElementType operator[](SizeType i) { return VectorFunctionModelElementType(this,i); }
    inline VectorFunctionModelType operator[](Range rng) { VectorFunctionModelType r=factory(*this).create_zeros(rng.size());
        for(SizeType i=0; i!=rng.size(); ++i) { r[i]=this->operator[](rng[i]); } return r; }

    inline DomainType const domain() const { return this->_ptr->domain(); }
    inline CodomainType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const { return this->_ptr->_range(); }
    inline Vector<ErrorType> const errors() const { return this->_ptr->_errors(); }
    inline ErrorType const error() const { return this->_ptr->_error(); }
    inline Void clobber() { this->_ptr->_clobber(); }
    inline Matrix<NumericType> const jacobian(const Vector<NumericType>& x) const;
//        Vector<Differential<NumericType>> dfx=this->_ptr->_evaluate(Differential<NumericType>::variables(1u,x));
//        return dfx.jacobian(); }

    inline Void restrict(const DomainType& d) { this->_ptr->restrict(d); }
  public:
    friend FunctionModelCreatorType factory(VectorFunctionModelType const& f) {
        FunctionModelFactoryType factory(f._ptr->_factory()); return FunctionModelCreatorType(f.domain(),factory); }
  public:
    friend inline ScalarFunctionModelType compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModelType& g) {
        return ScalarFunctionModelType(g._ptr->_compose(f)); }
    friend inline ScalarFunctionModelType compose(const ScalarMultivariateFunctionModelType& f, const VectorFunctionModelType& g) {
        return ScalarFunctionModelType(g._ptr->_compose(f)); }
    friend inline VectorFunctionModelType compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModelType& g) {
        return VectorFunctionModelType(g._ptr->_compose(f)); }
    friend inline VectorFunctionModelType compose(const VectorMultivariateFunctionModelType& f, const VectorFunctionModelType& g) {
        return VectorFunctionModelType(g._ptr->_compose(f)); }

    friend inline ScalarFunctionModelType unchecked_compose(const ScalarMultivariateFunctionModelType& f, const VectorFunctionModelType& g) {
        return ScalarFunctionModelType(g._ptr->_unchecked_compose(f)); }
    friend inline VectorFunctionModelType unchecked_compose(const VectorMultivariateFunctionModelType& f, const VectorFunctionModelType& g) {
        return VectorFunctionModelType(g._ptr->_unchecked_compose(f)); }

    friend inline VectorFunctionModelType operator+(const VectorFunctionModelType& f) {
        return VectorFunctionModelType(f._ptr->_clone()); }
    friend inline VectorFunctionModelType operator-(const VectorFunctionModelType& f) {
        VectorFunctionModelType r=f; for(SizeType i=0; i!=r.size(); ++i) { r[i]=-f[i]; } return r; }
    friend inline VectorFunctionModelType operator+(const VectorFunctionModelType& f1, const VectorFunctionModelType& f2) {
        VectorFunctionModelType r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+f2[i]; } return r; }
    friend inline VectorFunctionModelType operator-(const VectorFunctionModelType& f1, const VectorFunctionModelType& f2) {
        VectorFunctionModelType r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
    friend inline VectorFunctionModelType operator+(const VectorFunctionModelType& f1, const Vector<NumberModelType>& c2) {
        VectorFunctionModelType r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+c2[i]; } return r; }
    friend inline VectorFunctionModelType operator-(const VectorFunctionModelType& f1, const Vector<NumberModelType>& c2) {
        VectorFunctionModelType r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-c2[i]; } return r; }
    friend inline VectorFunctionModelType operator+(const Vector<NumberModelType>& c1, const VectorFunctionModelType& f2);
    friend inline VectorFunctionModelType operator-(const Vector<NumberModelType>& c1, const VectorFunctionModelType& f2);
    friend inline VectorFunctionModelType operator*(const NumberModelType& c1, const VectorFunctionModelType& f2) {
        VectorFunctionModelType r=f2; for(SizeType i=0; i!=r.size(); ++i) { r[i]=c1*f2[i]; } return r; }
    friend inline VectorFunctionModelType operator*(const VectorFunctionModelType& f1, const NumberModelType& c2) {
        VectorFunctionModelType r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]*c2; } return r; }
    friend inline VectorFunctionModelType operator/(const VectorFunctionModelType& f1, const NumberModelType& c2) {
        VectorFunctionModelType r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]/c2; } return r; }

    friend inline VectorFunctionModelType operator+(const VectorFunctionModelType& f1, const VectorMultivariateFunction<P>& f2) {
        return f1+factory(f1).create(f2); }
    friend inline VectorFunctionModelType operator-(const VectorFunctionModelType& f1, const VectorMultivariateFunction<P>& f2) {
        return f1-factory(f1).create(f2); }
    friend inline VectorFunctionModelType operator+(const VectorMultivariateFunction<P>& f1, const VectorFunctionModelType& f2) {
        return factory(f2).create(f1)+f2; }
    friend inline VectorFunctionModelType operator-(const VectorMultivariateFunction<P>& f1, const VectorFunctionModelType& f2) {
        return factory(f2).create(f1)-f2; }

  public:
    friend NormType norm(const FunctionModelType& f) {
        return f._ptr->_norm(); }
    friend FunctionModelType embed(const DomainType& d1, const FunctionModelType& f, const DomainType& d2) {
        return FunctionModelType(f._ptr->_embed(d1,d2)); }
    friend FunctionModelType embed(const DomainType& d, const FunctionModelType& f) {
        return embed(d,f,DomainType()); }
    friend FunctionModelType embed(const FunctionModelType& f, const BoxDomainType& d) {
        return embed(DomainType(),f,d); }
    friend FunctionModelType embed(const FunctionModelType& f, const IntervalDomainType& d) {
        return embed(f,DomainType(1,d)); }
    friend FunctionModelType restrict(const FunctionModelType& f, const DomainType& d) {
        FunctionModelInterfaceType* rptr=f._ptr->_clone(); rptr->restrict(d); return FunctionModelType(rptr); }
    friend FunctionModelType restriction(const FunctionModelType& f, const DomainType& d) {
        FunctionModelInterfaceType* rptr=f._ptr->_clone(); rptr->restrict(d); return FunctionModelType(rptr); }

    friend Vector<NumberModelType> unchecked_evaluate(const VectorFunctionModelType& f, const Vector<NumberModelType>& x) {
        return f._ptr->_unchecked_evaluate(x); }

    friend ScalarFunctionModelType unchecked_compose(const ScalarMultivariateFunctionType& f, const VectorFunctionModelType& g) {
        ScalarMultivariateFunctionModelInterfaceType const* fptr = dynamic_cast<ScalarMultivariateFunctionModelInterfaceType const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(ScalarMultivariateFunctionModelType(*fptr),g); } else { return compose(f,g); } }
    friend VectorFunctionModelType unchecked_compose(const VectorMultivariateFunctionType& f, const VectorFunctionModelType& g) {
        VectorMultivariateFunctionModelInterfaceType const* fptr = dynamic_cast<VectorMultivariateFunctionModelInterfaceType const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(VectorMultivariateFunctionModelType(*fptr),g); } else { return compose(f,g); } }

    //friend VectorFunctionModelType join(const ScalarFunctionModelType& f1, const ScalarFunctionModelType& f2);
    friend VectorFunctionModelType join(const ScalarFunctionModelType& f1, const VectorFunctionModelType& f2) {
        return join(VectorFunctionModelType(1u,f1),f2); }
    friend VectorFunctionModelType join(const VectorFunctionModelType& f1, const ScalarFunctionModelType& f2) {
        VectorFunctionModelType r=VectorFunctionModelType(f1._ptr->_clone()); r._ptr->_adjoin(f2); return r; }
    friend VectorFunctionModelType join(const VectorFunctionModelType& f1, const VectorFunctionModelType& f2) {
        return VectorFunctionModelType(f1._ptr->_join(f2)); }

    friend VectorFunctionModelType combine(const ScalarFunctionModelType& f1, const ScalarFunctionModelType& f2) {
        return VectorFunctionModelType(1,f1)._ptr->_combine(VectorFunctionModelType(1,f2)); };
    friend VectorFunctionModelType combine(const ScalarFunctionModelType& f1, const VectorFunctionModelType& f2) {
        return VectorFunctionModelType(1,f1)._ptr->_combine(f2); };
    friend VectorFunctionModelType combine(const VectorFunctionModelType& f1, const ScalarFunctionModelType& f2) {
        return VectorFunctionModelType(f1._ptr->_combine(VectorFunctionModelType(1,f2))); };
    friend VectorFunctionModelType combine(const VectorFunctionModelType& f1, const VectorFunctionModelType& f2) {
        return VectorFunctionModelType(f1._ptr->_combine(f2)); }

    friend inline ValidatedFunctionModelType refinement(const ValidatedFunctionModelType& f1, const ValidatedFunctionModelType& f2) {
        ARIADNE_ASSERT_MSG(f1.size()==f2.size(),"refinement(f1,f2): f1="<<f1<<", f2="<<f2<<")");
        ValidatedFunctionModelType r=+f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=refinement(f1[i],f2[i]); } return r; }

    friend VectorFunctionModelType antiderivative(const VectorFunctionModelType& f, SizeType j) {
        VectorFunctionModelType r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j); }
        return r;
    }
    friend VectorFunctionModelType antiderivative(const VectorFunctionModelType& f, SizeType j, const NumberModelType& c) {
        VectorFunctionModelType r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j,c); }
        return r;
    }

    // FIXME: Useful if NumberModelType and NumberType are not the same
//    friend VectorFunctionModelType antiderivative(const VectorFunctionModelType& f, SizeType j, const NumberType& c) {
//        return antiderivative(f,j,NumberModelType(c,f[0].value().precision())); }


    friend VectorFunctionModelType partial_evaluate(const VectorFunctionModelType& f, SizeType j, const NumberModelType& c) {
        return VectorFunctionModelType(f._ptr->_partial_evaluate(j,c)); }
    // FIXME: Useful if NumberModelType and NumberType are not the same
//    friend VectorFunctionModelType partial_evaluate(const VectorFunctionModel<P,D>& f, SizeType j, const NumberType& c) {
//        return VectorFunctionModelType(f._ptr->_partial_evaluate(j,c)); }

    friend OutputStream& operator<<(OutputStream& os, const VectorFunctionModelType& f) {
        return os <<  f.operator VectorMultivariateFunction<P>(); }

};

template<class P, class D> inline Number<P> evaluate(const ScalarFunctionModel<P,D>& f, const Vector<Number<P>>& x) {
    return f(x); }
template<class P, class D> inline Number<P> unchecked_evaluate(const ScalarFunctionModel<P,D>& f, const Vector<Number<P>>& x) {
    return f._ptr->_unchecked_evaluate(x); }
template<class P, class D> inline CanonicalNumericType<P,DP> unchecked_evaluate(const ScalarFunctionModel<P,D>& f, const Vector<CanonicalNumericType<P,DP>>& x);
template<class P, class D> inline CanonicalNumericType<P,MP> unchecked_evaluate(const ScalarFunctionModel<P,D>& f, const Vector<CanonicalNumericType<P,MP>>& x);

// FIXME: Implement for Multiple-Precision versions
template<class P> inline CanonicalNumericType<P,DP> unchecked_evaluate(const ScalarMultivariateFunction<P>& f, const Vector<CanonicalNumericType<P,DP>>& x) {
    ScalarMultivariateFunctionModelInterface<P> const* fptr = dynamic_cast<ScalarMultivariateFunctionModelInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(ScalarMultivariateFunctionModel<P>(*fptr),x); } else { return evaluate(f,x); } }

template<class P> inline Vector<CanonicalNumericType<P,DP>> unchecked_evaluate(const VectorMultivariateFunction<P>& f, const Vector<CanonicalNumericType<P,DP>>& x) {
    VectorMultivariateFunctionModelInterface<P> const* fptr = dynamic_cast<VectorMultivariateFunctionModelInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(VectorMultivariateFunctionModel<P>(*fptr),x); } else { return evaluate(f,x); } }

template<class P, class D> inline ScalarFunctionModel<P,D> unchecked_compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModel<P,D>& g) {
    ScalarFunctionModelInterface<P,D> const* fptr = dynamic_cast<ScalarFunctionModelInterface<P,D> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(ScalarFunctionModel<P,D>(*fptr),g); } else { return compose(f,g); } }
template<class P, class D> inline VectorFunctionModel<P,D> unchecked_compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModel<P,D>& g) {
    VectorFunctionModelInterface<P,D> const* fptr = dynamic_cast<VectorFunctionModelInterface<P,D> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(VectorFunctionModel<P,D>(*fptr),g); } else { return compose(f,g); } }

// FIXME: Should be unneeded
template<class D> ScalarFunctionModel<ValidatedTag,D> unchecked_compose(const ScalarMultivariateFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag,D>& g) {
    return VectorFunctionModel<ValidatedTag,D>(g._ptr->_unchecked_compose(f)); }
template<class D> VectorFunctionModel<ValidatedTag,D> unchecked_compose(const VectorMultivariateFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag,D>& g) {
    return VectorFunctionModel<ValidatedTag,D>(g._ptr->_unchecked_compose(f)); }
template<class D> ScalarFunctionModel<ValidatedTag,D> restrict(const ScalarFunctionModel<ValidatedTag,D>& f, const D& dom) {
    return ScalarFunctionModel<ValidatedTag,D>(f._ptr->_restriction(dom)); }
template<class D> VectorFunctionModel<ValidatedTag,D> restrict(const VectorFunctionModel<ValidatedTag,D>& f, const D& dom) {
    return VectorFunctionModel<ValidatedTag,D>(f._ptr->_restriction(dom)); }


template<class P, class D, class TR> VectorFunctionModelElement<P,D,TR>::operator const ScalarFunctionModel<P,D,TR> () const {
    return ScalarFunctionModel<P,D,TR>(_p->get(_i)); }





// Full output
template<class T> struct Representation { const T* pointer; Representation(const T& t) : pointer(&t) { } const T& reference() const { return *pointer; } };
template<class T> inline Representation<T> representation(const T& t) { return Representation<T>(t); }

template<class T> class HasRepr {
    template<class TT, class = decltype(declval<TT>().repr(declval<OutputStream&>()))> static True test(int);
    template<class TT> static False test(...);
  public:
    static const bool value = decltype(test<T>(1))::value;
};

template<class T, EnableIf<HasRepr<T>> =dummy> inline OutputStream& operator<<(OutputStream& os, const Representation<T>& obj) { obj.reference().repr(os); return os; }
template<class T, DisableIf<HasRepr<T>> =dummy> inline OutputStream& operator<<(OutputStream& os, const Representation<T>& obj) { return os << obj.reference(); }

} // namespace Ariadne

#endif
