/***************************************************************************
 *            function_model_mixin.hpp
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

/*! \file function_model_mixin.hpp
 *  \brief Mixin for concrete functions on bounded domains.
 */

#ifndef ARIADNE_FUNCTION_MODEL_MIXIN_HPP
#define ARIADNE_FUNCTION_MODEL_MIXIN_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../function/function_model_interface.hpp"
#include "../function/function_model.hpp"

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

template<class FM, class P, class D, class C, class TR=FunctionModelTraits<P>> class FunctionModelMixin;
template<class FM, class P, class D, class TR=FunctionModelTraits<P>> using ScalarMultivariateFunctionModelMixin = FunctionModelMixin<FM,P,D,IntervalDomainType,TR>;
template<class FM, class P, class D, class TR=FunctionModelTraits<P>> using VectorMultivariateFunctionModelMixin = FunctionModelMixin<FM,P,D,BoxDomainType,TR>;

template<class FCTRY, class P, class TR=FunctionModelTraits<P>> class FunctionModelFactoryMixin;

template<class FM, class P, class D, class TR> class FunctionModelMixin<FM,P,D,IntervalDomainType,TR>
    : public virtual ScalarFunctionModelInterface<P,D,TR>
    , public ScalarFunctionMixin<FM,P,D>
{
    using C=IntervalDomainType;
  public:
    using typename TR::CoefficientType;
    using typename TR::ErrorType;    
    using typename TR::NumericType;
    using typename TR::NormType;
    using typename TR::RangeType;    

    typedef typename TR::NumberModelType NumberModelType;
    typedef FunctionModel<P,D,C,TR> FunctionModelType;
    typedef FunctionModelInterface<P,D,C,TR> FunctionModelInterfaceType;
  public:
    FM apply(OperatorCode op) const;
  public:
    FunctionModelInterfaceType* _clone() const override {
        return new FM(static_cast<const FM&>(*this)); }
    CoefficientType const _value() const override { return static_cast<const FM&>(*this).value(); }
    CoefficientType const _gradient_value(SizeType i) const override { return static_cast<const FM&>(*this).gradient_value(i); }
    ErrorType const _error() const override { return static_cast<const FM&>(*this).error(); }
    Void _set_error(ErrorType const& e) override { static_cast<FM&>(*this).set_error(e); }
    Void _clobber() override { static_cast<FM&>(*this).clobber(); }

    NormType const _norm() const override {
        return norm(static_cast<const FM&>(*this)); }
    RangeType const _range() const override {
         return static_cast<RangeType>(static_cast<const FM&>(*this).range()); }

    FunctionModelInterfaceType* _antiderivative(SizeType j) const override {
        return new FM(antiderivative(static_cast<const FM&>(*this),j)); }
    FunctionModelInterfaceType* _antiderivative(SizeType j, NumberModelType c) const override {
        return new FM(antiderivative(static_cast<const FM&>(*this),j,c)); }
     FunctionModelInterfaceType* _restriction(const BoxDomainType& d) const override {
        return new FM(restriction(static_cast<const FM&>(*this),d)); }
    FunctionModelInterfaceType* _apply(OperatorCode op) const override {
        return new FM(this->apply(op)); }
    NumberModelType _unchecked_evaluate(const Vector<NumberModelType>& x) const override {
        return unchecked_evaluate(static_cast<const FM&>(*this),x); }
    FunctionModelInterfaceType* _partial_evaluate(SizeType j, const NumberModelType& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FM&>(*this),j,c)); }
    FunctionModelInterfaceType* _embed(const BoxDomainType& d1, const BoxDomainType& d2) const override {
        return new FM(embed(d1,static_cast<const FM&>(*this),d2)); }
    Boolean _refines(const FunctionModelInterfaceType& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return refines(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f)); }
    Boolean _inconsistent(const FunctionModelInterfaceType& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return inconsistent(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f)); }
    FunctionModelInterfaceType* _refinement(const FunctionModelInterfaceType& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return new FM(refinement(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    Void _iadd(const NumberModelType& c) override {
        static_cast<FM&>(*this)+=c; }
    Void _imul(const NumberModelType& c) override {
        static_cast<FM&>(*this)*=c; }
    Void _isma(const NumberModelType& c, const FunctionModelInterfaceType& f) override {
        static_cast<FM&>(*this)+=c*dynamic_cast<const FM&>(f); }
    Void _ifma(const FunctionModelInterfaceType& f1, const FunctionModelInterfaceType& f2) override {
        static_cast<FM&>(*this)+=dynamic_cast<const FM&>(f1)*dynamic_cast<const FM&>(f2); }
};

template<class FM, class P, class D, class TR> FM ScalarMultivariateFunctionModelMixin<FM,P,D,TR>::apply(OperatorCode op) const {
    const FM& f=static_cast<const FM&>(*this);
    switch(op) {
        case OperatorCode::NEG: return neg(f);
        case OperatorCode::REC: return rec(f);
        case OperatorCode::EXP: return exp(f);
        default: ARIADNE_FAIL_MSG("ScalarFunctionModel<P,D,TR>::apply(OperatorCode op): Operator op="<<op<<" not implemented\n");
    }
}


template<class FM, class P, class D, class TR> class FunctionModelMixin<FM,P,D,BoxDomainType,TR>
    : public virtual VectorFunctionModelInterface<P,D,TR>
    , public VectorFunctionMixin<FM,P,D>
{
    using C=BoxDomainType;
    using VFM=FM;
    using SFM=typename Element<VFM>::Type;
  public:
    typedef typename Element<FM>::Type ConcreteScalarFunctionModelType;

    using typename TR::CoefficientType;
    using typename TR::ErrorType;
    using typename TR::NumericType;
    using typename TR::NormType;
    typedef Box<typename TR::RangeType> RangeType;

    typedef typename TR::NumberModelType NumberModelType;
    typedef FunctionModelInterface<P,D,C,TR> FunctionModelInterfaceType;
    typedef FunctionModel<P,D,C,TR> FunctionModelType;

    typedef FunctionModelInterface<P,D,IntervalDomainType,TR> ScalarFunctionModelInterfaceType;
    typedef FunctionModelInterface<P,D,BoxDomainType,TR> VectorFunctionModelInterfaceType;

    typedef ScalarFunction<P,D> ScalarFunctionType;
    typedef ScalarMultivariateFunctionInterface<P> ScalarMultivariateFunctionInterfaceType;
    typedef VectorMultivariateFunctionInterface<P> VectorMultivariateFunctionInterfaceType;
    typedef ScalarMultivariateFunction<P> ScalarMultivariateFunctionType;
    typedef VectorMultivariateFunction<P> VectorMultivariateFunctionType;
  public:
    virtual FunctionModelInterfaceType* _clone() const override { return new FM(static_cast<const FM&>(*this)); }
    virtual Void _set(SizeType i, const ScalarFunctionModelInterfaceType& sf) override {
        if(!dynamic_cast<SFM const*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<FM&>(*this).FM::set(i,dynamic_cast<SFM const&>(sf)); }
    virtual FunctionModelInterfaceType* _derivative(SizeType j) const override {
        ARIADNE_NOT_IMPLEMENTED; }

    Vector<ErrorType> const _errors() const override { return static_cast<const FM&>(*this).errors(); }
    ErrorType const _error() const override { return static_cast<const FM&>(*this).error(); }
    Void _clobber() override { static_cast<FM&>(*this).clobber(); }

    RangeType const _range() const override {
         return static_cast<RangeType>(static_cast<const FM&>(*this).range()); }
    NormType const _norm() const override {
         return norm(static_cast<const FM&>(*this)); }

    FunctionModelInterfaceType* _embed(const BoxDomainType& d1, const BoxDomainType& d2) const override {
        return heap_copy(embed(d1,static_cast<const FM&>(*this),d2)); }
    FunctionModelInterfaceType* _restriction(const BoxDomainType& d) const override {
        return new FM(restriction(static_cast<const FM&>(*this),d)); }
    Void _adjoin(const ScalarFunctionModelInterfaceType& f) override {
        static_cast<FM&>(*this).FM::adjoin(dynamic_cast<SFM const&>(f)); }
    
    VectorFunctionModelInterfaceType* _join(const VectorFunctionModelInterfaceType& f) const override {
        return heap_copy(join(static_cast<const VFM&>(*this),dynamic_cast<const VFM&>(f))); }
    VectorFunctionModelInterfaceType* _combine(const VectorFunctionModelInterfaceType& f) const override {
        return heap_copy(combine(static_cast<const VFM&>(*this),dynamic_cast<const VFM&>(f))); }
    Vector<NumberModelType> _unchecked_evaluate(const Vector<NumberModelType>& x) const override {
        return unchecked_evaluate(static_cast<const FM&>(*this),x); }
    ScalarFunctionModelInterfaceType* _compose(const ScalarMultivariateFunctionInterfaceType& f) const override {
        return heap_copy(compose(f,static_cast<const FM&>(*this))); }
    VectorFunctionModelInterfaceType* _compose(const VectorMultivariateFunctionInterfaceType& f) const override {
        return heap_copy(compose(f,static_cast<const FM&>(*this))); }
    ScalarFunctionModelInterfaceType* _unchecked_compose(const ScalarMultivariateFunctionInterfaceType& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<SFM const&>(f),static_cast<VFM const&>(*this))); }
    VectorFunctionModelInterfaceType* _unchecked_compose(const VectorMultivariateFunctionInterfaceType& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const FM&>(f),static_cast<const FM&>(*this))); }
    FunctionModelInterfaceType* _partial_evaluate(SizeType j, const NumberModelType& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FM&>(*this),j,c)); }
};


template<class FCTRY, class P, class TR> class FunctionModelFactoryMixin
    : public FunctionModelFactoryInterface<P>
{
    typedef BoxDomainType VD;
    typedef IntervalDomainType SD;
    friend class FunctionModelFactory<P>;
  public:
    typedef VD VectorDomainType;
    typedef SD ScalarDomainType;

    typedef Number<P> NumberType;
    typedef typename TR::NumberModelType NumberModelType;
    
    typedef ScalarFunctionInterface<P,VD> ScalarMultivariateFunctionInterfaceType;
    typedef VectorFunctionInterface<P,VD> VectorMultivariateFunctionInterfaceType;
    
    typedef ScalarFunctionModelInterface<P,VD,TR> ScalarMultivariateFunctionModelInterfaceType;
    typedef VectorFunctionModelInterface<P,VD,TR> VectorMultivariateFunctionModelInterfaceType;

    typedef typename TR::NumericType NumericType;
    
  public:
    virtual FunctionModelFactoryInterface<P,TR>* clone() const override { return new FCTRY(this->upcast()); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << this->upcast(); }

/*
    NumberModelType create(const NumberType& number) const;
    ScalarMultivariateFunctionModelType create(const BoxDomainType& domain, const ScalarMultivariateFunctionInterfaceType& function) const;
    VectorMultivariateFunctionModelType create(const BoxDomainType& domain, const VectorMultivariateFunctionInterfaceType& function) const;
    ScalarMultivariateFunctionModelType create_zero(const BoxDomainType& domain) const;
    VectorMultivariateFunctionModelType create_zeros(SizeType result_size, const BoxDomainType& domain) const;
    ScalarMultivariateFunctionModelType create_constant(const BoxDomainType& domain, const NumberType& value) const;
    ScalarMultivariateFunctionModelType create_constant(const BoxDomainType& domain, const NumberModelType& value) const;
    VectorMultivariateFunctionModelType create_constants(const BoxDomainType& domain, const Vector<NumberType>& values) const;
    VectorMultivariateFunctionModelType create_constants(const BoxDomainType& domain, const Vector<NumberModelType>& values) const;
    ScalarMultivariateFunctionModelType create_coordinate(const BoxDomainType& domain, SizeType index) const;
    ScalarMultivariateFunctionModelType create_identity(const IntervalDomainType& domain) const;
    VectorMultivariateFunctionModelType create_identity(const BoxDomainType& domain) const;
    NumberModelType create_number(const NumberType& number) const;
*/
  private:
    template<class T> static inline T* heap_move(T&& t) { return new T(std::forward<T>(t)); }
    inline FCTRY const& upcast() const { return static_cast<FCTRY const&>(*this); }
    virtual NumberModelType _create(const NumberType& number) const override {
        return this->upcast().create(number); }
    virtual ScalarMultivariateFunctionModelInterfaceType* _create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterfaceType& function) const override {
        return heap_move(this->upcast().create(domain,function)); };
    virtual VectorMultivariateFunctionModelInterfaceType* _create(const VectorDomainType& domain, const VectorMultivariateFunctionInterfaceType& function) const override {
        return heap_move(this->upcast().create(domain,function)); };

    virtual ScalarMultivariateFunctionModelInterfaceType* _create_zero(const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_zero(domain)); };
    virtual ScalarMultivariateFunctionModelInterfaceType* _create_constant(const VectorDomainType& domain, const NumberType& value) const override {
        return heap_move(this->upcast().create_constant(domain,value)); };
    virtual ScalarMultivariateFunctionModelInterfaceType* _create_coordinate(const VectorDomainType& domain, SizeType j) const override {
        return heap_move(this->upcast().create_coordinate(domain,j)); };
    virtual VectorMultivariateFunctionModelInterfaceType* _create_zeros(SizeType n, const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_zeros(n,domain)); };
    virtual VectorMultivariateFunctionModelInterfaceType* _create_constants(const VectorDomainType& domain, const Vector<NumberType>& values) const override {
        return heap_move(this->upcast().create_constants(domain,values)); };
    virtual VectorMultivariateFunctionModelInterfaceType* _create_identity(const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_identity(domain)); };
};


} // namespace Ariadne

#endif
