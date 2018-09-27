/***************************************************************************
 *            function_model_interface.hpp
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

/*! \file function_model_interface.hpp
 *  \brief Interface for functions on bounded sets.
 */

#ifndef ARIADNE_FUNCTION_MODEL_INTERFACE_HPP
#define ARIADNE_FUNCTION_MODEL_INTERFACE_HPP

#include "../function/function.decl.hpp"
#include "../function/function_interface.hpp"

#include "../numeric/operators.hpp"

namespace Ariadne {

template<class P, class PR, class PRE> class FunctionModelFactoryInterface;
template<class P, class D, class PR, class PRE> class FunctionModelCreatorInterface;

template<class P, class D, class C, class PR, class PRE> class FunctionModelInterface;

template<class P, class D, class PR, class PRE> class FunctionModelInterface<P,D,IntervalDomainType,PR,PRE>
    : public virtual FunctionInterface<P,D,IntervalDomainType>
{
    static_assert(IsSame<D,IntervalDomainType>::value or IsSame<D,BoxDomainType>::value,"");
    typedef IntervalDomainType C;
  public:
    typedef D DomainType;
    typedef C CodomainType;
    
    typedef CanonicalCoefficientType<P,PR> CoefficientType;
    typedef CanonicalErrorType<P,PRE> ErrorType;
    typedef Interval<FloatUpperBound<PR>> RangeType;
    typedef FloatError<PR> NormType;

    typedef CanonicalNumericType<P,PR,PRE> NumberModelType;
    typedef ScalarFunctionModelInterface<P,D,PR,PRE> FunctionModelInterfaceType;
    typedef FunctionModelFactoryInterface<P,PR,PRE> FunctionModelFactoryInterfaceType;
  public:
    virtual CoefficientType const& value() const = 0;
    virtual CoefficientType const gradient_value(SizeType i) const = 0;
    virtual ErrorType const& error() const = 0;
    virtual Void set_error(const ErrorType& e) = 0;
    virtual Void clobber() = 0;

    virtual RangeType range() const = 0;
    virtual NormType const _norm() const = 0;

    virtual FunctionModelInterfaceType* _apply(OperatorCode op) const = 0;
    virtual NumberModelType _unchecked_evaluate(const Vector<NumberModelType>& x) const = 0;
    virtual FunctionModelInterfaceType* _partial_evaluate(SizeType j, const NumberModelType& c) const = 0;

    virtual FunctionModelFactoryInterfaceType* _factory() const = 0;
    virtual FunctionModelInterfaceType* _clone() const = 0;
    virtual FunctionModelInterfaceType* _create() const = 0;
    virtual FunctionModelInterfaceType* _embed(const DomainType& d1, const DomainType& d2) const = 0;
    virtual FunctionModelInterfaceType* _restriction(const DomainType& d) const = 0;

    virtual FunctionModelInterfaceType* _derivative(ElementIndexType<D> j) const = 0;
    virtual FunctionModelInterfaceType* _antiderivative(SizeType j) const = 0;
    virtual FunctionModelInterfaceType* _antiderivative(SizeType j, NumberModelType c) const = 0;

    virtual Boolean _refines(const FunctionModelInterfaceType& f) const = 0;
    virtual Boolean _inconsistent(const FunctionModelInterfaceType& f) const = 0;
    virtual FunctionModelInterfaceType* _refinement(const FunctionModelInterfaceType& f) const = 0;

    virtual Void _iadd(const NumberModelType& c) = 0;
    virtual Void _imul(const NumberModelType& c) = 0;
    virtual Void _isma(const NumberModelType& c, const FunctionModelInterfaceType& f) = 0;
    virtual Void _ifma(const FunctionModelInterfaceType& f1, const FunctionModelInterfaceType& f2) = 0;
};


template<class P, class D, class PR, class PRE> class FunctionModelInterface<P,D,BoxDomainType,PR,PRE>
    : public virtual FunctionInterface<P,D,BoxDomainType>
{
    static_assert(IsSame<D,IntervalDomainType>::value or IsSame<D,BoxDomainType>::value,"");
    typedef BoxDomainType C;
  public:
    typedef D DomainType;
    typedef C CodomainType;
    
    typedef CanonicalCoefficientType<P,PR> CoefficientType;
    typedef CanonicalErrorType<P,PRE> ErrorType;
    typedef Box<Interval<FloatUpperBound<PR>>> RangeType;
    typedef FloatError<PR> NormType;

    typedef CanonicalNumericType<P,PR,PRE> NumberModelType;
    typedef VectorFunctionModelInterface<P,D,PR,PRE> FunctionModelInterfaceType;
    typedef FunctionModelFactoryInterface<P,PR,PRE> FunctionModelFactoryInterfaceType;

    typedef ScalarFunctionModelInterface<P,D,PR,PRE> ScalarFunctionModelInterfaceType;
    typedef VectorFunctionModelInterface<P,D,PR,PRE> VectorFunctionModelInterfaceType;
  public:
    virtual RangeType const range() const = 0;
    virtual Vector<ErrorType> const errors() const = 0;
    virtual ErrorType const error() const = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual FunctionModelFactoryInterfaceType* _factory() const = 0;
    virtual FunctionModelInterfaceType* _clone() const = 0;
    virtual Void _set(SizeType, ScalarFunctionModelInterfaceType const&) = 0;
    virtual ScalarFunctionModelInterfaceType* _get(SizeType) const = 0;
    virtual FunctionModelInterfaceType* _embed(const DomainType& d1, const DomainType& d2) const = 0;
    virtual FunctionModelInterfaceType* _restriction(const DomainType& d) const = 0;
    virtual FunctionModelInterfaceType* _join(const FunctionModelInterfaceType& f2) const = 0;
    virtual FunctionModelInterfaceType* _combine(const FunctionModelInterfaceType& f2) const = 0;
    virtual Void _adjoin(const ScalarFunctionModelInterfaceType& f2) = 0;
    virtual Vector<NumberModelType> _unchecked_evaluate(const Vector<NumberModelType>& x) const = 0;
    virtual ScalarFunctionModelInterfaceType* _compose(const ScalarMultivariateFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterfaceType* _compose(const VectorMultivariateFunctionInterface<P>& f) const = 0;
    virtual ScalarFunctionModelInterfaceType* _unchecked_compose(const ScalarMultivariateFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterfaceType* _unchecked_compose(const VectorMultivariateFunctionInterface<P>& f) const = 0;
    virtual FunctionModelInterfaceType* _partial_evaluate(SizeType j, const NumberModelType& c) const = 0;
    virtual Void restrict(const DomainType& d) = 0;
};


template<class P, class PR, class PRE> class FunctionModelFactoryInterface
{
    typedef BoxDomainType VD;
    typedef IntervalDomainType SD;
    friend class FunctionModelFactory<P,PR,PRE>;
  public:
    typedef SD ScalarDomainType;
    typedef VD VectorDomainType;
    
    typedef Number<P> NumberType;
    typedef ScalarFunctionInterface<P,VD> ScalarMultivariateFunctionInterfaceType;
    typedef VectorFunctionInterface<P,VD> VectorMultivariateFunctionInterfaceType;
    typedef ScalarFunction<P,VD> ScalarMultivariateFunctionType;
    typedef VectorFunction<P,VD> VectorMultivariateFunctionType;

    typedef CanonicalNumericType<P,PR,PRE> NumberModelType;
    typedef ScalarFunctionModelInterface<P,VD,PR,PRE> ScalarMultivariateFunctionModelInterfaceType;
    typedef VectorFunctionModelInterface<P,VD,PR,PRE> VectorMultivariateFunctionModelInterfaceType;
    typedef ScalarFunctionModel<P,VD,PR,PRE> ScalarMultivariateFunctionModelType;
    typedef VectorFunctionModel<P,VD,PR,PRE> VectorMultivariateFunctionModelType;

    typedef FunctionModelFactoryInterface<P,PR,PRE> FunctionModelFactoryInterfaceType;
  public:
    virtual ~FunctionModelFactoryInterface() = default;
    virtual FunctionModelFactoryInterfaceType* clone() const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
    friend OutputStream& operator<<(OutputStream& os, FunctionModelFactoryInterfaceType const& factory) { factory._write(os); return os; }
  private:
    virtual NumberModelType _create(const NumberType& number) const = 0;
    virtual ScalarMultivariateFunctionModelInterfaceType* _create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterfaceType& function) const = 0;
    virtual VectorMultivariateFunctionModelInterfaceType* _create(const VectorDomainType& domain, const VectorMultivariateFunctionInterfaceType& function) const = 0;

    virtual ScalarMultivariateFunctionModelInterfaceType* _create_zero(const VectorDomainType& domain) const = 0;
    virtual ScalarMultivariateFunctionModelInterfaceType* _create_constant(const VectorDomainType& domain, const NumberType& value) const = 0;
    virtual ScalarMultivariateFunctionModelInterfaceType* _create_coordinate(const VectorDomainType& domain, SizeType index) const = 0;
    virtual VectorMultivariateFunctionModelInterfaceType* _create_zeros(SizeType result_size, const VectorDomainType& domain) const = 0;
    virtual VectorMultivariateFunctionModelInterfaceType* _create_constants(const VectorDomainType& domain, const Vector<NumberType>& values) const = 0;
    virtual VectorMultivariateFunctionModelInterfaceType* _create_identity(const VectorDomainType& domain) const = 0;
  public:
    NumberModelType create(const NumberType& number) const {
        return NumberModelType(this->_create(number)); }
    ScalarMultivariateFunctionModelType create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterfaceType& function) const {
        return ScalarMultivariateFunctionModelType(this->_create(domain,function)); }
    VectorMultivariateFunctionModelType create(const VectorDomainType& domain, const VectorMultivariateFunctionInterfaceType& function) const {
        return VectorMultivariateFunctionModelType(this->_create(domain,function)); }

    NumberModelType create_number(const NumberType& number) const {
        return NumberModelType(this->_create(number)); }

    ScalarMultivariateFunctionModelType create_zero(const VectorDomainType& domain) const {
        return ScalarMultivariateFunctionModelType(_create_zero(domain)); }
    ScalarMultivariateFunctionModelType create_constant(const VectorDomainType& domain, const NumberType& value) const {
        return ScalarMultivariateFunctionModelType(_create_constant(domain,value)); }
    ScalarMultivariateFunctionModelType create_constant(const VectorDomainType& domain, const NumberModelType& value) const {
        return ScalarMultivariateFunctionModelType(_create_constant(domain,NumberType(value))); }
    ScalarMultivariateFunctionModelType create_coordinate(const VectorDomainType& domain, SizeType index) const {
        return ScalarMultivariateFunctionModelType(_create_coordinate(domain,index)); }
    VectorMultivariateFunctionModelType create_zeros(SizeType result_size_, const VectorDomainType& domain) const {
        return VectorMultivariateFunctionModelType(_create_zeros(result_size_,domain)); }
    VectorMultivariateFunctionModelType create_constants(const VectorDomainType& domain, const Vector<NumberType>& values) const {
        return VectorMultivariateFunctionModelType(_create_constants(domain,values)); }
    VectorMultivariateFunctionModelType create_identity(const VectorDomainType& domain) const {
        return VectorMultivariateFunctionModelType(_create_identity(domain)); }

    ScalarMultivariateFunctionModelType create_identity(const ScalarDomainType& domain) const;
};


} // namespace Ariadne

#endif
