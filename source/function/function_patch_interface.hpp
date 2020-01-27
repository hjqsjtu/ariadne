/***************************************************************************
 *            function/function_patch_interface.hpp
 *
 *  Copyright  2011-20  Pieter Collins
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

/*! \file function/function_patch_interface.hpp
 *  \brief Interface for functions on bounded sets.
 */

#ifndef ARIADNE_FUNCTION_PATCH_INTERFACE_HPP
#define ARIADNE_FUNCTION_PATCH_INTERFACE_HPP

#include "../numeric/number.decl.hpp"
#include "../function/function.decl.hpp"

#include "../numeric/operators.hpp"
#include "../numeric/numeric.hpp"

#include "../algebra/algebra_interface.hpp"

#include "../function/domain.hpp"
#include "../function/function_interface.hpp"


namespace Ariadne {

typedef Interval<FloatDPUpperBound> IntervalRangeType;

template<class P, class D, class C> class FunctionPatch;

template<class P> class FunctionPatchFactoryInterface;
template<class P, class D> class FunctionPatchCreatorInterface;

template<class P, class D, class C> class FunctionPatchInterface;

template<class P, class D> using ScalarFunctionPatchInterface=FunctionPatchInterface<P,D,IntervalDomainType>;
template<class P, class D> using VectorFunctionPatchInterface=FunctionPatchInterface<P,D,BoxDomainType>;
template<class P, class C> using UnivariateFunctionPatchInterface=FunctionPatchInterface<P,IntervalDomainType,C>;
template<class P, class C> using MultivariateFunctionPatchInterface=FunctionPatchInterface<P,BoxDomainType,C>;
template<class P> using ScalarUnivariateFunctionPatchInterface=FunctionPatchInterface<P,IntervalDomainType,IntervalDomainType>;
template<class P> using ScalarMultivariateFunctionPatchInterface=FunctionPatchInterface<P,BoxDomainType,IntervalDomainType>;
template<class P> using VectorUnivariateFunctionPatchInterface=FunctionPatchInterface<P,IntervalDomainType,BoxDomainType>;
template<class P> using VectorMultivariateFunctionPatchInterface=FunctionPatchInterface<P,BoxDomainType,BoxDomainType>;

template<class P, class D, class C> class FunctionPatchAlgebraInterface;

template<class P, class D> class FunctionPatchAlgebraInterface<P,D,IntervalDomainType>
    : public virtual ElementaryAlgebraInterface<Number<P>>
{ };

template<class P, class D> class FunctionPatchAlgebraInterface<P,D,BoxDomainType>
{
    using C=BoxDomainType; using SC=IntervalDomainType;
  public:
    typedef typename ElementTraits<C>::IndexType ResultIndexType;
    virtual Void _set(ResultIndexType, ScalarFunctionPatchInterface<P,D> const&) = 0;
    virtual ScalarFunctionPatchInterface<P,D>* _get(ResultIndexType) const = 0;
    virtual Void _adjoin(const ScalarFunctionPatchInterface<P,D>& f2) = 0;
    virtual VectorFunctionPatchInterface<P,D>* _join(const VectorFunctionPatchInterface<P,D>& f2) const = 0;
    virtual VectorFunctionPatchInterface<P,D>* _combine(const VectorFunctionPatchInterface<P,D>& f2) const = 0;
};

template<class P, class D, class C> class FunctionPatchInterface
    : public virtual FunctionInterface<P,D,C>
    , public virtual FunctionPatchAlgebraInterface<P,D,C>
{
    static_assert(IsSame<D,IntervalDomainType>::value or IsSame<D,BoxDomainType>::value,"");
    static_assert(IsSame<C,IntervalDomainType>::value or IsSame<C,BoxDomainType>::value,"");
  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename ElementTraits<C>::BoundedRangeType RangeType;
    typedef PositiveValidatedUpperNumber NormType;
    typedef PositiveValidatedUpperNumber ErrorType;
    typedef ExactNumber CoefficientType;
    typedef typename ElementTraits<D>::IndexType ArgumentIndexType;

    template<class Y> using Argument = typename ElementTraits<DomainType>::template Type<Y>;
    template<class Y> using Result = typename ElementTraits<CodomainType>::template Type<Y>;
  public:
//    virtual Result<ErrorType> const _errors() const = 0;
//    virtual CoefficientType const _value() const = 0;
//    virtual CoefficientType const _gradient_value(ArgumentIndexType) const = 0;
//    virtual ErrorType const _error() const = 0;
    virtual NormType const _generic_norm() const = 0;
    virtual RangeType const _generic_range() const = 0;

    virtual FunctionPatchFactoryInterface<P>* _patch_factory() const { ARIADNE_NOT_IMPLEMENTED; }

    virtual DomainType const domain() const override = 0;
    virtual CodomainType const codomain() const override = 0;

    virtual Result<Number<P>> _unchecked_evaluate(const Argument<Number<P>>& x) const = 0;
    virtual FunctionPatchInterface<P,D,C>* _partial_evaluate(SizeType j, const Number<P>& c) const = 0;

    virtual FunctionPatchInterface<P,D,C>* _clone() const override = 0;
    virtual FunctionPatchInterface<P,D,C>* _create() const = 0;
    virtual FunctionPatchInterface<P,D,C>* _embed(const DomainType& d1, const DomainType& d2) const = 0;
    virtual FunctionPatchInterface<P,D,C>* _restriction(const DomainType& d) const = 0;

    virtual FunctionPatchInterface<P,D,C>* _derivative(ElementIndexType<D> j) const override = 0;
    virtual FunctionPatchInterface<P,D,C>* _antiderivative(ElementIndexType<D> j) const = 0;
    virtual FunctionPatchInterface<P,D,C>* _antiderivative(ElementIndexType<D> j, Number<P> c) const = 0;

/*
    virtual ScalarFunctionPatchInterface<P,D>* _compose(const ScalarMultivariateFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionPatchInterface<P,D>* _compose(const VectorMultivariateFunctionInterface<P>& f) const = 0;
    virtual ScalarFunctionPatchInterface<P,D>* _unchecked_compose(const ScalarMultivariateFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionPatchInterface<P,D>* _unchecked_compose(const VectorMultivariateFunctionInterface<P>& f) const = 0;
*/
    friend OutputStream& operator<<(OutputStream& os, FunctionPatchInterface<P,D,C> const& f) {
        return os << static_cast<FunctionInterface<P,D,C>const&>(f); }
};


template<class P> class FunctionPatchFactoryInterface
{
    typedef BoxDomainType VD;
    typedef IntervalDomainType SD;
  public:
    typedef SD ScalarDomainType;
    typedef VD VectorDomainType;
    virtual ~FunctionPatchFactoryInterface<P>() = default;
    virtual FunctionPatchFactoryInterface<P>* clone() const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
    friend OutputStream& operator<<(OutputStream& os, FunctionPatchFactoryInterface<P> const& factory) { factory._write(os); return os; }
  private: public:
    virtual ScalarFunctionPatchInterface<P,VD>* _create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const = 0;
    virtual VectorFunctionPatchInterface<P,VD>* _create(const VectorDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const = 0;

    virtual ScalarFunctionPatchInterface<P,VD>* _create_zero(const VectorDomainType& domain) const = 0;
    virtual ScalarFunctionPatchInterface<P,VD>* _create_constant(const VectorDomainType& domain, const Number<P>& value) const = 0;
    virtual ScalarFunctionPatchInterface<P,VD>* _create_coordinate(const VectorDomainType& domain, SizeType index) const = 0;
    virtual VectorFunctionPatchInterface<P,VD>* _create_zeros(SizeType rsize, const VectorDomainType& domain) const = 0;
    virtual VectorFunctionPatchInterface<P,VD>* _create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const = 0;
    virtual VectorFunctionPatchInterface<P,VD>* _create_projection(const VectorDomainType& domain, Range indices) const = 0;
    virtual VectorFunctionPatchInterface<P,VD>* _create_identity(const VectorDomainType& domain) const = 0;
  public:
    Number<P> create(const Number<P>& number) const {
        return Number<P>(this->_create(number)); }
    ScalarFunctionPatch<P,VD> create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const {
        return ScalarFunctionPatch<P,VD>(this->_create(domain,function)); }
    VectorFunctionPatch<P,VD> create(const VectorDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const {
        return VectorFunctionPatch<P,VD>(this->_create(domain,function)); }

    Number<P> create_number(const Number<P>& number) const {
        return Number<P>(this->_create(number)); }

    ScalarFunctionPatch<P,VD> create_zero(const VectorDomainType& domain) const {
        return ScalarFunctionPatch<P,VD>(_create_zero(domain)); }
    ScalarFunctionPatch<P,VD> create_constant(const VectorDomainType& domain, const Number<P>& value) const {
        return ScalarFunctionPatch<P,VD>(_create_constant(domain,value)); }
    ScalarFunctionPatch<P,VD> create_coordinate(const VectorDomainType& domain, SizeType index) const {
        return ScalarFunctionPatch<P,VD>(_create_coordinate(domain,index)); }
    VectorFunctionPatch<P,VD> create_zeros(SizeType rsize, const VectorDomainType& domain) const {
        return VectorFunctionPatch<P,VD>(_create_zeros(rsize,domain)); }
    VectorFunctionPatch<P,VD> create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const {
        return VectorFunctionPatch<P,VD>(_create_constants(domain,values)); }
    VectorFunctionPatch<P,VD> create_projection(const VectorDomainType& domain, Range indices) const {
        return VectorFunctionPatch<P,VD>(_create_projection(domain,indices)); }
    VectorFunctionPatch<P,VD> create_identity(const VectorDomainType& domain) const {
        return VectorFunctionPatch<P,VD>(_create_identity(domain)); }

    ScalarFunctionPatch<P,VD> create_identity(const ScalarDomainType& domain) const;
};


} // namespace Ariadne

#endif
