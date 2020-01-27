/***************************************************************************
 *            function/function_model_mixin.hpp
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

/*! \file function/function_model_mixin.hpp
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
#include "../algebra/algebra_mixin.hpp"

#include "../function/domain.hpp"

#include "../function/function_interface.hpp"
#include "../function/function_patch_interface.hpp"
#include "../function/function_mixin.hpp"
#include "../function/function.hpp"
#include "../function/function_patch_mixin.hpp"

namespace Ariadne {

template<class FM, class P, class D, class C, class PR=DoublePrecision, class PRE=PR> class FunctionModelMixin;
template<class FM, class P, class D, class PR=DoublePrecision, class PRE=PR> using ScalarMultivariateFunctionModelMixin = FunctionModelMixin<FM,P,D,IntervalDomainType,PR,PRE>;
template<class FM, class P, class D, class PR=DoublePrecision, class PRE=PR> using VectorMultivariateFunctionModelMixin = FunctionModelMixin<FM,P,D,BoxDomainType,PR,PRE>;

template<class FM, class P, class D, class PR, class PRE> class FunctionModelMixin<FM,P,D,IntervalDomainType,PR,PRE>
//    : public virtual ScalarFunctionModelInterface<P,D,PR,PRE>
//    , public ScalarFunctionMixin<FM,P,D>
//    , public ElementaryAlgebraMixin<FM,CanonicalNumericType<P,PR,PRE>>
//    , public ElementaryAlgebraMixin<FM,Number<P>>
    : public virtual ScalarFunctionModelInterface<P,D,PR,PRE>
    , public ScalarFunctionPatchMixin<FM,P,D>
    , public ElementaryAlgebraMixin<FM,CanonicalNumericType<P,PR,PRE>>
{
    using C=IntervalDomainType;
  public:
    typedef typename ElementTraits<D>::SizeType ArgumentSizeType;
    typedef typename ElementTraits<C>::SizeType ResultSizeType;
    typedef typename ElementTraits<D>::IndexType ArgumentIndexType;
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = ElementTraits<C>::template Type<X>;

    typedef ElementTraits<C>::BoundedRangeType GenericRangeType;
    typedef PositiveValidatedUpperNumber GenericNormType;
    typedef PositiveValidatedUpperNumber GenericErrorType;

    typedef Interval<FloatUpperBound<PR>> RangeType;
    typedef FloatError<PR> NormType;
    typedef CanonicalErrorType<P,PRE> ErrorType;

    typedef Number<P> Y;
    typedef CanonicalNumericType<P,PR,PRE> X;
  private: public:
    using ElementaryAlgebraMixin<FM,X>::_apply;
  public:

    ArgumentSizeType argument_size() const override {
        return this->domain().dimension(); }
    ResultSizeType result_size() const override {
        return this->codomain().dimension(); }

    virtual GenericRangeType const _generic_range() const override {
        return GenericRangeType(this->_range()); }
    virtual GenericNormType const _generic_norm() const override {
        return GenericNormType(this->_norm()); }


//    virtual Scalar<PositiveValidatedUpperNumber> const _errors() const override { return static_cast<FM const*>(this)->error(); }
//    virtual PositiveValidatedUpperNumber const _error() const override { return static_cast<FM const*>(this)->error(); }
//    virtual ExactNumber const _value() const override { return static_cast<FM const*>(this)->value(); }
//    virtual ExactNumber const _gradient_value(ArgumentIndexType k) const override { return this->gradient_value(k); }

    virtual Result<Number<P>> _unchecked_evaluate(Argument<Number<P>> const& x) const override {
        return this->_unchecked_evaluate(Argument<CanonicalNumericType<P,PR,PRE>>(x,static_cast<FM const&>(*this).precision())); }
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _partial_evaluate(SizeType j, Number<P> const& x) const override {
        return this->_partial_evaluate(j,CanonicalNumericType<P,PR,PRE>(x,static_cast<FM const&>(*this).precision())); }
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _antiderivative(ElementIndexType<D> j, Number<P> c) const override {
        return this->_antiderivative(j,CanonicalNumericType<P,PR,PRE>(c,static_cast<FM const&>(*this).precision())); }

    ScalarFunctionModelInterface<P,D,PR,PRE>* _clone() const override {
        return new FM(static_cast<const FM&>(*this)); }

    ScalarFunctionModelInterface<P,D,PR,PRE>* _create_copy() const override {
        return new FM(static_cast<const FM&>(*this)); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _create_zero() const override {
        return new FM(factory(static_cast<FM const&>(*this)).create_zero()); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _create_constant(CanonicalNumericType<P,PR,PRE> const& c) const override {
        return new FM(factory(static_cast<FM const&>(*this)).create_constant(c)); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _create_constant(Number<P> const& c) const override {
        return new FM(factory(static_cast<FM const&>(*this)).create_constant(c)); }

    RangeType const _range() const override {
         return static_cast<const FM&>(*this).range(); }
    ErrorType const _error() const override {
        return static_cast<const FM&>(*this).error(); }
    NormType const _norm() const override {
        return static_cast<const FM&>(*this).norm(); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _antiderivative(SizeType j) const override {
        return new FM(antiderivative(static_cast<const FM&>(*this),j)); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _antiderivative(SizeType j, CanonicalNumericType<P,PR,PRE> c) const override {
        return new FM(antiderivative(static_cast<const FM&>(*this),j,c)); }
     ScalarFunctionModelInterface<P,D,PR,PRE>* _restriction(const BoxDomainType& d) const override {
        return new FM(restriction(static_cast<const FM&>(*this),d)); }
    CanonicalNumericType<P,PR,PRE> _unchecked_evaluate(const Vector<CanonicalNumericType<P,PR,PRE>>& x) const override {
        return unchecked_evaluate(static_cast<const FM&>(*this),x); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _partial_evaluate(SizeType j, const CanonicalNumericType<P,PR,PRE>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FM&>(*this),j,c)); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _embed(const BoxDomainType& d1, const BoxDomainType& d2) const override {
        return new FM(embed(d1,static_cast<const FM&>(*this),d2)); }

    Boolean _refines(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return refines(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f)); }
    Boolean _inconsistent(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return inconsistent(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f)); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _refinement(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return new FM(refinement(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }

    OutputStream& _write(OutputStream& os) const override {
        return os << static_cast<FM const&>(*this); }
};


template<class FM, class P, class D, class PR, class PRE> class FunctionModelMixin<FM,P,D,BoxDomainType,PR,PRE>
    : public virtual VectorFunctionModelInterface<P,D,PR,PRE>
    , public VectorFunctionPatchMixin<FM,P,D>
{
    using C=BoxDomainType;
  public:
    typedef typename Element<FM>::Type ScalarMultivariateFunctionType;
    typedef Box<Interval<FloatUpperBound<PR>>> RangeType;
    typedef CanonicalErrorType<P,PRE> ErrorType;
    typedef FloatError<PR> NormType;
    typedef ElementTraits<C>::BoundedRangeType GenericRangeType;
    typedef PositiveValidatedUpperNumber GenericErrorType;
    typedef PositiveValidatedUpperNumber GenericNormType;
    typedef typename ElementTraits<D>::SizeType ArgumentSizeType;
    typedef typename ElementTraits<C>::SizeType ResultSizeType;
    template<class X> using Argument=typename ElementTraits<D>::template Type<X>;
    template<class X> using Result=ElementTraits<C>::template Type<X>;
  public:
    ArgumentSizeType argument_size() const override {
        return this->domain().dimension(); }
    ResultSizeType result_size() const override {
        return this->codomain().dimension(); }

    virtual GenericRangeType const _generic_range() const override { return GenericRangeType(this->_range()); }
    virtual GenericNormType const _generic_norm() const override { return GenericNormType(this->_norm()); }
//    virtual Vector<PositiveValidatedUpperNumber> const _errors() const override { return this->errors(); }
//    virtual PositiveValidatedUpperNumber const _error() const override { return static_cast<FM const*>(this)->error(); }
//    virtual ExactNumber const _value() const override { return static_cast<FM const*>(this)->value(); }
//    virtual ExactNumber const _gradient_value(ArgumentIndexType k) const override { return this->gradient_value(k); }
//    virtual PositiveValidatedUpperNumber const _norm() const override { return static_cast<FM const&>(*this).norm(); }
//    virtual BoxValidatedRangeType const _range() const override {
//        return BoxValidatedRangeType(static_cast<FM const&>(*this).range()); }

    virtual Result<Number<P>> _unchecked_evaluate(Argument<Number<P>> const& x) const override {
        return this->_unchecked_evaluate(Argument<CanonicalNumericType<P,PR,PRE>>(x,static_cast<FM const*>(this)->precision())); }
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _partial_evaluate(SizeType j, Number<P> const& x) const override {
        return this->_partial_evaluate(j, CanonicalNumericType<P,PR,PRE>(x,static_cast<FM const*>(this)->precision())); }
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _antiderivative(ElementIndexType<D> j) const override {
        return heap_copy(antiderivative(static_cast<FM const&>(*this),j)); }
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _antiderivative(ElementIndexType<D> j, Number<P> c) const override {
        return heap_copy(antiderivative(static_cast<FM const&>(*this),j,c)); }

    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _clone() const override { return new FM(static_cast<const FM&>(*this)); }

    ErrorType const _error() const override {
        return static_cast<const FM&>(*this).error(); }
    Vector<ErrorType> const _errors() const override {
        return static_cast<const FM&>(*this).errors(); }

    virtual Void _set(SizeType i, const ScalarFunctionPatchInterface<P,D>& sf) override {
        if(!dynamic_cast<const typename FM::ScalarMultivariateFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorMultivariateFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<FM&>(*this).FM::set(i,dynamic_cast<const ScalarMultivariateFunctionType&>(sf)); }
    Void _adjoin(const ScalarFunctionPatchInterface<P,D>& f) override {
        static_cast<FM&>(*this).FM::adjoin(dynamic_cast<const ScalarMultivariateFunctionType&>(f)); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _join(const VectorFunctionPatchInterface<P,D>& f) const override {
        return heap_copy(join(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _combine(const VectorFunctionPatchInterface<P,D>& f) const override {
        return heap_copy(combine(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }

    virtual Void _set(SizeType i, const ScalarFunctionModelInterface<P,D,PR,PRE>& sf) override {
        if(!dynamic_cast<const typename FM::ScalarMultivariateFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorMultivariateFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<FM&>(*this).FM::set(i,dynamic_cast<const ScalarMultivariateFunctionType&>(sf)); }
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _derivative(SizeType j) const override {
        ARIADNE_NOT_IMPLEMENTED; }
    RangeType const _range() const override {
         return static_cast<const FM&>(*this).range(); }
    NormType const _norm() const override {
         return static_cast<const FM&>(*this).norm(); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _embed(const BoxDomainType& d1, const BoxDomainType& d2) const override {
        return heap_copy(embed(d1,static_cast<const FM&>(*this),d2)); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _restriction(const BoxDomainType& d) const override {
        return new FM(restriction(static_cast<const FM&>(*this),d)); }
    Void _adjoin(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) override {
        static_cast<FM&>(*this).FM::adjoin(dynamic_cast<const ScalarMultivariateFunctionType&>(f)); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _join(const VectorFunctionModelInterface<P,D,PR,PRE>& f) const override {
        return heap_copy(join(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _combine(const VectorFunctionModelInterface<P,D,PR,PRE>& f) const override {
        return heap_copy(combine(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    Vector<CanonicalNumericType<P,PR,PRE>> _unchecked_evaluate(const Vector<CanonicalNumericType<P,PR,PRE>>& x) const override {
        return unchecked_evaluate(static_cast<const FM&>(*this),x); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _compose(const ScalarMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(compose(f,static_cast<const FM&>(*this))); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _compose(const VectorMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(compose(f,static_cast<const FM&>(*this))); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _unchecked_compose(const ScalarMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const ScalarMultivariateFunctionType&>(f),static_cast<const FM&>(*this))); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _unchecked_compose(const VectorMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const FM&>(f),static_cast<const FM&>(*this))); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _partial_evaluate(SizeType j, const CanonicalNumericType<P,PR,PRE>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FM&>(*this),j,c)); }
};


template<class FCTRY, class P, class PR, class PRE> class FunctionModelFactoryMixin
    : public FunctionModelFactoryInterface<P,PR,PRE>
{
    typedef BoxDomainType VD;
    typedef IntervalDomainType SD;
    friend class FunctionModelFactory<P,PR,PRE>;
  public:
    typedef VD VectorDomainType;
    typedef SD ScalarDomainType;

    virtual FunctionModelFactoryInterface<P,PR,PRE>* clone() const override { return new FCTRY(this->upcast()); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << this->upcast(); }
/*
    CanonicalNumericType<P,PR,PRE> create(const Number<P>& number) const;
    ScalarFunctionModel<P,D,PR,PRE> create(const BoxDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const;
    VectorFunctionModel<P,D,PR,PRE> create(const BoxDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const;
    ScalarFunctionModel<P,D,PR,PRE> create_zero(const BoxDomainType& domain) const;
    VectorFunctionModel<P,D,PR,PRE> create_zeros(SizeType result_size, const BoxDomainType& domain) const;
    ScalarFunctionModel<P,D,PR,PRE> create_constant(const BoxDomainType& domain, const Number<P>& value) const;
    ScalarFunctionModel<P,D,PR,PRE> create_constant(const BoxDomainType& domain, const CanonicalNumericType<P,PR,PRE>& value) const;
    VectorFunctionModel<P,D,PR,PRE> create_constants(const BoxDomainType& domain, const Vector<Number<P>>& values) const;
    VectorFunctionModel<P,D,PR,PRE> create_constants(const BoxDomainType& domain, const Vector<CanonicalNumericType<P,PR,PRE>>& values) const;
    ScalarFunctionModel<P,D,PR,PRE> create_coordinate(const BoxDomainType& domain, SizeType index) const;
    ScalarFunctionModel<P,D,PR,PRE> create_identity(const IntervalDomainType& domain) const;
    VectorFunctionModel<P,D,PR,PRE> create_identity(const BoxDomainType& domain) const;
    CanonicalNumericType<P,PR,PRE> create_number(const Number<P>& number) const;
*/
  private:
    template<class T> static inline T* heap_move(T&& t) { return new T(std::forward<T>(t)); }
    inline FCTRY const& upcast() const { return static_cast<FCTRY const&>(*this); }
    virtual CanonicalNumericType<P,PR,PRE> _create(const Number<P>& number) const override {
        return this->upcast().create(number); }
    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create(const VectorDomainType& domain, const ScalarFunctionInterface<P,VD>& function) const override {
        return heap_move(this->upcast().create(domain,function)); };
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create(const VectorDomainType& domain, const VectorFunctionInterface<P,VD>& function) const override {
        return heap_move(this->upcast().create(domain,function)); };

    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create_zero(const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_zero(domain)); };
    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create_constant(const VectorDomainType& domain, const Number<P>& value) const override {
        return heap_move(this->upcast().create_constant(domain,value)); };
    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create_coordinate(const VectorDomainType& domain, SizeType index) const override {
        return heap_move(this->upcast().create_coordinate(domain,index)); };
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create_zeros(SizeType rsize, const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_zeros(rsize,domain)); };
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const override {
        return heap_move(this->upcast().create_constants(domain,values)); };
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create_projection(const VectorDomainType& domain, Range indices) const override {
        return heap_move(this->upcast().create_projection(domain,indices)); };
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create_identity(const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_identity(domain)); };
};


} // namespace Ariadne

#endif
