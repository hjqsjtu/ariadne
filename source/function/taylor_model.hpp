/***************************************************************************
 *            taylor_model.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file taylor_model.hpp
 *  \brief ApproximateTag functions on a is_bounded domain with a sparse representation.
 */

#ifndef ARIADNE_TAYLOR_MODEL_HPP
#define ARIADNE_TAYLOR_MODEL_HPP

#include <map>

#include "../utility/macros.hpp"
#include "../utility/declarations.hpp"
#include "../utility/array.hpp"
#include "../utility/pointer.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/covector.hpp"
#include "../algebra/multi_index.hpp"
#include "../algebra/expansion.hpp"
#include "../algebra/sweeper.hpp"
#include "../algebra/operations.hpp"
#include "../algebra/evaluate.hpp"
#include "../function/domain.hpp"
#include "../function/scaling.hpp"

namespace Ariadne {

template<class T1, class T2> struct Product;

template<class P, class F, class I> class TaylorModel;
template<class P, class F> using UnivariateTaylorModel = TaylorModel<P,F,UniIndex>;
template<class P, class F> using MultivariateTaylorModel = TaylorModel<P,F,MultiIndex>;
template<class F> using ValidatedTaylorModel = TaylorModel<ValidatedTag,F,MultiIndex>;
template<class F> using ApproximateTaylorModel = TaylorModel<ApproximateTag,F,MultiIndex>;
using ValidatedTaylorModelDP = TaylorModel<ValidatedTag,FloatDP,MultiIndex>;
using ValidatedTaylorModelMP = TaylorModel<ValidatedTag,FloatMP,MultiIndex>;
using ApproximateTaylorModelDP = TaylorModel<ApproximateTag,FloatDP,MultiIndex>;
using ApproximateTaylorModelMP = TaylorModel<ApproximateTag,FloatMP,MultiIndex>;

template<class P, class F, class I> struct IsScalar< TaylorModel<P,F,I> > { static const Bool value = true; };

class IntersectionException;

class IntersectionException : public std::runtime_error {
  public:
    IntersectionException(const StringType& what) : std::runtime_error(what) { }
};

template<class P, class F, class I> struct ModelNumericTypedef;
template<class F, class I> struct ModelNumericTypedef<ValidatedTag,F,I> { typedef FloatBounds<typename F::PrecisionType> Type; };
template<class F, class I> struct ModelNumericTypedef<ApproximateTag,F,I> { typedef FloatApproximation<typename F::PrecisionType> Type; };
template<class P, class F, class I> using ModelNumericType = typename ModelNumericTypedef<P,F,I>::Type;

template<class P, class F, class I> struct AlgebraOperations<TaylorModel<P,F,I>>
    : NormedAlgebraOperations<TaylorModel<P,F,I>>
{
    typedef ModelNumericType<P,F,I> X;
    using NormedAlgebraOperations<TaylorModel<P,F,I>>::apply;
    static TaylorModel<P,F,I> apply(Pos,TaylorModel<P,F,I> tm);
    static TaylorModel<P,F,I> apply(Neg,TaylorModel<P,F,I> tm);
    static TaylorModel<P,F,I> apply(Add,TaylorModel<P,F,I> tm, X const& c);
    static TaylorModel<P,F,I> apply(Mul,TaylorModel<P,F,I> tm, X const& c);
    static TaylorModel<P,F,I> apply(Add,TaylorModel<P,F,I> const& tm1, TaylorModel<P,F,I> const& tm2);
    static TaylorModel<P,F,I> apply(Sub,TaylorModel<P,F,I> const& tm1, TaylorModel<P,F,I> const& tm2);
    static TaylorModel<P,F,I> apply(Mul,TaylorModel<P,F,I> const& tm1, TaylorModel<P,F,I> const& tm2);
    static TaylorModel<P,F,I> apply(Min,TaylorModel<P,F,I> const& tm1, TaylorModel<P,F,I> const& tm2);
    static TaylorModel<P,F,I> apply(Max,TaylorModel<P,F,I> const& tm1, TaylorModel<P,F,I> const& tm2);
    static TaylorModel<P,F,I> apply(Abs,TaylorModel<P,F,I> const& tm);
};



/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, ValidatedScalarTaylorFunctionModelDP, ValidatedVectorTaylorFunctionModelDP, TaylorConstrainedImageSet.
 */
template<class F, class I>
class TaylorModel<ValidatedTag,F,I>
    : public DispatchTranscendentalAlgebraOperations<TaylorModel<ValidatedTag,F,I>,CanonicalNumericType<ValidatedTag,typename F::PrecisionType>>
    , public DispatchOrderedAlgebraOperations<TaylorModel<ValidatedTag,F,I>,CanonicalNumericType<ValidatedTag,typename F::PrecisionType>>
    , public DispatchConcreteGenericAlgebraNumberOperations<TaylorModel<ValidatedTag,F,I>,CanonicalNumericType<ValidatedTag,typename F::PrecisionType>,ValidatedNumber>
{
    typedef typename F::PrecisionType PR;
    typedef typename F::PrecisionType PRE;
  public:
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;

    typedef F RawFloatType;
    typedef I IndexType;
    typedef FloatValue<PR> CoefficientType;
    typedef FloatError<PRE> ErrorType;
    typedef FloatError<PR> NormType;
    typedef ReverseLexicographicIndexLess ComparisonType;
    typedef SortedExpansion<IndexType,CoefficientType,ComparisonType> ExpansionType;
    typedef Sweeper<F> SweeperType;

    typedef typename ExpansionType::ArgumentSizeType ArgumentSizeType;
    typedef typename ExpansionType::ArgumentIndexType ArgumentIndexType;

    typedef typename IndexTraits<I>::UnitDomainType DomainType;
    typedef IntervalDomainType CodomainType;
    typedef Interval<FloatUpperBound<PR>> RangeType;

    template<class X> using Argument = typename IndexTraits<I>::template Argument<X>;
    template<class X> using Coargument = typename IndexTraits<I>::template Coargument<X>;

    //! \brief The computational paradigm.
    typedef ValidatedTag Paradigm;
    //! \brief The properties needed to define the TaylorModel calculus.
    typedef Sweeper<F> PropertiesType;

    //! \brief The type used for algebraic operations.
    typedef FloatBounds<PR> NumericType;
    typedef ValidatedNumber GenericNumericType;

    typedef FloatBounds<PR> ValidatedNumericType;
    typedef FloatApproximation<PR> ApproximateNumericType;

    typedef ScalarFunction<ValidatedTag,typename IndexTraits<I>::UnboundedDomainType> FunctionType;
    typedef ScalarFunction<ValidatedTag,typename IndexTraits<I>::UnboundedDomainType> ScalarFunctionType;
    typedef VectorFunction<ValidatedTag,typename IndexTraits<I>::UnboundedDomainType> VectorFunctionType;

    //! \brief The type used for the coefficients.
    typedef CoefficientType ValueType;

    //! \brief An Iterator through the (index,coefficient) pairs of the expansion.
    typedef typename ExpansionType::Iterator Iterator;
    //! \brief A constant Iterator through the (index,coefficient) pairs of the expansion.
    typedef typename ExpansionType::ConstIterator ConstIterator;
  private:
    ExpansionType _expansion;
    ErrorType _error;
    mutable SweeperType _sweeper;
  public:
    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    TaylorModel<ValidatedTag,F,I>();
    //! \brief Construct a TaylorModel in \a as arguments with the given accuracy control.
    TaylorModel<ValidatedTag,F,I>(ArgumentSizeType as, SweeperType swp);
    //! \brief Construct from a map giving the expansion, a constant giving the error, and an accuracy parameter.
    TaylorModel<ValidatedTag,F,I>(const Expansion<IndexType,double>& f, const double& e, SweeperType swp);
    TaylorModel<ValidatedTag,F,I>(const Expansion<IndexType,CoefficientType>& f, const ErrorType& e, SweeperType swp);
    TaylorModel<ValidatedTag,F,I>(const Expansion<IndexType,RawFloatType>& f, const RawFloatType& e, SweeperType swp);
    //! \brief Fast swap with another Taylor model.
    Void swap(TaylorModel<ValidatedTag,F,I>& tm);
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedTag,F,I> create() const;
    //! \brief The zero element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedTag,F,I> create_zero() const;
    //! \brief A constant element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedTag,F,I> create_constant(NumericType c) const;
    TaylorModel<ValidatedTag,F,I> create_constant(GenericNumericType) const;
    //! \brief The \a j<sup>th</sup> coordinate element of the algebra of Taylor models, with the same number of arguments and accuracy parameters.
    TaylorModel<ValidatedTag,F,I> create_coordinate(ArgumentIndexType j) const;
    //! \brief Set to zero.
    TaylorModel<ValidatedTag,F,I> create_ball(ErrorType e) const;
    //! \brief Set to zero.
    Void clear();

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<ValidatedTag,F,I>& operator=(const NumericType& c);
    TaylorModel<ValidatedTag,F,I>& operator=(const GenericNumericType& c);
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct the zero quantity in \a as independent variables.
    static TaylorModel<ValidatedTag,F,I> zero(ArgumentSizeType as, SweeperType swp) {
        TaylorModel<ValidatedTag,F,I> r(as,swp); return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedTag,F,I> constant(ArgumentSizeType as, const ValidatedNumericType& c, SweeperType swp) {
        TaylorModel<ValidatedTag,F,I> r(as,swp); r.set_value(CoefficientType(1,r.precision())); r*=c; return r; }
    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ValidatedTag,F,I> constant(ArgumentSizeType as, const ValidatedNumber& c, SweeperType swp) {
        TaylorModel<ValidatedTag,F,I> r(as,swp); r.set_value(CoefficientType(1,r.precision())); r*=c; return r; }
    //! \brief Construct the quantity with expansion \f$x_j\f$ in \a as independent variables.
    static TaylorModel<ValidatedTag,F,I> coordinate(ArgumentSizeType as, ArgumentIndexType j, SweeperType swp) {
        TaylorModel<ValidatedTag,F,I> r(as,swp); r.set_gradient(j,1); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<ValidatedTag,F,I> error(ArgumentSizeType as, ErrorType e, SweeperType swp) {
        TaylorModel<ValidatedTag,F,I> r(as,swp); r.set_error(e); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<ValidatedTag,F,I> ball(ArgumentSizeType as, ErrorType e, SweeperType swp) {
        TaylorModel<ValidatedTag,F,I> r(as,swp); r.set_error(e); return r; }
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a e
    static TaylorModel<ValidatedTag,F,I> unit_ball(ArgumentSizeType as, SweeperType swp) {
        TaylorModel<ValidatedTag,F,I> r(as,swp); r.set_error(1u); return r; }

    //! \brief Construct the quantity which scales the interval \a codom onto the unit interval.
    static TaylorModel<ValidatedTag,F,I> scaling(ArgumentSizeType as, ArgumentIndexType j, const IntervalDomainType& codom, SweeperType swp);

    //! \brief Return the vector of zero variables of size \a rs in \a as arguments.
    static Vector<TaylorModel<ValidatedTag,F,I>> zeros(ArgumentSizeType rs, ArgumentSizeType as, SweeperType swp);
    //! \brief Return the vector of constants with values \a c in \a as arguments.
    static Vector<TaylorModel<ValidatedTag,F,I>> constants(ArgumentSizeType as, const Vector<ValidatedNumericType>& c, SweeperType swp);
    //! \brief Return the vector of variables on the unit domain.
    static Vector<TaylorModel<ValidatedTag,F,I>> coordinates(ArgumentSizeType as, SweeperType swp);

    //! \brief Return the vector scaling the box \a codom onto the unit box.
    static TaylorModel<ValidatedTag,F,I> scalings(const IntervalDomainType& codom, SweeperType swp);
    static Vector<TaylorModel<ValidatedTag,F,I>> scalings(const Vector<IntervalDomainType>& codom, SweeperType swp);
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    friend Bool same(const TaylorModel<ValidatedTag,F,I>& tm1, const TaylorModel<ValidatedTag,F,I>& tm2) {
        return same(tm1._expansion, tm2._expansion) && same(tm1._error, tm2._error); }

    ValidatedLowerKleenean operator<(const TaylorModel<ValidatedTag,F,I>& sd) const {
        return (sd-*this)>0; }
    //! \brief Comparison with another Taylor model.
    ValidatedLowerKleenean operator>(const TaylorModel<ValidatedTag,F,I>& sd) const {
        return (*this-sd)>0; }

    //! \brief Comparison with a scalar.
    ValidatedLowerKleenean operator<(Int c) const {
        return this->range().upper()<c; }
    //! \brief Comparison with a scalar.
    ValidatedLowerKleenean operator>(Int c) const {
        return this->range().lower()>c; }
    //@}

    //@{
    /*! \name Data access */
    //! \brief The expansion.
    const ExpansionType& expansion() const { return this->_expansion; }
    //! \brief A reference to the expansion.
    ExpansionType& expansion() { return this->_expansion; }
    //! \brief The error of the expansion over the domain.
    const ErrorType& error() const { return this->_error; }
    //! \brief A reference to the error of the expansion over the domain.
    ErrorType& error() { return this->_error; }

    //! \brief The constant term in the expansion.
    const CoefficientType& value() const { return (*this)[IndexType::zero(this->argument_size())]; }
    //! \brief The coefficient of the gradient term \f$df/dx_j\f$.
    const CoefficientType& gradient_value(ArgumentIndexType j) const { return (*this)[IndexType::unit(this->argument_size(),j)]; }
    Coargument<CoefficientType> gradient_value() const;

    //! \brief The constant term in the expansion.
    CoefficientType average() const;
    //! \brief The radius of the smallest ball about a constant function containing the model.
    NormType radius() const;
    //! \brief An over-approximation to the supremum norm.
    NormType norm() const;
    //! \brief An over-approximation to the supremum norm.
    friend typename TaylorModel<ValidatedTag,F,I>::NormType norm(const TaylorModel<ValidatedTag,F,I>& tm) { return tm.norm(); }
    friend typename TaylorModel<ValidatedTag,F,I>::NormType mag(const TaylorModel<ValidatedTag,F,I>& tm) { return tm.norm(); }


    //! \brief A value \c e such that analytic functions are evaluated to a tolerance of \c e. Equal to the sweep threshold.
    F tolerance() const;


    //! \brief Set the constant term in the expansion.
    Void set_value(const CoefficientType& c) {
        this->_expansion.set(IndexType::zero(this->argument_size()),c); }
    //! \brief Set the coefficient of the term \f$df/dx_j\f$.
    Void set_gradient(ArgumentIndexType j, const CoefficientType& c) {
        this->_expansion.set(IndexType::unit(this->argument_size(),j),c); }
    Void set_gradient(ArgumentIndexType j,const Dyadic& c) {
        this->set_gradient(j,CoefficientType(c,this->precision())); }
     //! \brief Set the error of the expansion.
    Void set_error(const ErrorType& ne) { ARIADNE_ASSERT(ne.raw()>=0.0); this->_error=ne; }
    Void set_error(Nat m) { this->_error=m; }

    //! \brief The coefficient of the term in $x^a$.
    const CoefficientType& operator[](const IndexType& a) const { return this->_expansion[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    CoefficientType& operator[](const IndexType& a) { return this->_expansion.at(a); }

    //! \brief An Iterator to the first term in the expansion.
    Iterator begin() { return this->_expansion.begin(); }
    //! \brief A constant Iterator to the first term in the expansion.
    ConstIterator begin() const { return this->_expansion.begin(); }
    //! \brief An Iterator to the end of the expansion.
    Iterator end() { return this->_expansion.end(); }
    //! \brief A constant Iterator to the end of the expansion.
    ConstIterator end() const { return this->_expansion.end(); }

    //! \brief The number of variables in the argument of the quantity.
    ArgumentSizeType argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The maximum degree of terms in the expansion.
    DegreeType degree() const;
    //! \brief The number of nonzero terms in the expansion.
    SizeType number_of_terms() const { return this->_expansion.number_of_terms(); }
    SizeType number_of_nonzeros() const { return this->_expansion.number_of_nonzeros(); }
    //@}

    //@{
    /*! \name Function evaluation. */
    //! \brief The domain of the quantity, always given by \f$[-1,1]^{\mathrm{as}}\f$.
    UnitBox domain() const;
    //! \brief The codomain of the quantity.
    IntervalDomainType codomain() const;
    //! \brief An over-approximation to the range of the quantity.
    RangeType range() const;
    //! \brief Compute the gradient of the expansion with respect to the \a jth variable over the domain.
    RangeType gradient_range(ArgumentIndexType j) const;
    Coargument<RangeType> gradient_range() const;

    //! \brief Evaluate the quantity over the interval of points \a x.
    friend ApproximateNumericType evaluate(const TaylorModel<ValidatedTag,F,I>& f, const Argument<ApproximateNumericType>& x) {
        return TaylorModel<ValidatedTag,F,I>::_evaluate(f,x); }
    friend ValidatedNumericType evaluate(const TaylorModel<ValidatedTag,F,I>& f, const Argument<ValidatedNumericType>& x) {
        return TaylorModel<ValidatedTag,F,I>::_evaluate(f,x); }
    //! \brief Evaluate the gradient over the interval of points \a x.
    friend Coargument<FloatBounds<PR>> gradient(const TaylorModel<ValidatedTag,F,I>& f, const Argument<FloatBounds<PR>>& x) {
        return TaylorModel<ValidatedTag,F,I>::_gradient(f,x); }
    //! \brief Substite \a c for the \a k th variable.
    friend TaylorModel<ValidatedTag,F,I> partial_evaluate(const TaylorModel<ValidatedTag,F,I>& x, ArgumentSizeType k, ValidatedNumericType c) {
        return TaylorModel<ValidatedTag,F,I>::_partial_evaluate(x,k,c); }
    //! \relates TaylorModel<ValidatedTag,F Compose a vector of Taylor models with another.
    friend TaylorModel<ValidatedTag,F,I> compose(const Unscaling& uf, const TaylorModel<ValidatedTag,F,I>& tg) {
        return TaylorModel<ValidatedTag,F,I>::_compose(uf,tg); }
    //! \brief Evaluate the quantity over the interval of points \a x.
    friend TaylorModel<ValidatedTag,F,I> compose(const TaylorModel<ValidatedTag,F,UniIndex>& f, const TaylorModel<ValidatedTag,F,I>& g) {
        return TaylorModel<ValidatedTag,F,I>::_compose(f,g); }
    //! \brief Evaluate the quantity over the interval of points \a x.
    friend TaylorModel<ValidatedTag,F,I> compose(const TaylorModel<ValidatedTag,F,MultiIndex>& f, const Vector<TaylorModel<ValidatedTag,F,I>>& g) {
        return TaylorModel<ValidatedTag,F,I>::_compose(f,g); }
    //! \brief Scale the variable by post-composing with an affine map taking the interval ivl to the unit interval
    friend TaylorModel<ValidatedTag,F,I> compose(const TaylorModel<ValidatedTag,F,I>& tf, const VectorUnscaling& u);

    friend TaylorModel<ValidatedTag,F,I> evaluate(const TaylorModel<ValidatedTag,F,I>& f, const Vector<TaylorModel<ValidatedTag,F,I>>& g) { return compose(f,g); }
    template<class A> ArithmeticType<A,FloatValue<PR>> operator() (Vector<A> const&) const;
    //@}

    //@{
    /*! \name Inplace modifications. */
    //! \brief Scales the model by a function mapping \a dom into the unit interval.
    Void unscale(IntervalDomainType const& dom);
    //! \brief Compute the antiderivative (in place).
    Void antidifferentiate(ArgumentIndexType k);
    //! \brief Compute the weak derivative (in place).
    Void differentiate(ArgumentIndexType k);

    friend TaylorModel<ValidatedTag,F,I> unscale(TaylorModel<ValidatedTag,F,I> tm, IntervalDomainType const& dom) {
        tm.unscale(dom); return std::move(tm); }
    friend TaylorModel<ValidatedTag,F,I> antiderivative(TaylorModel<ValidatedTag,F,I> tm, ArgumentIndexType k) {
        tm.antidifferentiate(k); return std::move(tm); }
    friend TaylorModel<ValidatedTag,F,I> derivative(TaylorModel<ValidatedTag,F,I> tm, ArgumentIndexType k) {
        tm.differentiate(k); return std::move(tm); }
    //@}

    //@{
    /*! \name Operations on the domain. */
    //!\brief Split the Taylor model \a tm, subdividing along the independent variable \a k, taking the lower/middle/upper \a part.
    friend TaylorModel<ValidatedTag,F,I> split(const TaylorModel<ValidatedTag,F,I>& tm, ArgumentIndexType k, SplitPart part) {
        return TaylorModel<ValidatedTag,F,I>::_split(tm,k,part); }
    friend Pair<TaylorModel<ValidatedTag,F,I>,TaylorModel<ValidatedTag,F,I>> split(const TaylorModel<ValidatedTag,F,I>& tm, ArgumentSizeType k) {
        return make_pair(split(tm,k,SplitPart::LOWER),split(tm,k,SplitPart::UPPER)); }
    //! \relates TaylorModel<ValidatedTag,F,I> \brief Embed the model in a space of higher dimension
    friend TaylorModel<ValidatedTag,F,MultiIndex> embed(SizeType as1, const TaylorModel<ValidatedTag,F,I>& tm2, SizeType as3) {
        return TaylorModel<ValidatedTag,F,I>::_embed(as1,tm2,as3); }
    //@}

    //@{
    /*! \name ValidatedTag paradigm-based operations. */
    //! \brief Test if one model is disjoint from (is incompatible with) another.
    friend Bool consistent(const TaylorModel<ValidatedTag,F,I>& tm1, const TaylorModel<ValidatedTag,F,I>& tm2) {
        return TaylorModel<ValidatedTag,F,I>::_consistent(tm1,tm2); }
    //! \brief Test if one model is disjoint from (is incompatible with) another.
    friend Bool inconsistent(const TaylorModel<ValidatedTag,F,I>& tm1, const TaylorModel<ValidatedTag,F,I>& tm2) {
        return TaylorModel<ValidatedTag,F,I>::_inconsistent(tm1,tm2); }
    //! \brief Test if one model refines (is a subset of) another.
    friend Bool refines(const TaylorModel<ValidatedTag,F,I>& tm1, const TaylorModel<ValidatedTag,F,I>& tm2) {
        return TaylorModel<ValidatedTag,F,I>::_refines(tm1,tm2); }
    //! \brief An over-approximation to the common refinement of two Taylor models.
    //! Since the intersection of represented sets of functions cannot be represented
    //! exactly in the class of TaylorModels, truncation errors as well as roundoff errors
    //! may be present. In the absence of roundoff errors, the result is a subset of both
    //! arguments, and is guaranteed to contain any function contained in both arguments.
    friend TaylorModel<ValidatedTag,F,I> refinement(const TaylorModel<ValidatedTag,F,I>& tm1, const TaylorModel<ValidatedTag,F,I>& tm2) {
        return TaylorModel<ValidatedTag,F,I>::_refinement(tm1,tm2); }
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \relates TaylorModel<ValidatedTag,F,I> \brief Embed the model in a space of higher dimension, placing the error in the final variable.
    friend TaylorModel<ValidatedTag,F,I> embed_error(const TaylorModel<ValidatedTag,F,I>& tm) {
        return TaylorModel<ValidatedTag,F,I>::_embed_error(tm); }
    //! \relates TaylorModel<ValidatedTag,F,I> \brief Abstract away the given variables.
    //! For example, the model tm(x0,x1,x2,x3) becomes tm'(x0,x1)=tm(x0,[-1,+1],x1,[-1,+1]) on discarding x1 and x3.
    friend TaylorModel<ValidatedTag,F,I>  discard_variables(const TaylorModel<ValidatedTag,F,I>& tm, const Array<ArgumentSizeType>& variables) {
        return TaylorModel<ValidatedTag,F,I>::_discard_variables(tm,variables); }

    //! \brief Remove all terms based on the \a swp conditions.
    TaylorModel<ValidatedTag,F,I>& sweep(const SweeperType& swp);
    TaylorModel<ValidatedTag,F,I>& simplify(const PropertiesType& prp);

    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<ValidatedTag,F,I>& sweep();
    TaylorModel<ValidatedTag,F,I>& simplify();
    //! \brief Sorts keys.
    TaylorModel<ValidatedTag,F,I>& sort();
    //! \brief Remove terms with the same keys. Assumes sorted.
    TaylorModel<ValidatedTag,F,I>& unique();
    //! \brief Sort the terms in index order and combine terms with the same index.
    TaylorModel<ValidatedTag,F,I>& cleanup();

    //! \brief Set the error to zero.
    //! WARNING: This method does not preserve rigour of the model approximation.
    TaylorModel<ValidatedTag,F,I>& clobber();

    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief Specify a policy to use to remove low-impact terms.
    Void set_sweeper(SweeperType swp) { this->_sweeper=swp; }
    Void set_properties(PropertiesType prp) { this->_sweeper=prp; }
    //! \brief A shared pointer to an object using for removing low-impact terms.
    SweeperType sweeper() const { return this->_sweeper; }
    //! \brief The precision of the coefficients.
    PropertiesType properties() const { return this->sweeper(); }
    //! \brief The precision of the coefficients.
    PrecisionType precision() const { return this->sweeper().precision(); }
    //@}

    //@{
    /*! \name Inplace arithmetic operations. */
    //! \brief Add a constant numerical scalar \c r+=c .
    Void iadd(const NumericType& c);
    //! \brief Multiply by a numerical scalar \c r*=c .
    Void imul(const NumericType& c);
    //! \brief Scalar multiply and add \c r+=c*x .
    Void isma(const NumericType& c, const TaylorModel<ValidatedTag,F,I>& x);
    //! \brief Fused multiply and add \c r+=x1*x2 .
    Void ifma(const TaylorModel<ValidatedTag,F,I>& x1, const TaylorModel<ValidatedTag,F,I>& x2);

    //@}

    //@{
    /*! \name Order operators. */
    //! \brief The pointwise maximum.
    template<class FF> friend TaylorModel<ValidatedTag,FF,I> max(const TaylorModel<ValidatedTag,FF,I>& x, const TaylorModel<ValidatedTag,FF,I>& y);
    //! \brief The pointwise minimum.
    template<class FF> friend TaylorModel<ValidatedTag,FF,I> min(const TaylorModel<ValidatedTag,FF,I>& x, const TaylorModel<ValidatedTag,FF,I>& y);
    //! \brief The pointwise absolute value.
    //! \details If the range of \a x definitely does not include 0, returns +x or -x. Otherwise, uses a uniform polynomial approximation to abs.
    template<class FF> friend TaylorModel<ValidatedTag,FF,I> abs(const TaylorModel<ValidatedTag,FF,I>& x);
    //@}
    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<ValidatedTag,F,I>& x) { return x.str(os); }
    //@}

  public:
    OutputStream& str(OutputStream&) const;
    OutputStream& repr(OutputStream&) const;
  public: // FIXME: Should be private
    Void _set_error(const RawFloat<PR>& ne) { ARIADNE_ASSERT(ne>=0); this->_error=ErrorType(ne); }
    Void _append(IndexType const& a, CoefficientType const& v) { this->_expansion.append(a,v); }
    static Bool _consistent(const TaylorModel<ValidatedTag,F,I>& tm1, const TaylorModel<ValidatedTag,F,I>& tm2);
    static Bool _inconsistent(const TaylorModel<ValidatedTag,F,I>& tm1, const TaylorModel<ValidatedTag,F,I>& tm2);
    static Bool _refines(const TaylorModel<ValidatedTag,F,I>& tm1, const TaylorModel<ValidatedTag,F,I>& tm2);
    static TaylorModel<ValidatedTag,F,I> _refinement(const TaylorModel<ValidatedTag,F,I>& tm1, const TaylorModel<ValidatedTag,F,I>& tm2);
    static TaylorModel<ValidatedTag,F,I> _antiderivative(const TaylorModel<ValidatedTag,F,I>& tm, ArgumentSizeType k);
    static TaylorModel<ValidatedTag,F,I> _weak_derivative(const TaylorModel<ValidatedTag,F,I>& tm, ArgumentSizeType k);
    static TaylorModel<ValidatedTag,F,I> _embed_error(const TaylorModel<ValidatedTag,F,I>& tm);
    static TaylorModel<ValidatedTag,F,I>  _discard_variables(const TaylorModel<ValidatedTag,F,I>&, const Array<ArgumentSizeType>& variables);
    static TaylorModel<ValidatedTag,F,I> _split(const TaylorModel<ValidatedTag,F,I>& tm, ArgumentIndexType k, SplitPart part);
    static TaylorModel<ValidatedTag,F,MultiIndex> _embed(SizeType as1, const TaylorModel<ValidatedTag,F,I>& tm2, SizeType as3);
    static TaylorModel<ValidatedTag,F,I> _compose(TaylorModel<ValidatedTag,F,MultiIndex> const& tf, const Vector<TaylorModel<ValidatedTag,F,I>>& tg);
    static TaylorModel<ValidatedTag,F,I> _compose(Unscaling const& uf, TaylorModel<ValidatedTag,F,I> const& tg);
    static TaylorModel<ValidatedTag,F,I> _compose(TaylorModel<ValidatedTag,F,MultiIndex> const& tf, VectorUnscaling const& u, const Vector<TaylorModel<ValidatedTag,F,I>>& tg);
    static TaylorModel<ValidatedTag,F,I> _partial_evaluate(const TaylorModel<ValidatedTag,F,I>& x, ArgumentIndexType k, NumericType c);
    static Coargument<FloatBounds<PR>> _gradient(const TaylorModel<ValidatedTag,F,I>& x, Argument<FloatBounds<PR>> const& v);
    static FloatBounds<PR> _evaluate(const TaylorModel<ValidatedTag,F,I>& x, Argument<FloatBounds<PR>> const& v);
    static FloatApproximation<PR> _evaluate(const TaylorModel<ValidatedTag,F,I>& x, Argument<FloatApproximation<PR>> const& v);
};

// FIXME: Needed to dispatch gradient of ScaledFunctionPatch
Covector<FloatDPBounds> gradient(const MultivariateTaylorModel<ValidatedTag,FloatDP>& x, const Vector<FloatDPBounds>& v);
Covector<FloatMPBounds> gradient(const MultivariateTaylorModel<ValidatedTag,FloatMP>& x, const Vector<FloatMPBounds>& v);

// FIXME: Needed to dispatch gradient of ScaledFunctionPatch
template<class F, class I> template<class A> auto TaylorModel<ValidatedTag,F,I>::operator() (Vector<A> const& x) const -> ArithmeticType<A,FloatValue<PR>> {
    return horner_evaluate(this->expansion(),x)+FloatBounds<typename F::PrecisionType>(-this->error(),+this->error());
}


/*! \brief A class representing a power series expansion, scaled to the unit box, with an error term.
 *
 * See also Expansion, Polynomial, TaylorModel<ValidatedTag,F,I><IntervalDomainType>.
 */
template<class F, class I>
class TaylorModel<ApproximateTag,F,I>
    : public DispatchTranscendentalAlgebraOperations<TaylorModel<ApproximateTag,F,I>,FloatApproximation<PrecisionType<F>>>
    , public DispatchOrderedAlgebraOperations<TaylorModel<ApproximateTag,F,I>,FloatApproximation<PrecisionType<F>>>
    , public DispatchConcreteGenericAlgebraNumberOperations<TaylorModel<ApproximateTag,F,I>,FloatApproximation<PrecisionType<F>>,ApproximateNumber>
{
    typedef typename F::PrecisionType PR;
  public:
    typedef PR PrecisionType;
    typedef I IndexType;

    typedef FloatApproximation<PR> CoefficientType;
    typedef FloatApproximation<PR> ErrorType;
    typedef ReverseLexicographicIndexLess ComparisonType;
    typedef SortedExpansion<IndexType,CoefficientType,ComparisonType> ExpansionType;

    typedef typename ExpansionType::ArgumentSizeType ArgumentSizeType;
    typedef typename ExpansionType::ArgumentIndexType ArgumentIndexType;

    typedef typename IndexTraits<I>::UnitDomainType DomainType;
    typedef IntervalDomainType CodomainType;
    typedef Interval<FloatApproximation<PR>> RangeType;
    typedef FloatApproximation<PR> NormType;

    //! \brief The type used algebraic operations.
    typedef FloatApproximation<PR> NumericType;
    typedef ApproximateNumber GenericNumericType;
    typedef Sweeper<F> PropertiesType;
    //! \brief The type used for the coefficients.
    typedef CoefficientType ValueType;

    //! \brief The type used to determine the accuracy.
    typedef Sweeper<F> SweeperType;

    template<class X> using Argument = typename IndexTraits<I>::template Argument<X>;
    template<class X> using Coargument = typename IndexTraits<I>::template Coargument<X>;

    //! \brief An Iterator through the (index,coefficient) pairs of the expansion.
    typedef typename ExpansionType::Iterator Iterator;
    //! \brief A constant Iterator through the (index,coefficient) pairs of the expansion.
    typedef typename ExpansionType::ConstIterator ConstIterator;
  private:
    ExpansionType _expansion;
    mutable SweeperType _sweeper;
  private:
    static const CoefficientType _zero;

  public:
    //@{
    /*! \name Constructors and destructors. */
    explicit TaylorModel<ApproximateTag,F,I>();
    //! \brief Construct a TaylorModel<ApproximateTag,F,I> in \a as arguments with a default SweeperType.
    explicit TaylorModel<ApproximateTag,F,I>(ArgumentSizeType as);
    //! \brief Construct a TaylorModel<ApproximateTag,F,I> in \a as arguments with sweeper \a swp.
    TaylorModel<ApproximateTag,F,I>(ArgumentSizeType as, SweeperType swp);

    TaylorModel<ApproximateTag,F,I> create() const { return TaylorModel<ApproximateTag,F,I>(this->argument_size(),this->_sweeper); }
    TaylorModel<ApproximateTag,F,I> create_zero() const { return TaylorModel<ApproximateTag,F,I>(this->argument_size(),this->_sweeper); }
    TaylorModel<ApproximateTag,F,I> create_constant(NumericType) const;
    TaylorModel<ApproximateTag,F,I> create_constant(GenericNumericType) const;
    TaylorModel<ApproximateTag,F,I> create_variable(ArgumentSizeType i) const;
    TaylorModel<ApproximateTag,F,I> create_ball(ErrorType r) const;


    //! \brief Construct a constant quantity in \a as independent variables.
    static TaylorModel<ApproximateTag,F,I> constant(ArgumentSizeType as, const NumericType& c, SweeperType swp) {
        TaylorModel<ApproximateTag,F,I> r(as,swp); r[MultiIndex::zero(as)]=1; r*=c; return r; }
    //! \brief Construct the quantity with expansion \f$x_j\f$ in \a as independent variables.
    static TaylorModel<ApproximateTag,F,I> coordinate(ArgumentSizeType as, ArgumentSizeType j, SweeperType swp) {
        TaylorModel<ApproximateTag,F,I> r(as,swp); r[MultiIndex::unit(as,j)]=1; return r; }

    //! \brief Fast swap with another Taylor model.
    Void swap(TaylorModel<ApproximateTag,F,I>& other) { this->_expansion.swap(other._expansion); }
    //! \brief Set to zero.
    Void clear() { this->_expansion.clear(); }
    //! \brief A zero element with the same parameters.
    TaylorModel<ApproximateTag,F,I> null() const { return TaylorModel<ApproximateTag,F,I>(this->argument_size()); }
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to a constant, keeping the same number of arguments.
    TaylorModel<ApproximateTag,F,I>& operator=(const NumericType& c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),c); return *this; }
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    TaylorModel<ApproximateTag,F,I>& operator=(const GenericNumericType& c) { return (*this)=NumericType(c,this->precision()); }
    //@}

    //@{
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    friend Bool same(const TaylorModel<ApproximateTag,F,I>& tm1, const TaylorModel<ApproximateTag,F,I>& tm2) {
        return same(tm1._expansion, tm2._expansion); }
    //@}

    //@{
    /*! \name Data access */
    //! \brief The expansion.
    const ExpansionType& expansion() const { return this->_expansion; }
    //! \brief A reference to the expansion.
    ExpansionType& expansion() { return this->_expansion; }
    //! \brief The number of variables in the argument of the quantity.
    ArgumentSizeType argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The coefficient of the term in $x^a$.
    const CoefficientType& operator[](const MultiIndex& a) const { return this->_expansion[a]; }
    CoefficientType& operator[](const MultiIndex& a) { return this->_expansion.at(a); }
    //! \brief The constant term in the expansion.
    const CoefficientType& value() const { return (*this)[MultiIndex::zero(this->argument_size())]; }
    //! \brief The coefficient of the gradient term \f$df/dx_j\f$.
    const CoefficientType& gradient_value(ArgumentSizeType j) const { return (*this)[MultiIndex::unit(this->argument_size(),j)]; }
    Coargument<CoefficientType> gradient_value() const { Coargument<CoefficientType> r(this->argument_size());
        for(ArgumentIndexType j=0; j!=this->argument_size(); ++j) { r[j]=(*this)[MultiIndex::unit(this->argument_size(),j)]; } return r; }
    //! \brief The error of the expansion over the domain.
    Void error() const { }
    //@}

    //@{
    /*! \name Function evaluation. */
    //! \brief The domain of the quantity.
    UnitBox domain() const;
    //! \brief A coarse over-approximation to the range of the quantity.
    IntervalDomainType codomain() const;
    //! \brief An over-approximation to the range of the quantity.
    RangeType range() const;
    //! \brief Compute the gradient of the expansion with respect to the \a jth variable over the domain.
    RangeType gradient_range(ArgumentSizeType j) const;
    //@}

    //@{
    /*! \name Inplace modifications. */
    //! \brief Scale so that the old codomain maps into the unit interval.
    Void unscale(const IntervalDomainType& dom);
    //! \brief Compute the antiderivative (in place).
    Void antidifferentiate(ArgumentSizeType k);
    //! \brief Compute the derivative (in place).
    Void differentiate(ArgumentSizeType k);

    friend TaylorModel<ApproximateTag,F,I> unscale(TaylorModel<ApproximateTag,F,I> tm, const IntervalDomainType& dom) {
        tm.unscale(dom); return std::move(tm); }

    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \brief Specify a policy to use to remove low-impact terms.
    Void set_sweeper(SweeperType swp) { this->_sweeper=swp; }
    Void set_properties(PropertiesType prp) { this->set_sweeper(prp); }
    //! \brief A shared pointer to an object using for removing low-impact terms.
    SweeperType sweeper() const { return this->_sweeper; }
    PropertiesType properties() const { return this->_sweeper; }
    //! \brief The precision of the coefficients.
    PrecisionType precision() const { return this->sweeper().precision(); }
    //@}

    //@{
    /*! \name Simplification operations. */
    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    TaylorModel<ApproximateTag,F,I>& sweep();
    TaylorModel<ApproximateTag,F,I>& simplify();
    //! \brief Combine terms with the same index; requires sorted.
    TaylorModel<ApproximateTag,F,I>& unique();
    //! \brief Sort the terms in index order.
    TaylorModel<ApproximateTag,F,I>& sort();
    //@}

  public:
    //@{
    /*! \name Standard algebra interface. */
    //! \brief An approximation to the norm of the function.
    NormType norm() const;
    //! \brief An approximation to the average value of the function.
    CoefficientType average() const;
    //! \brief The tolerance to which analytic functions should be computed.
    F tolerance() const;
    //! \brief The radius of the ball containing the functions.
    NormType radius() const;
    //! \brief Write to an output stream.
    OutputStream& write(OutputStream&) const;
    //! \brief Inplace addition of a scalar constant.
    Void iadd(const FloatApproximation<PR>& c);
    //! \brief Inplace multiplication of a scalar constant.
    Void imul(const FloatApproximation<PR>& c);
    //! \brief Inplace addition of a scalar multiple of a Taylor model.
    Void isma(const FloatApproximation<PR>& c, const TaylorModel<ApproximateTag,F,I>& x);
    //! \brief Inplace addition of a product of Taylor models.
    Void ifma(const TaylorModel<ApproximateTag,F,I>& x1, const TaylorModel<ApproximateTag,F,I>& x2);

/*
    friend TaylorModel<ApproximateTag,F,I> operator-(TaylorModel<ApproximateTag,F,I> x) { x.imul(-1); return std::move(x); }
    friend TaylorModel<ApproximateTag,F,I>& operator+=(TaylorModel<ApproximateTag,F,I>& x, NumericType const& c) { x.iadd(c); return x; }
    friend TaylorModel<ApproximateTag,F,I>& operator*=(TaylorModel<ApproximateTag,F,I>& x, NumericType const& c) { x.imul(c); return x; }
    friend TaylorModel<ApproximateTag,F,I> operator+(TaylorModel<ApproximateTag,F,I> const& x1, TaylorModel<ApproximateTag,F,I> const& x2) {
        TaylorModel<ApproximateTag,F,I> r=x1; r.isma(+1,x2); return std::move(r); }
    friend TaylorModel<ApproximateTag,F,I> operator-(TaylorModel<ApproximateTag,F,I> const& x1, TaylorModel<ApproximateTag,F,I> const& x2) {
        TaylorModel<ApproximateTag,F,I> r=x1; r.isma(-1,x2); return std::move(r); }
    friend TaylorModel<ApproximateTag,F,I> operator*(TaylorModel<ApproximateTag,F,I> const& x1, TaylorModel<ApproximateTag,F,I> const& x2) {
        TaylorModel<ApproximateTag,F,I> r(x1.argument_size(),x1.sweeper()); r.ifma(x1,x2); return std::move(r); }

    template<class FF> friend TaylorModel<ApproximateTag,FF,I> max(const TaylorModel<ApproximateTag,FF,I>& x, const TaylorModel<ApproximateTag,FF,I>& y);
    template<class FF> friend TaylorModel<ApproximateTag,FF,I> min(const TaylorModel<ApproximateTag,FF,I>& x, const TaylorModel<ApproximateTag,FF,I>& y);
    template<class FF> friend TaylorModel<ApproximateTag,FF,I> abs(const TaylorModel<ApproximateTag,FF,I>& x);
*/
    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const TaylorModel<ApproximateTag,F,I>& x) {
        return x.str(os); }
    //@}

  public:
    OutputStream& str(OutputStream&) const;
    OutputStream& repr(OutputStream&) const;
};



template<class F, class I> Bool inconsistent(const Vector<TaylorModel<ValidatedTag,F,I>>& tm1, const Vector<TaylorModel<ValidatedTag,F,I>>& tm2) {
    ARIADNE_ASSERT(tm1.size()==tm2.size());
    for(SizeType i=0; i!=tm1.size(); ++i) { if(inconsistent(tm1[i],tm2[i])) { return true; } }
    return false;
}

template<class F, class I> Bool refines(const Vector<TaylorModel<ValidatedTag,F,I>>& tm1, const Vector<TaylorModel<ValidatedTag,F,I>>& tm2) {
    ARIADNE_ASSERT(tm1.size()==tm2.size());
    for(SizeType i=0; i!=tm1.size(); ++i) { if(not refines(tm1[i],tm2[i])) { return false; } }
    return true;
}

template<class F, class I> Vector<TaylorModel<ValidatedTag,F,I>> refinement(const Vector<TaylorModel<ValidatedTag,F,I>>& tm1, const Vector<TaylorModel<ValidatedTag,F,I>>& tm2) {
    ARIADNE_ASSERT(tm1.size()==tm2.size());
    Vector<TaylorModel<ValidatedTag,F,I>> r(tm1.size());
    for(SizeType i=0; i!=tm1.size(); ++i) { r[i]=refinement(tm1[i],tm2[i]); }
    return std::move(r);
}





template<class F> Vector<FloatBounds<PrecisionType<F>>> evaluate(const Vector<TaylorModel<ValidatedTag,F>>& tf, const Vector<FloatBounds<PrecisionType<F>>>& x) {
    Vector<FloatBounds<PrecisionType<F>>> r(tf.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=evaluate(tf[i],x); }
    return r;
}

template<class F> Vector<TaylorModel<ValidatedTag,F>> partial_evaluate(const Vector<TaylorModel<ValidatedTag,F>>& tf, SizeType k, const FloatBounds<PrecisionType<F>>& c) {
    Vector<TaylorModel<ValidatedTag,F>> r(tf.size(),ValidatedTaylorModel<F>::zero(tf.zero_element().argument_size()-1u,tf.zero_element().sweeper()));
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=partial_evaluate(tf[i],k,c); }
    return std::move(r);
}

template<class F> Vector<TaylorModel<ValidatedTag,F>> compose(const Vector<TaylorModel<ValidatedTag,F>>& tf, const Vector<TaylorModel<ValidatedTag,F>>& tg) {
    Vector<TaylorModel<ValidatedTag,F>> r(tf.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=compose(tf[i],tg); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> unscale(const Vector<TaylorModel<P,F>>& x, const Vector<IntervalDomainType>& dom) {
    Vector<TaylorModel<P,F>> r(x.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=unscale(x[i],dom[i]); }
    return std::move(r);
}

template<class P, class F, class I> Vector<TaylorModel<P,F,I>> compose(const VectorUnscaling& u, const Vector<TaylorModel<P,F,I>>& g) {
    ARIADNE_ASSERT_MSG(u.size()==g.size(),"u="<<u.domain()<<", g="<<g);
    Vector<TaylorModel<P,F,I>> r(u.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=compose(u[i],g[i]); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> compose(const Vector<TaylorModel<P,F>>& x,
                                                            const VectorUnscaling& u,
                                                            const Vector<TaylorModel<P,F>>& y) {
    return compose(x,compose(u,y));
}

template<class P, class F> Vector<TaylorModel<P,F>> antiderivative(const Vector<TaylorModel<P,F>>& x, SizeType k) {
    Vector<TaylorModel<P,F>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=antiderivative(x[i],k); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> derivative(const Vector<TaylorModel<P,F>>& x, SizeType k) {
    Vector<TaylorModel<P,F>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=derivative(x[i],k); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> combine(const Vector<TaylorModel<P,F>>& x1, const Vector<TaylorModel<P,F>>& x2) {
    return join(embed(0u,x1,x2.zero_element().argument_size()),embed(x1.zero_element().argument_size(),x2,0u));
}

template<class F> Vector<TaylorModel<ValidatedTag,F>> embed(SizeType as1, const Vector<TaylorModel<ValidatedTag,F>>& x2, SizeType as3) {
    Vector<TaylorModel<ValidatedTag,F>> r(x2.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=embed(as1,x2[i],as3); }
    return std::move(r);
}

template<class P, class F> Vector<TaylorModel<P,F>> split(const Vector<TaylorModel<P,F>>& x, SizeType j, SplitPart h) {
    Vector<TaylorModel<P,F>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=split(x[i],j,h); }
    return std::move(r);
}

template<class P, class F> Vector<typename TaylorModel<P,F>::RangeType> ranges(const Vector<TaylorModel<P,F>>& f) {
    Vector<typename TaylorModel<P,F>::RangeType> r(f.size()); for(SizeType i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return std::move(r);
}

template<class P, class F> Vector<typename TaylorModel<P,F>::ErrorType> errors(const Vector<TaylorModel<P,F>>& h) {
    Vector<typename TaylorModel<P,F>::ErrorType> e(h.size()); for(SizeType i=0; i!=h.size(); ++i) { e[i]=h[i].error(); } return e; }

template<class P, class F> Vector<typename TaylorModel<P,F>::NormType> norms(const Vector<TaylorModel<P,F>>& h) {
    Vector<NormType> r(h.size()); for(SizeType i=0; i!=h.size(); ++i) { r[i]=norm(h[i]); } return std::move(r); }

template<class P, class F> typename TaylorModel<P,F>::NormType norm(const Vector<TaylorModel<P,F>>& h) {
    typename TaylorModel<P,F>::NormType r=0u; for(SizeType i=0; i!=h.size(); ++i) { r=max(r,norm(h[i])); } return std::move(r);
}

template<class F> Matrix<ValidatedNumericType> jacobian(const Vector<TaylorModel<ValidatedTag,F>>& x, const Vector<ValidatedNumericType>& y);
template<class F> Matrix<ValidatedNumericType> jacobian(const Vector<TaylorModel<ValidatedTag,F>>& x, const Vector<ValidatedNumericType>& y, Array<SizeType>& p);
template<class F> Matrix<FloatDPValue> jacobian_value(const Vector<TaylorModel<ValidatedTag,F>>& x);
template<class F> Matrix<FloatDPValue> jacobian_value(const Vector<TaylorModel<ValidatedTag,F>>& x, const Array<SizeType>& p);
template<class F> Matrix<UpperIntervalType> jacobian_range(const Vector<TaylorModel<ValidatedTag,F>>& x);
template<class F> Matrix<UpperIntervalType> jacobian_range(const Vector<TaylorModel<ValidatedTag,F>>& x, const Array<SizeType>& p);




} // namespace Ariadne

#endif // ARIADNE_TAYLOR_MODEL_HPP
