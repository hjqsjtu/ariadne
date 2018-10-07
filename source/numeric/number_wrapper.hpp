/***************************************************************************
 *            numeric/number_wrapper.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/number_wrapper.hpp
 *  \brief
 */



#ifndef ARIADNE_NUMBER_WRAPPER_HPP
#define ARIADNE_NUMBER_WRAPPER_HPP

#include "../utility/module.hpp"
#include "../numeric/paradigm.hpp"

#include "number_interface.hpp"

#include "number.hpp"
#include "logical.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"
#include "float-user.hpp"

#include "../symbolic/templates.hpp"
#include "../numeric/operators.hpp"

namespace Ariadne {

/************ Number *********************************************************/

template<class... AWS> struct Aware;

template<class X, class I> class Mixin;
template<class X, class I> class Wrapper;

template<class X> using NumberMixin = Mixin<X,NumberInterface>;
template<class X> using NumberWrapper = Wrapper<X,NumberInterface>;

template<class X> using UpperNumberMixin = Mixin<X,UpperNumberInterface>;
template<class X> using UpperNumberWrapper = Wrapper<X,UpperNumberInterface>;


template<class X> inline X const* extract(NumberInterface const* y) {
     return dynamic_cast<NumberWrapper<X>const*>(y);
}

inline OutputStream& operator<<(OutputStream& os, NumberInterface const& y) { return y._write(os); }

inline OutputStream& operator<<(OutputStream& os, Comparison c) {
    return os << ( (c==Comparison::EQUAL) ? "EQUAL" : (c==Comparison::LESS) ? "LESS" : "GREATER" );
}
inline OutputStream& operator<<(OutputStream& os, Sign s) {
    return os << ( (s==Sign::ZERO) ? "ZERO" : (s==Sign::NEGATIVE) ? "NEGATIVE" : "POSITIVE" );
}

// FIXME: Should test for other potential infinities
inline Comparison cmp(NumberInterface const& y1, NumberInterface const& y2) {
    Comparison res;
    FloatDPValue const* x1=extract<FloatDPValue>(&y1);
    FloatDPValue const* x2=extract<FloatDPValue>(&y2);
    if(x1) {
        if(x2) { res= cmp(ExactDouble(x1->raw().get_d()),ExactDouble(x2->raw().get_d())); }
        else { res= cmp(ExactDouble(x1->raw().get_d()),y2._get_q()); }
    } else {
        if(x2) { res= cmp(y1._get_q(),ExactDouble(x2->raw().get_d())); }
        else { res= cmp(y1._get_q(),y2._get_q()); }
    }
    return res;
}



//------------ Dispatching code -----------------------

template<class X, class I, class OP, class J=I> struct UnaryOperationMixin;
template<class X, class I, class OP, class J=I> struct BinaryOperationMixin;
template<class X, class I, class OP, class Y> struct BinaryOperable;

template<class X, class OP> inline X const& _upcast(UnaryOperationMixin<X,NumberInterface,OP> const& y) { 
    return static_cast<Mixin<X,NumberInterface> const&>(y); }
template<class X, class OP> inline X const& _upcast(BinaryOperationMixin<X,NumberInterface,OP> const& y) { 
    return static_cast<Mixin<X,NumberInterface> const&>(y); }
template<class X, class OP, class Y> inline X const& _upcast(BinaryOperable<X,NumberInterface,OP,Y> const& y) { 
    return static_cast<Mixin<X,NumberInterface> const&>(y); }
    
template<class R> inline NumberInterface* _make_wrapper(R&& r, NumberInterface*) { return new NumberWrapper<R>(r); }


template<class X, class OP> inline X const& _upcast(UnaryOperationMixin<X,UpperNumberInterface,OP> const& y) { 
    return static_cast<UpperNumberMixin<X> const&>(y); }
template<class X, class OP> inline X const& _upcast(BinaryOperationMixin<X,UpperNumberInterface,OP> const& y) { 
    return static_cast<UpperNumberMixin<X> const&>(y); }
template<class X, class OP, class Y> inline X const& _upcast(BinaryOperable<X,UpperNumberInterface,OP,Y> const& y) { 
    return static_cast<UpperNumberMixin<X> const&>(y); }
template<class R> inline UpperNumberInterface* _make_wrapper(R&& r, UpperNumberInterface*) { return new UpperNumberWrapper<R>(r); }


template<class I, class R> inline I* _make_wrapper(R&& r) { I* p=nullptr; return _make_wrapper(std::forward<R>(r),p); }

//------------ Single-dispatching code -----------------------

template<class X, class I, class OP, class J> struct UnaryOperationMixin : public virtual J {
    virtual I* _apply(OP op) const final { return _make_wrapper<I>(op(_upcast<X>(*this))); }
};

//------------ Double-dispatching code -----------------------

template<class I, class OP, class Y> struct BinaryOperableInterface {
    virtual ~BinaryOperableInterface() = default;
    virtual I* _apply_left(OP op, Y const& y) const = 0;
    virtual I* _apply_right(OP op, Y const& y) const = 0;
};
template<class X, class I, class OP, class Y> struct BinaryOperable : virtual BinaryOperableInterface<I,OP,Y> {
    virtual I* _apply_left(OP op, Y const& other) const { return _make_wrapper<I>(op(_upcast<X>(*this),other)); }
    virtual I* _apply_right(OP op, Y const& other) const { return _make_wrapper<I>(op(other,_upcast<X>(*this))); }
};
template<class X, class I, class OP, class AW> struct BinaryOperableMixin;
template<class X, class I, class OP, class Y, class... YS> struct BinaryOperableMixin<X,I,OP,Aware<Y,YS...>>
    : BinaryOperable<X,I,OP,Y>, BinaryOperableMixin<X,I,OP,Aware<YS...>> { };
template<class X, class I, class OP> struct BinaryOperableMixin<X,I,OP,Aware<>> { };


template<class OP> inline NumberInterface* make_symbolic(OP op, NumberInterface const* yp1, NumberInterface const* yp2) {
    Handle<NumberInterface> y1(const_cast<NumberInterface*>(yp1)->shared_from_this());
    Handle<NumberInterface> y2(const_cast<NumberInterface*>(yp2)->shared_from_this());
    String yc1=yp1->_class_name(); String yc2=yp2->_class_name();    
    ARIADNE_THROW(DispatchException,op<<"(Number y1, Number y2) with y1="<<*yp1<<", y2="<<*yp2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
}

template<class OP> inline UpperNumberInterface* make_symbolic(OP op, UpperNumberInterface const* yp1, UpperNumberInterface const* yp2) {
    Handle<UpperNumberInterface> y1(const_cast<UpperNumberInterface*>(yp1)->shared_from_this());
    Handle<UpperNumberInterface> y2(const_cast<UpperNumberInterface*>(yp2)->shared_from_this());
    String yc1=yp1->_class_name(); String yc2=yp2->_class_name();    
    ARIADNE_THROW(DispatchException,op<<"(UpperNumber y1, UpperNumber y2) with y1="<<*yp1<<", y2="<<*yp2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
}



/*
template<class I, class OP> struct BinaryOperationInterface {
    virtual ~BinaryOperationInterface() = default;
    virtual I* _apply(OP op, I const* other) const = 0;
    virtual I* _rapply(OP op, I const* other) const = 0;
};
*/

template<class X, class I, class OP, class J> struct BinaryOperationMixin : public virtual J {
    virtual I* _apply(OP op, I const* other) const final;
    virtual I* _rapply(OP op, I const* other) const final;
};


template<class X, class I, class OP, class J> inline I* BinaryOperationMixin<X,I,OP,J>::_apply(OP op, I const* other) const {
    auto aware_other=dynamic_cast<BinaryOperableInterface<I,OP,X>const*>(other);
    if(aware_other) { X const& self = _upcast<X>(*this); return aware_other->_apply_right(op,self); }
    else { return other->_rapply(op,this); }
}
template<class X, class I, class OP, class J> inline I* BinaryOperationMixin<X,I,OP,J>::_rapply(OP op, I const* other) const {
    auto aware_other=dynamic_cast<BinaryOperableInterface<I,OP,X>const*>(other);
    if(aware_other) { X const& self = _upcast<X>(*this); return aware_other->_apply_left(op,self); }
    else { return make_symbolic(op,other,this); }
}

//------------ End dispatching code -----------------------

/*
template<class X, class I, class W> struct SameArithmeticMixin : public virtual I {
    X const& _cast(X const& self) { return self; }
    X const& _cast(I const& other) { return dynamic_cast<Wrapper<X,I>const&>(other); }
    I* _heap_move(X&& x) { return new W(x); }
    virtual I* _add(I const* other) const final { return _heap_move(add(_cast(*this),_cast(*other))); }
    virtual I* _sub(I const* other) const final { return _heap_move(sub(_cast(*this),_cast(*other))); }
    virtual I* _mul(I const* other) const final { return _heap_move(mul(_cast(*this),_cast(*other))); }
    virtual I* _div(I const* other) const final { return _heap_move(div(_cast(*this),_cast(*other))); }
    virtual I* _radd(I const* other) const final { return _heap_move(add(_cast(*other),_cast(*this))); }
    virtual I* _rsub(I const* other) const final { return _heap_move(sub(_cast(*other),_cast(*this))); }
    virtual I* _rmul(I const* other) const final { return _heap_move(mul(_cast(*other),_cast(*this))); }
    virtual I* _rdiv(I const* other) const final { return _heap_move(div(_cast(*other),_cast(*this))); }
};
*/

template<class X> class NumberGetterMixin : public virtual NumberInterface {
  public:
  //  operator X const& () const { return static_cast<Mixin<X,NumberInterface>const&>(*this); }
    static X const& _cast(NumberGetterMixin<X> const& self) { return static_cast<Mixin<X,NumberInterface>const&>(self); }
    static X& _cast(NumberGetterMixin<X>& self) { return static_cast<Mixin<X,NumberInterface>&>(self); }

    typedef Paradigm<X> P;
    friend class Number<P>;

    virtual NumberInterface* _copy() const override { return new NumberWrapper<X>(_cast(*this)); }
    virtual NumberInterface* _move() override { return new NumberWrapper<X>(std::move(_cast(*this))); }

    // FIXME: Proper comparisons for ExactNumber.
    virtual LogicalValue _apply(Equal, NumberInterface const* y) const override {
        if (this->_paradigm() == ParadigmCode::EXACT && y->_paradigm() == ParadigmCode::EXACT) {
            return LogicalValue( cmp(*this,*y)==Comparison::EQUAL ? LogicalValue::TRUE : LogicalValue::FALSE ); }
        if (this->_paradigm() == ParadigmCode::VALIDATED && y->_paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(OrderTag(),dp) == y->_get(OrderTag(),dp)); }
        else {
            return LogicalValue(this->_get(ApproximateTag(),dp) == y->_get(ApproximateTag(),dp)); }
    }
    virtual LogicalValue _apply(Less, NumberInterface const* y) const override {
        if (this->_paradigm() == ParadigmCode::EXACT && y->_paradigm() == ParadigmCode::EXACT) {
            return LogicalValue( cmp(*this,*y)==Comparison::LESS ? LogicalValue::TRUE : LogicalValue::FALSE ); }
        else if (this->_paradigm() == ParadigmCode::VALIDATED && y->_paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(OrderTag(),dp) < y->_get(OrderTag(),dp)); }
        else {
            return LogicalValue(this->_get(ApproximateTag(),dp) < y->_get(ApproximateTag(),dp));
        }
    }

    virtual Rational _get_q() const override {
        return this->_get_as<Rational>(); }

    virtual FloatDPBall _get(MetricTag,DoublePrecision pr,DoublePrecision pre) const override {
        return this->_get_as<FloatDPBall>(pr,pre); }
    virtual FloatDPBounds _get(OrderTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPBounds>(pr); }
    virtual FloatDPUpperBound _get(UpperTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPUpperBound>(pr); }
    virtual FloatDPLowerBound _get(LowerTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPLowerBound>(pr); }
    virtual FloatDPApproximation _get(ApproximateTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPApproximation>(pr); }
    virtual FloatMPDPBall _get(MetricTag,MultiplePrecision pr, DoublePrecision pre) const override {
        return this->_get_as<FloatMPDPBall>(pr,pre); }
    virtual FloatMPBall _get(MetricTag, MultiplePrecision pr, MultiplePrecision pre) const override {
        return this->_get_as<FloatMPBall>(pr,pre); }
    virtual FloatMPBounds _get(OrderTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPBounds>(pr); }
    virtual FloatMPUpperBound _get(UpperTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPUpperBound>(pr); }
    virtual FloatMPLowerBound _get(LowerTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPLowerBound>(pr); }
    virtual FloatMPApproximation _get(ApproximateTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPApproximation>(pr); }

    virtual ParadigmCode _paradigm() const override { return P::code(); }
    virtual String _class_name() const override { return class_name<X>(); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << _cast(*this); }

  private:
    template<class R, class... PRS, EnableIf<IsConstructible<R,X,PRS...>> = dummy>
        inline R _get_as(PRS... prs) const { return R(_cast(*this),prs...); }
    template<class R, class... PRS, DisableIf<IsConstructible<R,X,PRS...>> = dummy>
        inline R _get_as(PRS... prs) const { 
            std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>();
//            if constexpr(std::tuple_size<Tuple<PRS...>>::value==0) { std::cerr << " with precision " << std::get<0>(prs...); }
            //std::cerr<< " with precision " << pr << " and error precision " << pre << "\n"; 
            std::cerr<< "\n"; throw ParadigmError(); }
};

template<class X> struct DispatchingTraits { typedef Aware<X> AwareOfTypes; };
template<class X> using Awares = typename DispatchingTraits<X>::AwareOfTypes;

template<class X> class Mixin<X,NumberInterface>
    : public UnaryOperationMixin<X,NumberInterface,UnaryOperator>
    , public BinaryOperationMixin<X,NumberInterface,BinaryOperator>
//    , public BinaryOperationMixin<X,LogicalValue,LogicalOperator, NumberInterface>
    , public BinaryOperableMixin<X,NumberInterface,BinaryOperator,Awares<X>>
    , public NumberGetterMixin<X>
{
  public:
    operator X const& () const { return static_cast<Wrapper<X,NumberInterface>const&>(*this); }
    operator X& () { return static_cast<Wrapper<X,NumberInterface>&>(*this); }
};


template<class X> class Wrapper<X,NumberInterface>
    : public X, public Mixin<X,NumberInterface>
{
    static_assert(Not<IsSame<X,Handle<NumberInterface>>>::value,"X must be a concrete number, not a handle");
    static_assert(Not<IsSame<X,Number<Paradigm<X>>>>::value,"X must be a concrete number, not a generic number");
  public:
    Wrapper<X,NumberInterface>(const X& a) : X(a) { }
    Wrapper<X,NumberInterface>(X&& a) : X(std::forward<X>(a)) { }
};


template<class X> class UpperNumberGetterMixin : public virtual UpperNumberInterface {
    typedef Paradigm<X> P;
    static X const& _cast(UpperNumberGetterMixin<X> const& self) { return static_cast<UpperNumberMixin<X>const&>(self); }

  public:
    virtual UpperNumberInterface* _copy() const override { return new UpperNumberWrapper<X>(_cast(*this)); }
    virtual UpperNumberInterface* _move() override { return new UpperNumberWrapper<X>(std::move(_cast(*this))); }

    virtual FloatDPUpperBound _get(DoublePrecision pr) const override { return this->_get_as<FloatDPUpperBound>(pr); }
    virtual FloatMPUpperBound _get(MultiplePrecision pr) const override { return this->_get_as<FloatMPUpperBound>(pr); }

    virtual ParadigmCode _paradigm() const override { assert(false); }
    virtual String _class_name() const override { return class_name<X>(); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << _cast(*this); }

  private:
    template<class R, class... PRS, EnableIf<IsConstructible<R,X,PRS...>> = dummy>
        inline R _get_as(PRS... prs) const { return R(_cast(*this),prs...); }
    template<class R, DisableIf<IsConstructible<R,X>> = dummy>
        inline R _get_as() const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << "\n"; throw ParadigmError(); }
    template<class R, class PR, DisableIf<IsConstructible<R,X,PR>> = dummy>
        inline R _get_as(PR pr) const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << " with precision " << pr << "\n"; throw ParadigmError(); }
};

template<class X> class Mixin<X,UpperNumberInterface>
    : public UnaryOperationMixin<X,UpperNumberInterface,MonotoneUnaryOperator>
    , public BinaryOperationMixin<X,UpperNumberInterface,MonotoneBinaryOperator>
//    , public BinaryOperationMixin<X,LogicalValue,LogicalOperator, UpperNumberInterface>
    , public BinaryOperableMixin<X,UpperNumberInterface,MonotoneBinaryOperator,Awares<X>>
    , public UpperNumberGetterMixin<X>
{
  public:
    operator X const& () const { return static_cast<UpperNumberWrapper<X>const&>(*this); }
    operator X& () { return static_cast<UpperNumberWrapper<X>&>(*this); }
};

template<class X> class Wrapper<X,UpperNumberInterface>
    : public X, public Mixin<X,UpperNumberInterface>
{
    static_assert(Not<IsSame<X,Handle<UpperNumberInterface>>>::value,"X must be a concrete number, not a handle");
    static_assert(Not<IsSame<X,UpperNumber<Paradigm<X>>>>::value,"X must be a concrete number, not a generic number");
  public:
    Wrapper(const X& a) : X(a) { }
    Wrapper(X&& a) : X(std::forward<X>(a)) { }
};


} // namespace Ariadne

#endif /* ARIADNE_NUMBER_WRAPPER_HPP */
