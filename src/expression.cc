/***************************************************************************
 *            expression.cc
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
 
#include "expression.h"
#include "valuation.h"
#include "vector.h"

#include "numeric.h"
#include "taylor_model.h"
#include "differential.h"

#include "polynomial.h"
#include "affine.h"

#include "real.h"
#include "formula.h"
#include "function.h"

namespace Ariadne {


template<class T> inline
std::ostream& operator<<(std::ostream& os, const ExpressionInterface<T>& e) { return e.write(os); }


//! A constant, viewed as a function \f$c:\R^n\rightarrow\R\f$.
template<class R>
class ConstantExpression
    : public ExpressionInterface<R>
{
  public:
    ConstantExpression(const R& c) : _c(c) { }
    operator R () const { return _c; }
    R value() const { return _c; }
    virtual Operator type() const { return CNST; }
    virtual Set<String> arguments() const { return Set<String>(); }
    virtual ConstantExpression<R>* clone() const { return new ConstantExpression<R>(*this); }
    virtual std::ostream& write(std::ostream& os) const { return os << _c; }
  private:
    R _c;
};

//! A constant, viewed as a function \f$c:\R^n\rightarrow\R\f$.
template<>
class ConstantExpression<Real>
    : public ExpressionInterface<Real>
{
  public:
    ConstantExpression(const Float& c) : _c(c) { }
    ConstantExpression(const Interval& c) :  _c(c) { }
    operator Real() const { return _c; }
    Real value() const { return _c; }
    virtual Operator type() const { return CNST; }
    virtual Set<String> arguments() const { return Set<String>(); }
    virtual ConstantExpression<Real>* clone() const { return new ConstantExpression<Real>(*this); }
    virtual std::ostream& write(std::ostream& os) const {
        if(_c.lower()==_c.upper()) { os<<midpoint(_c); } else { os<<_c; } return os; }
  private:
    Real _c;
};


//! A projection onto a named variable.
template<class R>
class VariableExpression
    : public ExpressionInterface<R>
{
  public:
    explicit VariableExpression(const String& s) : _var(s) { }
    VariableExpression(const Variable<R>& v) : _var(v) { }
    const Variable<R>& variable() const { return this->_var; }
    String name() const { return this->_var.name(); }
    virtual Operator type() const { return VAR; }
    virtual VariableExpression<R>* clone() const { return new VariableExpression<R>(*this); }
    virtual Set<String> arguments() const { Set<String> r; r.insert(this->name()); return r; }
    virtual std::ostream& write(std::ostream& os) const { return os << this->name(); }
  private:
    Variable<R> _var;
};

template<class T> class CoordinateExpression;

//! A coordinate projection \f$\pi:\R^n\rightarrow\R\f$ given by \f$\pi(x)=x_j\f$.
template<>
class CoordinateExpression<Real>
    : public ExpressionInterface<Real>
{
    typedef unsigned int SizeType;
  public:
    explicit CoordinateExpression() : _as(0), _j(0) { }
    explicit CoordinateExpression(SizeType j) : _as(0), _j(j) { }
    explicit CoordinateExpression(unsigned int as, unsigned int j) : _as(as), _j(j) { }
    SizeType argument_size() const { return _as; }
    SizeType index() const { return _j; }
    SizeType coordinate() const { return _j; }
    String name() const { String s("x0"); s[1]+=_j; return s; }
    virtual Operator type() const { return IND; }
    virtual CoordinateExpression<Real>* clone() const { return new CoordinateExpression<Real>(*this); }
    virtual Set<String> arguments() const { Set<String> r; r.insert(this->name()); return r; }
    virtual std::ostream& write(std::ostream& os) const { return os << this->name(); }
  private:
    SizeType _as;
    SizeType _j;
};



//! A composed scalar function, based on a standard operator.
template<class R, class Op=Operator, class A=R> class UnaryExpression
    : public ExpressionInterface<R>
{
  public:
    UnaryExpression(const Op& op, const ExpressionInterface<A>& expr)
        : _op(op), _arg(expr.clone()) { }
    UnaryExpression(const Op& op, const ExpressionInterface<A>* expr)
        : _op(op), _arg(const_cast<ExpressionInterface<A>*>(expr)) { }
    UnaryExpression(const Op& op, shared_ptr< const ExpressionInterface<A> > expr)
        : _op(op), _arg(expr) { }
    virtual Operator type() const { return static_cast<Operator>(_op); }
    virtual UnaryExpression<R,Op,A>* clone() const { return new UnaryExpression<R,Op,A>(_op,_arg._ptr); }
    virtual Set<String> arguments() const { return this->_arg.arguments(); }
    virtual std::ostream& write(std::ostream& os) const;
  public:
    Op _op;
    Expression<A> _arg;
};

template<class R, class Op, class A> inline std::ostream& UnaryExpression<R,Op,A>::write(std::ostream& os) const {
    switch(_op) {
        case NEG: return os << '-' << _arg;
        case NOT: return os << '!' << _arg;
        default: return os << _op << "(" << _arg << ")";
    }
}


//! A composed scalar function, based on an arthmetic operator.
template<class R, class Op=Operator, class A1=R, class A2=A1> class BinaryExpression
    : public ExpressionInterface<R>
{
  public:
    BinaryExpression(Op op, const ExpressionInterface<A1>& expr1, const ExpressionInterface<A2>& expr2)
        : _op(op), _arg1(expr1.clone()), _arg2(expr2.clone()) { }
    BinaryExpression(Op op, const ExpressionInterface<A1>* expr1, const ExpressionInterface<A2>* expr2)
        : _op(op), _arg1(expr1), _arg2(expr2) { }
    BinaryExpression(Op op, shared_ptr< const ExpressionInterface<A1> > expr1, shared_ptr< const ExpressionInterface<A2> > expr2)
        : _op(op), _arg1(expr1), _arg2(expr2)  { }
    virtual Operator type() const { return static_cast<Operator>(_op); }
    virtual BinaryExpression<R,Op,A1,A2>* clone() const { return new BinaryExpression<R,Op,A1,A2>(_op,_arg1._ptr,_arg2._ptr); }
    virtual Set<String> arguments() const { return join(this->_arg1.arguments(),this->_arg2.arguments()); }
    virtual std::ostream& write(std::ostream& os) const {
        return os << "(" << _arg1 << symbol(_op) << _arg2 << ")"; }
/*
  private:
    template<class R, class A> inline
    void compute(R& r, const A& a) { r=Op()(_arg1->evaluate(a),_arg2->evaluate(a)); }
*/
  public:
    Op _op;
    Expression<A1> _arg1;
    Expression<A2> _arg2;
};

bool operator==(const Expression<Tribool>& e, bool v) {
    const ConstantExpression<Tribool>* expr = dynamic_cast<const ConstantExpression<Tribool>*>(e.ptr());
    return expr && expr->value()==v;
}


template<class R> Expression<R>::Expression(const R& c) : _ptr(new ConstantExpression<R>(c)) { }
template<class R> Expression<R>::Expression(const Variable<R>& v) : _ptr(new VariableExpression<R>(v)) { }

Expression<Real>::Expression(const double& c) : _ptr(new ConstantExpression<Real>(c)) { }
Expression<Real>::Expression(const Interval& c) : _ptr(new ConstantExpression<Real>(c)) { }
Expression<Real>::Expression(const Variable<Real>& v) : _ptr(new VariableExpression<Real>(v)) { }

template class Expression<Boolean>;
template class Expression<Tribool>;
template class Expression<String>;
template class Expression<Integer>;
template class Expression<Real>;






template<class R, class Op, class A> inline
Expression<R> make_expression(Op op, Expression<A> e) {
    return Expression<R>(new UnaryExpression<R,Op,A>(op,e._ptr)); }
template<class R, class Op, class A1, class A2> inline
Expression<R> make_expression(Op op, Expression<A1> e1, Expression<A2> e2) {
    return Expression<R>(new BinaryExpression<R,Op,A1,A2>(op,e1._ptr,e2._ptr)); }


Expression<Boolean> operator&&(Expression<Boolean> e1, Expression<Boolean> e2) {
    return make_expression<Boolean>(AND,e1,e2); }
Expression<Boolean> operator||(Expression<Boolean> e1, Expression<Boolean> e2) {
    return make_expression<Boolean>(OR,e1,e2); }
Expression<Boolean> operator!(Expression<Boolean> e) {
    return make_expression<Boolean>(NOT,e); }


Expression<tribool> operator&&(Expression<tribool> e1, Expression<tribool> e2) {
    return make_expression<tribool>(AND,e1,e2); }
Expression<tribool> operator||(Expression<tribool> e1, Expression<tribool> e2) {
    return make_expression<tribool>(OR,e1,e2); }
Expression<tribool> operator!(Expression<tribool> e) {
    return make_expression<tribool>(NOT,e); }


Expression<Boolean> operator==(Variable<String> v1, const String& s2) {
    return make_expression<Boolean>(EQ,Expression<String>(v1),Expression<String>(s2)); }


Expression<Boolean> operator==(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(EQ,e1,e2); }
Expression<Boolean> operator!=(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(NEQ,e1,e2); }
Expression<Boolean> operator>=(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(GEQ,e1,e2); }
Expression<Boolean> operator<=(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(LEQ,e1,e2); }
Expression<Boolean> operator>(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(GT,e1,e2); }
Expression<Boolean> operator<(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(LT,e1,e2); }



Expression<Integer> operator+(Expression<Integer> e) {
    return make_expression<Integer>(POS,e); }
Expression<Integer> operator-(Expression<Integer> e) {
    return make_expression<Integer>(NEG,e); }
Expression<Integer> operator+(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Integer>(ADD,e1,e2); }
Expression<Integer> operator-(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Integer>(SUB,e1,e2); }
Expression<Integer> operator*(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Integer>(MUL,e1,e2); }



Expression<Tribool> sgn(Expression<Real> e) {
    return make_expression<Tribool>(SGN,e); }

Expression<Tribool> operator==(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(EQ,e1,e2); }
Expression<Tribool> operator!=(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(NEQ,e1,e2); }
Expression<Tribool> operator>=(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(GEQ,e1,e2); }
Expression<Tribool> operator<=(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(LEQ,e1,e2); }
Expression<Tribool> operator>(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(GT,e1,e2); }
Expression<Tribool> operator<(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(LT,e1,e2); }


Expression<Real> operator+(Expression<Real> e) {
    return make_expression<Real>(POS,e); }
Expression<Real> operator-(Expression<Real> e) {
    return make_expression<Real>(NEG,e); }
Expression<Real> operator+(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(ADD,e1,e2); }
Expression<Real> operator-(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(SUB,e1,e2); }
Expression<Real> operator*(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(MUL,e1,e2); }
Expression<Real> operator/(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(DIV,e1,e2); }

//Expression<Real> pow(Expression<Real> e, int n) { }
//    return make_expression(POW,e,n); }

Expression<Real> neg(Expression<Real> e) {
    return make_expression<Real>(NEG,e); }
Expression<Real> rec(Expression<Real> e) {
    return make_expression<Real>(REC,e); }
Expression<Real> sqr(Expression<Real> e) {
    return make_expression<Real>(SQR,e); }
Expression<Real> sqrt(Expression<Real> e) {
    return make_expression<Real>(SQRT,e); }
Expression<Real> exp(Expression<Real> e) {
    return make_expression<Real>(EXP,e); }
Expression<Real> log(Expression<Real> e) {
    return make_expression<Real>(LOG,e); }
Expression<Real> sin(Expression<Real> e) {
    return make_expression<Real>(SIN,e); }
Expression<Real> cos(Expression<Real> e) {
    return make_expression<Real>(COS,e); }

Expression<Real> tan(Expression<Real> e) {
    return make_expression<Real>(TAN,e); }



inline void _set_constant(Float& r, const Interval& c) { r=midpoint(c); }
inline void _set_constant(Interval& r, const Interval& c) { r=c; }
inline void _set_constant(TaylorModel& r, const Interval& c) { r.clear(); r+=c; }
inline void _set_constant(Differential<Float>& r, const Interval& c) { r.clear(); r+=midpoint(c); }
inline void _set_constant(Differential<Interval>& r, const Interval& c) { r.clear(); r+=c; }

Boolean _compute(Comparison cmp, const String& s1, const String& s2) {
    switch(cmp) {
        case EQ:  return s1==s2;
        case NEQ: return s1!=s2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate comparison "<<cmp<<" on string arguments.");
    }
}

Boolean _compute(Comparison cmp, const Integer& z1, const Integer& z2) {
    switch(cmp) {
        case EQ:  return z1==z2;
        case NEQ: return z1!=z2;
        case LEQ: return z1<=z2;
        case GEQ: return z1>=z2;
        case LT:  return z1< z2;
        case GT:  return z1> z2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate comparison "<<cmp<<" on integer arguments.");
    }
}

template<class X> Tribool _compute(Comparison cmp, const X& x1, const X& x2) {
    switch(cmp) {
        case GT: case GEQ: return x1>x2;
        case LT: case LEQ: return x1<x2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate comparison "<<cmp<<" on real arguments.");
    }
}

Boolean _compute(Operator op, const Boolean& b) {
    switch(op) {
        case NOT: return !b;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

Boolean _compute(Operator op, const Boolean& b1, const Boolean& b2) {
    switch(op) {
        case AND: return b1 && b2;
        case OR: return b1 || b2;
        case XOR: return b1 ^ b2;
        case IMPL: return !b1 || b2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

Tribool _compute(Operator op, const Tribool& b) {
    switch(op) {
        case NOT: return !b;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

Tribool _compute(Operator op, const Tribool& b1, const Tribool& b2) {
    switch(op) {
        case AND: return b1 && b2;
        case OR: return b1 || b2;
        case XOR: return b1 ^ b2;
        case IMPL: return !b1 || b2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on a boolean arguments.");
    }
}

Integer _compute(Operator op, const Integer& x1, const Integer& x2) {
    switch(op) {
        case ADD: return x1+x2;
        case SUB: return x1-x2;
        case MUL: return x1*x2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on two integer arguments.");
    }
}

template<class X> X _compute(Operator op, const X& x1, const X& x2) {
    switch(op) {
        case ADD: return x1+x2;
        case SUB: return x1-x2;
        case MUL: return x1*x2;
        case DIV: return x1/x2;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on two real arguments.");
    }
}

Integer _compute(Operator op, const Integer& z) {
    switch(op) {
        case POS: return +z;
        case NEG: return -z;
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one integer argument.");
    }
}

template<class X>
X _compute(Operator op, const X& x) {
    switch(op) {
        case NEG: return -x;
        case REC: return 1/x;
        case EXP: return exp(x);
        case LOG: return log(x);
        case SIN: return sin(x);
        case COS: return cos(x);
        case TAN: return cos(x);
        default: ARIADNE_FAIL_MSG("Cannot evaluate operator "<<op<<" on one real argument.");
    }
}





Boolean evaluate(const Expression<Boolean>& e, const DiscreteValuation& x) {
    const ExpressionInterface<Boolean>* eptr=e.ptr();
    const BinaryExpression<Boolean>* bptr=dynamic_cast<const BinaryExpression<Boolean,Operator>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const UnaryExpression<Boolean>* uptr=dynamic_cast<const UnaryExpression<Boolean>*>(eptr);
    if(uptr) { return _compute(uptr->_op,evaluate(uptr->_arg,x)); }
    const ConstantExpression<Boolean>* cptr=dynamic_cast<const ConstantExpression<Boolean>*>(eptr);
    if(cptr) { return cptr->value(); }
    const BinaryExpression<Boolean,Comparison,String>* bsptr=dynamic_cast<const BinaryExpression<Boolean,Comparison,String>*>(eptr);
    if(bsptr) { return _compute(bsptr->_op,evaluate(bsptr->_arg1,x),evaluate(bsptr->_arg2,x)); }
    const BinaryExpression<Boolean,Comparison,Integer>* bzptr=dynamic_cast<const BinaryExpression<Boolean,Comparison,Integer>*>(eptr);
    if(bzptr) { return _compute(bsptr->_op,evaluate(bzptr->_arg1,x),evaluate(bzptr->_arg2,x)); }
    ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" to a Boolean using variables "<<x);
}

String evaluate(const Expression<String>& e, const DiscreteValuation& x) {
    const ExpressionInterface<String>* eptr=e.ptr();
    const ConstantExpression<String>* cptr=dynamic_cast<const ConstantExpression<String>*>(eptr);
    if(cptr) { return cptr->value(); }
    const VariableExpression<String>* vptr=dynamic_cast<const VariableExpression<String>*>(eptr);
    if(vptr) { return x[vptr->variable()]; }
    ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" to a String using variables "<<x);
}

Integer evaluate(const Expression<Integer>& e, const DiscreteValuation& x) {
    const ExpressionInterface<Integer>* eptr=e.ptr();
    const BinaryExpression<Integer>* bptr=dynamic_cast<const BinaryExpression<Integer>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const UnaryExpression<Integer>* uptr=dynamic_cast<const UnaryExpression<Integer>*>(eptr);
    if(uptr) { return _compute(uptr->_op,evaluate(uptr->_arg,x)); }
    const ConstantExpression<Integer>* cptr=dynamic_cast<const ConstantExpression<Integer>*>(eptr);
    if(cptr) { return cptr->value(); }
    const VariableExpression<Integer>* vptr=dynamic_cast<const VariableExpression<Integer>*>(eptr);
    if(vptr) { return x[vptr->variable()]; }
    ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" to an Integer using variables "<<x);
}


template<class X> Tribool evaluate(const Expression<Tribool>& e, const ContinuousValuation<X>& x) {
    const ExpressionInterface<Tribool>* eptr=e.ptr();
    const BinaryExpression<Tribool>* bptr=dynamic_cast<const BinaryExpression<Tribool>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const BinaryExpression<Tribool,Comparison,Real,Real>* brptr=dynamic_cast<const BinaryExpression<Tribool,Comparison,Real,Real>*>(eptr);
    if(brptr) { return _compute(bptr->_op,evaluate(brptr->_arg1,x),evaluate(brptr->_arg2,x)); }
    ARIADNE_FAIL_MSG("");
}

template<class X> X evaluate(const Expression<Real>& e, const ContinuousValuation<X>& x) {
    const ExpressionInterface<Real>* eptr=e.ptr();
    const BinaryExpression<Real>* bptr=dynamic_cast<const BinaryExpression<Real>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const UnaryExpression<Real>* uptr=dynamic_cast<const UnaryExpression<Real>*>(eptr);
    if(uptr) { return _compute(uptr->_op,evaluate(uptr->_arg,x)); }
    const ConstantExpression<Real>* cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
    if(cptr) { X r; _set_constant(r,cptr->value()); return r; }
    const VariableExpression<Real>* vptr=dynamic_cast<const VariableExpression<Real>*>(eptr);
    if(vptr) { return x[vptr->variable()]; }
    ARIADNE_FAIL_MSG("");
}

template<class X> Tribool evaluate(const Expression<Tribool>& e, const Vector<X>& x) {
    const ExpressionInterface<Tribool>* eptr=e.ptr();
    const BinaryExpression<Tribool>* bptr=dynamic_cast<const BinaryExpression<Tribool>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const BinaryExpression<Tribool,Comparison,Real,Real>* brptr=dynamic_cast<const BinaryExpression<Tribool,Comparison,Real,Real>*>(eptr);
    if(brptr) { return _compute(bptr->_op,evaluate(brptr->_arg1,x),evaluate(brptr->_arg2,x)); }
    ARIADNE_FAIL_MSG("");
}

template<class X> X evaluate(const Expression<Real>& e, const Vector<X>& x) {
    const ExpressionInterface<Real>* eptr=e.ptr();
    const BinaryExpression<Real>* bptr=dynamic_cast<const BinaryExpression<Real>*>(eptr);
    if(bptr) { return _compute(bptr->_op,evaluate(bptr->_arg1,x),evaluate(bptr->_arg2,x)); }
    const UnaryExpression<Real>* uptr=dynamic_cast<const UnaryExpression<Real>*>(eptr);
    if(uptr) { return _compute(uptr->_op,evaluate(uptr->_arg,x)); }
    const ConstantExpression<Real>* cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
    if(cptr) { X r; _set_constant(r,cptr->value()); return r; }
    const CoordinateExpression<Real>* iptr=dynamic_cast<const CoordinateExpression<Real>*>(eptr);
    if(iptr) { return x[iptr->index()]; }
    ARIADNE_FAIL_MSG("");
}

template Tribool evaluate(const Expression<Tribool>& e, const ContinuousValuation<Float>& x);
template Float evaluate(const Expression<Real>& e, const ContinuousValuation<Float>& x);





template<class X>
Polynomial<X> polynomial(const Expression<Real>& e, const Map<Variable<Real>,uint>& s)
{
    const ExpressionInterface<Real>* const eptr=e.ptr();
    const Operator op=eptr->type();

    const ConstantExpression<Real>* cptr;
    const VariableExpression<Real>* vptr;
    const CoordinateExpression<Real>* iptr;
    const UnaryExpression<Real,Operator>* uptr;
    const BinaryExpression<Real,Operator>* bptr;

    switch(op) {
        case CNST:
            cptr=static_cast<const ConstantExpression<Real>*>(eptr);
            return Polynomial<X>::constant(s.size(),cptr->value()); break;
        case VAR:
            vptr=static_cast<const VariableExpression<Real>*>(eptr);
            return Polynomial<X>::variable(s.size(),s[vptr->variable()]); break;
        case IND:
            iptr=static_cast<const CoordinateExpression<Real>*>(eptr);
            return Polynomial<X>::variable(s.size(),iptr->index()); break;
        case NEG:
            uptr=static_cast<const UnaryExpression<Real>*>(eptr);
            return -polynomial<X>(uptr->_arg,s); break;
        case ADD:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            return polynomial<X>(bptr->_arg1,s)+polynomial<X>(bptr->_arg2,s); break;
        case SUB:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            return polynomial<X>(bptr->_arg1,s)-polynomial<X>(bptr->_arg2,s); break;
        case MUL:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            return polynomial<X>(bptr->_arg1,s)*polynomial<X>(bptr->_arg2,s); break;
        case DIV:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            cptr=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg2.ptr());
            ARIADNE_ASSERT_MSG(cptr,"Cannot convert expression "<<e<<" to polynomial form.");
            return polynomial<X>(bptr->_arg1,s)/cptr->value(); break;
        default:
            ARIADNE_FAIL_MSG("Cannot convert expression "<<e<<" to polynomial form.");
    }
}


ScalarPolynomialFunction polynomial(const Expression<Real>& e, const Space<Real>& s) {
    return ScalarPolynomialFunction(polynomial<Interval>(e,s.indices()));
}

ScalarAffineFunction affine(const Expression<Real>& e, const Map<String,uint>& s) {

    const ExpressionInterface<Real>* eptr=e.ptr();
    Operator op=eptr->type();

    const ConstantExpression<Real>* cptr;
    const VariableExpression<Real>* vptr;
    const CoordinateExpression<Real>* iptr;
    const UnaryExpression<Real,Operator>* uptr;
    const BinaryExpression<Real,Operator>* bptr;
    const ConstantExpression<Real>* cptr1;
    const ConstantExpression<Real>* cptr2;

    switch(op) {
        case CNST:
            cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
            return ScalarAffineFunction::constant(s.size(),cptr->value()); break;
        case VAR:
            vptr=dynamic_cast<const VariableExpression<Real>*>(eptr);
            return ScalarAffineFunction::variable(s.size(),s[vptr->name()]); break;
        case IND:
            iptr=dynamic_cast<const CoordinateExpression<Real>*>(eptr);
            return ScalarAffineFunction::variable(s.size(),iptr->coordinate()); break;
        case NEG:
            uptr=static_cast<const UnaryExpression<Real>*>(eptr);
            return -affine(uptr->_arg,s); break;
        case ADD:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            return affine(bptr->_arg1,s)+affine(bptr->_arg2,s); break;
        case SUB:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            return affine(bptr->_arg1,s)-affine(bptr->_arg2,s); break;
        case DIV:
            bptr=dynamic_cast<const BinaryExpression<Real>*>(eptr);
            cptr=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg2.ptr());
            ARIADNE_ASSERT_MSG(cptr,"Cannot convert expression "<<e<<" to affine form.");
            return affine(bptr->_arg1,s)/cptr->value(); break;
        case MUL:
            bptr=static_cast<const BinaryExpression<Real>*>(eptr);
            cptr1=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg1.ptr());
            cptr2=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg2.ptr());
            ARIADNE_ASSERT_MSG(cptr1 || cptr2,"Cannot convert expression "<<e<<" to affine form.");
            if(cptr1) { return cptr1->value() * affine(bptr->_arg2,s); }
            else { return affine(bptr->_arg1,s) * cptr2->value(); }
            break;
        default:
            ARIADNE_FAIL_MSG("Cannot convert expression "<<e<<" to affine form.");
    }
}

template<class X>
Affine<X> affine(const Expression<Real>& e, const Array<String>& s) {
    Map<String,uint> vars;
    for(uint i=0; i!=s.size(); ++i) { vars.insert(s[i],i); }
    return affine<X>(e,vars);
}







} // namespace Ariadne