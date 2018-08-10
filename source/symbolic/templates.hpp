/***************************************************************************
 *            symbolic/templates.hpp
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

/*! \file symbolic/templates.hpp
 *  \brief
 */


#ifndef ARIADNE_EXPRESSION_TEMPLATES_HPP
#define ARIADNE_EXPRESSION_TEMPLATES_HPP

#include "../numeric/operators.hpp"

namespace Ariadne {

/************ ExpressionTemplate ************************************************/

class Effort;

template<class Op, class... Args> struct ExpressionTemplate;

template<class R, class A> struct ConstantFunction {
    R _r;
  public:
    ConstantFunction(R const& r) : _r(r) { }
    R operator() (A const& a) const { return _r; }
    friend OutputStream& operator<<(OutputStream& os, ConstantFunction<R,A> const& f) { return os << f._r; }
};

template<class C> struct ExpressionTemplate<Cnst,C> {
    Cnst _op; C _cnst;
    ExpressionTemplate(Cnst o, C c) : _op(o), _cnst(c) { }
    template<class T> operator T() const { return static_cast<T>(_cnst); }
    template<class... AS> decltype(auto) operator() (AS... vals) const { return _cnst; }
    template<class PR> decltype(auto) get(Effort const& eff, PR pr) const { return _cnst.compute(eff).get(pr); }
    friend OutputStream& operator<<(OutputStream& os, ExpressionTemplate expr) { return os << expr._cnst; }
};


template<class C, class A> struct ExpressionTemplate<Cnst,ConstantFunction<C,A>> {
    Cnst _op; ConstantFunction<C,A> _cnst;
    ExpressionTemplate(Cnst o, ConstantFunction<C,A> c) : _op(o), _cnst(c) { }
    template<class T> operator T() const { return static_cast<T>(_cnst); }
    template<class... AS> decltype(auto) operator() (AS... vals) const { return _cnst(vals...); }
    friend OutputStream& operator<<(OutputStream& os, ExpressionTemplate expr) { return os << expr._cnst; }
};

template<class I> struct ExpressionTemplate<Var,I> {
    Var _op; I _ind;
    ExpressionTemplate(Var o, I i) : _op(o), _ind(i) { }
    template<class AS> auto operator() (AS vals) const -> decltype(vals[_ind]) { return vals[_ind]; }
    template<class F> decltype(auto) apply(F const& f) const { return assert(false); return *this; }
    friend OutputStream& operator<<(OutputStream& os, ExpressionTemplate expr) { return os << expr._ind; }
};

struct Id { template<class A> A operator() (A const& a) { return a; } };

template<> struct ExpressionTemplate<Id> {
    Id _op;
    ExpressionTemplate(Id o) : _op(o) { }
    template<class A> auto operator() (A vals) const -> A { return vals; }
    friend OutputStream& operator<<(OutputStream& os, ExpressionTemplate expr) { return os << "id"; }
};

template<class O, class A> struct ExpressionTemplate<O,A> {
    O _op; A _arg;
    ExpressionTemplate(O o, A a) : _op(o), _arg(a) { }
    template<class T> operator T() const { return _op(static_cast<T>(_arg)); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg(vals...))) {
        return _op(_arg(vals...)); }
    template<class F> decltype(auto) apply(F const& f) const { return _op(f(_arg)); }
    template<class... PRS> decltype(auto) get(PRS const& ... prs) const { return _op(_arg.get(prs...)); }
    friend OutputStream& operator<<(OutputStream& os, ExpressionTemplate expr) {
        return os << expr._op.code() << "(" << expr._arg << ")"; }
};

template<class O, class A1, class A2> struct ExpressionTemplate<O,A1,A2> {
    O _op; A1 _arg1; A2 _arg2;
    ExpressionTemplate(O o, A1 a1, A2 a2) : _op(o), _arg1(a1), _arg2(a2) { }
    template<class T> operator T() const { return _op(static_cast<T>(_arg1),static_cast<T>(_arg2)); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg1(vals...),_arg2(vals...))) {
        return _op(_arg1(vals...),_arg2(vals...)); }
    template<class F> decltype(auto) apply(F const& f) const { return _op(f(_arg1),f(_arg2)); }
    template<class... PRS> decltype(auto) get(PRS const& ... prs) const { return _op(_arg1.get(prs...),_arg2.get(prs...)); }
    friend OutputStream& operator<<(OutputStream& os, ExpressionTemplate expr) {
        return os << expr._op.code() << "(" << expr._arg1 << "," << expr._arg2 << ")"; }
};

template<class A, class N> struct ExpressionTemplate<Pow,A,N> {
    Pow _op; A _arg; N _n;
    ExpressionTemplate(Pow o, A a, N n) : _op(o), _arg(a), _n(n) { }
    template<class T> operator T() const { return _op(static_cast<T>(_arg),_n); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg(vals...),_n)) {
        return _op(_arg(vals...),_n); }
    template<class F> decltype(auto) apply(F const& f) const { return _op(f(_arg),_n); }
    friend OutputStream& operator<<(OutputStream& os, ExpressionTemplate expr) {
        return os << expr._op.code() << "(" << expr._arg << "," << expr._n << ")"; }
};

template<class O, class A1, class A2, class A3> struct ExpressionTemplate<O,A1,A2,A3> {
    O _op; A1 _arg1; A2 _arg2; A3 _arg3;
    ExpressionTemplate(O o, A1 a1, A2 a2, A3 a3) : _op(o), _arg1(a1), _arg3(a3) { }
    template<class T> operator T() const { return _op(static_cast<T>(_arg1),static_cast<T>(_arg2),static_cast<T>(_arg3)); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_arg1(vals...),_arg2(vals...),_arg3(vals...))) {
        return _op(_arg1(vals...),_arg2(vals...)); }
    template<class F> decltype(auto) apply(F const& f) const { return _op(f(_arg1),f(_arg2),f(_arg3)); }
    friend OutputStream& operator<<(OutputStream& os, ExpressionTemplate expr) {
        return os << expr._op.code() << "(" << expr._arg1 << "," << expr._arg2 << "," << expr._arg3 << ")"; }
};

template<class O, class... Args> ExpressionTemplate<O,Args...> make_expression_template(O o, Args... as) {
    return ExpressionTemplate<O,Args...>(o, as...); }

template<class O, class... Args> struct TemporaryExpression;

template<class O, class A> struct TemporaryExpression<O,A> {
    O _op; A const& _arg;
    TemporaryExpression(O o, A const& a) : _op(o), _arg(a) { }
};

template<class O, class A1, class A2> struct TemporaryExpression<O,A1,A2> {
    O _op; A1 const& _arg1; A2 const& _arg2;
    TemporaryExpression(O o, A1 const& a1, A2 const& a2) : _op(o), _arg1(a1), _arg2(a2) { }
};

template<class O, class A1, class A2, class A3> struct TemporaryExpression<O,A1,A2,A3> {
    O _op; A1 const& _arg1; A2 const& _arg2; A3 const& _arg3;
    TemporaryExpression(O o, A1 const& a1, A2 const& a2, A3 const& a3) : _op(o), _arg1(a1), _arg3(a3) { }
};


} // namespace Ariadne

#endif
