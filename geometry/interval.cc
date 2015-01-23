/***************************************************************************
 *            interval.cc
 *
 *  Copyright 2008-10  Pieter Collins
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

#include "utility/standard.h"
#include "config.h"

#include <iostream>
#include <iomanip>
#include <cassert>
#include "utility/container.h"

#include "utility/macros.h"
#include "utility/exceptions.h"
#include "numeric/integer.h"
#include "numeric/float.h"
#include "numeric/rational.h"
#include "numeric/float-exact.h"
#include "numeric/decimal.h"
#include "numeric/dyadic.h"
#include "geometry/interval.h"


namespace Ariadne {


const UpperInterval pi_ivl=ExactInterval(pi_down,pi_up);

Nat ExactInterval::output_precision = 6;

ExactInterval widen(ExactInterval x)
{
    rounding_mode_t rm=get_rounding_mode();
    const double& xl=internal_cast<const double&>(x.lower().raw());
    const double& xu=internal_cast<const double&>(x.upper().raw());
    const double m=std::numeric_limits<float>::min();
    set_rounding_upward();
    volatile double wu=xu+m;
    volatile double mwl=-xl+m;
    volatile double wl=-mwl;
    set_rounding_mode(rm);
    assert(wl<xl); assert(wu>xu);
    return ExactInterval(wl,wu);
}

ExactInterval narrow(ExactInterval x)
{
    rounding_mode_t rm=get_rounding_mode();
    const double& xl=internal_cast<const double&>(x.lower().raw());
    const double& xu=internal_cast<const double&>(x.upper().raw());
    const double m=std::numeric_limits<float>::min();
    set_rounding_upward();
    volatile double mnu=-xu+m;
    volatile double nu=-mnu;
    volatile double nl=xl+m;
    set_rounding_mode(rm);
    assert(xl<nl); assert(nu<xu);
    return ExactInterval(nl,nu);
}

ExactInterval trunc(ExactInterval x)
{

    rounding_mode_t rm=get_rounding_mode();
    const double& xl=internal_cast<const double&>(x.lower().raw());
    const double& xu=internal_cast<const double&>(x.upper().raw());
    // Use machine epsilon instead of minimum to move away from zero
    const float fm=std::numeric_limits<float>::epsilon();
    volatile float tu=xu;
    if(tu<xu) { set_rounding_upward(); tu+=fm; }
    volatile float tl=xl;
    if(tl>xl) { set_rounding_downward(); tl-=fm; }
    set_rounding_mode(rm);
    assert(tl<=xl); assert(tu>=xu);
    return ExactInterval(double(tl),double(tu));
}

ExactInterval trunc(ExactInterval x, Nat n)
{
    ExactInterval e=ExactInterval(ExactFloat(std::pow(2.0,52-(Int)n)));
    UpperInterval y=x+e;
    UpperInterval r=y-e;
    return ExactInterval(r.lower_raw(),r.upper_raw());
}



ExactInterval::ExactInterval(const Float& x) : ExactInterval(x,x) { }

ExactInterval::ExactInterval(const Dyadic& b) : ExactInterval(b.operator Rational()) { }

ExactInterval::ExactInterval(const Decimal& d) : ExactInterval(d.operator Rational()) { }


#ifdef HAVE_GMPXX_H

ExactInterval::ExactInterval(const Integer& z) : ExactInterval(Rational(z)) {
}

ExactInterval::ExactInterval(const Rational& q) : ExactInterval(q,q) {
}

ExactInterval::ExactInterval(const Rational& ql, const Rational& qu) : l(ql.get_d()), u(qu.get_d())  {
    static const double min_dbl=std::numeric_limits<double>::min();
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(downward);
    while(Rational(l)>ql) { l=sub_rnd(l,min_dbl); }
    set_rounding_mode(upward);
    while(Rational(u)<qu) { u=add_rnd(u,min_dbl); }
    set_rounding_mode(rounding_mode);
}

ExactInterval::ExactInterval(const Real& lower, const Real& upper)
    : l(ValidatedFloat(lower).lower().raw()), u(ValidatedFloat(upper).upper().raw()) {
}

ExactInterval& ExactInterval::operator=(const Rational& q) {
    return *this = ExactInterval(q);
}


#endif // HAVE_GMPXX_H


OutputStream&
operator<<(OutputStream& os, const ExactInterval& ivl)
{
    //if(ivl.lower()==ivl.upper().raw()) { return os << "{" << std::setprecision(ExactInterval::output_precision) << ivl.lower().raw().get_d() << ; }
    rounding_mode_t rnd=get_rounding_mode();
    os << '{';
    set_rounding_downward();
    os << std::showpoint << std::setprecision(ExactInterval::output_precision) << ivl.lower().raw().get_d();
    os << ':';
    set_rounding_upward();
    os << std::showpoint << std::setprecision(ExactInterval::output_precision) << ivl.upper().raw().get_d();
    set_rounding_mode(rnd);
    os << '}';
    return os;

}

/*
OutputStream&
operator<<(OutputStream& os, const ExactInterval& ivl)
{
    return os << '[' << ivl.l << ':' << ivl.u << ']';
}
*/

/*
OutputStream&
operator<<(OutputStream& os, const ExactInterval& ivl)
{
    if(ivl.lower().raw()==ivl.upper().raw()) {
        return os << std::setprecision(18) << ivl.lower().raw();
    }

    StringStream iss,uss;
    iss << std::setprecision(18) << ivl.lower().raw();
    uss << std::setprecision(18) << ivl.upper().raw();

    StringType lstr,ustr;
    iss >> lstr; uss >> ustr;

    // Test if one endpoint is an integer and the other is not
    // If this is the case, append ".0" to to integer value
    if( (lstr.find('.')==lstr.size()) xor (ustr.find('.')==lstr.size()) ) {
        if(lstr.find('.')==lstr.size()) {
            lstr+=".0";
        } else {
            ustr+=".0";
        }
    }

    // Write common head
    Nat i;
    for(i=0; (i<std::min(lstr.size(),ustr.size()) && lstr[i]==ustr[i]); ++i) {
        os << lstr[i];
    }

    os << "[";
    if(i==lstr.size()) {
        os << "0";
    }
    for(Nat li=i; li != lstr.size(); ++li) {
        os << lstr[li];
    }
    os << ":";
    if(i==ustr.size()) {
        os << "0";
    }
    for(Nat ui=i; ui != ustr.size(); ++ui) {
        os << ustr[ui];
    }
    os << "]";
    return os;

}
*/
InputStream&
operator>>(InputStream& is, ExactInterval& ivl)
{
    Float l,u;
    char cl,cm,cr;
    is >> cl >> l >> cm >> u >> cr;
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    ivl.set(l,u);
    return is;
}



} // namespace Ariadne

