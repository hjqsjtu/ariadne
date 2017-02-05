/***************************************************************************
 *            numeric/twoexp.h
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/twoexp.h
 *  \brief
 */

#ifndef ARIADNE_TWOEXP_H
#define ARIADNE_TWOEXP_H

#include <cmath>
#include "float.decl.h"

namespace Ariadne {

/************ TwoExp ***************************************************/

//! \ingroup NumericModule
//! \brief A class representing a number of the form  \c 2<sup>\it n</sup> for some \it n.
//! Useful since floating-point numbers can be exactly multiplied and divided by powers of \c 2.
class TwoExp {
    Int _n;
  public:
    explicit TwoExp(Int n) : _n(n) { }
    Int exponent() const { return this->_n; }
    // NOTE: Use std::pow(2.0,_n) not (1<<_n) since latter does not handle very large exponents
    Float64 get_raw(Precision64 pr) const;
    FloatMP get_raw(PrecisionMP pr) const;
    double get_d() const { return std::pow(2.0,this->_n); }
    operator Float64Value () const;
    operator Float64Error () const;
    operator Float64Ball () const;
    operator Float64Bounds () const;
    friend Float64Value operator+(TwoExp);
    friend Float64Value operator-(TwoExp);
};
inline TwoExp two_exp(Int n) { return TwoExp(n); }

/*
inline Float64 TwoExp::get_raw(Precision64 pr) const {
    return Float64(this->get_d());
}
inline FloatMP TwoExp::get_raw(PrecisionMP pr) const {
    return FloatMP(this->get_d(),pr);
}
*/

} // namespace Ariadne

#endif