/***************************************************************************
 *            multifunction.hpp
 *
 *  Copyright 2008-12  Pieter Collins
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

/*! \file multifunction.hpp
 *  \brief Multivalued functions
 */

#ifndef ARIADNE_MULTIFUNCTION_HPP
#define ARIADNE_MULTIFUNCTION_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "utility/declarations.hpp"
#include "utility/macros.hpp"
#include "utility/pointer.hpp"
#include "utility/container.hpp"
#include "utility/metaprogramming.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/differential.hpp"

#include "function/function.hpp"
#include "function/function_model.hpp"

#include "geometry/interval.hpp"
#include "geometry/box.hpp"
#include "geometry/function_set.hpp"

namespace Ariadne {


template<class P> class GenericPoint {
    class Interface {
        virtual ~Interface() = 0;
        virtual Vector<Interval<Number<P>>> _bounds() const;
        GenericPoint<P> _apply(VectorFunction<P> const& f) const;
    };
    SharedPointer<const Interface> _ptr;
  public:
    GenericPoint(SharedPointer<const Interface> ptr) : _ptr(ptr) { }
    Vector<Interval<Number<P>>> bounds() { return this->_ptr->_bounds(); }
    Interval<Number<P>> operator[](SizeType i) { return this->bounds()[i]; }
    friend GenericPoint<P> apply(VectorFunction<P> const& f, GenericPoint<P> const& pt) { return pt._ptr->_apply(f); }
};

template<class P> class ImageGenericPoint : public GenericPoint<P>::Interface {
    VectorFunctionModel<P> _h;
    ImageGenericPoint(VectorFunctionModel<P> const& h) : _h(h) { }
    friend ImageGenericPoint<P> apply(VectorFunction<P> const& f, GenericPoint<P> const& pt) {
        return ImageGenericPoint<P>(f(pt._h)); }
};

class LocatedSet;
class ValidatedLocatedSet;

template<class P> struct LocatedSetTypedef;
template<class P> using LocatedSetType = typename LocatedSetTypedef<P>::Type;
template<> struct LocatedSetTypedef<ValidatedTag> { typedef ValidatedLocatedSet Type; };
template<> struct LocatedSetTypedef<EffectiveTag> { typedef LocatedSet Type; };

template<class P> class Multifunction;
template<class P> class MultifunctionInterface {
    friend class Multifunction<P>;
    virtual ~MultifunctionInterface<P>() = default;
    virtual LocatedSetType<P> _call(Vector<Number<P>> const& x) const;
};

template<class P> class Multifunction {
    SharedPointer<MultifunctionInterface<P>> _ptr;
    typedef BoxDomain D;
    typedef Number<P> Y;
    typedef LocatedSetType<P> S;

    typedef D DomainType;

    LocatedSetType<P> operator() (Vector<Number<P>> const& x) const { return this->_ptr->_call(x); }
    friend LocatedSetType<P> apply (Multifunction<P> const& f, LocatedSetType<P> const& x);
};


using ValidatedImageSet = ValidatedConstrainedImageSet;

template<class P> class FunctionModelMultifunction {
    VectorFunctionModel<P> _f;
    SizeType _ne;
    BoxDomain domain();
    BoxDomain error_domain();
    ValidatedImageSet operator() (Vector<ValidatedNumber> const& x) {
        return ValidatedImageSet(compose(_f,_f.create_constant(error_domain(),x),_f.create_identity(error_domain()))); }
};

} // namespace Ariadne

#endif
