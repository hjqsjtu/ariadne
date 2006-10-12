/***************************************************************************
 *            python/export_polyhedron.cc
 *
 *  21 October 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "real_typedef.h"

#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;


template<typename R>
void export_polyhedron() 
{
  typedef Polyhedron<R> RPolyhedron;
  typedef Rectangle<R> RRectangle;
  typedef Polytope<R> RPolytope;

  def("disjoint", (bool(*)(const RPolyhedron&, const RPolyhedron&))(&disjoint));
  def("interiors_intersect", (bool(*)(const RPolyhedron&, const RPolyhedron&))(&interiors_intersect));
  def("inner_subset", (bool(*)(const RPolyhedron&, const RPolyhedron&))(&inner_subset));
  def("subset", (bool(*)(const RPolyhedron&, const RPolyhedron&))(&subset));
  def("convex_hull", (RPolyhedron(*)(const RPolyhedron&, const RPolyhedron&))(&convex_hull));

  class_<RPolyhedron>("Polyhedron",init<size_type>())
    .def(init<RPolyhedron>())
    .def(init<RRectangle>())
    .def("dimension", &RPolyhedron::dimension)
    .def(self_ns::str(self))
  ;
  
}

template void export_polyhedron<MPFloat>();
