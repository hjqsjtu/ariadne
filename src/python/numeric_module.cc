/***************************************************************************
 *            python/numeric_module.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include <boost/python.hpp>

#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"

using namespace Ariadne::Numeric;

void export_numeric();
template<class R> void export_function();
template<class R> void export_interval();

BOOST_PYTHON_MODULE(numeric)
{
  export_numeric();
  export_function<MPFloat>();
  export_interval<Float64>();
  export_interval<MPFloat>();
}
