/***************************************************************************
 *            macros/assert.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_MACROS_ASSERT_H
#define ARIADNE_MACROS_ASSERT_H

#include <iostream>
#include <cassert>
#include <exception>

#include "macros/throw.h"

/*! \brief Evaluates \a expression in a boolean context and checks if the result is \a true. */
#define ARIADNE_ASSERT(expression) \
{ \
  bool result = (expression); \
  if(!result) { \
    ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<__PRETTY_FUNCTION__,"Assertion `" << #expression << "' failed.\n"); \
  } \
} \

#endif // ARIADNE_MACROS_ASSERT_H
