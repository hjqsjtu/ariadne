/***************************************************************************
 *            test_misc.cpp
 *
 *  Copyright 2008--17 Pieter Collins
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

#include <iostream>
#include "config.h"
#include "test.hpp"

#include "function/c1_taylor_function.hpp"

namespace Ariadne { }

using namespace Ariadne;

#define ARIADNE_PRINT(expr) { std::cout << #expr << ": " << (expr) << "\n"; }

Int main() {
    C1TaylorSeries x=C1TaylorSeries::coordinate();
    //C1TaylorSeries e=C1TaylorSeries::uniform_ball();
    C1TaylorSeries e0=C1TaylorSeries::uniform_ball();
    C1TaylorSeries e1=C1TaylorSeries::derivative_ball();
    auto y=x;

    ARIADNE_PRINT(x);
    ARIADNE_PRINT(x*x);
    ARIADNE_PRINT(2*x);
    ARIADNE_PRINT(2*(x*x));
    ARIADNE_PRINT((2*x)*x);
    auto f=2*x*x-1+e0/8+e1/8;
    auto g=(y/2+1)*y+1;
    auto h=compose(g,f);
    ARIADNE_PRINT(f);
    ARIADNE_PRINT(g);
    ARIADNE_PRINT(h);
    ARIADNE_PRINT(f*f);
    ARIADNE_PRINT(f*f/2);
}
