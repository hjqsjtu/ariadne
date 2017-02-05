/***************************************************************************
 *            test_constraint_solver.cc
 *
 *  Copyright  2010  Pieter Collins
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
#include <fstream>

#include "config.h"
#include "test.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "function/function.h"
#include "function/constraint.h"
#include "function/procedure.h"
#include "solvers/constraint_solver.h"
#include "geometry/box.h"


using namespace std;
using namespace Ariadne;

class TestConstraintSolver
{
    Nat verbosity;
  public:
    TestConstraintSolver(Nat v) : verbosity(v) { }

    Void test() {
        ARIADNE_TEST_CALL(test_empty_reduce_inequality());
        ARIADNE_TEST_CALL(test_empty_reduce_equality());
        ARIADNE_TEST_CALL(test_empty_reduce_mixed());
        ARIADNE_TEST_CALL(test_empty_hull_reduce());
        ARIADNE_TEST_CALL(test_empty_box_reduce());
        ARIADNE_TEST_CALL(test_hull_reduce());
        ARIADNE_TEST_CALL(test_box_reduce());
        ARIADNE_TEST_CALL(test_monotone_reduce());
        ARIADNE_TEST_CALL(test_split());
    }

    Void test_empty_reduce_inequality() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,1.0},{0.0,1.0}};
        List<EffectiveConstraint> c = {4<=2*x[0]+x[1]};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_reduce_equality() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,1.0},{0.0,1.0}};
        List<EffectiveConstraint> c = {2*x[0]+x[1]==4};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_reduce_mixed() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,0.25},{0.0, 2.0}};
        List<EffectiveConstraint> c = {x[1]<=1,x[0]+x[1]==2};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_hull_reduce() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,0.25},{0.0,2.0}};
        List<EffectiveConstraint> c = {x[1]<=1, x[0]+x[1]==2};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[1]));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_box_reduce() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,0.25},{0.0, 2.0}};
        List<EffectiveConstraint> c = {x[1]<=1,x[0]+x[1]==2};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[0],0));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[1],0));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[0],1));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[1],1));
        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[1]));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_hull_reduce() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,2.0},{0.0,2.0}};
        List<EffectiveConstraint> c = {-2<=2*x[0]+x[1]<=1};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0,0.5},{0.0,1.0}}));
    }

    Void test_box_reduce() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,2.0},{0.0,2.0}};
        EffectiveConstraint c = (-2<=2*x[0]+x[1]<=1);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,0));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0,0.75},{0.0,2.0}}));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,1));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0,0.75},{0.0,1.25}}));
    }


    Void test_monotone_reduce() {
        List<EffectiveScalarFunction> x=EffectiveScalarFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,2.0},{0.0,2.0}};
        EffectiveConstraint c = (-2<=2*x[0]+x[1]<=1);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,0));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0,0.75},{0.0,2.0}}));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,1));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0,0.75},{0.0,1.25}}));
    }

    Void test_split() {
        ARIADNE_TEST_WARN("test_split: Not implemented");
    }
};

Int main(Int argc, const char* argv[]) {
    TestConstraintSolver(get_verbosity(argc,argv)).test();
    return ARIADNE_TEST_FAILURES;
}
