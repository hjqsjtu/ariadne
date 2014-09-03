/***************************************************************************
 *            test_solvers.cc
 *
 *  Copyright  2008-10  Pieter Collins
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
#include <sstream>
#include <string>

#include "config.h"

#include "solver.h"
#include "function.h"
#include "taylor_function.h"
#include "vector.h"
#include "expression.h"
#include "space.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

class TestSolver
{
  private:
    std::unique_ptr<SolverInterface> solver;
    std::string solver_class_name;
  public:
    TestSolver(const SolverInterface& s,const char* n)
        : solver(s.clone()), solver_class_name(n) { }

    int test() {
        ARIADNE_TEST_PRINT(*solver);
        ARIADNE_TEST_CALL(test_solve());
        ARIADNE_TEST_CALL(test_implicit());
        ARIADNE_TEST_CALL(test_scalar_implicit());
        return 0;
    }

    void test_solve() {
        RealScalarFunction x=RealScalarFunction::coordinate(1,0);
        IntervalVector d({Interval(0.0,1.0)});
        RealVectorFunction f({(x*x+1)*x-1});
        IntervalVector p=solver->solve(f,d);
        ARIADNE_TEST_BINARY_PREDICATE(subset,p[0],Interval(0.6823,0.6824));
    }

    void test_implicit() {
        //TaylorModelAccuracy::set_default_sweep_threshold(1e-12);

        RealScalarFunction aa=RealScalarFunction::coordinate(1,0);
        RealScalarFunction a=RealScalarFunction::coordinate(2,0);
        RealScalarFunction x=RealScalarFunction::coordinate(2,1);
        RealScalarFunction bb;
        IntervalVector p,r;
        RealVectorFunction f;
        IntervalVectorFunctionModel h;
        RealVectorFunction e;

        // Test solution of x-a=0. This should be very easy to solve.
        p=IntervalVector({Interval(-0.25,0.25)});
        r=IntervalVector({Interval(-2.0,2.0)});
        f=RealVectorFunction({x-a});
        ARIADNE_TEST_PRINT(f);
        h=solver->implicit(f,p,r);
        ARIADNE_TEST_PRINT(h);
        e=RealVectorFunction(1u,aa);
        ARIADNE_TEST_PRINT(e);
        ARIADNE_TEST_COMPARE(norm((h-e).range()),<,1e-8);

        // Test solution of 4x^2+x-4-a=0 on [0.875,1.125]. There is a unique solution with positive derivative.
        p=IntervalVector({Interval(0.875,1.125)});
        r=IntervalVector({Interval(0.25,1.25)});
        f=RealVectorFunction({(x*x+1)*x-a});
        ARIADNE_TEST_PRINT(f);
        h=solver->implicit(f,p,r);
        ARIADNE_TEST_PRINT(h);
        bb=RealScalarFunction(aa-numeric_cast<Real>(p[0].midpoint()))/numeric_cast<Real>(p[0].radius());
        Decimal a0(0.682328), a1(0.0521547), a2(-0.0023232), a3(0.000147778);
        e=RealVectorFunction( { a0+bb*(a1+bb*(a2+bb*a3)) } );
        ARIADNE_TEST_PRINT(e);
        ARIADNE_TEST_COMPARE(norm((h-e).range()),<,1e-4);

        // Test solution of 4x^2+x-4-a=0 on [-0.25,0.25]. There is a unique solution with positive derivative.
        p=IntervalVector({Interval(-0.25,0.25)});
        r=IntervalVector({Interval(0.25,2.0)});
        f=RealVectorFunction({4*x+x*x-a-4});
        ARIADNE_TEST_PRINT(f);
        h=solver->implicit(f,p,r);
        ARIADNE_TEST_PRINT(h);
        bb=RealScalarFunction(aa-numeric_cast<Real>(p[0].midpoint()))/numeric_cast<Real>(p[0].radius());
        Decimal c0(0.828427), c1(0.0441942), c2(-0.000345267), c3(0.00000539468);
        e=RealVectorFunction( { c0+bb*(c1+bb*(c2+bb*c3)) } );
        ARIADNE_TEST_PRINT(e);
        ARIADNE_TEST_COMPARE(norm((h-e).range()),<,1e-4);

        // Test solution of x-2*a=0 on [-1,+1], taking values in [-1,+1].
        // There is at most one solution, but this lies partially outside the range.
        // Should obtain PartialSolutionException
        p=IntervalVector({Interval(-1,1)});
        r=IntervalVector({Interval(-1,1)});
        f=RealVectorFunction({x-2*a});
        ARIADNE_TEST_PRINT(f);
        try {
            h=solver->implicit(f,p,r);
            ARIADNE_TEST_NOTIFY(solver_class_name<<" silently returns partially defined vector implicit function.");
        }
        catch(SolverException) {
            ARIADNE_TEST_THROWS(solver->implicit(f,p,r),SolverException);
            ARIADNE_TEST_NOTIFY(solver_class_name<<" throws error on partially defined vector implicit function.");
        }

    }

    void test_scalar_implicit() {
        //TaylorModelAccuracy::set_default_sweep_threshold(1e-12);

        RealScalarFunction aa=RealScalarFunction::coordinate(1,0);
        RealScalarFunction a=RealScalarFunction::coordinate(2,0);
        RealScalarFunction x=RealScalarFunction::coordinate(2,1);
        // Test solution of x-2*a=0 on [-1,+1], taking values in [-1,+1]. There is at most one solution.
        // Uses scalar implicit
        IntervalVector p; Interval r;
        RealScalarFunction e,f,s; // s is unscaling functions
        IntervalScalarFunctionModel h;

        ARIADNE_TEST_PRINT(*solver);

        // Test solution of 4x^2+x-4-a=0 on [0.875,1.125]. There is a unique solution with positive derivative.
        p=IntervalVector({Interval(0.875,1.125)});
        r=Interval(0.25,1.25);
        f=RealScalarFunction((x*x+1)*x-a);
        ARIADNE_TEST_PRINT(f);
        h=solver->implicit(f,p,r);
        ARIADNE_TEST_PRINT(h);
        s=RealScalarFunction(aa-numeric_cast<Real>(p[0].midpoint()))/numeric_cast<Real>(p[0].radius());
        Decimal a0(0.682328), a1(0.0521547), a2(-0.0023232), a3(0.000147778);
        e=RealScalarFunction( a0+s*(a1+s*(a2+s*a3)) );
        ARIADNE_TEST_PRINT(e);
        ARIADNE_TEST_COMPARE(mag((h-e).range()),<,1e-4);

        // Test solution of x-2*a=0 on [-1,+1], taking values in [-1,+1].
        // There is at most one solution, but this lies partially outside the range.
        // Should obtain PartialSolutionException
        p=IntervalVector({Interval(-1,1)});
        r=Interval(-1,1);
        f=RealScalarFunction(x-2*a);
        ARIADNE_TEST_PRINT(f);
        IntervalScalarFunctionModel g=ScalarTaylorFunction(join(p,r),f,ThresholdSweeper(1e-12));
        ARIADNE_TEST_PRINT(g);
        try {
            h=solver->implicit(g,p,r);
            ARIADNE_TEST_NOTIFY(solver_class_name<<" silently returns partially defined scalar implicit function.");
        }
        catch(SolverException) {
            ARIADNE_TEST_THROWS(solver->implicit(g,p,r),SolverException);
            ARIADNE_TEST_NOTIFY(solver_class_name<<" throws error on partially defined scalar implicit function.");
        }

    }

};

#include "differential.h"

int main(int argc, const char **argv) {
/*
    IntervalVector D={{-1,+1},{-1,+1}};
    VectorTaylorFunction x=VectorTaylorFunction::identity(D,ThresholdSweeper(1e-10));
    ScalarTaylorFunction f=2*x[0]-x[1];
    std::cerr<<"D="<<D<<"\n";
    Interval D0=D[0];
    std::cerr<<"f="<<representation(f)<<"\n";
    Vector<Differential<Interval>> dx=Differential<Interval>::variables(2,2,1u,D);
    Differential<Interval> dx0=dx[0];
    std::cerr<<"dx="<<dx<<"\n";
    std::cerr<<"unscale(dx[0],D[0])="<<unscale(dx0,D[0])<<"\n";
    Interval c(add_ivl(D0.lower()/2,D0.upper()/2));
    Interval r(sub_ivl(D0.upper()/2,D0.lower()/2));
    std::cerr<<"(dx[0]-c)/r="<<((dx0-c)/r)<<" c="<<c<<" r="<<r<<"\n";
    std::cerr<<"unscale(dx,D)="<<unscale(dx,D)<<"\n";
    std::cerr<<"f(dx)="<<f.evaluate(dx)<<"\n";
    return 0;
*/
    int verbosity=get_verbosity(argc,argv);

    IntervalNewtonSolver interval_newton_solver(maximum_error=1e-5,maximum_number_of_steps=12);
    interval_newton_solver.verbosity=verbosity;
    TestSolver(interval_newton_solver,"IntervalNewtonSolver").test();

    KrawczykSolver krawczyk_solver(maximum_error=1e-5,maximum_number_of_steps=12);
    ARIADNE_TEST_PRINT(krawczyk_solver.function_factory());
    krawczyk_solver.verbosity=verbosity;
    TestSolver(krawczyk_solver,"KrawczykSolver").test();

    FactoredKrawczykSolver factored_krawczyk_solver(maximum_error=1e-5,maximum_number_of_steps=12);
    factored_krawczyk_solver.verbosity=verbosity;
    TestSolver(factored_krawczyk_solver,"FactoredKrawczykSolver").test();

    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
