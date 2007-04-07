/***************************************************************************
 *            test.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*!\file test.h 
 * \brief Macros for test suite. 
 */

#include <iostream>
#include <exception>

int ARIADNE_TEST_FAILURES=0;

// This needs to be a function since we do not want to evaluate the result twice,
// and can't store it in a variable since we don't know it's type.
template<class R, class ER> 
bool 
ariadne_check(std::ostream& os, const R& r, const ER& er) {
  os << r; return (r==er);
}

/*! \brief Catches an exception and writes a diagnostic to standard output and standard error. */
#define ARIADNE_CATCH(message) \
  catch(const std::exception& e) {          \
    ++ARIADNE_TEST_FAILURES;                                            \
    std::cout << "exception: \"" << e.what() << "\"\n" << std::endl;;                     \
    std::cerr << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ << ": " << message << " throwed \"" << e.what() << "\"." << std::endl; \
  }                                                                     \
  catch(...) { \
    ++ARIADNE_TEST_FAILURES;                \
    std::cout << "unknown exception\n" << std::endl; \
    std::cerr << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ << ": " << message << " throwed an unknown exception." << std::endl; \
  } \


/*! \brief Tries to execute \a statement, writing the statement to standard output. Writes a diagnostic report to standard error if an exception is thrown. <br> <b>Important:</b> Use the ARIADNE_CONSTRUCT() macro if \a statement declares a variable and calls a constructor. */
#define ARIADNE_TRY(statement) \
{ \
  std::cout << #statement << ": " << std::flush; \
  try { \
    statement;                      \
    std::cout << "ok\n" << std::endl; \
  } \
  ARIADNE_CATCH("Statement `" << #statement << "'")       \
}\


/*! \brief Tries to evaluate \a expression, writing the expression and the result to standard ouput. Writes a diagnostic report to standard error if an exception is thrown., prints the result and gives a diagnostic if an exception is thrown. */
#define ARIADNE_EVALUATE(expression) \
{ \
  std::cout << #expression << ": " << std::flush; \
  try { \
    std::cout << (expression) << "\n" << std::endl;      \
  } \
  ARIADNE_CATCH("Expression `" << #expression << "'")       \
}\

/*! \brief Evaluates \a expression in a boolean context and checks if the result is \a true. */
#define ARIADNE_TEST_ASSERT(expression)               \
{ \
  std::cout << #expression << ": " << std::flush; \
  bool result = (expression);                     \
  if(result) { \
    std::cout << "true\n" << std::endl; \
  } else { \
    ++ARIADNE_TEST_FAILURES; \
    std::cout << "false" << std::endl; \
    std::cerr << __FILE__ << ":" << __LINE__ << ": " << __FUNCTION__ << ": Assertion `" << #expression << "' failed.\n" << std::endl; \
  } \
}\



/*! \brief Evaluates \a expression and checks if the result is equal to \a expected. */
#define ARIADNE_CHECK(expression,expected) \
{ \
  std::cout << #expression << ": " << std::flush; \
  bool ok = ariadne_check(std::cout,expression,expected); \
  if(ok) { \
    std::cout << "\n" << std::endl; \
  } else { \
    ++ARIADNE_TEST_FAILURES; \
    std::cout << " (expected: " << #expected << ")\n" << std::endl;        \
    std::cerr << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ << ": Check `" << #expression << "==" << #expected << "' failed." << std::endl; \
  } \
}\


/*! \brief Constructs object \a variable of type \a class from \a expression. */
#define ARIADNE_CONSTRUCT(class,variable,expression) \
{ \
  std::cout << #class << " " << #variable << "" << #expression << ": " << std::flush; \
  try { \
    class variable expression;                                         \
    std::cout << #variable << "==" << variable << "\n" << std::endl;     \
  }                                       \
  ARIADNE_CATCH("Constructor `" << #class << "" << #variable << "" << #expression << "'") \
} \
class variable expression; \


/*! \brief Assigns object \a variable from \a expression. */
#define ARIADNE_ASSIGN(variable, expression) \
{ \
  std::cout << #variable << "=" << #expression << ": " << std::flush; \
  try { \
    variable=(expression);                      \
    std::cout << variable << "\n" << std::endl; \
  }                                       \
  ARIADNE_CATCH("Assignment `" << #variable << "=" << #expression << "'")   \
} \

