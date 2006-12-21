/***************************************************************************
 *            rectangle.h
 *
 *  Mon 2 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file rectangle.h
 *  \brief Rectangles and cuboids.
 */

#ifndef _ARIADNE_RECTANGLE_H
#define _ARIADNE_RECTANGLE_H

#include <iosfwd>

#include "../declarations.h"

#include "../base/array.h"
#include "../base/iterator.h"
#include "../base/tribool.h"
#include "../exceptions.h"

#include "../numeric/arithmetic.h"
#include "../numeric/function.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"

#include "../geometry/point.h"
#include "../geometry/rectangle_expression.h"

namespace Ariadne {
  namespace Geometry {

    template<> 
    inline bool is_a<Rectangle,Rectangle>() { return true; }
    template<> 
    inline bool is_a<Rectangle,Parallelotope>() { return true; }
    template<> 
    inline bool is_a<Rectangle,Zonotope>() { return true; }
    template<> 
    inline bool is_a<Rectangle,Polyhedron>() { return true; }

    template<class R> class RectangleVerticesIterator;
    
    /*! \ingroup BasicSet
     *  \brief A cuboid of arbitrary dimension.
     * 
     * The most important geomtric object in %Ariadne, the %Rectangle class 
     * describes a rectangle, cuboid or hypercuboid in Euclidean space.
     * Rectangles are closed under most major geometric operations, and geometric
     * predicates and operations are simple to compute.
     * Rectangular sets require little data to describe, and are used as cells in grid-based and lattice based sets,
     * making them ubiquitous in storage-representations.
     * Further, rectangles are easily and exactly convertible to the other major polyhedral
     * set representations, namely, Zonotope,  Polytope and Polyhedron.
     * 
     * Rectangles are decribed by the lower and upper bounds in each of their
     * dimensions, as accessed by the lower_bound(dimension_type) const and upper_bound(dimension_type) const methods. 
     * 
     * Rectangles are by default ordered by the lexicographic order on lower-left corner.
     * If these are equal, then points are ordered by the lexicographic order 
     * on the upper-right corner.
     *
     * A rectangle is described by the literal "[a0,b0]x[a1,b1]x...".
     * 
     * Rectangles are always bounded, but may be empty. A zero-dimensional rectangle is considered to be non-empty.
     *
     * \b Storage: A %Rectangle of dimension \a d is specified by \a 2d real 
     * number giving the lower and upper bounds in each coordinate.
     * The lower bound in the \a i th coordinate is the \a 2i th element, and the 
     * upper bound is the \a 2i+1 st element.
     */
    template<class R>
    class Rectangle 
      : public RectangleExpression< Rectangle<R> >
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the rectangle. */
      typedef Point<R> state_type;
      /*! \brief An iterator to the vertices of the rectangle. */
      typedef RectangleVerticesIterator<R> vertices_iterator;
      typedef RectangleVerticesIterator<R> vertices_const_iterator;
     public:
      //@{
      //! \name Constructors
      /*! \brief Construct an empty rectangle with dimension \a d. */
      explicit Rectangle(size_type d=0);
      
      /*! \brief Construct a rectangle from a range of interval values. */
      template<class ForwardIterator> explicit Rectangle(ForwardIterator b, ForwardIterator e);
      
      /*! \brief Construct from an array of intervals. */
      explicit Rectangle(const array< Interval<R> >& a);
      
      /*! \brief Construct from a std::vector of intervals. */
      explicit Rectangle(const std::vector< Interval<R> >& v);

      /*! \brief Construct a degenerate rectangle from a single point. */
      explicit Rectangle(const Point<R>& pt);
      
      /*! \brief Construct a rectangle from an interval point. */
      explicit Rectangle(const Point< Interval<R> >& pt);;
      
      /*! \brief Construct from two corners. */
      explicit Rectangle(const Point<R>& pt1, const Point<R>& pt2);
      
      /*! \brief Construct from a string literal. */
      explicit Rectangle(const std::string& s);
      
      /*! \brief Construct from an interval vector. */
      explicit Rectangle(const LinearAlgebra::Vector< Interval<R> >& iv);
      
      /*! \brief Convert from a rectangle expression. */
      template<class E> Rectangle(const RectangleExpression<E>& r);
    
    
      /*! \brief Copy constructor. */
      Rectangle(const Rectangle<R>& r);
    
      /*! \brief Copy assignment operator. */
      Rectangle<R>& operator=(const Rectangle<R>& r);

      /*! \brief Assign from a rectangle expression. */
      template<class E> Rectangle<R>& operator=(const RectangleExpression<E>& r);
    
      //@}
      
      
      //@{
      //! \name Conversion operators
      /*! \brief Convert to an interval point. */
      operator Point< Interval<R> >() const;
      //@}
      
      
      //{@
      //! \name Comparison operators
      /*! \brief The equality operator */
      bool operator==(const Rectangle<R>& A) const;
      
      /*! \brief The inequality operator */
      bool operator!=(const Rectangle<R>& A) const;
      //@}
      

      //@{
      //! \name Data access
      /*! \brief Returns the projection onto the \a i th coordinate. */
      Interval<R>& operator[] (dimension_type i);
     
      /*! \brief The lower bound of the \a i th coordinate */
      const R& lower_bound(dimension_type i) const;
      
      /*! \brief A reference to the lower bound of the \a i th coordinate */
      R& lower_bound(dimension_type i);
      
      /*! \brief The upper bound of the \a i th coordinate */
      const R& upper_bound(dimension_type i) const;
      
      /*! \brief A reference to the upper bound of the \a i th coordinate */
      R& upper_bound(dimension_type i);
      
      /*! \brief The projection onto the \a i th coordinate. */
      const Interval<R>& operator[] (dimension_type i) const;
      
      /*! \brief The interval of values in the \a i th coordinate. */
      const Interval<R>& interval(dimension_type i) const;
      
      /*! \brief The lower corner. */
      Point<R> lower_corner() const;
      
      /*! \brief The upper corner. */
      Point<R> upper_corner() const;
      
      /*! \brief The set of position vectors of the rectangle. */
      LinearAlgebra::Vector< Interval<R> > position_vectors() const;
      //@}
      
      
      //@{ 
      //! \name Modifying operations
      /*! \brief Makes the rectangle empty. */
      void clear();
      
      /*! \brief Sets the \a i th interval. */
      void set_interval(dimension_type i, Interval<R> x);
      
      /*! \brief Sets the lower bound of the \a i th coordinate to \a r. */
      void set_lower_bound(dimension_type i, const R& l);
      
      /*! \brief Sets the upper bound of the \a i th coordinate to \a u. */
      void set_upper_bound(dimension_type i, const R& u);

      /*! \brief Expand the %Rectangle by \a delta in each direction. */
      Rectangle<R>& expand_by(const R& delta);

      /*! \brief Return a copy of the %Rectangle expanded by \a delta in each direction. */
      Rectangle<R> expand(const R& delta) const;
      //@}
      
      
      //@{
      //! \name Rectangle geometric operations
      /*! \brief The dimension of the Euclidean space the rectangle lies in. */
      size_type dimension() const;
      
      /*! \brief True if the rectangle is empty. A zero-dimensional rectangle is considered empty. */
      tribool empty() const;
      
      /*! \brief The centre. */
      Point<R> centre() const;
      
      /*! \brief The radius in the sup norm. */
      R radius() const;
      
      /*! An approximation to the volume. */
      R volume() const;
      
      /*! \brief Compute a quadrant of the Rectangle determined by \a q.
       *  \a q is a binary word such that the ith bit of q is 0 if the lower half
       *  of the rectangle in the ith coordinate is used, and 1 if the upper
       *  half is used.
       */
      Rectangle<R> quadrant(const Combinatoric::BinaryWord& q) const;
      /*! \brief Subdivide into smaller pieces. */
      ListSet<R,Geometry::Rectangle> subdivide() const;
      
      /*! \brief The vertices of the rectangle. */
      PointList<R> vertices() const;

      /*! \brief An iterator to the first vertex of the rectangle. */
      vertices_const_iterator vertices_begin() const;
      /*! \brief An iterator to the end vertex of the rectangle. */
      vertices_const_iterator vertices_end() const;

      /*! \brief The number of vertices. */
      size_type number_of_vertices() const;
        
      /*! \brief The \a i th vertex. */
      Point<R> vertex(size_type i) const;
        
      /*! \brief Tests if \a point is included into a rectangle. */
      tribool contains(const Point<R>& pt) const;
      
      /*! \brief Checks for boundedness. */
      tribool bounded() const;
      
      /*! \brief A rectangle containing the given rectangle; returns a copy. */
      Rectangle bounding_box() const;
      //@}
      
#ifdef DOXYGEN
      //@{ 
      //! \name Binary geometric predicates
      /*! \brief Set equality operator. */
      friend tribool equal(const Rectangle<R>& A, const Rectangle<R>& B) const;
      /*! \brief Tests disjointness with \a r. */
      friend tribool disjoint(const Rectangle<R>& A, const Rectangle<R>& B) const;
      /*! \brief Tests if the rectangle is a subset of another rectangle \a r. */
      friend tribool subset(const Rectangle<R>& A, const Rectangle<R>& B) const;
      //@}

      //@{ 
      //! \name Binary geometric operations
      /*! \brief The intersection of \a A and \a B. */
      friend Rectangle<R> intersection(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The closure of the intersection of the interiors of \a A and \a B. */
      friend Rectangle<R> regular_intersection(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The smallest rectangle containing \a A and \a B. */
      friend Rectangle<R> rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The componentwise sum of rectangles \a A and \a B. */
      friend Rectangle<R> minkowski_sum(const Rectangle<R>& A, const Rectangle<R>& B); 
      /*! \brief The componentwise difference of rectangles \a A and \a B. */
      friend Rectangle<R> minkowski_difference(const Rectangle<R>& A, const Rectangle<R>& B); 
      
      /*! \brief The difference between two rectangles. */
      friend LinearAlgebra::Vector< Interval<R> > operator-(const Rectangle<R>& A, const Rectangle& B);
      /*! \brief Adds a vector to a rectangle. */
      friend Rectangle<R> operator+(const Rectangle<R>& r, const LinearAlgebra::Vector<R>& v);
      /*! \brief Adds an interval vector to a rectangle. */
      friend Rectangle<R> operator+(const Rectangle<R>& r, const LinearAlgebra::Vector< Interval<R> >& v);
      /*! \brief Subtracts a vector from a rectangle. */
      friend Rectangle<R> operator-(const Rectangle<R>& r, const LinearAlgebra::Vector<R>& v);
      /*! \brief Subtracts an interval vector from a rectangle. */
      friend Rectangle<R> operator-(const Rectangle<R>& r, const LinearAlgebra::Vector< Interval<R> >& v);
      //@}
#endif
      
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
     private:
      array<R> _bounds;
    };
    
    
    
    
    template<class R>
    class Rectangle< Interval<R> > 
      : public RectangleExpression< Rectangle< Interval<R> > >
    {
      typedef Interval<R> I;
     public:
      explicit Rectangle(dimension_type d=0);
      template<class E> Rectangle(const RectangleExpression<E>& e);
      template<class E> Rectangle< Interval<R> >& operator=(const RectangleExpression<E>& e);
      dimension_type dimension() const;
      const Interval<R>& lower_bound(const dimension_type& i) const;
      const Interval<R>& upper_bound(const dimension_type& i) const;
      void set_lower_bound(const dimension_type& i, const Interval<R>& x);
      void set_upper_bound(const dimension_type& i, const Interval<R>& x);
     private:
      template<class RE> void assign(const RE& re);
     private:
      array<I> _bounds;
    };
    
    
    
    
    template<class R>
    class RectangleVerticesIterator 
      : public boost::iterator_facade<RectangleVerticesIterator<R>,
                                      Point<R>,
                                      boost::forward_traversal_tag,
                                      Point<R> const&,
                                      Point<R> const*
                                     >
    {
     public:
      RectangleVerticesIterator(const Rectangle<R>& r, const bool end);
      bool equal(const RectangleVerticesIterator<R>& other) const;
      const Point<R>& dereference() const;
      void increment();
     private:
      const Rectangle<R>* _r; 
      long unsigned int _i; 
      bool _parity; 
      Point<R> _pt;
    };
    
    
    
    
    
    template<class R> tribool equal(const Rectangle<R>& A, const Rectangle<R>& B);
      
    template<class R> tribool disjoint(const Rectangle<R>& A, const Rectangle<R>& B);
      
    template<class R> tribool subset(const Rectangle<R>& A, const Rectangle<R>& B);
    
    template<class R> Rectangle<R> closed_intersection(const Rectangle<R>& A, const Rectangle<R>& B);
      
    template<class R> Rectangle<R> open_intersection(const Rectangle<R>& A, const Rectangle<R>& B);
      
    template<class R> Rectangle<R> rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B);
      
    template<class R1, class R2>
    Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> 
    minkowski_sum(const Rectangle<R1>& A, const Rectangle<R2>& B);

    template<class R1, class R2> 
    Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> 
    minkowski_difference(const Rectangle<R1>& A, const Rectangle<R2>& B);
    
    
    
    template<class R> Rectangle<R> over_approximation(const Rectangle< Numeric::Interval<R> >& ir);
      
    template<class R> Rectangle<R> under_approximation(const Rectangle< Numeric::Interval<R> >& ir);
    
    
    
    template<class R>
    LinearAlgebra::Vector< Interval<R> > 
    operator-(const Geometry::Rectangle<R>& r1, 
              const Geometry::Rectangle<R>& r2);
    
    template<class R, class E>
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::VectorExpression<E>& v);
    
    template<class R>
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v);
    
    template<class R>
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector< Interval<R> >& v);
    
    template<class R>
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v);
    
    template<class R>
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector< Interval<R> >& v);
    
    template<class R>
    Geometry::Rectangle<R> 
    scale(const Geometry::Rectangle<R>& r, const R& scale_factor);
    
    
    
  }
}

#include "rectangle.inline.h"

#endif /* _ARIADNE_RECTANGLE_H */
