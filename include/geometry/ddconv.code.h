/***************************************************************************
 *            ddconv.code.h
 *
 *  Copyright 2006  Alberto Casagrande, Pieter Collins
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
 
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/matrix.code.h"
#include "../base/stlio.h"

namespace Ariadne {
namespace Geometry {

int verbosity=0;
  
  
/*! \param argument is a list of constraints (or generators).
 *  \param result is an output parameter storing a list of generators (or constraints).
 */  
template<class R>
void
ddconv(std::vector< LinearAlgebra::Vector<R> >&  result,
       const std::vector< LinearAlgebra::Vector<R> >&  argument)
{  
  
  // 'dimension' is the dimension of the state space
  dimension_type dimension=argument[0].size();
  
  // The number of lines (as opposed to rays) in 'generators'
  size_type number_of_lines;

  // 'constraints' is a list of vectors $a$ representing the constraint $a^Tx\geq 0$
  const std::vector< LinearAlgebra::Vector<R> >& constraints=argument;

  // 'generators' is a list of vectors $r$ representing the ray $x=\lambda r, \lambda\geq0$
  std::vector< LinearAlgebra::Vector<R> >& generators=result;
  

  // temporary storage
  std::vector< LinearAlgebra::Vector<R> > new_generators;
  
  // Constants
  R zero=0;
  R one=1;
  R minus_one=0;
  
  //std::cout << "C=" << constraints << "\n\n";
  
  LinearAlgebra::Vector<R> v(dimension);
  LinearAlgebra::Vector<R> w(dimension);
  
  // Initialize 'lines' to include lines parallel to coordinate directions.
  number_of_lines=dimension;
  for(dimension_type i=0; i!=dimension; ++i) {
    v=LinearAlgebra::Vector<R>::zero(dimension);
    v(i)=one;
    generators.push_back(v);
  }
  
  // Add constraints sequentially and compute new generators
  for(size_type k=0; k!=constraints.size(); ++k) {
    if(verbosity > 1) { std::cerr << "G=" << generators << "\n"; }
    if(verbosity > 1) { std::cerr << "c=" << constraints[k] << "\n"; }
    
    // If any line does not saturate the new constraint, add a ray containing 
    // the first such line to the generators. If more lines do not saturate the
    // constraint, add a multiple of the first so that the rest do.
    
    // Find non-saturating line to use a pivot
    size_type imax=number_of_lines;
    R max=zero;
    for(size_type i=0; i!=number_of_lines; ++i) {
      R dot=abs(inner_product(constraints[k],generators[i]));
      if(dot>max) { 
        max=dot;
        imax=i;
      }
    }
    
    if(imax!=number_of_lines) { // Test if non-saturating line found
      // Reduce number of lines
      number_of_lines-=1;

      // Swap pivot line
      if(imax!=number_of_lines) {
        std::swap(generators[imax],generators[number_of_lines]);
      }
      
      // Alias pivot element
      LinearAlgebra::Vector<R>& p=generators[number_of_lines];

      // Ensure 'p' satisfies new constraint
      if(inner_product(constraints[k],p)<zero) {
        p=-p;
      }
      
      // Add linear multiple of 'p' to remaining lines to ensure they saturate new constraint
      R pdot=max;
      for(size_type i=0; i!=number_of_lines; ++i) {
        R dot=inner_product(constraints[k],generators[i]);
        generators[i]-=(dot/pdot)*p;
      }
     
      // Add linear multiple of 'p' to generating rays to ensure they saturate new constraint
      for(size_type i=number_of_lines+1; i!=generators.size(); ++i) {
        R dot=inner_product(constraints[k],generators[i]);
        generators[i]-=(dot/pdot)*p;
      }
     
      //std::cout << "G=" << generators << "\n"; 
    
    } else { // No pivot line found, so deal with new constraint normally
  
      LinearAlgebra::Matrix<int> saturation_matrix(k+1,generators.size());
      for(size_type i=0; i!=k+1; ++i) {
        for(size_type j=0; j!=generators.size(); ++j) {
          R dot=inner_product(constraints[i],generators[j]);
          if(dot<zero) {
            saturation_matrix(i,j)=-1;
          } else if(dot>zero) {
            saturation_matrix(i,j)=+1;
          } else {
            saturation_matrix(i,j)=0;
          }
        }
      }
      //std::cout << "S=" << saturation_matrix  << "\n\n";
      //std::vector< LinearAlgebra::Vector<R> > generators;
    
      size_type number_of_generators=generators.size();
      for(size_type j=number_of_lines; j!=number_of_generators; ++j) {
        if(saturation_matrix(k,j)==-1) {
          bool found_adjacent=false;
          v=generators[j];
          for(size_type js=number_of_lines; js!=number_of_generators; ++js) {
            if(saturation_matrix(k,js)==1) {
              int asum=0;
              for(size_type i=0; i!=k; ++i) {
                if(saturation_matrix(i,j)!=saturation_matrix(i,js)) {
                  ++asum;
                }
              }
              if(asum<=2) {
                R dot=inner_product(constraints[k],v);
                R dots=inner_product(constraints[k],generators[js]);
                w=v-(dot/dots)*generators[js];
                if(found_adjacent) {
                  generators.push_back(w);
                } else {
                  found_adjacent=true;
                  generators[j]=w;
                }
              }
            }
          }
        }
      }
    }
  }
  
  // Normalize generators
  for(size_type i=0; i!=generators.size(); ++i) {
    if(generators[i](dimension-1)!=zero) {
      generators[i]/=abs(generators[i](dimension-1));
    }
  }
  
  //std::cout << "Result:\n  C=" << constraints << "\n  G=" << generators << std::endl;
  LinearAlgebra::Matrix<int> saturation_matrix(constraints.size(),generators.size());
  for(size_type i=0; i!=constraints.size(); ++i) {
    for(size_type j=0; j!=generators.size(); ++j) {
      R dot=inner_product(constraints[i],generators[j]);
      if(dot<zero) {
        saturation_matrix(i,j)=-1;
      } else if(dot>zero) {
        saturation_matrix(i,j)=+1;
      } else {
        saturation_matrix(i,j)=0;
      }
    }
  }
  //std::cout << "  S=" << saturation_matrix  << "\n\n";
    
}

}
}
