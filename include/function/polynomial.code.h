/***************************************************************************
 *            polynomial.code.h
 *
 *  Copyright  2007  Alberto Casagrande
 *  pieter.collins@cwi.nl
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

#include "polynomial.h"
#include "exceptions.h"

#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../base/stlio.h"
#include "../base/array.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../function/position_index.h"
#include "../function/sorted_index.h"
#include "../function/multi_index.h"

#include "../output/texstream.h"

namespace Ariadne {

template<class R>
Function::Polynomial<R>::Polynomial() 
  : _result_size(0), 
    _argument_size(0), 
    _degree(0),
    _data() 
{
}


template<class R>
Function::Polynomial<R>::Polynomial(const std::string& s) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Function::Polynomial<R>::Polynomial(const size_type& rs, const size_type& as, const size_type& d)
  : _result_size(rs), 
    _argument_size(as), 
    _degree(d),
    _data(rs*Numeric::choose(d+as,as),static_cast<R>(0)) 
{
}


template<class R>
void
Function::Polynomial<R>::resize(const size_type& rs, const size_type& as, const size_type& d)
{
  this->_result_size=rs; 
  this->_argument_size=as; 
  this->_degree=d;
  this->_data.resize(rs*Numeric::choose(d+as,as));
}


template<class R>
Function::Polynomial<R>
Function::Polynomial<R>::zero(const size_type& rs, const size_type& as)  
{
  return Polynomial<R>(rs,as,0);
}


template<class R>
Function::Polynomial<R>
Function::Polynomial<R>::one(const size_type& as)  
{
  R o=1;
  return Polynomial<R>(1u,as,0,&o);
}


template<class R>
Function::Polynomial<R>
Function::Polynomial<R>::constant(const size_type& as, const R& c)  
{
  return Polynomial<R>(1u,as,0,&c);
}


template<class R>
bool
Function::Polynomial<R>::operator==(const Polynomial<R>& p2) const
{
  const Polynomial<R>& p1=*this;
  if(p1._result_size!=p2._result_size || p1._argument_size!=p2._argument_size) {
    return false;
  }
  if(p1._degree==p2._degree) {
    return p1._data==p2._data;
  } else {
    size_type smin=std::min(p1._data.size(),p2._data.size());
    size_type smax=std::max(p1._data.size(),p2._data.size());
    for(size_type i=0; i!=smin; ++i) {
      if(p1._data[i]!=p2._data[i]) {
        return false;
      }
    } 
    if(p1._data.size()>p2._data.size()) {
      for(size_type i=smin; i!=smax; ++i) {
        if(p1._data[i]!=0) {
          return false;
        }
      } 
    } else {
      for(size_type i=smin; i!=smax; ++i) {
        if(p2._data[i]!=0) {
          return false;
        }
      } 
    }
    return true;
  }
}


template<class R>
bool
Function::Polynomial<R>::operator!=(const Polynomial<R>& p2) const
{
  return !(*this==p2);
}


template<class R>
size_type 
Function::Polynomial<R>::argument_size() const
{ 
  return this->_argument_size; 
}

template<class R>
size_type 
Function::Polynomial<R>::result_size() const 
{ 
  return this->_result_size;
}

template<class R>
size_type 
Function::Polynomial<R>::degree() const 
{
  return this->_degree; 
}
      
template<class R>
size_type 
Function::Polynomial<R>::smoothness() const 
{ 
  return (size_type) -1; 
}
      
template<class R>
const array<R>&
Function::Polynomial<R>::data() const 
{ 
  return this->_data;
}
      

template<class R>
void
Function::Polynomial<R>::set(const size_type& i, const MultiIndex& j, const R& x) 
{ 
  assert(i<this->result_size());
  assert(j.degree()<=this->degree());
  this->_data[this->result_size()*j.position()+i]=x;
}
      
template<class R>
R&
Function::Polynomial<R>::at(const size_type& i, const MultiIndex& j)
{ 
  assert(i<this->result_size());
  assert(j.degree()<=this->degree());
  return this->_data[this->result_size()*j.position()+i];
}
      
template<class R>
const R&
Function::Polynomial<R>::get(const size_type& i, const MultiIndex& j) const 
{ 
  assert(i<this->result_size());
  assert(j.degree()<=this->degree());
  return this->_data[this->result_size()*j.position()+i];
}
      


template<class R>
Function::Polynomial<R> 
Function::Polynomial<R>::component(const size_type& i) const
{
  Polynomial<R> result(1u,this->argument_size(),this->degree());
  R* rptr=result._data.begin();
  R* eptr=result._data.end();
  const R* aptr=this->_data.begin()+i;
  const size_type& inc=this->result_size();
  while(rptr!=eptr) {
    *rptr=*aptr;
    rptr+=1;
    aptr+=inc;
  }
  return result;
}



template<class R>
Function::Polynomial<typename Function::Polynomial<R>::I> 
Function::Polynomial<R>::truncate(const size_type& degree, const size_type& smoothness, const Geometry::Rectangle<R>& domain) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
LinearAlgebra::Vector<typename Function::Polynomial<R>::F> 
Function::Polynomial<R>::evaluate(const LinearAlgebra::Vector<F>& x) const
{
  if(this->argument_size()!=x.size()) {
    ARIADNE_THROW(IncompatibleSizes,"Polynomial::evaluate(Vector)","Incompatible argument size");
  }

  // TODO: Make this more efficient
  LinearAlgebra::Vector<F> result(this->result_size());
  for(MultiIndex j(this->argument_size()); j.degree()<=this->degree(); ++j) {
    F xa=1;
    for(size_type k=0; k!=j.number_of_variables(); ++k) {
      xa*=Numeric::pow(x[k],j[k]);
    }
    for(size_type i=0; i!=this->result_size(); ++i) {
      result[i]+=this->get(i,j)*xa;
    }
  }
  return result;
}



template<class R0, class R1, class R2>
void
Function::add(Polynomial<R0>& p0, const Polynomial<R1>& p1, const Polynomial<R2>& p2)
{
  if(p1.result_size()!=p2.result_size()) {
    ARIADNE_THROW(IncompatibleSizes,"add(Polynomial,Polynomial)","Incompatible result sizes");
  }
  if(p1.argument_size()!=p2.argument_size()) {
    ARIADNE_THROW(IncompatibleSizes,"add(Polynomial,Polynomial)","Incompatible argument sizes");
  }
  size_type rs=p1.result_size();
  size_type as=p1.argument_size();
  size_type d1=p1.degree();
  size_type d2=p2.degree();
  size_type d=std::max(d1,d2);

  p0.resize(rs,as,d);
  size_type kmin=std::min(p1._data.size(),p2._data.size());
  size_type kmax=p0._data.size();
  for(size_type i=0; i!=kmin; ++i) {
    p0._data[i]=p1._data[i]+p2._data[i];
  }
  if(d1>d2) {
    for(size_type i=kmin; i!=kmax; ++i) {
      p0._data[i]=p1._data[i];
    }
  } else {
    for(size_type i=kmin; i!=kmax; ++i) {
      p0._data[i]=p2._data[i];
    }
  }
}


template<class R0, class R1, class R2>
void
Function::sub(Polynomial<R0>& p0, const Polynomial<R1>& p1, const Polynomial<R2>& p2)
{
  if(p1.result_size()!=p2.result_size()) {
    ARIADNE_THROW(IncompatibleSizes,"add(Polynomial,Polynomial)","Incompatible result sizes");
  }
  if(p1.argument_size()!=p2.argument_size()) {
    ARIADNE_THROW(IncompatibleSizes,"add(Polynomial,Polynomial)","Incompatible argument sizes");
  }
  size_type rs=p1.result_size();
  size_type as=p1.argument_size();
  size_type d1=p1.degree();
  size_type d2=p2.degree();
  size_type d=std::max(d1,d2);

  p0.resize(rs,as,d);
  size_type kmin=std::min(p1._data.size(),p2._data.size());
  size_type kmax=p0._data.size();
  for(size_type i=0; i!=kmin; ++i) {
    p0._data[i]=p1._data[i]-p2._data[i];
  }
  if(d1>d2) {
    for(size_type i=kmin; i!=kmax; ++i) {
      p0._data[i]=p1._data[i];
    }
  } else {
    for(size_type i=kmin; i!=kmax; ++i) {
      p0._data[i]=-p2._data[i];
    }
  }
}


template<class R0, class R1, class R2>
void
Function::mul(Polynomial<R0>& p0, const Polynomial<R1>& p1, const Polynomial<R2>& p2)
{
  if(p1.result_size()!=1u) {
    ARIADNE_THROW(IncompatibleSizes,"mul(Polynomial,Polynomial)","p1.result_size()="<<p1.result_size());
  }
  if(p2.result_size()!=1u) {
    ARIADNE_THROW(IncompatibleSizes,"mul(Polynomial,Polynomial)","p2.result_size()="<<p2.result_size());
  }
  if(p1.argument_size()!=p2.argument_size()) {
    ARIADNE_THROW(IncompatibleSizes,"add(Polynomial p1,Polynomial p2)","p1.result_size()="<<p1.result_size()<<", p2.result_size()="<<p2.result_size());
  }

  size_type rs=1u;
  size_type as=p1.argument_size();
  size_type d1=p1.degree();
  size_type d2=p2.degree();
  size_type d=d1+d2;

  p0.resize(rs,as,d);
  MultiIndex j0(as);
  for(MultiIndex j1(as); j1.degree()<=d1; ++j1) {
    const R1& x1=p1.get(0u,j1);
    for(MultiIndex j2(as); j2.degree()<=d2; ++j2) {
      const R2& x2=p2.get(0u,j2);
      j0=j1+j2;
      p0.at(0u,j0)+=x1*x2;
    }
  }
}


template<class R0, class R1, class R2>
void
Function::div(Polynomial<R0>& p0, const Polynomial<R1>& p1, const Polynomial<R2>& p2)
{
  Polynomial<R0> p3;
  Function::recip(p3,p2);
  Function::mul(p0,p1,p3);
}


template<class R0,class R1>
void
Function::pow(Polynomial<R0>& p0, const Polynomial<R1>& p1, const unsigned int& n)
{
  assert(p1.result_size()==1);

  if(n==1) {
    p0=p1;
    return;
  }

  R0 one=1;
  p0=Polynomial<R0>::one(p1.argument_size());
  if(n==0) {
    return;
  }

  Polynomial<R0> tmp(p1);
  for(uint i=1; i<=n; i*=2) {
    if(i&n) {
      p0=tmp*p0;
    }
    tmp=tmp*tmp;
  }
}



template<class R0,class R1>
void
Function::scale(Polynomial<R0>& p0, const R1& x1)
{
  for(size_type i=0; i!=p0._data.size(); ++i) {
    p0._data[i]*=x1;
  }
}

template<class R0,class R1>
void
Function::recip(Polynomial<R0>& p0, const Polynomial<R1>& p1)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R0, class R1, class R2>
void
Function::compose(Polynomial<R0>& p0, const Polynomial<R1>& p1, const Polynomial<R2>& p2)
{
  // TODO: Improve this algorithm as it's critical!!
  if(p1.argument_size()!=p2.result_size()) {
    ARIADNE_THROW(IncompatibleSizes,"compose(Polynomial p1,Polynomial p2)","p1.argument_size()="<<p1.argument_size()<<", p2.result_size()="<<p2.result_size());
  }

  if(p1.degree()==0) { 
    p0=static_cast< Polynomial<R0> >(p1); 
    return;
  }

  p0.resize(p1.result_size(),p2.argument_size(),p1.degree()*p2.degree());
  
  Polynomial<R0>* all_powers=new Polynomial<R0>[p2.result_size()*(p1.degree()+1)];
  Polynomial<R0>* powers[p2.result_size()];
  for(size_type i=0; i!=p2.result_size(); ++i) {
    powers[i]=all_powers+i*(p1.degree()+1);
  }
  for(size_type i=0; i!=p2.result_size(); ++i) {
    powers[i][0]=Polynomial<R0>::one(p2.argument_size());
    powers[i][1]=p2.component(i);
    if(p1.degree()>=2) {
      powers[i][2]=Function::pow(powers[i][1],2);
    }
    for(size_type j=3; j<=p1.degree(); ++j) {
      powers[i][j]=powers[i][2]*powers[i][j-2];
    }
  }
  
  Polynomial<R0>* results=new Polynomial<R0>[p1.result_size()];
  for(size_type i=0; i!=p1.result_size(); ++i) {
    results[i]=Polynomial<R0>::zero(1u,p2.argument_size());
  }

  for(size_type i=0; i!=p1.result_size(); ++i) {
    for(MultiIndex j(p1.argument_size()); j.degree()<=p1.degree(); ++j) {
      Polynomial<R0> t=Polynomial<R0>::constant(p2.argument_size(),p1.get(i,j));
      for(size_type k=0; k!=p1.argument_size(); ++k) {
        t=t*powers[k][j[k]];
      }
      results[i]=results[i]+t;
    }
  }
  
  for(size_type i=0; i!=p0.result_size(); ++i) {
    for(size_type j=0; j!=p0.data().size()/p0.result_size(); ++j) {
      p0._data[i+j*p0.result_size()]=results[i].data()[j];
    }
  }
  
  delete[] results;
  delete[] all_powers;
}


template<class R>
LinearAlgebra::Matrix<typename Function::Polynomial<R>::F> 
Function::Polynomial<R>::jacobian(const LinearAlgebra::Vector<F>& s) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
std::ostream&
Function::Polynomial<R>::write(std::ostream& os) const 
{
  os << "Polynomial(\n";
  for(uint i=0; i!=this->result_size(); ++i) {
    os << "  " << i << ": " << std::flush;
    bool first=true;
    for(MultiIndex j(this->argument_size()); j.degree()<=this->degree(); ++j) {
      if(first) { first=false; } else { os << ", "; } 
      os << j << ":" << j.position() << ":"  << std::flush; os << this->get(i,j) << std::flush;
    }
    os << ";\n";
  }
  os << ")";
  return os;
}


template<class R>
std::ostream&
Function::operator<<(std::ostream& os, const Polynomial<R>& p)
{
  return p.write(os);
}


template<class R>
Output::texstream&
Output::operator<<(Output::texstream& texs, const Function::Polynomial<R>& p)
{
  using namespace Function;
  texs << "%Polynomial\n";
  texs << "\\ensuremath{\n";
  texs << "\\left( \\begin{array}{c}\n";
  char var='x';
  for(size_type i=0; i!=p.result_size(); ++i) {
    bool first = true;
    if(i!=0) { texs << "\\\\"; }
    for(MultiIndex j(p.argument_size()); j.degree()<=p.degree(); ++j) {
      const R& a=p.get(i,j);
      if(a!=0) {
        if(first) { first=false; }
        else { if(a>0) { texs << '+'; } }
        if(a==1) { if(j.degree()==0) { texs << a; } }
        else if(a==-1) { if(j.degree()==0) { texs << a; } else { texs << '-'; } }
        else { texs << a << ' '; }
        for(size_type k=0; k!=p.argument_size(); ++k) {
          if(j[k]!=0) {
            texs << var << "_{ " << k+1 << "}";
            if(j[k]!=1) {
              texs << "^{" << j[k] << "}";
            }
          }
          texs << " ";
        }
      } 
    }
    texs << "\n"; 
  }
  texs << "\\end{array}\\right)\n}\n";
  return texs;
}


template<class R>
void
Function::Polynomial<R>::instantiate()
{
  typedef typename Numeric::traits<R>::arithmetic_type I;
  R* x=0;
  Polynomial<R>* p=0;
  std::ostream* os = 0;
  Output::texstream* texs = 0;

  Function::operator+(*p,*p);
  Function::operator-(*p,*p);
  Function::operator*(*p,*p);
  Function::operator/(*p,*p);
  Function::operator*(*p,*x);
  Function::operator/(*p,*x);
  Function::pow(*p,0u);
  Function::compose(*p,*p);
  *os << *p;
  *texs << *p;
}




}
