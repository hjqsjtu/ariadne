/***************************************************************************
 *            taylor_function.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include "functional.h"
#include "config.h"

#include <iostream>
#include <iomanip>

#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "polynomial.h"
#include "differential.h"
#include "taylor_model.h"

#include "function.h"
#include "taylor_function.h"

#include "../src/function_mixin.tcc"
#include "taylor_function.h"

#define VOLATILE ;

namespace Ariadne {

typedef unsigned int uint;

static double TAYLOR_FUNCTION_WRITING_ACCURACY = 1e-8;

void _set_scaling(ScalarTaylorFunction& x, const Interval& ivl, uint j)
{
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(upward);
    const Float& l=ivl.lower();
    const Float& u=ivl.upper();
    VOLATILE Float pc=u; pc+=l;
    VOLATILE Float nc=-u; nc-=l;
    VOLATILE Float pg=u; pg-=l;
    VOLATILE Float ng=l; ng-=u;
    x.error()=(pc+nc+pg+ng)/4;
    set_rounding_mode(to_nearest);
    MultiIndex a(x.argument_size());
    x.expansion().append(a,(l+u)/2);
    ++a[j];
    x.expansion().append(a,(l+u)/2);
    set_rounding_mode(rounding_mode);
}



ScalarFunctionModel<ValidatedTag>& ScalarFunctionModel<ValidatedTag>::operator=(const ScalarTaylorFunction& f) {
    this->_ptr=clone_on_copy_ptr< ScalarFunctionModelInterface<ValidatedTag> >(new ScalarTaylorFunction(f)); return *this;
}


ScalarTaylorFunction::ScalarTaylorFunction()
    : _domain(), _model()
{ }

ScalarTaylorFunction::ScalarTaylorFunction(const Box& d, Sweeper swp)
    : _domain(d), _model(d.size(),swp)
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const Box& d, const Expansion<Float>& p, const Float& e, const Sweeper& swp)
    : _domain(d), _model(p,e,swp)
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const Box& d, const ValidatedTaylorModel& m)
    : _domain(d), _model(m)
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const Box& d, const ValidatedScalarFunction& f, Sweeper swp)
    : _domain(d), _model(f.argument_size(),swp)
{
    ARIADNE_ASSERT_MSG(d.size()==f.argument_size(),"d="<<d<<" f="<<f);
    Vector<ValidatedTaylorModel> x=ValidatedTaylorModel::scalings(d,swp);
    this->_model=f.evaluate(x);
    this->_model.sweep();
}

ScalarTaylorFunction::ScalarTaylorFunction(const ValidatedScalarFunctionModel& f) {
     ARIADNE_ASSERT_MSG(dynamic_cast<const ScalarTaylorFunction*>(f._ptr.operator->())," f="<<f);
     *this=dynamic_cast<const ScalarTaylorFunction&>(*f._ptr);
}

ScalarTaylorFunction& ScalarTaylorFunction::operator=(const ValidatedScalarFunctionModel& f)
{
    return (*this)=ScalarTaylorFunction(this->domain(),f,this->sweeper());
}



ScalarTaylorFunction ScalarTaylorFunction::zero(const Box& d, Sweeper swp)
{
    return ScalarTaylorFunction(d,ValidatedTaylorModel::zero(d.size(),swp));
}

ScalarTaylorFunction ScalarTaylorFunction::constant(const Box& d, const ValidatedNumberType& c, Sweeper swp)
{
    return ScalarTaylorFunction(d,ValidatedTaylorModel::constant(d.size(),c,swp));
}

ScalarTaylorFunction ScalarTaylorFunction::identity(const Interval& ivl, Sweeper swp)
{
    Vector<Interval> d(1,ivl);
    return ScalarTaylorFunction(d,ValidatedTaylorModel::scaling(d.size(),0u,d[0u],swp));
}

ScalarTaylorFunction ScalarTaylorFunction::coordinate(const Box& d, uint j, Sweeper swp)
{
    ARIADNE_ASSERT(j<d.size());
    return ScalarTaylorFunction(d,ValidatedTaylorModel::scaling(d.size(),j,d[j],swp));
}

ScalarTaylorFunction ScalarTaylorFunction::variable(const Box& d, uint j, Sweeper swp)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::variable","Use ScalarTaylorFunction::coordinate instead");
    ARIADNE_ASSERT(j<d.size());
    return ScalarTaylorFunction(d,ValidatedTaylorModel::scaling(d.size(),j,d[j],swp));
}


Vector<ScalarTaylorFunction> ScalarTaylorFunction::constants(const Box& d, const Vector<ValidatedNumberType>& c, Sweeper swp)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::constants","Use VectorTaylorFunction::constant instead");
    Vector<ScalarTaylorFunction> x(c.size(),ScalarTaylorFunction(d,swp));
    for(uint i=0; i!=c.size(); ++i) {
        x[i]=c[i];
    }
    return x;
}

Vector<ScalarTaylorFunction> ScalarTaylorFunction::variables(const Box& d, Sweeper swp)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::variables","Use VectorTaylorFunction::identity instead");
    return variables(d,0u,d.size(),swp);
}

Vector<ScalarTaylorFunction> ScalarTaylorFunction::variables(const Box& d, uint imin, uint imax, Sweeper swp)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::variables","Use VectorTaylorFunction::projection instead");
    ARIADNE_ASSERT(imin<=imax);
    ARIADNE_ASSERT(imax<=d.size());

    Vector<ScalarTaylorFunction> x(imax-imin);
    for(uint i=imin; i!=imax; ++i) {
        x[i-imin]=ScalarTaylorFunction::variable(d,i,swp);
    }
    return x;
}


ScalarTaylorFunction* ScalarTaylorFunction::_clone() const
{
    return new ScalarTaylorFunction(*this);
}

ScalarTaylorFunction* ScalarTaylorFunction::_create() const
{
    return new ScalarTaylorFunction(this->domain(),this->_model.sweeper());
}

VectorFunctionModelInterface<ValidatedTag>* ScalarTaylorFunction::_create_vector(Nat i) const
{
    return new VectorTaylorFunction(i,this->domain(),this->_model.sweeper());
}



Sweeper ScalarTaylorFunction::sweeper() const
{
    return this->_model.sweeper();
}


Polynomial<ValidatedNumberType> polynomial(const ValidatedTaylorModel& tm);

Polynomial<ValidatedNumberType>
ScalarTaylorFunction::polynomial() const
{
    Polynomial<ValidatedNumberType> p(this->argument_size());
    p=Ariadne::polynomial(this->model());

    Vector<Polynomial<ValidatedNumberType> > s(this->argument_size());
    for(uint j=0; j!=this->argument_size(); ++j) {
        if(this->domain()[j].radius()<=0) {
            ARIADNE_ASSERT(this->domain()[j].radius()==0);
            s[j]=Polynomial<ValidatedNumberType>::constant(this->argument_size(),this->domain()[j]);
        } else {
            s[j]=Ariadne::polynomial(ValidatedTaylorModel::unscaling(this->argument_size(),j,this->domain()[j],this->sweeper()));
        }
    }

    return compose(p,s);
}

ValidatedScalarFunction
ScalarTaylorFunction::function() const
{
    return new ScalarTaylorFunction(*this);
}


bool ScalarTaylorFunction::operator==(const ScalarTaylorFunction& tv) const
{
    return this->_domain==tv._domain && this->_model==tv._model;
}



ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const ScalarTaylorFunction& y) {
    ARIADNE_ASSERT(subset(x.domain(),y.domain()));
    if(x.domain()==y.domain()) { x._model+=y._model; }
    else { x._model+=restrict(y,x.domain())._model; }
    return x;
}

ScalarTaylorFunction& operator-=(ScalarTaylorFunction& x, const ScalarTaylorFunction& y) {
    ARIADNE_ASSERT(subset(x.domain(),y.domain()));
    if(x.domain()==y.domain()) { x._model-=y._model; }
    else { x._model-=restrict(y,x.domain())._model; }
    return x;
}


ScalarTaylorFunction& operator+=(ScalarTaylorFunction& f, const ValidatedNumberType& c) {
    f._model+=c;
    return f;
}

ScalarTaylorFunction& operator-=(ScalarTaylorFunction& f, const ValidatedNumberType& c) {
    f._model-=c;
    return f;
}

ScalarTaylorFunction& operator*=(ScalarTaylorFunction& f, const ValidatedNumberType& c) {
    f._model*=c;
    return f;
}

ScalarTaylorFunction& operator/=(ScalarTaylorFunction& f, const ValidatedNumberType& c) {
    f._model/=c;
    return f;
}

ScalarTaylorFunction operator+(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,x._model);
}

ScalarTaylorFunction operator-(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,-x._model);
}


ScalarTaylorFunction operator+(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model+x2._model); }
    else {
        Box domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model+restrict(x2,domain)._model);}
}

ScalarTaylorFunction operator-(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model-x2._model); }
    else {
        Box domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model-restrict(x2,domain)._model);}
}

ScalarTaylorFunction operator*(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model*x2._model); }
    else {
        Box domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model*restrict(x2,domain)._model);}
}

ScalarTaylorFunction operator/(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model/x2._model); }
    else {
        Box domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model/restrict(x2,domain)._model);}
}





ScalarTaylorFunction max(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,max(x1._model,x2._model)); }
    else {
        Box domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,max(restrict(x1,domain)._model,restrict(x2,domain)._model));}
}

ScalarTaylorFunction min(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,min(x1._model,x2._model)); }
    else {
        Box domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,min(restrict(x1,domain)._model,restrict(x2,domain)._model));}
}



ScalarTaylorFunction* ScalarTaylorFunction::_derivative(uint j) const
{
    return new ScalarTaylorFunction(Ariadne::derivative(*this,j));
}


Float
ScalarTaylorFunction::evaluate(const Vector<Float>& x) const
{
    const ScalarTaylorFunction& f=*this;
    if(!contains(f.domain(),x)) {
        ARIADNE_THROW(DomainException,"f.evaluate(x) with f="<<f<<", x="<<x,"x is not an element of f.domain()="<<f.domain());
    }
    Vector<Float> sx=Ariadne::unscale(x,f._domain);
    return Ariadne::evaluate(this->_model.expansion(),sx);
}

ValidatedNumberType
ScalarTaylorFunction::evaluate(const Vector<ValidatedNumberType>& x) const
{
    return Ariadne::evaluate(*this,x);
}

ValidatedNumberType
ScalarTaylorFunction::operator()(const Vector<ValidatedNumberType>& x) const
{
    return Ariadne::evaluate(*this,x);
}


void ScalarTaylorFunction::restrict(const Box& d) {
    ARIADNE_ASSERT(subset(d,this->domain()));
    const Box& od=this->domain();
    for(uint j=0; j!=d.size(); ++j) {
        if(od[j]!=d[j]) { *this=Ariadne::restrict(*this,j,d[j]); }
    }
}

ScalarTaylorFunction restrict(const ScalarTaylorFunction& tv, const Box& d) {
    ARIADNE_ASSERT(subset(d,tv.domain()));
    const Box& od=tv.domain();
    ScalarTaylorFunction r=tv;
    for(uint j=0; j!=d.size(); ++j) {
        if(od[j]!=d[j]) { r=restrict(r,j,d[j]); }
    }
    return r;
}

ScalarTaylorFunction extend(const ScalarTaylorFunction& tv, const Box& d) {
    const Box& domain=tv.domain();
    ARIADNE_ASSERT(subset(domain,d));
    for(uint i=0; i!=d.size(); ++i) {
        ARIADNE_ASSERT(domain[i]==d[i] || domain[i].lower()==domain[i].upper());
    }
    return ScalarTaylorFunction(d,tv._model);
}

ValidatedNumberType
evaluate(const ScalarTaylorFunction& f, const Vector<ValidatedNumberType>& x) {
    if(!subset(x,f.domain())) {
        ARIADNE_THROW(DomainException,std::setprecision(17)<<"evaluate(f,x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
    }
    return unchecked_evaluate(f,x);
}

ValidatedNumberType
unchecked_evaluate(const ScalarTaylorFunction& f, const Vector<ValidatedNumberType>& x) {
    return evaluate(f.model(),unscale(x,f.domain()));
}


ScalarTaylorFunction
compose(const EffectiveScalarFunction& g, const VectorTaylorFunction& f)
{
    return ScalarTaylorFunction(f.domain(),g.evaluate(f.models()));
}

ScalarTaylorFunction
compose(const ValidatedScalarFunction& g, const VectorTaylorFunction& f)
{
    return ScalarTaylorFunction(f.domain(),g.evaluate(f.models()));
}

ScalarTaylorFunction
compose(const ValidatedScalarFunctionInterface& g, const VectorTaylorFunction& f)
{
    return ScalarTaylorFunction(f.domain(),g.evaluate(f.models()));
}

ScalarTaylorFunction
compose(const ScalarTaylorFunction& g, const VectorTaylorFunction& f)
{
    if(!subset(f.codomain(),g.domain())) {
        ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
    }
    return unchecked_compose(g,f);
}

ScalarTaylorFunction
unchecked_compose(const ScalarTaylorFunction& g, const VectorTaylorFunction& f)
{
    return ScalarTaylorFunction(f.domain(),compose(g.model(),unscale(f.models(),g.domain())));
}



ScalarTaylorFunction
partial_evaluate(const ScalarTaylorFunction& tf, uint k, const ExactNumberType& c)
{
    return partial_evaluate(tf,k,ValidatedNumberType(c));
}

ScalarTaylorFunction
partial_evaluate(const ScalarTaylorFunction& te, uint k, const ValidatedNumberType& c)
{
    // Scale c to domain
    const uint as=te.argument_size();
    ARIADNE_ASSERT(k<as);
    const Box& domain=te.domain();
    const Interval& dk=domain[k];
    ValidatedNumberType sc=(c-med_ivl(dk))/rad_ivl(dk);

    Box new_domain(as-1);
    for(uint i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
    for(uint i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

    ValidatedTaylorModel new_model=partial_evaluate(te.model(),k,sc);

    return ScalarTaylorFunction(new_domain,new_model);
}


// To scale from a model on [a,b] to a model on [c,d], use scale factor s=(d-c)/(b-a)
// and translation t=((c+d)-(a+b))/(b-a)
// Because we are scaling the model on [-1,+1], this is not the same as
// the mapping taking [a,b] to [c,d]
ScalarTaylorFunction partial_restrict(const ScalarTaylorFunction& tv, uint k, const Interval& new_ivl) {
    ARIADNE_ASSERT(k<tv.argument_size())
    const Interval& old_ivl=tv.domain()[k];
    ARIADNE_ASSERT(subset(new_ivl,old_ivl));
    if(new_ivl==old_ivl) { return tv; }
    Float a=old_ivl.lower(); Float b=old_ivl.upper();
    Float c=new_ivl.lower(); Float d=new_ivl.upper();
    if(a==b) { ARIADNE_ASSERT( a<b || (a==b && c==d) ); return tv; }
    ValidatedNumberType s=sub_ivl(d,c)/sub_ivl(b,a);
    // ValidatedNumberType t=(mul_ivl(b,c)-mul_ivl(a,d))/sub_ivl(b,a);  // WRONG!!
    ValidatedNumberType t=(add_ivl(c,d)-add_ivl(a,b))/sub_ivl(b,a);
    Box new_dom=tv.domain();
    new_dom[k]=new_ivl;
    return ScalarTaylorFunction(new_dom,preaffine(tv.model(),k,s,t));
}


ScalarTaylorFunction restrict(const ScalarTaylorFunction& tv, uint k, const Interval& new_ivl) {
    return partial_restrict(tv,k,new_ivl);
}

ScalarTaylorFunction antiderivative(const ScalarTaylorFunction& f, uint k, Float c)
{
    ARIADNE_ASSERT(k<f.argument_size());
    ARIADNE_ASSERT(contains(f.domain()[k],c));

    ScalarTaylorFunction g = antiderivative(f,k);
    VectorTaylorFunction h = VectorTaylorFunction::identity(g.domain(),g.sweeper());
    h[k] = ScalarTaylorFunction::constant(g.domain(),c,g.sweeper());

    return g-compose(g,h);
}





pair<ScalarTaylorFunction,ScalarTaylorFunction>
split(const ScalarTaylorFunction& tv, uint j)
{
    pair<ValidatedTaylorModel,ValidatedTaylorModel> models=split(tv.model(),j);
    pair<Box,Box> subdomains=split(tv.domain(),j);
    return make_pair(ScalarTaylorFunction(subdomains.first,models.first),
                     ScalarTaylorFunction(subdomains.second,models.second));

}

bool refines(const ScalarTaylorFunction& tv1, const ScalarTaylorFunction& tv2)
{
    if(tv1.domain()==tv2.domain()) { return refines(tv1.model(),tv2.model()); }
    if(subset(tv2.domain(),tv1.domain())) { return refines(restrict(tv1,tv2.domain()).model(),tv2.model()); }
    else { return false; }
}

bool disjoint(const ScalarTaylorFunction& tv1, const ScalarTaylorFunction& tv2)
{
    if(tv1.domain()==tv2.domain()) {
        return disjoint(tv1.model(),tv2.model());
    } else {
        Box domain=intersection(tv1.domain(),tv2.domain());
        return disjoint(restrict(tv1,domain).model(),restrict(tv2,domain).model());
    }
}

ScalarTaylorFunction intersection(const ScalarTaylorFunction& tv1, const ScalarTaylorFunction& tv2)
{
    ARIADNE_ASSERT(tv1.domain()==tv2.domain());
    return ScalarTaylorFunction(tv1.domain(),intersection(tv1.model(),tv2.model()));
}

Float norm(const ScalarTaylorFunction& f) {
    return norm(f.model());
}

Float distance(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2) {
    return norm(f1-f2);
}

Float distance(const ScalarTaylorFunction& f1, const EffectiveScalarFunction& f2) {
    return distance(f1,ScalarTaylorFunction(f1.domain(),f2,f1.sweeper()));
}

Vector<ScalarTaylorFunction> compose(const Vector<ScalarTaylorFunction>& x, const Vector<ScalarTaylorFunction>& y) {
    ARIADNE_NOT_IMPLEMENTED; }
ScalarTaylorFunction compose(const ScalarTaylorFunction& x, const Vector<ScalarTaylorFunction>& y) {
    ARIADNE_NOT_IMPLEMENTED; }

Vector<ScalarTaylorFunction>
prod(const Matrix<ValidatedNumberType>& A,
     const Vector<ScalarTaylorFunction>& x)
{
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(A.column_size()==x.size());
    for(uint i=0; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x.zero_element().argument_size()); }

    Vector<ScalarTaylorFunction> r(A.row_size(),x.zero_element());
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=A.column_size(); ++j) {
            //r[i]+=A[i][j]*x[j];
            const ValidatedNumberType& Aij=A[i][j]; const ScalarTaylorFunction& xj=x[j]; ScalarTaylorFunction& ri=r[i]; ri+=Aij*xj;
        }
    }
    return r;
}

Matrix<Interval>
jacobian(const Vector<ScalarTaylorFunction>& tv, const Box& x);

ScalarTaylorFunction
midpoint(const ScalarTaylorFunction& f)
{
    ValidatedTaylorModel tm=f.model();
    tm.set_error(0.0);
    return ScalarTaylorFunction(f.domain(),tm);
}


ScalarTaylorFunction
compose(const EffectiveScalarFunction& f, const Vector<ScalarTaylorFunction>& g)
{
    ARIADNE_ASSERT(f.argument_size()==g.size());
    for(uint i=0; i!=g.size(); ++i) {
        ARIADNE_ASSERT(g.zero_element().domain()==g[i].domain());
    }

    Box gdomain=g.zero_element().domain();
    Vector<ValidatedTaylorModel> gmodels(g.size());
    for(uint i=0; i!=g.size(); ++i) { gmodels[i]=g[i].model(); }

    return ScalarTaylorFunction(gdomain,f.evaluate(gmodels));
}

std::ostream& operator<<(std::ostream& os, const Representation<Float>& flt_repr)
{
    const Float& flt=*flt_repr.pointer;
    int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "Float(" << flt << ")";
    os.precision(precision); os.flags(flags);
    return os;
}

std::ostream& operator<<(std::ostream& os, const Representation<Interval>& ivl_repr)
{
    const Interval& ivl=*ivl_repr.pointer;
    int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "Interval("<<ivl.lower()<<","<<ivl.upper()<<")";
    os.precision(precision); os.flags(flags);
    return os;
}

std::ostream& operator<<(std::ostream& os, const Representation< Expansion<Float> >& exp_repr)
{
    const Expansion<Float>& exp=*exp_repr.pointer;
    int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "Expansion<Float>(" << exp.argument_size() << "," << exp.number_of_nonzeros();
    for(Expansion<Float>::const_iterator iter=exp.begin(); iter!=exp.end(); ++iter) {
        for(uint j=0; j!=iter->key().size(); ++j) {
            os << "," << uint(iter->key()[j]);
        }
        os << "," << iter->data();
    }
    os << ")";
    os.precision(precision); os.flags(flags);
    return os;
}

template<class X> std::ostream& operator<<(std::ostream& os, const Representation< Vector<X> >& vec_repr)
{
    const Vector<X>& vec=*vec_repr.pointer;
    ARIADNE_ASSERT(vec.size()!=0);
    os << "(";
    for(uint i=0; i!=vec.size(); ++i) {
        if(i!=0) { os << ","; }
        os << representation(vec[i]);
    }
    os << ")";
    return os;
}
template std::ostream& operator<< <Float> (std::ostream&, const Representation< Vector<Float> >&);
template std::ostream& operator<< <Interval> (std::ostream&, const Representation< Vector<Interval> >&);

std::ostream& operator<<(std::ostream& os, const Representation< Box >& box_repr)
{
    const Vector<Interval>& vec=*box_repr.pointer;
    ARIADNE_ASSERT(vec.size()!=0);
    os << "(";
    for(uint i=0; i!=vec.size(); ++i) {
        if(i!=0) { os << ","; }
        os << representation(vec[i]);
    }
    os << ")";
    return os;
}
std::ostream& operator<<(std::ostream& os, const Representation< Sweeper >& swp_repr)
{
    const Sweeper& swp=*swp_repr.pointer;
    int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << swp;
    os.precision(precision); os.flags(flags);
    return os;
}

template<class T> std::ostream& operator<<(std::ostream& os, const Representation< List<T> >& lst_repr)
{
    const List<T>& lst=*lst_repr.pointer;
    ARIADNE_ASSERT(lst.size()!=0);
    os << "(";
    for(uint i=0; i!=lst.size(); ++i) {
        if(i!=0) { os << ","; }
        os << representation(lst[i]);
    }
    os << ")";
    return os;
}


std::ostream&
ScalarTaylorFunction::write(std::ostream& os) const
{
    os << midpoint(this->polynomial());
    if(this->error()>0.0) { os << "+/-" << this->error(); }
    return os;
}

std::ostream&
ScalarTaylorFunction::repr(std::ostream& os) const
{
    return os << "ScalarTaylorFunction(" << representation(this->domain()) << ", " << representation(this->model().expansion())
              << "," << representation(this->error())<<","<<this->sweeper()<<")";
}

std::ostream&
operator<<(std::ostream& os, const ScalarTaylorFunction& tf)
{
    return tf.write(os);
}

std::ostream& operator<<(std::ostream& os, const Representation<ScalarTaylorFunction>& tf)
{
    return tf.pointer->repr(os);
}

/*
std::ostream& operator<<(std::ostream& os, const Representation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& function=*frepr.pointer;
    ScalarTaylorFunction truncated_function=function;
    truncated_function.set_error(0.0);
    truncated_function.sweep(ThresholdSweeper(TAYLOR_FUNCTION_WRITING_ACCURACY));

    os << midpoint(truncated_function.polynomial());
    if(truncated_function.error()>0.0) { os << "+/-" << truncated_function.error(); }
    if(function.error()>0.0) { os << "+/-" << function.error(); }
    // TODO: Use Unicode +/- literal when this becomes avaialable in C++0x
    return os;
}
*/
/*
std::ostream& operator<<(std::ostream& os, const ModelsRepresentation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& f=*frepr.pointer;
    Float truncatation_error = 0.0;
    os << "<"<<f.domain()<<"\n";
    for(ValidatedTaylorModel::const_iterator iter=f.begin(); iter!=f.end(); ++iter) {
        if(abs(iter->data())>frepr.threshold) { truncatation_error+=abs(iter->data()); }
        else { os << iter->key() << ":" << iter->data() << ","; }
    }
    os << "+/-" << truncatation_error << "+/-" << f.error();
    return os;
}
*/

std::ostream& operator<<(std::ostream& os, const ModelsRepresentation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& f=*frepr.pointer;
    ScalarTaylorFunction tf=f;
    tf.clobber();
    tf.sweep(ThresholdSweeper(frepr.threshold));
    os << "("<<tf.model()<<"+/-"<<f.error();
    return os;
}

std::ostream& operator<<(std::ostream& os, const PolynomialRepresentation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& function=*frepr.pointer;
    ScalarTaylorFunction truncated_function = function;
    truncated_function.clobber();
    truncated_function.sweep(ThresholdSweeper(frepr.threshold));
    Float truncatation_error = truncated_function.error();
    truncated_function.clobber();
    Polynomial<Float> polynomial_function = midpoint(polynomial(truncated_function));
    if(frepr.names.empty()) { os << polynomial_function; }
    else { os << named_argument_repr(polynomial_function,frepr.names); }
    os << "+/-" << truncatation_error << "+/-" << function.error();
    return os;
}





bool
check(const Vector<ScalarTaylorFunction>& tv)
{
    for(uint i=0; i!=tv.size(); ++i) {
        if(tv.zero_element().domain()!=tv[i].domain()) { return false; }
    }
    return true;
}

Vector< Expansion<Float> >
expansion(const Vector<ScalarTaylorFunction>& x)
{
    Vector< Expansion<Float> > r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].expansion();
    }
    return r;
}

Vector<Float>
error(const Vector<ScalarTaylorFunction>& x)
{
    Vector<Float> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].error();
    }
    return r;
}

Vector<Float>
value(const Vector<ScalarTaylorFunction>& x)
{
    Vector<Float> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

Vector<Interval>
ranges(const Vector<ScalarTaylorFunction>& x)
{
    Vector<Interval> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].range();
    }
    return r;
}


Vector<ValidatedNumberType>
evaluate(const Vector<ScalarTaylorFunction>& tv, const Vector<ValidatedNumberType>& x)
{
    Vector<ValidatedNumberType> r(tv.size());
    for(uint i=0; i!=tv.size(); ++i) {
        r[i]=evaluate(tv[i],x);
    }
    return r;
}

Matrix<ValidatedNumberType>
jacobian(const Vector<ScalarTaylorFunction>& tv, const Vector<ValidatedNumberType>& x)
{
    ARIADNE_ASSERT(check(tv));
    const Vector<ValidatedNumberType>& dom=tv.zero_element().domain();
    const uint n=dom.size();
    Vector< Differential<ValidatedNumberType> > s(n,n,1u);
    for(uint j=0; j!=n; ++j) {
        Interval dj=dom[j];
        s[j].set_value((x[j]-med_ivl(dj))/rad_ivl(dj));
        s[j].set_gradient(j,rec(rad_ivl(dj)));
    }
    Vector< Expansion<Float> > p=expansion(tv);
    Vector< Differential<ValidatedNumberType> > d=evaluate(p,s);
    return d.jacobian();
}





VectorTaylorFunction::VectorTaylorFunction()
    : _domain(), _models()
{
}

VectorTaylorFunction::VectorTaylorFunction(uint k)
    : _domain(), _models(k)
{
}

VectorTaylorFunction::VectorTaylorFunction(uint m, const Box& d, Sweeper swp)
    : _domain(d), _models(m,ValidatedTaylorModel(d.size(),swp))
{
}

VectorTaylorFunction::VectorTaylorFunction(uint k, const ScalarTaylorFunction& f)
    : _domain(f.domain()), _models(k,f.model())
{
}

VectorTaylorFunction::VectorTaylorFunction(const ValidatedVectorFunctionModel& f)
    : _domain(), _models()
{
    ARIADNE_ASSERT(dynamic_cast<const VectorTaylorFunction*>(&f.reference()));
    *this = dynamic_cast<const VectorTaylorFunction&>(f.reference());
}

VectorTaylorFunction& VectorTaylorFunction::operator=(const ValidatedVectorFunctionModel& f)
{
    ARIADNE_ASSERT(dynamic_cast<const VectorTaylorFunction*>(&f.reference()));
    *this = dynamic_cast<const VectorTaylorFunction&>(f.reference());
    return *this;
}


VectorTaylorFunction substitute(const VectorTaylorFunction& f, uint k, const ScalarTaylorFunction& h) {
    ARIADNE_ASSERT_MSG(f.argument_size()==h.argument_size()+1u,"f="<<f<<", k="<<k<<", h="<<h);
    const Interval& dk=f.domain()[k];
    ValidatedTaylorModel unscaled_model=unscale(h.model(),dk);
    return VectorTaylorFunction(h.domain(),substitute(f.models(),k,unscaled_model));
}

ScalarTaylorFunction substitute(const ScalarTaylorFunction& f, uint k, const ScalarTaylorFunction& h) {
    const Interval& dk=f.domain()[k];
    ValidatedTaylorModel unscaled_model=unscale(h.model(),dk);
    return ScalarTaylorFunction(h.domain(),substitute(f.model(),k,unscaled_model));
}


VectorTaylorFunction::VectorTaylorFunction(const Box& d,
                                           const Vector<ValidatedTaylorModel>& f)
    : _domain(d), _models(f)
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT_MSG(d.size()==f[i].argument_size(),"d="<<d<<", f="<<f);
    }
}

VectorTaylorFunction::VectorTaylorFunction(const Box& d,
                                           const Vector< Expansion<Float> >& f,
                                           Sweeper swp)
    : _domain(d), _models(f.size())
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
        _models[i]=ValidatedTaylorModel(f[i],0.0,swp);
    }
}

VectorTaylorFunction::VectorTaylorFunction(const Box& d,
                                           const Vector< Expansion<Float> >& f,
                                           const Vector<Float>& e,
                                           Sweeper swp)
    : _domain(d), _models(f.size(),ValidatedTaylorModel(d.size(),swp))
{
    ARIADNE_ASSERT(f.size()==e.size());
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
        _models[i]=ValidatedTaylorModel(f[i],e[i],swp);
    }
}



VectorTaylorFunction::VectorTaylorFunction(const Box& d,
                                           const ValidatedVectorFunction& f,
                                           const Sweeper& swp)
    : _domain(d), _models(f.result_size())
{
    //ARIADNE_ASSERT_MSG(f.result_size()>0, "d="<<d<<", f="<<f<<", swp="<<swp);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<ValidatedTaylorModel> x=ValidatedTaylorModel::scalings(d,swp);
    this->_models=f.evaluate(x);
    this->sweep();
}

/*
VectorTaylorFunction::VectorTaylorFunction(const Box& d,
                                           const EffectiveVectorFunction& f,
                                           const Sweeper& swp)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<ValidatedTaylorModel> x=ValidatedTaylorModel::scalings(d,swp);
    for(uint i=0; i!=x.size(); ++i) { x[i].set_sweeper(swp); }
    this->_models=f.evaluate(x);
    this->sweep();
}
*/

VectorTaylorFunction::VectorTaylorFunction(const Vector<ScalarTaylorFunction>& v)
    : _domain(), _models(v.size())
{
    ARIADNE_ASSERT(v.size()>0);
    for(uint i=0; i!=v.size(); ++i) { ARIADNE_ASSERT(v[i].domain()==v.zero_element().domain()); }
    this->_domain=v.zero_element().domain();
    for(uint i=0; i!=v.size(); ++i) {
        this->_models[i]=v[i].model();
    }
}

VectorTaylorFunction::VectorTaylorFunction(const List<ScalarTaylorFunction>& v)
    : _domain(), _models(v.size())
{
    ARIADNE_ASSERT(v.size()>0);
    for(uint i=1; i!=v.size(); ++i) { ARIADNE_ASSERT(v[i].domain()==v[0].domain()); }
    this->_domain=v[0].domain();
    for(uint i=0; i!=v.size(); ++i) {
        this->_models[i]=v[i].model();
    }
}

VectorTaylorFunction::VectorTaylorFunction(std::initializer_list<ScalarTaylorFunction> lst)
    : _domain(), _models(lst.size())
{
    *this=VectorTaylorFunction(List<ScalarTaylorFunction>(lst));
}


VectorTaylorFunction* VectorTaylorFunction::_clone() const
{
    return new VectorTaylorFunction(*this);
}

VectorTaylorFunction* VectorTaylorFunction::_create() const
{
    return new VectorTaylorFunction(this->result_size(), ScalarTaylorFunction(this->domain(),this->sweeper()));
}

ScalarTaylorFunction* VectorTaylorFunction::_create_zero() const
{
    return new ScalarTaylorFunction(this->domain(),this->sweeper());
}

VectorTaylorFunction* VectorTaylorFunction::_create_identity() const
{
    Sweeper sweeper=this->sweeper();
    VectorTaylorFunction* result = new VectorTaylorFunction(this->domain().size(), ScalarTaylorFunction(this->domain(),sweeper));
    for(uint i=0; i!=result->size(); ++i) { (*result)[i]=ScalarTaylorFunction::coordinate(this->domain(),i,sweeper); }
    return result;
}

Void VectorTaylorFunction::adjoin(const ScalarTaylorFunction& sf)
{
    ARIADNE_ASSERT_MSG(sf.domain()==this->domain(),"sf="<<sf);
    this->_models=join(this->_models,sf.model());
}







VectorTaylorFunction
VectorTaylorFunction::constant(const Box& d, const Vector<Float>& c, Sweeper swp)
{
    return VectorTaylorFunction(d,ValidatedTaylorModel::constants(d.size(),c,swp));
}

VectorTaylorFunction
VectorTaylorFunction::constant(const Box& d, const Vector<ValidatedNumberType>& c, Sweeper swp)
{
    return VectorTaylorFunction(d,ValidatedTaylorModel::constants(d.size(),c,swp));
}

VectorTaylorFunction
VectorTaylorFunction::identity(const Box& d, Sweeper swp)
{
    return VectorTaylorFunction(d,ValidatedTaylorModel::scalings(d,swp));
}

VectorTaylorFunction
VectorTaylorFunction::projection(const Box& d, uint imin, uint imax, Sweeper swp)
{
    return VectorTaylorFunction(ScalarTaylorFunction::variables(d,imin,imax,swp));
}


Polynomial<ValidatedNumberType> polynomial(const ValidatedTaylorModel& tm) {
    return Polynomial<ValidatedNumberType>(tm.expansion())+ValidatedNumberType(-tm.error(),+tm.error());
}

Vector< Polynomial<ValidatedNumberType> >
VectorTaylorFunction::polynomials() const
{
    Vector<Polynomial<ValidatedNumberType> > p(this->result_size(),Polynomial<ValidatedNumberType>(this->argument_size()));
    for(uint i=0; i!=this->result_size(); ++i) {
        p[i]=static_cast<ScalarTaylorFunction>((*this)[i]).polynomial();
    }
    return p;
}

Vector< Expansion<Float> > const
VectorTaylorFunction::expansions() const
{
    Vector< Expansion<Float> > e(this->result_size(),Expansion<Float>(this->argument_size()));
    for(uint i=0; i!=this->result_size(); ++i) {
        e[i]=this->models()[i].expansion();
    }
    return e;
}

Vector<Float> const
VectorTaylorFunction::errors() const
{
    Vector<Float> e(this->result_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        e[i]=this->models()[i].error();
    }
    return e;
}

Float const
VectorTaylorFunction::error() const
{
    Float e=0.0;
    for(uint i=0; i!=this->result_size(); ++i) {
        e=max(e,this->models()[i].error());
    }
    return e;
}

ValidatedVectorFunction
VectorTaylorFunction::function() const
{
    return new VectorTaylorFunction(*this);
}

bool
VectorTaylorFunction::operator==(const VectorTaylorFunction& tm) const
{
    return this->_models==tm._models;
}



bool
VectorTaylorFunction::operator!=(const VectorTaylorFunction& p2) const
{
    return !(*this==p2);
}



Sweeper
VectorTaylorFunction::sweeper() const
{
    ARIADNE_ASSERT(this->size()>0); return this->_models[0].sweeper();
}


void
VectorTaylorFunction::set_sweeper(Sweeper swp)
{
    for(uint i=0; i!=this->result_size(); ++i) {
        this->_models[i].set_sweeper(swp);
    }
}

const Box&
VectorTaylorFunction::domain() const
{
    return this->_domain;
}

const Box
VectorTaylorFunction::codomain() const
{
    Box result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].range();
    }
    return result;
}


const Box
VectorTaylorFunction::range() const
{
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].range();
    }
    return result;
}


const Vector<Float>
VectorTaylorFunction::centre() const
{
    Vector<Float> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].value();
    }
    return result;
}


const Vector<ValidatedTaylorModel>&
VectorTaylorFunction::models() const
{
    return this->_models;
}

const ValidatedTaylorModel&
VectorTaylorFunction::model(uint i) const
{
    return this->_models[i];
}

ValidatedTaylorModel&
VectorTaylorFunction::model(uint i)
{
    return this->_models[i];
}




uint
VectorTaylorFunction::argument_size() const
{
    return this->_domain.size();
}


uint
VectorTaylorFunction::result_size() const
{
    return this->_models.size();
}


ScalarTaylorFunction
VectorTaylorFunction::operator[](uint i) const
{
    return this->get(i);
}

VectorTaylorFunctionElementReference
VectorTaylorFunction::operator[](uint i)
{
    return VectorTaylorFunctionElementReference(*this,i);
}

ScalarTaylorFunction
VectorTaylorFunction::get(uint i) const
{
    return ScalarTaylorFunction(this->_domain,this->_models[i]);
}

void
VectorTaylorFunction::set(uint i, const ScalarTaylorFunction& e)
{
    ARIADNE_ASSERT_MSG(this->size()>i,"Cannot set "<<i<<"th element of VectorTaylorFunction "<<(*this));
    if(this->domain().size()!=0) {
        ARIADNE_ASSERT_MSG(e.domain()==this->domain(),"Domain of "<<e<<" conflicts with existing domain "<<this->domain());
    } else {
        this->_domain=e.domain();
    }
    this->_models[i]=e.model();
}














template<class T>
void
VectorTaylorFunction::_compute(Vector<T>& r, const Vector<T>& a) const
{
    typedef typename T::NumericType X;
    const VectorTaylorFunction& f=*this;
    Vector<T> sx=Ariadne::unscale(a,f._domain);
    for(uint i=0; i!=r.size(); ++i) {
        T ri=Ariadne::evaluate(this->_models[i].expansion(),sx);
        X e=convert_error<X>(this->_models[i].error());
        r[i]=ri+e;
    }
}




VectorTaylorFunction&
VectorTaylorFunction::sweep()
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].sweep();
    }
    return *this;
}

VectorTaylorFunction&
VectorTaylorFunction::sweep(const SweeperInterface& sweeper)
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].sweep(sweeper);
    }
    return *this;
}


Void
VectorTaylorFunction::clobber()
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].clobber();
    }
}



Vector<Float>
VectorTaylorFunction::evaluate(const Vector<Float>& x) const
{
    const VectorTaylorFunction& f=*this;
    if(!contains(f.domain(),x)) {
        ARIADNE_THROW(DomainException,"f.evaluate(x) with f="<<f<<", x="<<x,"x is not an element of f.domain()="<<f.domain());
    }
    Vector<Float> sx=Ariadne::unscale(x,f._domain);
    Vector<Float> r(this->result_size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=Ariadne::evaluate(this->_models[i].expansion(),sx);
    }
    return r;
}


Vector<ValidatedNumberType>
VectorTaylorFunction::evaluate(const Vector<ValidatedNumberType>& x) const
{
    const VectorTaylorFunction& f=*this;
    if(!subset(x,f.domain())) {
        ARIADNE_THROW(DomainException,"f.evaluate(x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
    }
    Vector<ValidatedNumberType> sx=Ariadne::unscale(x,f._domain);
    return Ariadne::evaluate(f._models,sx);
}

Vector<ValidatedNumberType>
VectorTaylorFunction::operator()(const Vector<ValidatedNumberType>& x) const
{
    return this->evaluate(x);
}

Matrix<ValidatedNumberType>
VectorTaylorFunction::jacobian(const Vector<ValidatedNumberType>& x) const
{
    Matrix<ValidatedNumberType> J=Ariadne::jacobian(this->_models,unscale(x,this->_domain));
    for(uint j=0; j!=J.column_size(); ++j) {
        ValidatedNumberType rad=rad_ivl(this->_domain[j]);
        for(uint i=0; i!=J.row_size(); ++i) {
            J[i][j]/=rad;
        }
    }
    return J;
}

Void
VectorTaylorFunction::restrict(const Box& x)
{
    *this=Ariadne::restrict(*this,x);
}


VectorTaylorFunction
join(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    ARIADNE_ASSERT_MSG(f1.domain()==f2.domain(),"f1="<<f1<<", f2="<<f2);
    return VectorTaylorFunction(f1.domain(),join(f1.models(),f2.model()));
}

VectorTaylorFunction
join(const VectorTaylorFunction& f, const VectorTaylorFunction& g)
{
    ARIADNE_ASSERT(f.domain()==g.domain());
    return VectorTaylorFunction(f.domain(),join(f.models(),g.models()));
}

VectorTaylorFunction
join(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    ARIADNE_ASSERT(f1.domain()==f2.domain());
    return VectorTaylorFunction(f1.domain(),join(f1.model(),f2.model()));
}

VectorTaylorFunction
join(const ScalarTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    ARIADNE_ASSERT(f1.domain()==f2.domain());
    return VectorTaylorFunction(f1.domain(),join(f1.model(),f2.models()));
}

VectorTaylorFunction
combine(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    return VectorTaylorFunction(join(f1.domain(),f2.domain()),combine(f1.model(),f2.model()));
}

VectorTaylorFunction
combine(const ScalarTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    return VectorTaylorFunction(join(f1.domain(),f2.domain()),combine(f1.model(),f2.models()));
}

VectorTaylorFunction
combine(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    return VectorTaylorFunction(join(f1.domain(),f2.domain()),combine(f1.models(),f2.model()));
}

VectorTaylorFunction
combine(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    return VectorTaylorFunction(join(f1.domain(),f2.domain()),combine(f1.models(),f2.models()));
}


VectorTaylorFunction
embed(const VectorTaylorFunction& f, const Interval& d)
{
    return embed(Box(),f,Box(1u,d));
}

VectorTaylorFunction
embed(const VectorTaylorFunction& f, const Box& d)
{
    return embed(Box(),f,d);
}

VectorTaylorFunction
embed(const Box& d, const VectorTaylorFunction& f)
{
    return embed(d,f,Box());
}

VectorTaylorFunction
embed(const Box& d1, const VectorTaylorFunction& f, const Box& d2)
{
    return VectorTaylorFunction(join(d1,f.domain(),d2),embed(embed(d1.size(),f.models()),d2.size()));
}

VectorTaylorFunction
restrict(const VectorTaylorFunction& f, const Box& d)
{
    ARIADNE_ASSERT_MSG(subset(d,f.domain()),"Cannot restrict "<<f<<" to non-sub-domain "<<d);
    if(d==f.domain()) { return f; }
    VectorTaylorFunction r(f.result_size(),d,f.sweeper());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,restrict(f[i],d));
    }
    return r;
}

std::pair<VectorTaylorFunction,VectorTaylorFunction>
split(const VectorTaylorFunction& tf, uint j)
{
    pair< Vector<ValidatedTaylorModel>,Vector<ValidatedTaylorModel> > models=split(tf.models(),j);
    pair<Box,Box> subdomains=split(tf.domain(),j);
    return make_pair(VectorTaylorFunction(subdomains.first,models.first),
                     VectorTaylorFunction(subdomains.second,models.second));

}

bool
refines(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(uint i=0; i!=f1.result_size(); ++i) {
        if(!refines(f1[i],f2[i])) { return false; }
    }
    return true;
}

bool
disjoint(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(uint i=0; i!=f1.result_size(); ++i) {
        if(disjoint(f1[i],f2[i])) { return true; }
    }
    return false;
}

VectorTaylorFunction
intersection(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    VectorTaylorFunction r(f1.result_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r[i]=intersection(f1[i],f2[i]);
    }
    return r;
}

VectorTaylorFunction&
operator+=(VectorTaylorFunction& f, const VectorTaylorFunction& g)
{
    ARIADNE_ASSERT(f.result_size()==g.result_size());
    ARIADNE_ASSERT(subset(f.domain(),g.domain()));
    ARIADNE_ASSERT(f.domain()==g.domain());
    f._models+=g._models;
    return f;
}

VectorTaylorFunction&
operator-=(VectorTaylorFunction& f, const VectorTaylorFunction& g)
{
    ARIADNE_ASSERT(f.result_size()==g.result_size());
    ARIADNE_ASSERT(subset(f.domain(),g.domain()));
    ARIADNE_ASSERT(f.domain()==g.domain());
    f._models+=g._models;
    return f;
}

VectorTaylorFunction&
operator+=(VectorTaylorFunction& f, const Vector<ValidatedNumberType>& c)
{
    ARIADNE_ASSERT(f.result_size()==c.size());
    f._models+=c;
    return f;
}

VectorTaylorFunction&
operator-=(VectorTaylorFunction& f, const Vector<ValidatedNumberType>& c)
{
    ARIADNE_ASSERT(f.result_size()==c.size());
    f._models-=c;
    return f;
}

VectorTaylorFunction&
operator*=(VectorTaylorFunction& f, const ValidatedNumberType& c)
{
    f._models*=c;
    return f;
}

VectorTaylorFunction&
operator/=(VectorTaylorFunction& f, const ValidatedNumberType& c)
{
    f._models/=c;
    return f;
}


VectorTaylorFunction
operator+(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    ARIADNE_ASSERT_MSG(!empty(intersection(f1.domain(),f2.domain())),
                       "operator+(VectorTaylorFunction f1, VectorTaylorFunction f2) with f1="<<f1<<" f2="<<f2<<
                       ": domains are disjoint");
    if(f1.domain()==f2.domain()) {
        return VectorTaylorFunction(f1.domain(),Vector<ValidatedTaylorModel>(f1.models()+f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator+(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}


VectorTaylorFunction
operator-(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    ARIADNE_ASSERT(!empty(intersection(f1.domain(),f2.domain())));
    if(f1.domain()==f2.domain()) {
        return VectorTaylorFunction(f1.domain(),Vector<ValidatedTaylorModel>(f1.models()-f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator-(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}

VectorTaylorFunction
operator*(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    ARIADNE_ASSERT(!empty(intersection(f1.domain(),f2.domain())));
    if(f1.domain()==f2.domain()) {
        return VectorTaylorFunction(f1.domain(),Vector<ValidatedTaylorModel>(f1.models()*f2.model()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator*(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}

VectorTaylorFunction
operator*(const ScalarTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    ARIADNE_ASSERT(!empty(intersection(f1.domain(),f2.domain())));
    if(f1.domain()==f2.domain()) {
        return VectorTaylorFunction(f1.domain(),Vector<ValidatedTaylorModel>(f1.model()*f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator*(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}



VectorTaylorFunction
operator-(const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(-f.models()));
}

VectorTaylorFunction
operator*(const ValidatedNumberType& c, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(f.models()*c));
}

VectorTaylorFunction
operator*(const VectorTaylorFunction& f, const ValidatedNumberType& c)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(f.models()*c));
}

VectorTaylorFunction
operator/(const VectorTaylorFunction& f, const ValidatedNumberType& c)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(f.models()/c));
}

VectorTaylorFunction
operator+(const VectorTaylorFunction& f, const Vector<ValidatedNumberType>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(f.models()+c));
}

VectorTaylorFunction
operator-(const VectorTaylorFunction& f, const Vector<ValidatedNumberType>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(f.models()-c));
}

VectorTaylorFunction
operator*(const Matrix<Float>& A, const VectorTaylorFunction& f)
{
    ARIADNE_PRECONDITION(A.column_size()==f.size());
    Vector<ValidatedTaylorModel> models(A.row_size(),ValidatedTaylorModel(f.argument_size(),f.sweeper()));
    for(size_t i=0; i!=A.row_size(); ++i) {
        for(size_t j=0; j!=A.column_size(); ++j) {
            models[i] += A.get(i,j) * f.model(j);
        }
    }
    return VectorTaylorFunction(f.domain(),models);
}

VectorTaylorFunction
operator*(const Matrix<ValidatedNumberType>& A, const VectorTaylorFunction& f)
{
    ARIADNE_PRECONDITION(A.column_size()==f.size());
    Vector<ValidatedTaylorModel> models(A.row_size(),ValidatedTaylorModel(f.argument_size(),f.sweeper()));
    for(size_t i=0; i!=A.row_size(); ++i) {
        for(size_t j=0; j!=A.column_size(); ++j) {
            models[i] += A.get(i,j) * f.model(j);
        }
    }
    return VectorTaylorFunction(f.domain(),models);
}


VectorTaylorFunction
operator-(const VectorTaylorFunction& f1, const EffectiveVectorFunction& f2) {
    return f1 - VectorTaylorFunction(f1.domain(),f2,f1.sweeper());
}





VectorTaylorFunction
partial_evaluate(const VectorTaylorFunction& tf, uint k, const ExactNumberType& c)
{
    return partial_evaluate(tf,k,ValidatedNumberType(c));
}

VectorTaylorFunction
partial_evaluate(const VectorTaylorFunction& tf, uint k, const ValidatedNumberType& c)
{
    // Scale c to domain
    const uint as=tf.argument_size();
    ARIADNE_ASSERT(k<as);
    const Vector<Interval>& domain=tf.domain();
    const Interval& dk=domain[k];
    ValidatedNumberType sc=(c-med_ivl(dk))/rad_ivl(dk);

    Vector<Interval> new_domain(as-1);
    for(uint i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
    for(uint i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

    Vector<ValidatedTaylorModel> new_models=partial_evaluate(tf.models(),k,sc);

    return VectorTaylorFunction(new_domain,new_models);
}


VectorTaylorFunction
partial_restrict(const VectorTaylorFunction& tf, uint k, const Interval& d)
{
    VectorTaylorFunction r(tf.result_size(),tf.domain(),tf.sweeper());
    for(uint i=0; i!=tf.result_size(); ++i) {
        r[i]=partial_restrict(tf[i],k,d);
    }
    return r;
}

VectorTaylorFunction
restrict(const VectorTaylorFunction& tf, uint k, const Interval& d)
{
    return partial_restrict(tf,k,d);
}


Vector<ValidatedNumberType>
evaluate(const VectorTaylorFunction& f, const Vector<ValidatedNumberType>& x) {
    if(!subset(x,f.domain())) {
        ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
    }
    return unchecked_evaluate(f,x);
}

Vector<ValidatedNumberType>
unchecked_evaluate(const VectorTaylorFunction& f, const Vector<ValidatedNumberType>& x) {
    return evaluate(f.models(),unscale(x,f.domain()));
}

VectorTaylorFunction
compose(const EffectiveVectorFunction& g, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),g.evaluate(f.models()));
}

VectorTaylorFunction
compose(const ValidatedVectorFunction& g, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),g.evaluate(f.models()));
}

VectorTaylorFunction
compose(const VectorTaylorFunction& g, const VectorTaylorFunction& f)
{
    if(!subset(f.codomain(),g.domain())) {
        ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
    }
    return unchecked_compose(g,f);
}


VectorTaylorFunction
unchecked_compose(const VectorTaylorFunction& g, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),compose(g.models(),unscale(f.models(),g.domain())));
}



VectorTaylorFunction
derivative(const VectorTaylorFunction& f, uint k)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    ValidatedNumberType fdomkrad=rad_ivl(f.domain()[k].lower(),f.domain()[k].upper());
    VectorTaylorFunction g=f;
    for(uint i=0; i!=g.size(); ++i) {
        g[i]=derivative(f[i],k);
    }
    return g;
}

VectorTaylorFunction
antiderivative(const VectorTaylorFunction& f, uint k)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    ValidatedNumberType fdomkrad=rad_ivl(f.domain()[k].lower(),f.domain()[k].upper());
    VectorTaylorFunction g=f;
    for(uint i=0; i!=g.size(); ++i) {
        g._models[i].antidifferentiate(k);
        g._models[i]*=fdomkrad;
    }
    return g;
}

VectorTaylorFunction
antiderivative(const VectorTaylorFunction& f, uint k, ExactNumberType c)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    ValidatedNumberType fdomkrad=rad_ivl(f.domain()[k].lower(),f.domain()[k].upper());
    VectorTaylorFunction g=f;
    for(uint i=0; i!=g.size(); ++i) {
        g[i]=antiderivative(f[i],k,c);
    }
    return g;
}






Float norm(const VectorTaylorFunction& f) {
    Float res=0.0;
    for(uint i=0; i!=f.result_size(); ++i) {
        res=max(res,norm(f[i]));
    }
    return res;
}

Float distance(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2) {
    return norm(f1-f2);
}

Float distance(const VectorTaylorFunction& f1, const EffectiveVectorFunction& f2) {
    return distance(f1,VectorTaylorFunction(f1.domain(),f2,f1.sweeper()));
}


std::ostream&
VectorTaylorFunction::write(std::ostream& os) const
{
    Vector< Polynomial<Float> > p=midpoint(this->polynomials());
    Vector<Float> e=this->errors();

    os << "[";
    for(uint i=0; i!=this->result_size(); ++i) {
        if(i!=0) { os << ", "; }
        os << p[i];
        if(e[i]!=0) { os << "+/-" << e[i]; }
    }
    os << "]";
    return os;
}

std::ostream&
VectorTaylorFunction::repr(std::ostream& os) const
{
    return os << "VectorTaylorFunction("
              << representation(this->domain()) << ", " << representation(this->expansions()) << ", "
              << representation(this->domain()) << ", " << representation(this->sweeper()) << ")";
}

std::ostream& operator<<(std::ostream& os, const Representation<VectorTaylorFunction>& repr)
{
    return repr.pointer->repr(os);
}

std::ostream& operator<<(std::ostream& os, const PolynomialRepresentation<VectorTaylorFunction>& repr)
{
    const VectorTaylorFunction& function = *repr.pointer;
    os << "[";
    for(uint i=0; i!=function.result_size(); ++i) {
        if(i!=0) { os << ","; }
        os << polynomial_repr(function[i],repr.threshold,repr.names);
    }
    return os << "]";
}

std::ostream& operator<<(std::ostream& os, const PolynomialRepresentation< List<ScalarTaylorFunction> >& repr)
{
    const List<ScalarTaylorFunction>& functions = *repr.pointer;
    os << "[";
    for(uint i=0; i!=functions.size(); ++i) {
        if(i!=0) { os << ","; }
        os << polynomial_repr(functions[i],repr.threshold);
    }
    return os << "]";
}



std::ostream&
operator<<(std::ostream& os, const VectorTaylorFunction& p)
{
    return p.write(os);
}

Polynomial<ValidatedNumberType> polynomial(const ScalarTaylorFunction& tfn) {
    return Polynomial<ValidatedNumberType>(tfn.polynomial());
}

Vector< Polynomial<ValidatedNumberType> > polynomial(const VectorTaylorFunction& tfn) {
    return tfn.polynomials();
}

List< Polynomial<ValidatedNumberType> > polynomials(const List<ScalarTaylorFunction>& tfns) {
    List< Polynomial<ValidatedNumberType> > result;
    for(uint i=0; i!=tfns.size(); ++i) {
        result.append(polynomial(tfns[i]));
    }
    return result;
}


ScalarTaylorFunction
TaylorFunctionFactory::create(const Box& domain, const ValidatedScalarFunctionInterface& function) const
{
    return ScalarTaylorFunction(domain,function,this->_sweeper);
}

ScalarTaylorFunction*
TaylorFunctionFactory::_create(const Box& domain, const ValidatedScalarFunctionInterface& function) const
{
    return new ScalarTaylorFunction(domain,function,this->_sweeper);
}

VectorTaylorFunction
TaylorFunctionFactory::create(const Box& domain, const ValidatedVectorFunctionInterface& function) const
{
    return VectorTaylorFunction(domain,function,this->_sweeper);
}

VectorTaylorFunction*
TaylorFunctionFactory::_create(const Box& domain, const ValidatedVectorFunctionInterface& function) const
{
    return new VectorTaylorFunction(domain,function,this->_sweeper);
}

ScalarTaylorFunction
TaylorFunctionFactory::create_zero(const Box& domain) const
{
    return ScalarTaylorFunction::zero(domain,this->_sweeper);
}

ScalarTaylorFunction
TaylorFunctionFactory::create_constant(const Box& domain, ValidatedNumberType value) const
{
    return ScalarTaylorFunction::constant(domain,value,this->_sweeper);
}

ScalarTaylorFunction
TaylorFunctionFactory::create_coordinate(const Box& domain, uint k) const
{
    return ScalarTaylorFunction::coordinate(domain,k,this->_sweeper);
}

ScalarTaylorFunction
TaylorFunctionFactory::create_identity(const Interval& domain) const
{
    return ScalarTaylorFunction::identity(domain,this->_sweeper);
}

VectorTaylorFunction
TaylorFunctionFactory::create_identity(const Box& domain) const
{
    return VectorTaylorFunction::identity(domain,this->_sweeper);
}


FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory() {
    return new TaylorFunctionFactory(Sweeper());
}

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory(double sweep_threshold) {
    return new TaylorFunctionFactory(ThresholdSweeper(sweep_threshold));
}

/*
latexstream&
operator<<(Output::latexstream& texs, const VectorTaylorFunction& p)
{
    using namespace Function;
    texs << "%VectorTaylorFunction\n";
    texs << "\\ensuremath{\n";
    texs << "\\left( \\begin{Array}{c}\n";
    char var='x';
    for(uint i=0; i!=p.result_size(); ++i) {
        bool first = true;
        if(i!=0) { texs << "\\\\"; }
        for(MultiIndex j(p.argument_size()); j.degree()<=p.order(); ++j) {
            const ValidatedNumberType& a=p.centre_derivatives()[i][j];
            if(a!=0) {
                if(first) { first=false; }
                else { if(a>0) { texs << '+'; } }
                if(a==1) { if(j.degree()==0) { texs << a; } }
                else if(a==-1) { if(j.degree()==0) { texs << a; } else { texs << '-'; } }
                else { texs << a << ' '; }
                for(uint k=0; k!=p.argument_size(); ++k) {
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
    texs << "\\end{Array}\\right)\n}\n";
    return texs;
}
*/


} // namespace Ariadne
