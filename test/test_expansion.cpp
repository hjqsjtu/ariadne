/***************************************************************************
 *            test_expansion.cc
 *
 *  Copyright 2009--17  Pieter Collins
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
#include <vector>
#include "config.h"
#include "numeric/numeric.h"
#include "algebra/expansion.h"

#include "test.h"
using namespace std;
using namespace Ariadne;

template<class K, class V> struct MapValue {
    typedef K key_type;
    K first; V second;
    MapValue(const K& k, const V& v) : first(k), second(v) { }
    const K& key() const { return first; }
    const V& data() const { return second; }
};


template<class F> class TestExpansion
{
    typedef PrecisionType<F> PR;
    typedef MultiIndex MI;
    typedef Expansion<F> E;
    typedef typename Expansion<F>::Iterator Iterator;
    typedef typename Expansion<F>::ConstIterator ConstIterator;
  public:
    PR prec; F zero;
    GradedLess graded_less;
    LexicographicLess lexicographic_less;
    ReverseLexicographicLess reverse_lexicographic_less;
  public:
    TestExpansion(F const& z);
    Void test();
  private:
    Void test_working();
    Void test_concept();
    Void test_iterator_concept();
    Void test_data_access();
    Void test_equality();
    Void test_cleanup();
    Void test_constructors();
    //Void test_indexing();
    Void test_find();
    Void test_embed();
};


template<class F> TestExpansion<F>::TestExpansion(F const& z)
    : prec(z.precision()), zero(z)
{
}

template<class F> Void TestExpansion<F>::test()
{
    ARIADNE_TEST_CALL(test_working());
    ARIADNE_TEST_CALL(test_data_access());
    ARIADNE_TEST_CALL(test_equality());
    ARIADNE_TEST_CALL(test_cleanup());
    ARIADNE_TEST_CALL(test_constructors());
    //ARIADNE_TEST_CALL(test_indexing());
    ARIADNE_TEST_CALL(test_find());
    ARIADNE_TEST_CALL(test_embed());
}


template<class F> Void TestExpansion<F>::test_working()
{
    Expansion<F> e=Expansion<F>(3,zero);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_PRINT(e.size());
    // Append values
    ARIADNE_TEST_EXECUTE(e.append({0,0,0},2.0));
    ARIADNE_TEST_EXECUTE(e.append({1,0,0},3.0));
    ARIADNE_TEST_EXECUTE(e.append({0,1,0},5.0));
    ARIADNE_TEST_EXECUTE(e.append({1,0,1},7.0));
    ARIADNE_TEST_PRINT(e);

    ARIADNE_TEST_EQUAL(e.find({0,0,0}),e.begin());
    //assert(false);
}


template<class F> Void TestExpansion<F>::test_concept()
{
    PR prec=zero.precision();
    F x(5,prec);
    SizeType as(3);

    Expansion<F> e(as,prec);
    const Expansion<F> ce(as,prec);

    e=Expansion<F>(as,zero);
    e=Expansion<F>(ce);

    //e=Expansion<F>(3,1, {0.0, 0.0,0.0,0.0}, prec);
    //e=Expansion<F>(3,1, {1, 2,3,5.0}, prec);
    e=Expansion<F>({ {{0,0},1}, {{1,0,0},2}, {{0,1,0},3}, {{0,0,1},5.0} }, prec);
    e=Expansion<F>({ {{0,0},1}, {{1,0,0},2}, {{0,1,0},3}, {{0,0,1},5.0} }, prec);

    MultiIndex a(as);
    e.reserve(2u);
    e.set(a,x);
    e.prepend(a,x);
    e.append(a,x);
    e.append_sum(a,a,x);
    e.clear();

    e.index_sort(GradedLess());
    e.index_sort(LexicographicLess());
    e.index_sort(GradedIndexLess());
    e.sort(ReverseLexicographicIndexLess());

    x=ce[a];

    ce.number_of_terms();
    ce.argument_size();

    e.erase(e.begin());

    ce.check();
}

template<class F> Void TestExpansion<F>::test_iterator_concept()
{
    MultiIndex a(3);
    Expansion<F> e(3,zero);
    const Expansion<F> cp(3,zero);

    typename Expansion<F>::Iterator iter=e.begin(); iter=e.end(); iter=e.find(a);
    typename Expansion<F>::ConstIterator citer=e.begin(); citer=e.end(); citer=e.find(a);
    citer=e.begin(); citer=cp.end(); citer=cp.find(a);

    typename Expansion<F>::ValueType val=*iter;
    typename Expansion<F>::Reference ref=*iter;
    typename Expansion<F>::ConstReference ncref=*iter;

    typename Expansion<F>::ValueType cval=*citer;
    typename Expansion<F>::ConstReference cref=*citer;

    Bool res;

    ++iter; --iter;
    ++citer; --citer;

    res=(iter==iter); res=(iter!=iter); res=(citer==citer); res=(citer!=citer);
    res=(citer==iter); res=(citer!=iter); res=(iter==citer); res=(iter!=citer);

    ref=cref; ref=ncref;
}

// Test dereferencing of iterators
template<class F> Void TestExpansion<F>::test_data_access()
{
    Expansion<F> e(3,prec);

    // Append values
    ARIADNE_TEST_EXECUTE(e.append({0,0,0},2.0));
    ARIADNE_TEST_EXECUTE(e.append({1,0,0},3.0));
    ARIADNE_TEST_EXECUTE(e.append({0,1,0},5.0));
    ARIADNE_TEST_EXECUTE(e.append({1,0,1},7.0));

    // Test Iterator difference
    ARIADNE_TEST_PRINT(e.begin());
    ARIADNE_TEST_PRINT(e.end());
    ARIADNE_TEST_EQUAL(e.begin()+e.number_of_terms(),e.end());
    ARIADNE_TEST_EQUAL(e.end()-e.number_of_terms(),e.begin());
    ARIADNE_TEST_EQUAL(e.end()-e.begin(),PointerDifferenceType(e.number_of_terms()));

    // Test derefencing of iterators
    ARIADNE_TEST_PRINT(e.begin());
    ARIADNE_TEST_PRINT(*e.begin());
    ARIADNE_TEST_EQUAL(e.begin()->key(),MultiIndex({0,0,0}));
    ARIADNE_TEST_EQUAL(e.begin()->data(),2.0);

    Expansion<F> const& ce=e;
    ARIADNE_TEST_PRINT(ce.begin());
    ARIADNE_TEST_PRINT(*ce.begin());
    ARIADNE_TEST_EQUAL(ce.begin()->key(),MultiIndex({0,0,0}));
    ARIADNE_TEST_EQUAL(ce.begin()->data(),2.0);


    // The behaviour of iterators is rather odd and not what might be expected
    // A MultiIndex reference assigned to by iter->key() changes its value
    // when the Iterator is incremented, but a F reference does not.
    // This behaviour should be changed in future versions if technologically
    // feasible.
    Iterator iter=e.begin();
    const MultiIndex& aref=iter->key();
    const F& xref=iter->data();
    F x1=iter->data();
    ++iter;
    MultiIndex a2=iter->key();
    ARIADNE_TEST_ASSERT(a2==aref);
    ARIADNE_TEST_ASSERT(x1==xref);

    // Test finding of values of iterators
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(iter=e.find(MultiIndex({1,0,0})));
    ARIADNE_TEST_EQUAL(iter->key(),MultiIndex({1,0,0}));
    ARIADNE_TEST_EQUAL(iter->data(),3.0);



    // Test hand-coded swap of values
    ARIADNE_TEST_CONSTRUCT(Iterator,iter1,(e.begin()+1));
    ARIADNE_TEST_CONSTRUCT(Iterator,iter2,(e.begin()+3));

    // Perform swap
    ARIADNE_TEST_CONSTRUCT(typename Expansion<F>::ValueType,tmp,(*iter2));
    ARIADNE_TEST_ASSERT(tmp.key()==MultiIndex({1,0,1}));
    ARIADNE_TEST_ASSERT(tmp.data()==7.0);
    ARIADNE_TEST_EXECUTE(*iter2=*iter1);
    ARIADNE_TEST_ASSERT(iter2->key()==MultiIndex({1,0,0}));
    ARIADNE_TEST_ASSERT(iter2->data()==3.0);
    ARIADNE_TEST_EXECUTE(*iter1=tmp);
    ARIADNE_TEST_ASSERT(iter1->key()==MultiIndex({1,0,1}));
    ARIADNE_TEST_ASSERT(iter1->data()==7.0);


}

template<class F> Void TestExpansion<F>::test_equality()
{
    MultiIndex a(2);
    MultiIndex b(2); ++b;
    Expansion<F> e1(2,zero),e2(2,zero);
    e1.append(a,1.0); e1.append(b,2.0);
    e2.append(a,1.0); e2.append(b,3.0);
    ARIADNE_TEST_COMPARE(e1,!=,e2);
    e2.clear(); e2.append(a,1.0); e2.append(b,2.0);
    ARIADNE_TEST_EQUAL(e1,e2);
    e1.clear(); e1.append(b,2.0);
    e2.clear(); e2.append(a,0.0); e2.append(b,2.0);
    if(!(e1==e2)) { ARIADNE_TEST_NOTIFY("Expansion<F> objects differing by explicit zeros are considered nonequal."); }
    e1.clear(); e1.append(a,-0.0);
    e1.clear(); e1.append(a,+0.0);
    if(!(e1==e2)) { ARIADNE_TEST_NOTIFY("Expansion<F> objects differing by +0 versus -0 coefficients are considered nonequal."); }
    e1.clear(); e1.append(a,1.0); e1.append(b,2.0);
    e2.clear(); e2.append(b,2.0); e2.append(a,1.0);
    if(!(e1==e2)) { ARIADNE_TEST_NOTIFY("Expansion<F> objects differing by order of set operators are considered nonequal."); }
}

template<class F> Void TestExpansion<F>::test_cleanup()
{
    // Test to see if the cleanup/sort operations work.
    // Since these are used in the constructors, we can't use the main constructors to test this
    MultiIndex a(3);
    MultiIndex b(3); ++b;
    ARIADNE_TEST_PRINT(a);
    ARIADNE_TEST_PRINT(b);

    Expansion<F> e(3,zero);
    ARIADNE_TEST_PRINT(e);
    for(Nat i=0; i!=2; ++i) {
        if(i%2) { e.append(a,1/(1.+i)); ++b; ++b; a=b; ++b; } else { e.append(b,1/(1.+i));}
        ARIADNE_TEST_PRINT(e);
    }

    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(e.index_sort(graded_less));
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(e.combine_terms());
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(e.remove_zeros());
    ARIADNE_TEST_PRINT(e);

}

template<class F> Void TestExpansion<F>::test_constructors()
{
    // Empty initialiser list causes failure
    ARIADNE_TEST_FAIL(Expansion<F> e0({}));

    // Empty expansion
    ARIADNE_TEST_CONSTRUCT(Expansion<F>,e1,(3,prec));
    // Expansion with all entries; useful for checking ordering of indices
    //ARIADNE_TEST_CONSTRUCT(Expansion<F>,e2,(3,4, {1., 2.,3.,4., 5.,6.,7.,8.,9.,10.,
    //    11.,12.,13.,14.,15.,16.,17.,18.,19.,20., 21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35}));

    // Dense expansion
    ARIADNE_TEST_CONSTRUCT(Expansion<F>,e3,({{{2,0,0},5.},{{1,1,0},2.},{{1,0,1},0.},{{0,2,0},0.},{{0,1,1},3.},{{0,2,0},0.}}, prec));
    //ARIADNE_TEST_CONSTRUCT(Expansion<F>,e3,(3,2, {0., 0.,0.,0., 5.,2.,0.,0.,3.,0.}));
    ARIADNE_TEST_PRINT(e3);
    ARIADNE_TEST_COMPARE(e3.find(MultiIndex({2,0,0})),!=,e3.end());
    ARIADNE_TEST_EQUAL(e3[MultiIndex({2,0,0})],5.0);

    // Sparse expansion with unordered indiced
    Expansion<F> pp3({{{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0}}, prec);
    ARIADNE_TEST_CONSTRUCT(Expansion<F>,p3,({{{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0}}, prec));
    ARIADNE_TEST_EQUAL(p3[MultiIndex({1,2})],5.0);
    ARIADNE_TEST_EQUAL(p3[MultiIndex({0,0})],2.0);

    // Unordered indices
    ARIADNE_TEST_COMPARE(E({{{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0}}, prec),!=,E({{{0,0},2.0}, {{1,0},3.0}, {{0,1},11.0}, {{3,0},7.0}, {{1,2},5.0}}, prec));
    // Repeated indices; do not sum in expansion class
    ARIADNE_TEST_COMPARE(E({{{1,0},2.0}, {{1,0},3.0}, {{0,2},7.0}}, prec),!=,E({{{1,0},5.0}, {{0,2},7.0}}, prec));

    // Regression tests for expansions with only preces
    ARIADNE_TEST_CONSTRUCT(Expansion<F>,pr2,({{{},0.0}},prec));
    ARIADNE_TEST_CONSTRUCT(Expansion<F>,pr1,({{{3,0,0},0.0}},prec));

    // Regression tests for higher-order expansions
    ARIADNE_TEST_CONSTRUCT(Expansion<F>,ho1,({{{0,1,0,0,0},2.0}, {{0,1,0,0,1},3.0}, {{2,0,1,0,0},5.0}, {{0,0,0,0,0},7.0}},prec));
}



template<class F> Void TestExpansion<F>::test_find()
{
    Expansion<F> e({ {{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0} }, prec);
    MultiIndex a(2);
    a[0]=1; a[1]=2;
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_PRINT(e.find(a)-e.begin());
    ARIADNE_TEST_COMPARE(e.find(a),!=,e.end());
    ARIADNE_TEST_EQUAL(e.find(a)->key(),a);
    ARIADNE_TEST_EQUAL(e.find(a)->data(),5.0);
    a[1]=1;
    ARIADNE_TEST_EQUAL(e.find(a),e.end());

}

template<class F> Void TestExpansion<F>::test_embed()
{
    ARIADNE_TEST_CONSTRUCT(Expansion<F>,e,({ {{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0} }));
    ARIADNE_TEST_EQUAL(embed(0,e,2),Expansion<F>({ {{1,2,0,0},5.0}, {{0,0,0,0},2.0}, {{1,0,0,0},3.0}, {{3,0,0,0},7.0}, {{0,1,0,0},11.0} }));
    ARIADNE_TEST_EQUAL(embed(1,e,0),Expansion<F>({ {{0,1,2},5.0}, {{0,0,0},2.0}, {{0,1,0},3.0}, {{0,3,0},7.0}, {{0,0,1},11.0} }));
    ARIADNE_TEST_EQUAL(embed(1,e,2),Expansion<F>({ {{0,1,2,0,0},5.0}, {{0,0,0,0,0},2.0}, {{0,1,0,0,0},3.0}, {{0,3,0,0,0},7.0}, {{0,0,1,0,0},11.0} }));

}

Int main() {
    Float64 zero_64{Precision64()};
    FloatMP zero_mp{PrecisionMP(128)};
    TestExpansion<Float64>(zero_64).test();
    TestExpansion<FloatMP>(zero_mp).test();
    return ARIADNE_TEST_FAILURES;
}