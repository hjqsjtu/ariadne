/***************************************************************************
 *            test_expansion.cc
 *
 *  Copyright 2009  Pieter Collins
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


class TestExpansion
{
    typedef MultiIndex MI;
    typedef Expansion<Float64> E;
  public:
    Precision64 prec64;
    GradedLess graded_less;
    LexicographicLess lexicographic_less;
    ReverseLexicographicLess reverse_lexicographic_less;
  public:
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


Void TestExpansion::test()
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


Void TestExpansion::test_working()
{
    Expansion<Float64> e(3,prec64);
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


Void TestExpansion::test_concept()
{
    Precision64 pr;
    Float64 z(pr);
    Float64 x(5,pr);
    SizeType as(3);

    Expansion<Float64> e(as,pr);
    const Expansion<Float64> ce(as,pr);

    e=Expansion<Float64>(as,z);
    e=Expansion<Float64>(ce);

    //e=Expansion<Float64>(3,1, {0.0, 0.0,0.0,0.0}, pr);
    //e=Expansion<Float64>(3,1, {1, 2,3,5.0}, pr);
    e=Expansion<Float64>({ {{0,0},1}, {{1,0,0},2}, {{0,1,0},3}, {{0,0,1},5.0} }, pr);
    e=Expansion<Float64>({ {{0,0},1}, {{1,0,0},2}, {{0,1,0},3}, {{0,0,1},5.0} }, pr);

    MultiIndex a(as);
    e.reserve(2u);
    e.set(a,x);
    e.prepend(a,x);
    e.append(a,x);
    e.append_sum(a,a,x);
    e.clear();

    e.index_sort(GradedLess());
    e.sort(ReverseLexicographicIndexLess());

    x=ce[a];

    ce.number_of_terms();
    ce.argument_size();

    e.erase(e.begin());

    ce.check();
}

Void TestExpansion::test_iterator_concept()
{
    MI a(3);
    Expansion<Float64> e(3,prec64);
    const Expansion<Float64> cp(3,prec64);

    Expansion<Float64>::Iterator iter=e.begin(); iter=e.end(); iter=e.find(a);
    Expansion<Float64>::ConstIterator citer=e.begin(); citer=e.end(); citer=e.find(a);
    citer=e.begin(); citer=cp.end(); citer=cp.find(a);

    Expansion<Float64>::ValueType val=*iter;
    Expansion<Float64>::Reference ref=*iter;
    Expansion<Float64>::ConstReference ncref=*iter;

    Expansion<Float64>::ValueType cval=*citer;
    Expansion<Float64>::ConstReference cref=*citer;

    Bool res;

    ++iter; --iter;
    ++citer; --citer;

    res=(iter==iter); res=(iter!=iter); res=(citer==citer); res=(citer!=citer);
    res=(citer==iter); res=(citer!=iter); res=(iter==citer); res=(iter!=citer);

    ref=cref; ref=ncref;
}

// Test dereferencing of iterators
Void TestExpansion::test_data_access()
{
    Expansion<Float64> e(3,Precision64());

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



    // The behaviour of iterators is rather odd and not what might be expected
    // A MultiIndex reference assigned to by iter->key() changes its value
    // when the Iterator is incremented, but a Float64 reference does not.
    // This behaviour should be changed in future versions if technologically
    // feasible.
    Expansion<Float64>::Iterator iter=e.begin();
    const MultiIndex& aref=iter->key();
    const Float64& xref=iter->data();
    Float64 x1=iter->data();
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
    ARIADNE_TEST_CONSTRUCT(Expansion<Float64>::Iterator,iter1,(e.begin()+1));
    ARIADNE_TEST_CONSTRUCT(Expansion<Float64>::Iterator,iter2,(e.begin()+3));

    // Perform swap
    ARIADNE_TEST_CONSTRUCT(Expansion<Float64>::ValueType,tmp,(*iter2));
    ARIADNE_TEST_ASSERT(tmp.key()==MultiIndex({1,0,1}));
    ARIADNE_TEST_ASSERT(tmp.data()==7.0);
    ARIADNE_TEST_EXECUTE(*iter2=*iter1);
    ARIADNE_TEST_ASSERT(iter2->key()==MultiIndex({1,0,0}));
    ARIADNE_TEST_ASSERT(iter2->data()==3.0);
    ARIADNE_TEST_EXECUTE(*iter1=tmp);
    ARIADNE_TEST_ASSERT(iter1->key()==MultiIndex({1,0,1}));
    ARIADNE_TEST_ASSERT(iter1->data()==7.0);


}

Void TestExpansion::test_equality()
{
    MI a(2);
    MI b(2); ++b;
    Expansion<Float64> e1(2,prec64),e2(2,prec64);
    e1.append(a,1.0); e1.append(b,2.0);
    e2.append(a,1.0); e2.append(b,3.0);
    ARIADNE_TEST_COMPARE(e1,!=,e2);
    e2.clear(); e2.append(a,1.0); e2.append(b,2.0);
    ARIADNE_TEST_EQUAL(e1,e2);
    e1.clear(); e1.append(b,2.0);
    e2.clear(); e2.append(a,0.0); e2.append(b,2.0);
    if(!(e1==e2)) { ARIADNE_TEST_NOTIFY("Expansion<Float64> objects differing by explicit zeros are considered nonequal."); }
    e1.clear(); e1.append(a,-0.0);
    e1.clear(); e1.append(a,+0.0);
    if(!(e1==e2)) { ARIADNE_TEST_NOTIFY("Expansion<Float64> objects differing by +0 versus -0 coefficients are considered nonequal."); }
    e1.clear(); e1.append(a,1.0); e1.append(b,2.0);
    e2.clear(); e2.append(b,2.0); e2.append(a,1.0);
    if(!(e1==e2)) { ARIADNE_TEST_NOTIFY("Expansion<Float64> objects differing by order of set operators are considered nonequal."); }
}


Void TestExpansion::test_cleanup()
{
    // Test to see if the cleanup/sort operations work.
    // Since these are used in the constructors, we can't use the main constructors to test this
    MI a(3);
    MI b(3); ++b;
    ARIADNE_TEST_PRINT(a);
    ARIADNE_TEST_PRINT(b);

    Expansion<Float64> e(3,prec64);
    ARIADNE_TEST_PRINT(e);
    for(Nat i=0; i!=2; ++i) {
        if(i%2) { e.append(a,1/(1.+i)); ++b; ++b; a=b; ++b; } else { e.append(b,1/(1.+i));}
        ARIADNE_TEST_PRINT(e);
    }

    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(e.index_sort(graded_less));
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EXECUTE(e.remove_zeros());
    ARIADNE_TEST_PRINT(e);

}

Void TestExpansion::test_constructors()
{
    // Empty expansion
    ARIADNE_TEST_CONSTRUCT(Expansion<Float64>,e1,(3,prec64));
    // Expansion with all entries; useful for checking ordering of indices
    //ARIADNE_TEST_CONSTRUCT(Expansion<Float64>,e2,(3,4, {1., 2.,3.,4., 5.,6.,7.,8.,9.,10.,
    //    11.,12.,13.,14.,15.,16.,17.,18.,19.,20., 21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35}));

    // Dense expansion
    ARIADNE_TEST_CONSTRUCT(Expansion<Float64>,e3,({{{2,0,0},5.},{{1,1,0},2.},{{1,0,1},0.},{{0,2,0},0.},{{0,1,1},3.},{{0,2,0},0.}}, prec64));
    //ARIADNE_TEST_CONSTRUCT(Expansion<Float64>,e3,(3,2, {0., 0.,0.,0., 5.,2.,0.,0.,3.,0.}));
    ARIADNE_TEST_PRINT(e3);
    ARIADNE_TEST_COMPARE(e3.find(MultiIndex({2,0,0})),!=,e3.end());
    ARIADNE_TEST_EQUAL(e3[MultiIndex({2,0,0})],5.0);

    // Sparse expansion with unordered indiced
    Expansion<Float64> pp3({{{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0}}, prec64);
    ARIADNE_TEST_CONSTRUCT(Expansion<Float64>,p3,({{{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0}}, prec64));
    ARIADNE_TEST_EQUAL(p3[MultiIndex({1,2})],5.0);
    ARIADNE_TEST_EQUAL(p3[MultiIndex({0,0})],2.0);

    // Unordered indices
    ARIADNE_TEST_COMPARE(E({{{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0}}, prec64),!=,E({{{0,0},2.0}, {{1,0},3.0}, {{0,1},11.0}, {{3,0},7.0}, {{1,2},5.0}}, prec64));
    // Repeated indices; do not sum in expansion class
    ARIADNE_TEST_COMPARE(E({{{1,0},2.0}, {{1,0},3.0}, {{0,2},7.0}}, prec64),!=,E({{{1,0},5.0}, {{0,2},7.0}}, prec64));

    // Regression tests for expansions with only zeroes
    ARIADNE_TEST_CONSTRUCT(Expansion<Float64>,pr2,({{{},0.0}},prec64));
    ARIADNE_TEST_CONSTRUCT(Expansion<Float64>,pr1,({{{3,0,0},0.0}},prec64));

    // Regression tests for higher-order expansions
    ARIADNE_TEST_CONSTRUCT(Expansion<Float64>,ho1,({{{0,1,0,0,0},2.0}, {{0,1,0,0,1},3.0}, {{2,0,1,0,0},5.0}, {{0,0,0,0,0},7.0}},prec64));
}




/* Not needed since we cannot look up by index without order
Void TestExpansion::test_indexing()
{
    Expansion<Float64> e(3,4, 0,0,0,2.0,  1,0,0,3.0, 1,0,1,5.0, 2,1,0,7.0);
    const Expansion<Float64>& pc=e;
    ARIADNE_TEST_EQUAL(e[MI(3, 1,0,0)],3.0);

    e[MI(3, 1,0,0)]-=0.5;
    ARIADNE_TEST_EQUAL(e[MI(3, 1,0,0)],2.5);

    e[MI(3, 1,1,0)]=11.0;
    ARIADNE_TEST_EQUAL(e[MI(3, 1,1,0)],11.0);

    e.clear();
    e[MI(3, 0,0,0)]=2.0;
    e[MI(3, 0,1,0)]=3.0;
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EQUAL(e.number_of_terms(),2);
    ARIADNE_TEST_EQUAL(e[MI(3, 0,0,0)],2.0);
    ARIADNE_TEST_EQUAL(e[MI(3, 0,1,0)],3.0);
    ARIADNE_TEST_EXECUTE(e[MI(3, 1,0,0)]=5.0);
    ARIADNE_TEST_EQUAL(e[MI(3, 1,0,0)],5.0);
    ARIADNE_TEST_EQUAL(e[MI(3, 0,1,0)],3.0);
    ARIADNE_TEST_EXECUTE(e[MI(3, 0,0,1)]=7.0);
    ARIADNE_TEST_EQUAL(e[MI(3, 0,0,1)],7.0);

    // Test insert at beginning
    e.clear();
    e[MI(3, 0,1,0)]=2.0;
    e[MI(3, 0,0,1)]=3.0;
    e[MI(3, 2,1,0)]=5.0;
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EQUAL(pc[MI(3, 0,0,0)],0.0);
    ARIADNE_TEST_EQUAL(e.number_of_terms(),3);
    ARIADNE_TEST_EXECUTE(e[MI(3, 0,0,0)]=7);
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_EQUAL(e.number_of_terms(),4);
    Expansion<Float64>::ConstIterator iter=e.begin();
    ARIADNE_TEST_EQUAL(iter->key(),MI(3, 0,1,0));
    ARIADNE_TEST_EQUAL(iter->data(),2.0);
    ++iter;
    ARIADNE_TEST_EQUAL(iter->key(),MI(3, 0,0,1));
    ARIADNE_TEST_EQUAL(iter->data(),3.0);
    ++iter;
    ARIADNE_TEST_EQUAL(iter->key(),MI(3, 2,1,0));
    ARIADNE_TEST_EQUAL(iter->data(),5.0);
    ++iter;
    ARIADNE_TEST_EQUAL(iter->key(),MI(3, 0,0,0));
    ARIADNE_TEST_EQUAL(iter->data(),7.0);
}
*/


Void TestExpansion::test_find()
{
    Expansion<Float64> e({ {{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0} });
    MI a(2);
    a[0]=1; a[1]=2;
    ARIADNE_TEST_PRINT(e);
    ARIADNE_TEST_PRINT(e.find(a)-e.begin());
    ARIADNE_TEST_COMPARE(e.find(a),!=,e.end());
    ARIADNE_TEST_EQUAL(e.find(a)->key(),a);
    ARIADNE_TEST_EQUAL(e.find(a)->data(),5.0);
    a[1]=1;
    ARIADNE_TEST_EQUAL(e.find(a),e.end());

}

Void TestExpansion::test_embed()
{
    ARIADNE_TEST_CONSTRUCT(Expansion<Float64>,e,({ {{1,2},5.0}, {{0,0},2.0}, {{1,0},3.0}, {{3,0},7.0}, {{0,1},11.0} }));
    ARIADNE_TEST_EQUAL(embed(0,e,2),Expansion<Float64>({ {{1,2,0,0},5.0}, {{0,0,0,0},2.0}, {{1,0,0,0},3.0}, {{3,0,0,0},7.0}, {{0,1,0,0},11.0} }));
    ARIADNE_TEST_EQUAL(embed(1,e,0),Expansion<Float64>({ {{0,1,2},5.0}, {{0,0,0},2.0}, {{0,1,0},3.0}, {{0,3,0},7.0}, {{0,0,1},11.0} }));
    ARIADNE_TEST_EQUAL(embed(1,e,2),Expansion<Float64>({ {{0,1,2,0,0},5.0}, {{0,0,0,0,0},2.0}, {{0,1,0,0,0},3.0}, {{0,3,0,0,0},7.0}, {{0,0,1,0,0},11.0} }));

}

Int main() {
    TestExpansion().test();
    return ARIADNE_TEST_FAILURES;
}
