/***************************************************************************
 *            grid_submodule.cc
 *
 *  Copyright 2008--17  Pieter Collins
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


#include "boost_python.h"
#include "utilities.h"

#include <boost/python.hpp>

#include "numeric/numeric.h"
#include "geometry/function_set.h"
#include "geometry/grid_set.h"
#include "geometry/list_set.h"

using namespace boost::python;
using namespace Ariadne;



// For some reason, Boost::Python doesn't like returning a const GridCell& from an Iterator dereference. We therefore
// wrap the usual GridTreeConstIterator in a class which returns a value upon dereferenceing, and use the range
// Iterator functionality.
// I have tried several approaches to the return value, and none of them seem to work.
class PyGridTreeSetConstIterator : public GridTreeConstIterator {
  public:
    typedef GridCell value_type;
    typedef const GridCell reference;
    typedef const GridCell* pointer;

    PyGridTreeSetConstIterator(GridTreeConstIterator iter) : GridTreeConstIterator(iter) { }
    const GridCell operator*() const { return this->GridTreeConstIterator::operator*(); }
};

PyGridTreeSetConstIterator py_cells_begin(const GridTreeSubset& gts) {
    return PyGridTreeSetConstIterator(gts.begin()); }
PyGridTreeSetConstIterator py_cells_end(const GridTreeSubset& gts) {
    return PyGridTreeSetConstIterator(gts.end()); }


// Iterator through the boxes of a GridTreeSet set
class PyGridTreeSetBoxConstIterator
: public boost::iterator_facade< GridTreeConstIterator, ExactBoxType const, boost::forward_traversal_tag > {
  public:
    //typedef boost::forward_traversal_tag iterator_category;
    typedef ExactBoxType value_type;
    typedef ExactBoxType value;
    typedef SizeType difference_type;
    typedef const ExactBoxType reference;
    typedef const ExactBoxType* pointer;

    PyGridTreeSetBoxConstIterator(PyGridTreeSetConstIterator iter) : _iter(iter) { }
    const ExactBoxType operator*() const { return this->_iter->box(); }
    PyGridTreeSetBoxConstIterator& operator++() { ++this->_iter; return *this; }
    PyGridTreeSetBoxConstIterator operator++(Int) { PyGridTreeSetBoxConstIterator ret(*this); ++this->_iter; return ret; }
    Bool operator==(const PyGridTreeSetBoxConstIterator& other) { return this->_iter==other._iter; }
  private:
    GridTreeSet::ConstIterator _iter;
};

PyGridTreeSetBoxConstIterator py_boxes_begin(const GridTreeSubset& gts) {
    return PyGridTreeSetBoxConstIterator(gts.begin()); }
PyGridTreeSetBoxConstIterator py_boxes_end(const GridTreeSubset& gts) {
    return PyGridTreeSetBoxConstIterator(gts.end()); }







Void export_grid()
{
    typedef Vector<Float64> RVector;
    typedef Vector<ExactIntervalType> IVector;

    class_< Grid > grid_class("Grid",no_init);
    grid_class.def(init<Nat>());
    grid_class.def(init<Nat,Float64>());
    grid_class.def(init< Vector<Float64>, Vector<Float64> >());
    grid_class.def("dimension", &Grid::dimension);
    grid_class.def("origin", &Grid::origin, return_value_policy<copy_const_reference>());
    grid_class.def("lengths", &Grid::lengths, return_value_policy<copy_const_reference>());
    grid_class.def(self_ns::str(self));
}


Void export_grid_cell()
{

    class_<GridCell> grid_cell_class("GridCell",no_init);
    grid_cell_class.def("dimension", &GridCell::dimension);
    grid_cell_class.def("depth", &GridCell::depth);
    grid_cell_class.def("split", (Pair<GridCell,GridCell>(GridCell::*)()const) &GridCell::split);
    grid_cell_class.def("split", (GridCell(GridCell::*)(Bool)const) &GridCell::split);
    grid_cell_class.def("box", &GridCell::box, return_value_policy<copy_const_reference>());
    grid_cell_class.def(self_ns::str(self));

    def("smallest_enclosing_primary_cell", &GridCell::smallest_enclosing_primary_cell);

    to_python< std::pair<GridCell,GridCell> >();
}


Void export_grid_tree_set() {

    class_<GridTreeSubset> grid_tree_subset_class("GridTreeSubset",no_init);

    class_<GridTreeSet, bases<GridTreeSubset,DrawableInterface> > grid_tree_set_class("GridTreeSet",init<GridTreeSet>());
    grid_tree_set_class.def(init<Nat>());
    grid_tree_set_class.def(init<Grid>());
    grid_tree_set_class.def("bounding_box", &GridTreeSet::bounding_box);
    grid_tree_set_class.def("is_empty", &GridTreeSet::is_empty);
    grid_tree_set_class.def("size", &GridTreeSet::size);
    grid_tree_set_class.def("dimension", &GridTreeSet::dimension);
    grid_tree_set_class.def("clear", &GridTreeSet::clear);
    grid_tree_set_class.def("mince", &GridTreeSet::mince);
    grid_tree_set_class.def("recombine", &GridTreeSet::recombine);
    grid_tree_set_class.def("grid", &GridTreeSet::grid,return_value_policy<copy_const_reference>());
    grid_tree_set_class.def("measure", &GridTreeSet::measure);
    grid_tree_set_class.def("adjoin", (Void(GridTreeSet::*)(const GridCell&))(&GridTreeSet::adjoin));
    grid_tree_set_class.def("adjoin", (Void(GridTreeSet::*)(const GridTreeSubset&))(&GridTreeSet::adjoin));
    grid_tree_set_class.def("restrict", (Void(GridTreeSet::*)(const GridTreeSubset&))(&GridTreeSet::restrict));
    grid_tree_set_class.def("remove", (Void(GridTreeSet::*)(const GridTreeSubset&))(&GridTreeSet::remove));
    grid_tree_set_class.def("adjoin_over_approximation", (Void(GridTreeSet::*)(const ExactBoxType&,const Nat)) &GridTreeSet::adjoin_over_approximation);
    grid_tree_set_class.def("adjoin_outer_approximation", (Void(GridTreeSet::*)(const CompactSetInterface&,const Nat)) &GridTreeSet::adjoin_outer_approximation);
    grid_tree_set_class.def("adjoin_inner_approximation", (Void(GridTreeSet::*)(const OpenSetInterface&,const Nat,const Nat)) &GridTreeSet::adjoin_inner_approximation);
    grid_tree_set_class.def("__len__", &GridTreeSet::size);
    //grid_tree_set_class.def("__iter__", boost::python::Iterator<const GridTreeSet>());
    grid_tree_set_class.def("__iter__", range(&py_cells_begin,&py_cells_end));
    grid_tree_set_class.def("cells", range(&py_cells_begin,&py_cells_end));
    grid_tree_set_class.def("boxes", range(&py_boxes_begin,&py_boxes_end));
    grid_tree_set_class.def(self_ns::str(self));

    def("union",(GridTreeSet(*)(const GridTreeSubset&,const GridTreeSubset&))(&join));
    def("difference",(GridTreeSet(*)(const GridTreeSubset&,const GridTreeSubset&))(&difference));
    def("intersection",(GridTreeSet(*)(const GridTreeSubset&,const GridTreeSubset&))(&intersection));
    def("intersect",(Bool(*)(const GridTreeSubset&,const GridTreeSubset&))(&intersect));
    def("subset",(Bool(*)(const GridTreeSubset&,const GridTreeSubset&))(&subset));
    def("intersect",(Bool(*)(const GridCell&,const GridTreeSubset&))(&intersect));
    def("subset",(Bool(*)(const GridCell&,const GridTreeSubset&))(&subset));

    def("outer_approximation",(GridTreeSet(*)(const CompactSetInterface&,const Grid&,const Nat)) &outer_approximation);
    def("inner_approximation",(GridTreeSet(*)(const OpenSetInterface&,const Grid&,const Nat,const Nat)) &inner_approximation);

}


Void storage_submodule()
{
    export_grid();
    export_grid_cell();
    export_grid_tree_set();
}
