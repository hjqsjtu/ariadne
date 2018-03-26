/***************************************************************************
 *            valve-permissive.hpp
 *
 *  Copyright  2017 Luca Geretti
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

#include <cstdarg>
#include "ariadne.hpp"

using namespace Ariadne;

AtomicHybridAutomaton getController2()
{

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant hmin("hmin",5.75_decimal);
    RealConstant hmax("hmax",7.75_decimal);
    RealConstant delta("delta",0.002_decimal);

    // Declare the shared system variables
    RealVariable height("height2");

    // Declare the events we use
    DiscreteEvent e_can_open("open2");
    DiscreteEvent e_can_close("close2");
    DiscreteEvent e_must_open("mustopen2");
    DiscreteEvent e_must_close("mustclose2");

    AtomicHybridAutomaton controller("controller2");

    // Declare the values the valve can variable can have
    AtomicDiscreteLocation rising("rising2");
    AtomicDiscreteLocation falling("falling2");

    controller.new_mode(rising,List<RealAssignment>());
    controller.new_mode(falling,List<RealAssignment>());

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    controller.new_invariant(falling,height<=hmin-delta,e_must_open);
    controller.new_invariant(rising,height>=hmax+delta,e_must_close);

    controller.new_transition(falling,e_can_open,rising,height<=hmin+delta,permissive);
    controller.new_transition(rising,e_can_close,falling,height>=hmax-delta,permissive);

    return controller;
}
