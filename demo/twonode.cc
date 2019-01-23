/* Continuous Time Bayesian Network Reasoning and Learning Engine
 * Copyright (C) 2009 The Regents of the University of California
 *
 * see docs/AUTHORS for contributor list
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "ctbn.h"
#include "markovdyn.h"
#include "bn.h"
#include "multirv.h"
#include <iostream>

using namespace std;

using namespace ctbn;

int main(int argc, char **argv) {
	Context nullc; // set of no variables
	Context var0,var1; // two variable sets, each with a single variable:
	var0.AddVar(0,2); // add variable 0 with 2 states
	var1.AddVar(1,2); // add variable 1 with 2 states

	// info for variable 0:
	MultiRV init0(var0,nullc); // the initial distribution, no parents
	vectr p0(2); // the initial distribution's parameters
	p0[0] = 0.5; // this variable's initial distribution is uniform
	p0[1] = 0.5;
	init0[0].SetDist(p0); // set the initial distribution ([0] b/c it has
	// no parents and so there is only 1 distribution)

	MarkovDyn dyn0(var0,var1); // the CIM has "var0" as its variables 
	// with "var1" variables as parents
	// when var1=0, Q = [-1 1 ; 2 -2]
	dyn0(0)->Intensity()[0][0] = -1;
	dyn0(0)->Intensity()[0][1] = 1;
	dyn0(0)->Intensity()[1][0] = 2;
	dyn0(0)->Intensity()[1][1] = -2;
	// when var1=1, Q = [-10 10 ; 20 -20]
	dyn0(1)->Intensity()[0][0] = -10;
	dyn0(1)->Intensity()[0][1] = 10;
	dyn0(1)->Intensity()[1][0] = 20;
	dyn0(1)->Intensity()[1][1] = -20;

	// info for variable 1:
	MultiRV init1(var1,var0); // the initial distribution, 0 as parent
	vectr p1(2); 
	p1[0] = 0.9; // when var0=0, var1 tends to start in state 0
	p1[1] = 0.1;
	init1[0].SetDist(p1); // set the distribution
	p1[0] = 0.1; // when var0=1, var1 tends to start in state 1
	p1[1] = 0.9;
	init1[1].SetDist(p1); // set the distribution

	MarkovDyn dyn1(var1,var0); // the CIM has "var0" as its variables 
	// with "var1" variables as parents
	// when var0=0, Q = [-5 5; 6 -6]
	dyn1(0)->Intensity()[0][0] = -5;
	dyn1(0)->Intensity()[0][1] = 5;
	dyn1(0)->Intensity()[1][0] = 6;
	dyn1(0)->Intensity()[1][1] = -6;
	// when var0=1, Q = [-7 7; 8 -8]
	dyn1(1)->Intensity()[0][0] = -7;
	dyn1(1)->Intensity()[0][1] = 7;
	dyn1(1)->Intensity()[1][0] = 8;
	dyn1(1)->Intensity()[1][1] = -8;

	// now place into CTBN:
	CTBNDyn ctbndyn;
	ctbndyn.AddNode(dyn0.Clone());
	ctbndyn.AddNode(dyn1.Clone());
	BN bn;
	bn.AddNode(init0.Clone());
	bn.AddNode(init1.Clone());

	CTBN ctbn(bn.Clone(),ctbndyn.Clone());

	ctbn.Save(cout);
}
