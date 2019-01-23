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

// Usage: makechain [options] > <output file>
//
// This code creates a chain structured network and saves
// it into a file. 
//
// The intensities of leaving each state for the head node
// are the same, which can be set using the parameter "headrate".
// When leaving the current state, the head node has equal
// probability to jump to other states.
//
// For the other nodes, each node stays in its current state 
// if it matches its parents and otherwise transitions 
// to its parent’s state with a high probability (set using
// parameter "chaserate".) 

/* Command line parameters:
 *
 * Set the number of nodes
 * -Dnumofnodes <integer> (default: 3)
 *
 * Set the number of states for each node
 * -Dnumofstates <integer> (default: 2)
 *
 * Set the intensity of leaving each state for the head node.
 * -Dheadrate <double> (default: 1.0)
 *
 * Set the intensity of transitioning to the state that is the same 
 * as the current state of the parent for nodes other than the head node
 * -Dchaserate <double> (default: 1.0)
 *
 * Set the intensity of transitioning to other states 
 * -Drandrate <double> (default: 0.1)
 */

#include "markovdyn.h"
#include "markov.h"
#include "streamextra.h"
#include "context.h"
#include <iostream>
#include "multirv.h"
#include "trajectory.h"
#include "params.h"
#include "bn.h"
#include "ctbndyn.h"
#include <iomanip>

using namespace std;
using namespace ctbn;

int main(int argc, char **argv) {
	InitParams(argc, argv);
	cout << setprecision(4);

	// Set various parameters including the number of nodes, states per node,
	// initial distribution values, and Q matrix values.
	string paramval;
	int nnode = ParamInt("numofnodes", 3);
	int nstate = ParamInt("numofstates", 2);
	double headrate = ParamDouble("headrate", 1.0);
	double chaserate  = ParamDouble("chaserate", 1.0);
	double randrate  = ParamDouble("randrate", 0.1 );

	// Generate the head node
	Context nullc;
	Context headc;
	headc.AddVar(0, nstate);
	MarkovDyn Phead(headc, nullc);
	MultiRV P0head(headc, nullc);
	vectr p0(nstate, 0.0);
	for (int i=0; i<nstate; i++) {
		for (int j=0; j<nstate; j++)
			Phead(0)->Intensity()[i][j] = (i==j ? -headrate : headrate/(nstate-1));
		p0[i] = (i==0 ? 1.0 : 0.0);
	}
	P0head[0].SetDist(p0);
	CTBNDyn ctbndyn;
	ctbndyn.AddNode(Phead.Clone());
	BN bn;
	bn.AddNode(P0head.Clone());

	// Generate the rest of the nodes
	for (int i=1; i<nnode; i++) {
		Context pcontext;
		Context nodecontext;
		pcontext.AddVar(i-1, nstate);
		nodecontext.AddVar(i, nstate);
		MarkovDyn node(nodecontext, pcontext);
		for (int j=0; j<nstate; j++)
			for (int k=0; k<nstate; k++)
				for (int l=0; l<nstate; l++)
					node(j)->Intensity()[k][l] = (j==k ?
							((k==l) ? -randrate*(nstate-1) : randrate)
							: ((k==l) ? -chaserate-randrate*(nstate-2) :
								(l==j ? chaserate : randrate)));
		ctbndyn.AddNode(node.Clone());

		Context initp;
		initp.AddVar(i, nstate);
		MultiRV P0node(initp, nullc);
		P0node[0].SetDist(p0);
		bn.AddNode(P0node.Clone());
	}

	// Create the final object to represent the chain network
	Markov m(bn.Clone(), ctbndyn.Clone());

	m.Save(cout);
}
