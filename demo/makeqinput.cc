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
#include "markovdyn.h"
#include "streamextra.h"
#include "context.h"
#include "multirv.h"
#include "params.h"
#include "bn.h"
#include "ctbndyn.h"
#include "em.h"
#include "utils.h"
#include <iostream>
#include "markov.h"

using namespace std;
using namespace ctbn;


// this file is used to create an example CTBN model and a partial trajectory sampled from the model
// the model, the evidence, and the true results of queries will be saved to a file "queryinput.data"
// "queryinput.data" will be feeded as an input to:
// exactquery.cc, importancequery.cc and gibbsquery.cc
// the above three .cc file will do the queries, and compare the results with the true answers in "queryinput.data"
// the queries will be specified in these three files

// In case you want to change the example model, you can modify this file and rewrite the output to "queryinput.data"
// and then do the exactquery.cc/importancequery.cc/gibbsquery.cc 

// Usage:
// no inputs to run: ./makequeryinput >output_file
int main (int argc, char *argv) {
	randomizer.Reset(1);
	// build up the context
	Context X,Y,Z;
	Context Null;
	X.AddVar(0,2);
	Y.AddVar(1,2);
	Z.AddVar(2,3);

	// build up the dynamics
	MarkovDyn PX(X,Null); // process of X
	MarkovDyn PY(Y,X); // process of Y given X
	MarkovDyn PZ(Z,Null);
	
	PX(0)->Intensity()[0][0] = -1;
	PX(0)->Intensity()[0][1] = 1;
	PX(0)->Intensity()[1][0] = 2;
	PX(0)->Intensity()[1][1] = -2;


	PY(0)->Intensity()[0][0] = -3;
	PY(0)->Intensity()[0][1] = 3;
	PY(0)->Intensity()[1][0] = 4;
	PY(0)->Intensity()[1][1] = -4;

	PY(1)->Intensity()[0][0] = -5;
	PY(1)->Intensity()[0][1] = 5;
	PY(1)->Intensity()[1][0] = 6;
	PY(1)->Intensity()[1][1] = -6;


	PZ(0)->Intensity()[0][0] = -6;
	PZ(0)->Intensity()[0][1] = 3;
	PZ(0)->Intensity()[0][2] = 3;
	PZ(0)->Intensity()[1][1] = -12;
	PZ(0)->Intensity()[1][0] = 6;
	PZ(0)->Intensity()[1][2] = 6;
	PZ(0)->Intensity()[2][2] = -18;
	PZ(0)->Intensity()[2][0] = 9;
	PZ(0)->Intensity()[2][1] = 9;

	// build up the start distribution
	MultiRV P0X(X,Null);
	MultiRV P0Y(Y,X);
	MultiRV P0Z(Z,Null);

	vectr p0(2,0.0);
	p0[0] = 0.5;  p0[1] = 0.5;
	P0X[0].SetDist(p0);

	p0[0] = 0.8; p0[1] = 0.2;
	P0Y[0].SetDist(p0);

	p0[0] = 0.3; p0[1] = 0.7;
	P0Y[1].SetDist(p0);

	vectr pz(3,0.0);
	pz[0] = 0.2; pz[1] = 0.3; pz[2] = 0.5;
	P0Z[0].SetDist(pz);

	// set up CTBN:
	CTBNDyn ctbndyn;
	ctbndyn.AddNode(PX.Clone());
	ctbndyn.AddNode(PY.Clone());
	ctbndyn.AddNode(PZ.Clone());

	BN bn;
	bn.AddNode(P0X.Clone());
	bn.AddNode(P0Y.Clone());
	bn.AddNode(P0Z.Clone());

	Markov ctbn(bn.Clone(),ctbndyn.Clone());

	Context context = ctbndyn.Domain() + ctbndyn.CondDomain();

	// sample a partial trajectory as the input evidence
	double begintime = ParamDouble("BeginTime", 0.0);
	double endtime = ParamDouble("EndTime", 20.0);
	Trajectory evid;
	evid.SetBeginTime(begintime);
	evid.SetEndTime(endtime);
	ctbn.Sample(evid);
	RemoveNodesInformation(evid, context, ctbndyn.NumofNodes(), ctbndyn.NumofNodes()*4 , 0.2);

	// output the model
	ctbn.Save(cout);
	cout << endl;
	evid.Save(cout);
	cout << endl;
	// output the expected query results (exact):
	cout << 8.1538 << " " << 42.5293 << endl;
	return 0;
}
