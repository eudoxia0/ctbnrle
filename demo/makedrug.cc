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

// This code creates the drug effect network as first seen in 
// Continuous Time Bayesian Networks by Nodelman et al., 2002 and
// save it to stdout 
//
// Usage: makedrug [options] > <filename>
// options:
// 
//  Set whether to save only the dynamics (as opposed to the full CTBN):
//  -Dsavedyn <integer> (default: 0 [=false])
//
//  Set whether to save as a pointer (which can be loaded into any superclass
//  -Dsaveptr <integer> (default: 0 [=false])


#include "bn.h"
#include "em.h"
#include "context.h"
#include "ctbn.h"
#include "ctbndyn.h"
#include "factoredmatrix.h"
#include "factoredvector.h"
#include "markov.h"
#include "markovdyn.h"
#include "multirv.h"
#include "params.h"
#include "streamextra.h"
#include "trajectory.h"
#include "utils.h"

#include <iomanip>
#include <iostream>


#include <iomanip>
#include <iostream>

using namespace std;
using namespace ctbn;


int main(int argc, char **argv) {

	InitParams(argc, argv);
	bool savedyn = ParamInt("savedyn",false);
	bool saveasptr = ParamInt("saveptr",false);
	
	cout << setprecision(4);

	//*******************************************************

	// Create contexts for each variable
	Context nullc;
	Context eatc, fullc, hungryc, uptakec, concentrationc, barometerc, painc, drowsyc;
	eatc.AddVar(0, 2);
	fullc.AddVar(1,3);
	hungryc.AddVar(2,2);
	uptakec.AddVar(3, 2);
	concentrationc.AddVar(4,3);
	barometerc.AddVar(5,3);
	painc.AddVar(6, 2);
	drowsyc.AddVar(7,2);

	// Define the connections between variables. For a variable with multiple 
	// parents, the context used is the union of the two parent contexts
	MarkovDyn eat(eatc, hungryc);
	MarkovDyn full(fullc, eatc);
	MarkovDyn hungry(hungryc, fullc);
	MarkovDyn uptake(uptakec, nullc);
	Context full_and_uptakec(fullc, uptakec);
	MarkovDyn concentration(concentrationc, full_and_uptakec);
	MarkovDyn barometer(barometerc, nullc);
	MarkovDyn pain(painc, Context(barometerc, concentrationc));
	MarkovDyn drowsy(drowsyc, concentrationc);

	// Generate starting distributions for each variable
	MultiRV eatrv(eatc, nullc);
	MultiRV fullrv(fullc, nullc);
	MultiRV hungryrv(hungryc, nullc);
	MultiRV uptakerv(uptakec, nullc);
	MultiRV concentrationrv(concentrationc, nullc);
	MultiRV barometerrv(barometerc, nullc);
	MultiRV painrv(painc, nullc);
	MultiRV drowsyrv(drowsyc, nullc);

	// The following fills the conditional intensity matrices for each 
	// variable and the initial distributions (multinomial distributions).

	//------eat--------
	eat(0)->Intensity()[0][0] = -0.01;
	eat(0)->Intensity()[0][1] = 0.01;
	eat(0)->Intensity()[1][0] = 10.0;
	eat(0)->Intensity()[1][1] = -10.0;
	eat(1)->Intensity()[0][0] = -2.0;
	eat(1)->Intensity()[0][1] = 2.0;
	eat(1)->Intensity()[1][0] = 1.0;
	eat(1)->Intensity()[1][1] = -1.0;
	vectr peat(2, 0.0);
	peat[0] = 0.1;
	peat[1] = 0.9;
	eatrv[0].SetDist(peat);

	//------full--------
	full(0)->Intensity()[0][0] = -0.02;
	full(0)->Intensity()[0][1] = 0.01;
	full(0)->Intensity()[0][2] = 0.01;
	full(0)->Intensity()[1][0] = 0.3;
	full(0)->Intensity()[1][1] = -0.31;
	full(0)->Intensity()[1][2] = 0.01;
	full(0)->Intensity()[2][0] = 0.01;
	full(0)->Intensity()[2][1] = 1.0;
	full(0)->Intensity()[2][2] = -1.01;

	full(1)->Intensity()[0][0] = -6.01;
	full(1)->Intensity()[0][1] = 6.0;
	full(1)->Intensity()[0][2] = 0.01;
	full(1)->Intensity()[1][0] = 0.01;
	full(1)->Intensity()[1][1] = -3.01;
	full(1)->Intensity()[1][2] =  3.0;
	full(1)->Intensity()[2][0] = 0.01;
	full(1)->Intensity()[2][1] = 0.01;
	full(1)->Intensity()[2][2] = -0.02;

	vectr pfull(3, 0.0);
	pfull[0] = 0.7;
	pfull[1] = 0.2;
	pfull[2] = 0.1;
	fullrv[0].SetDist(pfull);

	//------hungry--------
	hungry(0)->Intensity()[0][0] = -2.0;
	hungry(0)->Intensity()[0][1] = 2.0;
	hungry(0)->Intensity()[1][0] = 0.01;
	hungry(0)->Intensity()[1][1] = -0.01;
	hungry(1)->Intensity()[0][0] = -0.01;
	hungry(1)->Intensity()[0][1] = 0.01;
	hungry(1)->Intensity()[1][0] = 0.01;
	hungry(1)->Intensity()[1][1] = -0.01;
	hungry(2)->Intensity()[0][0] = -0.01;
	hungry(2)->Intensity()[0][1] = 0.01;
	hungry(2)->Intensity()[1][0] = 10.0;
	hungry(2)->Intensity()[1][1] = -10.0;

	vectr phungry(2, 0.0);
	phungry[0] = 0.6;
	phungry[1] = 0.4;	
	hungryrv[0].SetDist(phungry);

	//------uptake--------
	uptake(0)->Intensity()[0][0] = -0.0;
	uptake(0)->Intensity()[0][1] = 0.0;
	uptake(0)->Intensity()[1][0] = 0.5;
	uptake(0)->Intensity()[1][1] = -0.5;

	vectr puptake(2, 0.0);
	puptake[0] = 0.45;	
	puptake[1] = 0.55;
	uptakerv[0].SetDist(puptake);


	// For this variable there is more than one parent.
	// The CIM is filled by creating an Instantiation object that represents 
	// the current intensity matrix we are entering values for. The 
	// MarkovDyn objects can be then indexed with Instantiation objects.

	//------concentration------
	Instantiation i(Context(fullc, uptakec));

	// Set uptake to state 0; full to state 0
	i.SetVal(3,0); i.SetVal(1,0);
	concentration(i)->Intensity()[0][0] = -0.02;
	concentration(i)->Intensity()[0][1] = 0.01;
	concentration(i)->Intensity()[0][2] = 0.01;
	concentration(i)->Intensity()[1][0] = 0.25;
	concentration(i)->Intensity()[1][1] = -0.26;
	concentration(i)->Intensity()[1][2] = 0.01;
	concentration(i)->Intensity()[2][0] = 0.01;
	concentration(i)->Intensity()[2][1] = 0.5;
	concentration(i)->Intensity()[2][2] = -0.51;
	// Set uptake to state 1; full to state 0
	// (We repeat this process for each parent instantiation.)
	i.SetVal(3,1); i.SetVal(1,0);
	concentration(i)->Intensity()[0][0] = -2.01;
	concentration(i)->Intensity()[0][1] = 2.0;
	concentration(i)->Intensity()[0][2] = 0.01;
	concentration(i)->Intensity()[1][0] = 0.01;
	concentration(i)->Intensity()[1][1] = -1.01;
	concentration(i)->Intensity()[1][2] = 1.0;
	concentration(i)->Intensity()[2][0] = 0.01;
	concentration(i)->Intensity()[2][1] = 0.01;
	concentration(i)->Intensity()[2][2] = -0.02;
	i.SetVal(3,0); i.SetVal(1,1);
	concentration(i)->Intensity()[0][0] = -0.02;
	concentration(i)->Intensity()[0][1] = 0.01;
	concentration(i)->Intensity()[0][2] = 0.01;
	concentration(i)->Intensity()[1][0] = 0.25;
	concentration(i)->Intensity()[1][1] = -0.26;
	concentration(i)->Intensity()[1][2] = 0.01;
	concentration(i)->Intensity()[2][0] = 0.01;
	concentration(i)->Intensity()[2][1] = 0.5;
	concentration(i)->Intensity()[2][2] = -0.51;
	i.SetVal(3,1); i.SetVal(1,1);
	concentration(i)->Intensity()[0][0] = -1.01;
	concentration(i)->Intensity()[0][1] = 1.0;
	concentration(i)->Intensity()[0][2] = 0.01;
	concentration(i)->Intensity()[1][0] = 0.01;
	concentration(i)->Intensity()[1][1] = -2.01;
	concentration(i)->Intensity()[1][2] = 2.0;
	concentration(i)->Intensity()[2][0] = 0.01;
	concentration(i)->Intensity()[2][1] = 0.01;
	concentration(i)->Intensity()[2][2] = -0.02;

	i.SetVal(3,0); i.SetVal(1,2);
	concentration(i)->Intensity()[0][0] = -0.02;
	concentration(i)->Intensity()[0][1] = 0.01;
	concentration(i)->Intensity()[0][2] = 0.01;
	concentration(i)->Intensity()[1][0] = 0.25;
	concentration(i)->Intensity()[1][1] = -0.26;
	concentration(i)->Intensity()[1][2] = 0.01;
	concentration(i)->Intensity()[2][0] = 0.01;
	concentration(i)->Intensity()[2][1] = 0.5;
	concentration(i)->Intensity()[2][2] = -0.51;

	i.SetVal(3,1); i.SetVal(1,2);
	concentration(i)->Intensity()[0][0] = -0.51;
	concentration(i)->Intensity()[0][1] = 0.5;
	concentration(i)->Intensity()[0][2] = 0.01;
	concentration(i)->Intensity()[1][0] = 0.01;
	concentration(i)->Intensity()[1][1] = -4.01;
	concentration(i)->Intensity()[1][2] = 4.0;
	concentration(i)->Intensity()[2][0] = 0.01;
	concentration(i)->Intensity()[2][1] = 0.01;
	concentration(i)->Intensity()[2][2] = -0.02;

	vectr pconcentration(3, 0.0);
	pconcentration[0] = 0.7;//0.1;
	pconcentration[1] = 0.2;//0.2;
	pconcentration[2] = 0.1;//0.7;
	concentrationrv[0].SetDist(pconcentration);

	//------barometer------
	barometer(0)->Intensity()[0][0] = -0.21;
	barometer(0)->Intensity()[0][1] = 0.2;
	barometer(0)->Intensity()[0][2] = 0.01;
	barometer(0)->Intensity()[1][0] = 0.05;
	barometer(0)->Intensity()[1][1] = -0.1;
	barometer(0)->Intensity()[1][2] = 0.05;
	barometer(0)->Intensity()[2][0] = 0.01;
	barometer(0)->Intensity()[2][1] = 0.2;
	barometer(0)->Intensity()[2][2] = -0.21;

	vectr pbarometer(3, 0.0);
	pbarometer[0] = 0.3;
	pbarometer[1] = 0.7;
	barometerrv[0].SetDist(pbarometer);

	//------pain------
	Instantiation j(Context(barometerc, concentrationc));
	j.SetVal(4,0); j.SetVal(5,0);
	pain(j)->Intensity()[0][0] = -6.0;
	pain(j)->Intensity()[0][1] = 6.0;
	pain(j)->Intensity()[1][0] = 0.1;
	pain(j)->Intensity()[1][1] = -0.1;

	j.SetVal(4,1); j.SetVal(5,0);
	pain(j)->Intensity()[0][0] = -1.0;
	pain(j)->Intensity()[0][1] = 1.0;
	pain(j)->Intensity()[1][0] = 0.1;
	pain(j)->Intensity()[1][1] = -0.1;
	j.SetVal(4,2); j.SetVal(5,0);
	pain(j)->Intensity()[0][0] = -6.0;
	pain(j)->Intensity()[0][1] = 6.0;
	pain(j)->Intensity()[1][0] = 0.1;
	pain(j)->Intensity()[1][1] = -0.1;
	j.SetVal(4,0); j.SetVal(5,1);
	pain(j)->Intensity()[0][0] = -1.0;
	pain(j)->Intensity()[0][1] = 1.0;
	pain(j)->Intensity()[1][0] = 0.3;
	pain(j)->Intensity()[1][1] = -0.3;
	j.SetVal(4,1); j.SetVal(5,1);
	pain(j)->Intensity()[0][0] = -0.3;
	pain(j)->Intensity()[0][1] = 0.3;
	pain(j)->Intensity()[1][0] = 0.3;
	pain(j)->Intensity()[1][1] = -0.3;
	j.SetVal(4,2); j.SetVal(5,1);
	pain(j)->Intensity()[0][0] = -1.0;
	pain(j)->Intensity()[0][1] = 1.0;
	pain(j)->Intensity()[1][0] = 0.3;
	pain(j)->Intensity()[1][1] = -0.3;
	j.SetVal(4,0); j.SetVal(5,2);
	pain(j)->Intensity()[0][0] = -0.01;
	pain(j)->Intensity()[0][1] = 0.01;
	pain(j)->Intensity()[1][0] = 2.0;
	pain(j)->Intensity()[1][1] = -2.0;
	j.SetVal(4,1); j.SetVal(5,2);
	pain(j)->Intensity()[0][0] = -0.01;
	pain(j)->Intensity()[0][1] = 0.01;
	pain(j)->Intensity()[1][0] = 2.0;
	pain(j)->Intensity()[1][1] = -2.0;
	j.SetVal(4,2); j.SetVal(5,2);
	pain(j)->Intensity()[0][0] = -0.01;
	pain(j)->Intensity()[0][1] = 0.01;
	pain(j)->Intensity()[1][0] = 2.0;
	pain(j)->Intensity()[1][1] = -2.0;

	vectr ppain(2, 0.0);
	ppain[0] = 0.1;
	ppain[1] = 0.9;
	painrv[0].SetDist(ppain);

	//------drowsy------
	drowsy(0)->Intensity()[0][0] = -0.17;
	drowsy(0)->Intensity()[0][1] = 0.17;
	drowsy(0)->Intensity()[1][0] = 0.5;
	drowsy(0)->Intensity()[1][1] = -0.5;
	drowsy(1)->Intensity()[0][0] = -0.33;
	drowsy(1)->Intensity()[0][1] = 0.33;
	drowsy(1)->Intensity()[1][0] = 0.33;
	drowsy(1)->Intensity()[1][1] = -0.33;
	drowsy(2)->Intensity()[0][0] = -1.0;
	drowsy(2)->Intensity()[0][1] = 1.0;
	drowsy(2)->Intensity()[1][0] = 0.1;
	drowsy(2)->Intensity()[1][1] = -0.1;

	vectr pdrowsy(2, 0.0);
	pdrowsy[0] = 0.1;
	pdrowsy[1] = 0.9;
	drowsyrv[0].SetDist(pdrowsy);


	// The dynamics of a CTBN are defined by the dynamics of each conditional
	// Markov process. Here it is the MarkovDyn objects created above.
	CTBNDyn ctbndyn;
	ctbndyn.AddNode(eat.Clone());
	ctbndyn.AddNode(full.Clone());
	ctbndyn.AddNode(hungry.Clone());

	BN bn;
	bn.AddNode(eatrv.Clone());
	bn.AddNode(fullrv.Clone());
	bn.AddNode(hungryrv.Clone());

	//(ctbndyn.JointMatrix()).niceprint(cout);

	ctbndyn.AddNode(uptake.Clone());
	ctbndyn.AddNode(concentration.Clone());
	ctbndyn.AddNode(barometer.Clone());
	ctbndyn.AddNode(pain.Clone());
	ctbndyn.AddNode(drowsy.Clone());

	// The initial distribution of a CTBN is defined by the underlying 
	// Bayesian network containing a multinomial distribution for each
	// variable.

	bn.AddNode(uptakerv.Clone());
	bn.AddNode(concentrationrv.Clone());
	bn.AddNode(barometerrv.Clone());
	bn.AddNode(painrv.Clone());
	bn.AddNode(drowsyrv.Clone());

	/*
	Dynamics *dyn = &ctbndyn;
	FactoredMatrix apprm(dyn);
	RV *rvbn = &bn;
	FactoredVector apprv(rvbn);
	apprv.Print();
	//cout << apprv.GetDistInst(0,0) << endl;
	FactoredVector bbb = apprv;
	//cout << bbb.GetDistInst(0,0) << endl;
	//apprm.Mult(1, bbb);
	
	Context context = ctbndyn.Domain() + ctbndyn.CondDomain();
	Instantiation x(context, -1);
	int querynode = ParamInt("QueryNode", 5);
	int querytimeval = ParamInt("QueryTimeVal", 1);
	x.SetVal(querynode, querytimeval);

	querynode = ParamInt("QueryNode", 1);
	querytimeval = ParamInt("QueryTimeVal", 0);
	x.SetVal(querynode, querytimeval);

	querynode = ParamInt("QueryNode", 2);
	querytimeval = ParamInt("QueryTimeVal", 0);
	x.SetVal(querynode, querytimeval);

	querynode = ParamInt("QueryNode",3);
	querytimeval = ParamInt("QueryTimeVal", 0);
	x.SetVal(querynode, querytimeval);

	querynode = ParamInt("QueryNode", 4);
	querytimeval = ParamInt("QueryTimeVal", 0);
	x.SetVal(querynode, querytimeval);

	RVSimple *v = &apprv;
	apprm.Cond(0.0, 1.0, x);
	apprm.Mult(v);
	 */
	 
	// Construct the final object that represents the CTBN.
	Markov m(bn.Clone(), ctbndyn.Clone());
	Context context = ctbndyn.Domain() + ctbndyn.CondDomain();
	
	// sample a partial trajectory as the input evidence
	double begintime = ParamDouble("BeginTime", 0.0);
	double endtime = ParamDouble("EndTime", 10.0);
	Trajectory evid;
	evid.SetBeginTime(begintime);
	evid.SetEndTime(endtime);
	
	evid.AddPointEvidence(0, 0.5, 1);
	evid.AddTransition(1, 2.5, 2);
	evid.AddTransition(6, 5.0, 1);
	//m.Sample(evid);
	//RemoveNodesInformation(evid, context, ctbndyn.NumofNodes(), ctbndyn.NumofNodes()*6 , 0.3);
	
	m.Save(cout);
	cout << endl << evid << endl;

	// Output the CTBN to the cout strem
/*	if (savedyn) {
		//if (saveasptr) ctbndyn.SavePtr(cout);
		//else ctbndyn.Save(cout);
		ctbndyn.Save(cout);
	} else {
		//if (saveasptr) m.SavePtr(cout);
		//else m.Save(cout);
		m.Save(cout);
	}
*/
}
