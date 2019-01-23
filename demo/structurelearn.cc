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

// This code generates the drug effect network as first seen in 
// Continuous Time Bayesian Networks by Nodelman et al., 2002.
// To demonstrate structure search on complete data, trajectories
// are sampled from the network and used to learn the structure.

/* Command line parameters:
 *
 * Set the number of trajectories
 * -Dnumtraj <integer> (default: 1000)
 *
 * Set the length of each trajectory
 * -Dtlen <double> (default: 5.0)
 *
 * Set the alpha prior (prior belief for number of transitions)
 * -Dalpha <double> (defualt: 1.0)
 *
 * Set the tau prior (prior belief for amount of time spent)
 * -Dtau <double> (defualt: 1.0)
 */

#include "markovdyn.h"
#include "markov.h"
#include "streamextra.h"
#include "context.h"
#include "multirv.h"
#include "trajectory.h"
#include "params.h"
#include "ctbndyn.h"
#include "bn.h"
#include "em.h"
#include "ctbn.h"
#include "nonexpsuffstatsquery.h"
#include "famscore.h"
#include "grapheditsearch.h"
#include <iomanip>
#include <iostream>

using namespace std;
using namespace ctbn;

int main(int argc, char **argv) {
	InitParams(argc, argv);
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
	MultiRV fullrv(fullc, eatc);
	MultiRV hungryrv(hungryc, eatc);
	MultiRV uptakerv(uptakec, nullc);
	MultiRV concentrationrv(concentrationc, uptakec+hungryc);
	MultiRV barometerrv(barometerc, nullc);
	MultiRV painrv(painc, barometerc);
	MultiRV drowsyrv(drowsyc, concentrationc+painc);

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
	peat[0] = 0.3;
	peat[1] = 0.7;
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
	pfull[0] = 0.5;
	pfull[1] = 0.25;
	pfull[2] = 0.25;
	fullrv[0].SetDist(pfull);
	pfull[0] = 0.2;
	pfull[1] = 0.25;
	pfull[2] = 0.55;
	fullrv[1].SetDist(pfull);
	

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
	phungry[0] = 0.35;
	phungry[1] = 0.65;
	hungryrv[0].SetDist(phungry);
	phungry[0] = 0.45;
	phungry[1] = 0.55;
	hungryrv[1].SetDist(phungry);

	//------uptake--------
	uptake(0)->Intensity()[0][0] = -0.0;
	uptake(0)->Intensity()[0][1] = 0.0;
	uptake(0)->Intensity()[1][0] = 0.5;
	uptake(0)->Intensity()[1][1] = -0.5;

	vectr puptake(2, 0.0);
	puptake[1] = 1;
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
	pconcentration[0] = 0.1;
	pconcentration[1] = 0.3;
	pconcentration[2] = 0.6;
	concentrationrv[0].SetDist(pconcentration);
	pconcentration[0] = 0.1;
	pconcentration[1] = 0.4;
	pconcentration[2] = 0.5;
	concentrationrv[1].SetDist(pconcentration);
	pconcentration[0] = 0.4;
	pconcentration[1] = 0.2;
	pconcentration[2] = 0.4;
	concentrationrv[2].SetDist(pconcentration);
	pconcentration[0] = 0.5;
	pconcentration[1] = 0.4;
	pconcentration[2] = 0.1;
	concentrationrv[3].SetDist(pconcentration);

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
	pbarometer[1] = 1;
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
	ppain[0] = 0.3;
	ppain[1] = 0.7;
	painrv[1].SetDist(ppain);
	ppain[0] = 0.9;
	ppain[1] = 0.1;
	painrv[2].SetDist(ppain);
	
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
	pdrowsy[0] = 0.8;
	pdrowsy[1] = 0.2;
	drowsyrv[0].SetDist(pdrowsy);
	pdrowsy[0] = 0.7;
	pdrowsy[1] = 0.3;
	drowsyrv[1].SetDist(pdrowsy);
	pdrowsy[0] = 0.3;
	pdrowsy[1] = 0.7;
	drowsyrv[2].SetDist(pdrowsy);
	pdrowsy[0] = 0.4;
	pdrowsy[1] = 0.6;
	drowsyrv[3].SetDist(pdrowsy);
	pdrowsy[0] = 0.2;
	pdrowsy[1] = 0.8;
	drowsyrv[4].SetDist(pdrowsy);
	pdrowsy[0] = 0.8;
	pdrowsy[1] = 0.2;
	drowsyrv[5].SetDist(pdrowsy);
	

	// The dynamics of a CTBN are defined by the dynamics of each conditional
	// Markov process. Here it is the MarkovDyn objects created above.
	CTBNDyn ctbndyn;
	ctbndyn.AddNode(eat.Clone());
	ctbndyn.AddNode(full.Clone());
	ctbndyn.AddNode(hungry.Clone());
	ctbndyn.AddNode(uptake.Clone());
	ctbndyn.AddNode(concentration.Clone());
	ctbndyn.AddNode(barometer.Clone());
	ctbndyn.AddNode(pain.Clone());
	ctbndyn.AddNode(drowsy.Clone());

	// The initial distribution of a CTBN is defined by the underlying 
	// Bayesian network containing a multinomial distribution for each
	// variable.
	BN bn;
	bn.AddNode(eatrv.Clone());
	bn.AddNode(fullrv.Clone());
	bn.AddNode(hungryrv.Clone());
	bn.AddNode(uptakerv.Clone());
	bn.AddNode(concentrationrv.Clone());
	bn.AddNode(barometerrv.Clone());
	bn.AddNode(painrv.Clone());
	bn.AddNode(drowsyrv.Clone());

	// Construct the final object that represents the CTBN.
	Markov m(bn.Clone(), ctbndyn.Clone());

	cout << "Demo of learning a CTBN structure" << endl 
		<< "================================="<< endl;
	// Load the actual structure into a structure object
	Structure sA;
	ctbndyn.GetStructure(sA);
	cout << "Actual Structure:\n";
	sA.Print(cout);
	cout << endl << endl;

	// Set the number and length of trajectories to sample
	int numTrajs = ParamInt("numtraj",1000);
	double tlen = ParamDouble("tlen",5.0);

	// Sample the trajectories
	cout << "Sampling trajectories...\n\n";
	vector<Trajectory> trVec;
	for(int i=0; i<numTrajs; i++) {
		Trajectory tr;
		tr.SetBeginTime(0);
		tr.SetEndTime(tlen);
		m.Sample(tr);
		trVec.push_back(tr);
	}

	// Create an object for returning sufficient statistics
	// for use with the scoring object
	NonExpSuffStatsQuery* ssQuery = new NonExpSuffStatsQuery();
	ssQuery->SetData(&trVec);
	ssQuery->SetDynamics(&ctbndyn);

	// Create the scoring object by specifying priors and the
	// sufficient statistics object to use for scoring
	double alphaPrior = ParamDouble("alpha",1.0);
	double tauPrior = ParamDouble("tau",1.0);
	CTBNFamScore* cfs = new CTBNFamScore(alphaPrior,tauPrior,ssQuery);

	// Create the search object by providing the context to search
	// over and the scoring function to use
	GraphEditSearch ges(ctbndyn.Domain(),cfs);

	cout << "Learning the structure...\n\n";

	// Calling LearnStructure() returns a structure object with the
	// learned structure
	Structure s = ges.LearnStructure();

	cout << "Learned Structure:\n";
	s.Print(cout);
	cout << endl;

	// Compared the learned structure against the actual structure
	// with Hamming distance
	cout << "Hamming Distance: " << s.HammingDist(sA) << endl;
	
	delete ssQuery;
	delete cfs;

	cout << endl << endl;
	cout << "Demo of learning a BN structure" << endl 
		<< "================================="<< endl;
	// Load the actual structure into a structure object
	Structure sABN;
	bn.GetStructure(sABN);
	cout << "Actual Structure:\n";
	sABN.Print(cout);
	cout << endl << endl;


	// Sample the trajectories (point evidence for BNs)
	cout << "Sampling...\n\n";
	vector<Trajectory> trVecBN;
	for(int i=0; i<numTrajs*2; i++) {
		Trajectory trBN;
		trBN.SetBeginTime(0);
		trBN.SetEndTime(1e-10);
		Instantiation x0;
		bn.Sample(x0);
		trBN.AddTransition(x0,trBN.TimeBegin());
		trVecBN.push_back(trBN);
	}

	// Create an object for returning sufficient statistics
	// for use with the scoring object
	NonExpSuffStatsQuery* ssQueryBN = new NonExpSuffStatsQuery();
	ssQueryBN->SetData(&trVecBN);
	ssQueryBN->SetInitRV(&bn);

	// Create the scoring object by specifying priors and the
	// sufficient statistics object to use for scoring
	BNFamScore* bfs = new BNFamScore(alphaPrior,ssQueryBN);

	// Create the search object by providing the context to search
	// over and the scoring function to use
	// We need to specify true as the 3rd parameter here to restrict
	// to structures without cycles
	GraphEditSearch gesBN(bn.Domain(),bfs,true);

	cout << "Learning the structure...\n\n";

	// Calling LearnStructure() returns a structure object with the
	// learned structure
	Structure sBN = gesBN.LearnStructure();

	cout << "Learned Structure:\n";
	sBN.Print(cout);
	cout << endl;

	// Compared the learned structure against the actual structure
	// by checking the difference between the number of parameters
	cout << "Difference in number of parameters: " 
		<< abs(sABN.GetNumParams(bn.Domain()) -
			sBN.GetNumParams(bn.Domain())) << endl;

	delete ssQueryBN;
	delete bfs;
}

