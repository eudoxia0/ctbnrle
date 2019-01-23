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
// To demonstrate structure search on incomplete data, trajectories
// are sampled from the network and evidence is removed. SEM is then 
// used (with an underlying inference method) to learn the structure
// and parameters of the model.

/* Command line parameters:
 *
 * Set the number of trajectories
 * -Dnumtraj <integer> (default: 1000)
 *
 * Set the length of each trajectory
 * -Dtlen <double> (default: 5.0)
 *
 * Set the alpha prior (prior belief for number of transitions)
 * -Dalpha <double> (default: 1.0)
 *
 * Set the tau prior (prior belief for amount of time spent)
 * -Dtau <double> (default: 1.0)
 *
 * Set the amount of loss in the trajectory data
 * -Dloss <double> (default: 0.25)
 */

#include "markovdyn.h"
#include "markov.h"
#include "streamextra.h"
#include "context.h"
#include "multirv.h"
#include "trajectory.h"
#include "samplinginf.h"
#include "expmethod.h"
#include "importancesampler.h"
#include "params.h"
#include "ctbndyn.h"
#include "bn.h"
#include "em.h"
#include "ctbn.h"
#include "nonexpsuffstatsquery.h"
#include "famscore.h"
#include "grapheditsearch.h"
#include "utils.h"
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
	eatc.AddVar(0,2);
	fullc.AddVar(1,3);
	hungryc.AddVar(2,2);
	uptakec.AddVar(3,2);
	concentrationc.AddVar(4,3);
	barometerc.AddVar(5,3);
	painc.AddVar(6,2);
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

	//Building independent structure
	MultiRV eatrvi(eatc, nullc);
	MultiRV fullrvi(fullc, nullc);
	MultiRV hungryrvi(hungryc, nullc);
	MultiRV uptakervi(uptakec, nullc);
	MultiRV concentrationrvi(concentrationc, nullc);
	MultiRV barometerrvi(barometerc, nullc);
	MultiRV painrvi(painc, nullc);
	MultiRV drowsyrvi(drowsyc, nullc);

	BN bnindep;
	bnindep.AddNode(eatrvi.Clone());
	bnindep.AddNode(fullrvi.Clone());
	bnindep.AddNode(hungryrvi.Clone());
	bnindep.AddNode(uptakervi.Clone());
	bnindep.AddNode(concentrationrvi.Clone());
	bnindep.AddNode(barometerrvi.Clone());
	bnindep.AddNode(painrvi.Clone());
	bnindep.AddNode(drowsyrvi.Clone());
	
	
	MarkovDyn eatindep(eatc, nullc);
	MarkovDyn fullindep(fullc, nullc);
	MarkovDyn hungryindep(hungryc, nullc);
	MarkovDyn uptakeindep(uptakec, nullc);
	MarkovDyn concentrationindep(concentrationc, nullc);
	MarkovDyn barometerindep(barometerc, nullc);
	MarkovDyn painindep(painc, nullc);
	MarkovDyn drowsyindep(drowsyc, nullc);
	
	CTBNDyn ctbndynindep;
	ctbndynindep.AddNode(eatindep.Clone());
	ctbndynindep.AddNode(fullindep.Clone());
	ctbndynindep.AddNode(hungryindep.Clone());
	ctbndynindep.AddNode(uptakeindep.Clone());
	ctbndynindep.AddNode(concentrationindep.Clone());
	ctbndynindep.AddNode(barometerindep.Clone());
	ctbndynindep.AddNode(painindep.Clone());
	ctbndynindep.AddNode(drowsyindep.Clone());

	Markov mindep(bnindep.Clone(), ctbndynindep.Clone());

	cout << "Demo of Structural EM" << endl 
		<< "================================="<< endl;
	// Load the actual structure into a structure object
	Structure sABN;
	Structure sACTBN;
	bn.GetStructure(sABN);
	ctbndyn.GetStructure(sACTBN);
	cout << "Actual Structure (BN):\n";
	sABN.Print(cout);
	cout << endl << endl;

	cout << "Actual Structure (CTBN):\n";
	sABN.Print(cout);
	cout << endl << endl;

	// Set the number and length of trajectories to sample
	int numTrajs = ParamInt("numtraj",500);
	double tlen = ParamDouble("tlen",5.0);
	double loss = ParamDouble("loss",0.25);
	double factor = loss/0.25;

	const Context& con = ctbndyn.Domain();
	int numNodes = ctbndyn.NumofNodes();

	// Sample the trajectories and remove evidence randomly 
	// in 0.25 time step segments. Also removes starting evidence 
	// from half of the variables at random from each trajectory.
	cout << "Sampling trajectories...\n\n";
	vector<Trajectory> trVec;
	for(int k=0; k<numTrajs; k++) {
		Trajectory tr; tr.SetEndTime(tlen);
		m.Sample(tr);
		int nit = (int)round(numNodes*tlen*factor);
		RemoveNodesInformation(tr, con, numNodes, nit, 0.25/tlen);
		for(int i=0; i<numNodes/2; i++) {
			int var = randomizer.RandInt(numNodes);
			RemoveInformation(tr, con, var, 0.00, 0.25);
		}
		trVec.push_back(tr);
	}

	// Use the independent structure as the starting point for SEM
	Process *sem_ctbn = mindep.Clone();

	// Initialize the parameters to some random values.
	dynamic_cast<Markov *>(sem_ctbn)->Scramble();

	// Structure classes to load the learned structures into
	Structure sLBN;
	Structure sLCTBN;

	// Use importance sampling as the inference method for SEM
	ExpMethod method;
	SamplingInf inf;
	Sampler *sampler;
	sampler = new ImportanceSampler(sem_ctbn, &trVec[0], &method);
	inf.SetSampler(sampler);

	cout << "Running SEM...\n";

	// Pass in data, starting process, and inference method to run SEM
	SEM(trVec,sem_ctbn,&inf);

	Markov *learned = dynamic_cast <Markov *>(sem_ctbn);

	// Get the learned Bayesian network
	const BN *learnedBN = dynamic_cast <const BN *>(
						learned->GetStartDist());
	learnedBN->GetStructure(sLBN);
	cout << "Learned structure (BN):\n";
	sLBN.Print(cout); 
	cout << endl;

	// Compute difference in number of parameters to measure accuracy
	int paramd = abs(sABN.GetNumParams(bn.Domain()) -
				sLBN.GetNumParams(learnedBN->Domain()));
	cout << "Param Distance: " << paramd << endl << endl;

	// Get the learned CTBN 
	const CTBNDyn* learnedCTBN = dynamic_cast <const CTBNDyn *>(
						learned->GetDynamics());
	learnedCTBN->GetStructure(sLCTBN);
	cout << "Learned structure (CTBN):\n";
	sLCTBN.Print(cout); 
	cout << endl;

	// Compute Hamming distance to measure accuracy	
	int hamming = sACTBN.HammingDist(sLCTBN);
	cout << "Hamming Distance: " << hamming << endl << endl; 
}
