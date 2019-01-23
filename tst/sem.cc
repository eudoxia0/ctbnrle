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
#include "markov.h"
#include "streamextra.h"
#include "context.h"
#include "multirv.h"
#include "trajectory.h"
#include "exactmarkovinf.h"
#include "params.h"
#include "ctbndyn.h"
#include "bn.h"
#include "em.h"
#include "structure.h"
#include "suffstatsquery.h"
#include "nonexpsuffstatsquery.h"
#include "brutestructuresearch.h"
#include "grapheditsearch.h"
#include "utils.h"
#include "samplinginf.h"
#include "expmethod.h"
#include "importancesampler.h"

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
using namespace ctbn;

namespace {
     int testing_result = 0;
     double kEpsilon = 0.1;
}

#define EXPECT(cond) do { \
        if (!(cond)) {                                                  \
               std::cerr << "Test failure at " << __FILE__ << ":" << __LINE__ \
                      << " in " << __func__ \
                      << ", failed condition is " #cond << std::endl; \
            testing_result = 1; \
        } \
    } while (0)



int main(int argc, char **argv) {
  SetParam("EMLLHImprov","0.05");
  SetParam("EMItt","10");
  SetParam("VerboseEM","0");
  InitParams(argc, argv);
  randomizer.Reset(11324432); // could get unlucky, but this one works (as do many others)

  Context nullc;
  Context eatc, fullc, hungryc, uptakec, concentrationc, barometerc, painc;//, drowsyc;
  eatc.AddVar(0,2);
  fullc.AddVar(1,2);
  hungryc.AddVar(2,2);
  uptakec.AddVar(3,2);
  concentrationc.AddVar(4,2);
  barometerc.AddVar(5,2);
  painc.AddVar(6,2);
//  drowsyc.AddVar(7,2);

  MarkovDyn eat(eatc, hungryc);
  MarkovDyn full(fullc, eatc);
  MarkovDyn hungry(hungryc, fullc);
  MarkovDyn uptake(uptakec, nullc);
  Context full_and_uptakec(fullc, uptakec);
  MarkovDyn concentration(concentrationc, full_and_uptakec);
  MarkovDyn barometer(barometerc, nullc);
  MarkovDyn pain(painc, Context(barometerc, concentrationc));
//  MarkovDyn drowsy(drowsyc, concentrationc);

	// Generate starting distributions for each variable
	MultiRV eatrv(eatc, nullc);
	MultiRV fullrv(fullc, eatc);
	MultiRV hungryrv(hungryc, eatc);
	MultiRV uptakerv(uptakec, nullc);
	MultiRV concentrationrv(concentrationc, uptakec+hungryc);
	MultiRV barometerrv(barometerc, nullc);
	MultiRV painrv(painc, barometerc);
//	MultiRV drowsyrv(drowsyc, concentrationc+painc);

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
  full(0)->Intensity()[0][1] = 0.02;
      full(0)->Intensity()[1][0] = 0.3;
      full(0)->Intensity()[1][1] = -0.3;

  full(1)->Intensity()[0][0] = -6.0;
  full(1)->Intensity()[0][1] = 6.0;
      full(1)->Intensity()[1][0] = 3.01;
      full(1)->Intensity()[1][1] = -3.01;

  vectr pfull(2, 0.0);
	pfull[0] = 0.75;
	pfull[1] = 0.25;
	fullrv[0].SetDist(pfull);
	pfull[0] = 0.55;
	pfull[1] = 0.45;
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
  puptake[0] = 0.3;
  puptake[1] = 0.7;
  uptakerv[0].SetDist(puptake);

  //------concentration------
  Instantiation i(Context(fullc, uptakec));
  i.SetVal(3,0); i.SetVal(1,0);
  concentration(i)->Intensity()[0][0] = -0.02;
  concentration(i)->Intensity()[0][1] = 0.02;
      concentration(i)->Intensity()[1][0] = 0.26;
      concentration(i)->Intensity()[1][1] = -0.26;
  i.SetVal(3,1); i.SetVal(1,0);
  concentration(i)->Intensity()[0][0] = -2.01;
  concentration(i)->Intensity()[0][1] = 2.01;
      concentration(i)->Intensity()[1][0] = 1.01;
      concentration(i)->Intensity()[1][1] = -1.01;
  i.SetVal(3,0); i.SetVal(1,1);
  concentration(i)->Intensity()[0][0] = -0.01;
  concentration(i)->Intensity()[0][1] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.25;
      concentration(i)->Intensity()[1][1] = -0.25;
  i.SetVal(3,1); i.SetVal(1,1);
  concentration(i)->Intensity()[0][0] = -1.01;
  concentration(i)->Intensity()[0][1] = 1.01;
      concentration(i)->Intensity()[1][0] = 2.01;
      concentration(i)->Intensity()[1][1] = -2.01;

 	vectr pconcentration(2, 0.0);
	pconcentration[0] = 0.7;
	pconcentration[1] = 0.3;
	concentrationrv[0].SetDist(pconcentration);
	pconcentration[0] = 0.6;
	pconcentration[1] = 0.4;
	concentrationrv[1].SetDist(pconcentration);
	pconcentration[0] = 0.8;
	pconcentration[1] = 0.2;
	concentrationrv[2].SetDist(pconcentration);
	pconcentration[0] = 0.5;
	pconcentration[1] = 0.5;
	concentrationrv[3].SetDist(pconcentration);

  //------barometer------
  barometer(0)->Intensity()[0][0] = -0.2;
  barometer(0)->Intensity()[0][1] = 0.2;
      barometer(0)->Intensity()[1][0] = 0.1;
      barometer(0)->Intensity()[1][1] = -0.1;

  vectr pbarometer(2, 0.0);
  pbarometer[0] = 0.2;
  pbarometer[0] = 0.8;
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


	vectr ppain(2, 0.0);
	ppain[0] = 0.1;
	ppain[1] = 0.9;
	painrv[0].SetDist(ppain);
	ppain[0] = 0.3;
	ppain[1] = 0.7;
	painrv[1].SetDist(ppain);

/*
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
*/	

  CTBNDyn ctbndyn;
  ctbndyn.AddNode(eat.Clone());
  ctbndyn.AddNode(full.Clone());
  ctbndyn.AddNode(hungry.Clone());
  ctbndyn.AddNode(uptake.Clone());
  ctbndyn.AddNode(concentration.Clone());
  ctbndyn.AddNode(barometer.Clone());
  ctbndyn.AddNode(pain.Clone());
//  ctbndyn.AddNode(drowsy.Clone());

  BN bn;
  bn.AddNode(eatrv.Clone());
  bn.AddNode(fullrv.Clone());
  bn.AddNode(hungryrv.Clone());
  bn.AddNode(uptakerv.Clone());
  bn.AddNode(concentrationrv.Clone());
  bn.AddNode(barometerrv.Clone());
  bn.AddNode(painrv.Clone());
//  bn.AddNode(drowsyrv.Clone());

  Markov m(bn.Clone(), ctbndyn.Clone());

	
	//Building independent structure
	MultiRV eatrvi(eatc, nullc);
	MultiRV fullrvi(fullc, nullc);
	MultiRV hungryrvi(hungryc, nullc);
	MultiRV uptakervi(uptakec, nullc);
	MultiRV concentrationrvi(concentrationc, nullc);
	MultiRV barometerrvi(barometerc, nullc);
	MultiRV painrvi(painc, nullc);
//	MultiRV drowsyrvi(drowsyc, nullc);

	BN bnindep;
	bnindep.AddNode(eatrvi.Clone());
	bnindep.AddNode(fullrvi.Clone());
	bnindep.AddNode(hungryrvi.Clone());
	bnindep.AddNode(uptakervi.Clone());
	bnindep.AddNode(concentrationrvi.Clone());
	bnindep.AddNode(barometerrvi.Clone());
	bnindep.AddNode(painrvi.Clone());
//	bnindep.AddNode(drowsyrvi.Clone());
	
	
	MarkovDyn eatindep(eatc, nullc);
	MarkovDyn fullindep(fullc, nullc);
	MarkovDyn hungryindep(hungryc, nullc);
	MarkovDyn uptakeindep(uptakec, nullc);
	MarkovDyn concentrationindep(concentrationc, nullc);
	MarkovDyn barometerindep(barometerc, nullc);
	MarkovDyn painindep(painc, nullc);
//	MarkovDyn drowsyindep(drowsyc, nullc);
	
	CTBNDyn ctbndynindep;
	ctbndynindep.AddNode(eatindep.Clone());
	ctbndynindep.AddNode(fullindep.Clone());
	ctbndynindep.AddNode(hungryindep.Clone());
	ctbndynindep.AddNode(uptakeindep.Clone());
	ctbndynindep.AddNode(concentrationindep.Clone());
	ctbndynindep.AddNode(barometerindep.Clone());
	ctbndynindep.AddNode(painindep.Clone());
//	ctbndynindep.AddNode(drowsyindep.Clone());



	Markov mindep(bnindep.Clone(), ctbndynindep.Clone());

//	mindep.Scramble();
	//=============
	//STRUCTURE SEARCH TESTING
	//=============
	
	int NUM_TRAJ = ParamInt("numtraj", 500); 
	double TRAJ_LEN = ParamDouble("tlen", 10.0);
	double loss = ParamDouble("loss", 0.25);

	double factor = loss/0.25;

	const Context& con = ctbndyn.Domain();
	int numNodes = ctbndyn.NumofNodes();

	//Generate trajectories for training
	vector <Trajectory> trVec;
	if (ParamInt("TrajFromInput",0)) {
		cin >> trVec;
	} else {
		for(int i=0; i<NUM_TRAJ; i++) {
			Trajectory tr; tr.SetEndTime(TRAJ_LEN);
			m.Sample(tr);
			int nit = (int)round(numNodes*TRAJ_LEN*factor);
			RemoveNodesInformation(tr, con, numNodes, nit, 0.25/TRAJ_LEN);
			trVec.push_back(tr);
		}
	}
	if (ParamInt("TrajToOutput",0)) {
		cout << trVec << endl;
		exit(0);
	}
	Process *sem_ctbn = mindep.Clone();
	dynamic_cast<Markov *>(sem_ctbn)->Scramble();

	Structure sABN;
	Structure sACTBN;

	Structure sLBN;
	Structure sLCTBN;

	bn.GetStructure(sABN);

	ctbndyn.GetStructure(sACTBN);

	ExactMarkovInf inf;
		
	SEM(trVec,sem_ctbn,&inf);

	Markov* learned = dynamic_cast <Markov *>(sem_ctbn);

	const CTBNDyn* learnedCTBN = dynamic_cast <const CTBNDyn *>(learned->GetDynamics());
	learnedCTBN->GetStructure(sLCTBN);
	sACTBN.Print(cout); cout << endl;
	sLCTBN.Print(cout); cout << endl;

	EXPECT(sACTBN.HammingDist(sLCTBN)==0);
	if (testing_result==0)
		cout << "PASS sem unittest." << endl;
	return testing_result!=0;
}
