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
#include <iostream>
#include "multirv.h"
#include "trajectory.h"
#include "params.h"
#include "ctbndyn.h"
#include "bn.h"
#include "em.h"
#include "structure.h"
#include "suffstatsquery.h"
#include "nonexpsuffstatsquery.h"
#include "brutestructuresearch.h"
#include "grapheditsearch.h"
#include "famscore.h"

#include <vector>
#include <fstream>

using namespace std;
using namespace ctbn;

//const int TRIALS(10);
double PRIOR_TRANS(1);
double PRIOR_TRANSBN(1);
double PRIOR_TIME(1);
int NUM_TRAJ(2000);
double TRAJ_LEN(10);

double BNFUZZ(0.2);

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



int main(int argc, char **argv)
{
  InitParams(argc, argv);


  Context nullc;
  Context eatc, fullc, hungryc, uptakec, concentrationc, barometerc, painc, drowsyc;
  eatc.AddVar(120,2);
  fullc.AddVar(241,3);
  hungryc.AddVar(312,2);
  uptakec.AddVar(463,2);
  concentrationc.AddVar(594,3);
  barometerc.AddVar(685,3);
  painc.AddVar(776,2);
  drowsyc.AddVar(837,2);

  MarkovDyn eat(eatc, hungryc);
  MarkovDyn full(fullc, eatc);
  MarkovDyn hungry(hungryc, fullc);
  MarkovDyn uptake(uptakec, nullc);
  Context full_and_uptakec(fullc, uptakec);
  MarkovDyn concentration(concentrationc, full_and_uptakec);
  MarkovDyn barometer(barometerc, nullc);
  MarkovDyn pain(painc, Context(barometerc, concentrationc));
  MarkovDyn drowsy(drowsyc, concentrationc);

  MultiRV eatrv(eatc, nullc);
  MultiRV fullrv(fullc, nullc);
  MultiRV hungryrv(hungryc, nullc);
  MultiRV uptakerv(uptakec, nullc);
  MultiRV concentrationrv(concentrationc, nullc);
  MultiRV barometerrv(barometerc, nullc);
  MultiRV painrv(painc, nullc);
  MultiRV drowsyrv(drowsyc, nullc);
  //------eat--------
  eat(0)->Intensity()[0][0] = -0.01;
  eat(0)->Intensity()[0][1] = 0.01;
      eat(0)->Intensity()[1][0] = 10.0;
      eat(0)->Intensity()[1][1] = -10.0;
  eat(1)->Intensity()[0][0] = -2.0;
  eat(1)->Intensity()[0][1] = 2.0;
      eat(1)->Intensity()[1][0] = 1.0;
      eat(1)->Intensity()[1][1] = -1.0;
  vectr peat(2, BNFUZZ);
  peat[0] = 1-BNFUZZ;
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

  vectr pfull(3, BNFUZZ);
  pfull[0] = 1-BNFUZZ*2;
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

  vectr phungry(2, BNFUZZ);
  phungry[0] = 1-BNFUZZ;
  hungryrv[0].SetDist(phungry);

  //------uptake--------
  uptake(0)->Intensity()[0][0] = -0.0;
  uptake(0)->Intensity()[0][1] = 0.0;
      uptake(0)->Intensity()[1][0] = 0.5;
      uptake(0)->Intensity()[1][1] = -0.5;

  vectr puptake(2, BNFUZZ);
  puptake[1] = 1-BNFUZZ;
  uptakerv[0].SetDist(puptake);

  //------concentration------
  Instantiation i(Context(fullc, uptakec));
  i.SetVal(463,0); i.SetVal(241,0);
  concentration(i)->Intensity()[0][0] = -0.02;
  concentration(i)->Intensity()[0][1] = 0.01;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.25;
      concentration(i)->Intensity()[1][1] = -0.26;
      concentration(i)->Intensity()[1][2] = 0.01;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.5;
          concentration(i)->Intensity()[2][2] = -0.51;
  i.SetVal(463,1); i.SetVal(241,0);
  concentration(i)->Intensity()[0][0] = -2.01;
  concentration(i)->Intensity()[0][1] = 2.0;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.01;
      concentration(i)->Intensity()[1][1] = -1.01;
      concentration(i)->Intensity()[1][2] = 1.0;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.01;
          concentration(i)->Intensity()[2][2] = -0.02;
  i.SetVal(463,0); i.SetVal(241,1);
  concentration(i)->Intensity()[0][0] = -0.02;
  concentration(i)->Intensity()[0][1] = 0.01;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.25;
      concentration(i)->Intensity()[1][1] = -0.26;
      concentration(i)->Intensity()[1][2] = 0.01;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.5;
          concentration(i)->Intensity()[2][2] = -0.51;
  i.SetVal(463,1); i.SetVal(241,1);
  concentration(i)->Intensity()[0][0] = -1.01;
  concentration(i)->Intensity()[0][1] = 1.0;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.01;
      concentration(i)->Intensity()[1][1] = -2.01;
      concentration(i)->Intensity()[1][2] = 2.0;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.01;
          concentration(i)->Intensity()[2][2] = -0.02;

  i.SetVal(463,0); i.SetVal(241,2);
  concentration(i)->Intensity()[0][0] = -0.02;
  concentration(i)->Intensity()[0][1] = 0.01;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.25;
      concentration(i)->Intensity()[1][1] = -0.26;
      concentration(i)->Intensity()[1][2] = 0.01;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.5;
          concentration(i)->Intensity()[2][2] = -0.51;

  i.SetVal(463,1); i.SetVal(241,2);
  concentration(i)->Intensity()[0][0] = -0.51;
  concentration(i)->Intensity()[0][1] = 0.5;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.01;
      concentration(i)->Intensity()[1][1] = -4.01;
      concentration(i)->Intensity()[1][2] = 4.0;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.01;
          concentration(i)->Intensity()[2][2] = -0.02;

  vectr pconcentration(3, BNFUZZ);
  pconcentration[0] = 1-2*BNFUZZ;
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

  vectr pbarometer(3, BNFUZZ);
  pbarometer[1] = 1-2*BNFUZZ;
  barometerrv[0].SetDist(pbarometer);

  //------pain------
  Instantiation j(Context(barometerc, concentrationc));
  j.SetVal(594,0); j.SetVal(685,0);
  pain(j)->Intensity()[0][0] = -6.0;
  pain(j)->Intensity()[0][1] = 6.0;
      pain(j)->Intensity()[1][0] = 0.1;
      pain(j)->Intensity()[1][1] = -0.1;

  j.SetVal(594,1); j.SetVal(685,0);
  pain(j)->Intensity()[0][0] = -1.0;
  pain(j)->Intensity()[0][1] = 1.0;
      pain(j)->Intensity()[1][0] = 0.1;
      pain(j)->Intensity()[1][1] = -0.1;
  j.SetVal(594,2); j.SetVal(685,0);
  pain(j)->Intensity()[0][0] = -6.0;
  pain(j)->Intensity()[0][1] = 6.0;
      pain(j)->Intensity()[1][0] = 0.1;
      pain(j)->Intensity()[1][1] = -0.1;
  j.SetVal(594,0); j.SetVal(685,1);
  pain(j)->Intensity()[0][0] = -1.0;
  pain(j)->Intensity()[0][1] = 1.0;
      pain(j)->Intensity()[1][0] = 0.3;
      pain(j)->Intensity()[1][1] = -0.3;
  j.SetVal(594,1); j.SetVal(685,1);
  pain(j)->Intensity()[0][0] = -0.3;
  pain(j)->Intensity()[0][1] = 0.3;
      pain(j)->Intensity()[1][0] = 0.3;
      pain(j)->Intensity()[1][1] = -0.3;
  j.SetVal(594,2); j.SetVal(685,1);
  pain(j)->Intensity()[0][0] = -1.0;
  pain(j)->Intensity()[0][1] = 1.0;
      pain(j)->Intensity()[1][0] = 0.3;
      pain(j)->Intensity()[1][1] = -0.3;
  j.SetVal(594,0); j.SetVal(685,2);
  pain(j)->Intensity()[0][0] = -0.01;
  pain(j)->Intensity()[0][1] = 0.01;
      pain(j)->Intensity()[1][0] = 2.0;
      pain(j)->Intensity()[1][1] = -2.0;
  j.SetVal(594,1); j.SetVal(685,2);
  pain(j)->Intensity()[0][0] = -0.01;
  pain(j)->Intensity()[0][1] = 0.01;
      pain(j)->Intensity()[1][0] = 2.0;
      pain(j)->Intensity()[1][1] = -2.0;
  j.SetVal(594,2); j.SetVal(685,2);
  pain(j)->Intensity()[0][0] = -0.01;
  pain(j)->Intensity()[0][1] = 0.01;
      pain(j)->Intensity()[1][0] = 2.0;
      pain(j)->Intensity()[1][1] = -2.0;

  vectr ppain(2, BNFUZZ);
  ppain[1] = 1-BNFUZZ;
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

  vectr pdrowsy(2, BNFUZZ);
  pdrowsy[0] = 1-BNFUZZ;
  drowsyrv[0].SetDist(pdrowsy);


  CTBNDyn ctbndyn;
  ctbndyn.AddNode(hungry.Clone());
  ctbndyn.AddNode(eat.Clone());
  ctbndyn.AddNode(drowsy.Clone());
  ctbndyn.AddNode(full.Clone());
  ctbndyn.AddNode(uptake.Clone());
  ctbndyn.AddNode(pain.Clone());
  ctbndyn.AddNode(barometer.Clone());
  ctbndyn.AddNode(concentration.Clone());

  BN bn;
  bn.AddNode(hungryrv.Clone());
  bn.AddNode(eatrv.Clone());
  bn.AddNode(drowsyrv.Clone());
  bn.AddNode(fullrv.Clone());
  bn.AddNode(uptakerv.Clone());
  bn.AddNode(painrv.Clone());
  bn.AddNode(barometerrv.Clone());
  bn.AddNode(concentrationrv.Clone());

  Markov m(bn.Clone(), ctbndyn.Clone());

	//=============
	//STRUCTURE SEARCH TESTING
	//=============
	
	//Generate trajectories
	vector <Trajectory> trVec;
	vector <double> w(NUM_TRAJ,1);
	for(int i=0; i<NUM_TRAJ; i++) {
		Trajectory tr; 
		tr.SetBeginTime(0);
		tr.SetEndTime(TRAJ_LEN);
		m.Sample(tr);
		trVec.push_back(tr);
	}

	Structure sA;
	ctbndyn.GetStructure(sA);

	Structure s;
	NonExpSuffStatsQuery* ssQuery = new NonExpSuffStatsQuery();
	ssQuery->SetData(&trVec);
	CTBNFamScore* cfs = new CTBNFamScore(PRIOR_TRANS,PRIOR_TIME,ssQuery);
	BruteStructureSearch bss(ctbndyn.Domain(),cfs,3);
	GraphEditSearch ges(ctbndyn.Domain(),cfs);

	s = ges.LearnStructure();
	EXPECT(sA.HammingDist(s)==0);

	s = bss.LearnStructure();
	EXPECT(sA.HammingDist(s)==0);

	delete cfs;

	BNFamScore* bfs = new BNFamScore(PRIOR_TRANSBN,ssQuery);
	GraphEditSearch gesBN(bn.Domain(),bfs,true);
	s = gesBN.LearnStructure();
	Structure sBN;
	bn.GetStructure(sBN); // should be empty (no parents)
	EXPECT(sBN.HammingDist(s)==0);
	
	delete bfs;

	delete ssQuery;
	if (testing_result==0)
		cout << "PASS structsearch unittest." << endl;
	return testing_result!=0;

}
