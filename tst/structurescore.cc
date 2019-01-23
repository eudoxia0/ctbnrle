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
#include "ctbn.h"

#include <vector>
#include <fstream>

using namespace std;
using namespace ctbn;

const int TRIALS(50);
const double NUM_TRAJ(40);
const double TRAJ_LEN(2);
double PRIOR_TRANS(1);
double PRIOR_TIME(1);

void score(const Markov& m, const vector <Trajectory>& trVec,
		double &results);

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
  InitParams(argc, argv);

  Context nullc;
  Context eatc, fullc, hungryc, uptakec, concentrationc, barometerc, painc, drowsyc;
  eatc.AddVar(20,2);
  fullc.AddVar(41,3);
  hungryc.AddVar(12,2);
  uptakec.AddVar(63,2);
  concentrationc.AddVar(94,3);
  barometerc.AddVar(85,3);
  painc.AddVar(76,2);
  drowsyc.AddVar(37,2);

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
  vectr peat(2, 0.0);
  peat[0] = 1;
  peat[1] = 0;
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
  pfull[0] = 1;
  //pfull[1] = 0.25;
  //pfull[2] = 0.25;
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
  phungry[0] = 1;
  hungryrv[0].SetDist(phungry);

  //------uptake--------
  uptake(0)->Intensity()[0][0] = -0.0;
  uptake(0)->Intensity()[0][1] = 0.0;
      uptake(0)->Intensity()[1][0] = 0.5;
      uptake(0)->Intensity()[1][1] = -0.5;

  vectr puptake(2, 0.0);
  puptake[1] = 1;
  uptakerv[0].SetDist(puptake);

  //------concentration------
  Instantiation i(Context(fullc, uptakec));
  i.SetVal(63,0); i.SetVal(41,0);
  concentration(i)->Intensity()[0][0] = -0.02;
  concentration(i)->Intensity()[0][1] = 0.01;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.25;
      concentration(i)->Intensity()[1][1] = -0.26;
      concentration(i)->Intensity()[1][2] = 0.01;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.5;
          concentration(i)->Intensity()[2][2] = -0.51;
  i.SetVal(63,1); i.SetVal(41,0);
  concentration(i)->Intensity()[0][0] = -2.01;
  concentration(i)->Intensity()[0][1] = 2.0;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.01;
      concentration(i)->Intensity()[1][1] = -1.01;
      concentration(i)->Intensity()[1][2] = 1.0;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.01;
          concentration(i)->Intensity()[2][2] = -0.02;
  i.SetVal(63,0); i.SetVal(41,1);
  concentration(i)->Intensity()[0][0] = -0.02;
  concentration(i)->Intensity()[0][1] = 0.01;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.25;
      concentration(i)->Intensity()[1][1] = -0.26;
      concentration(i)->Intensity()[1][2] = 0.01;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.5;
          concentration(i)->Intensity()[2][2] = -0.51;
  i.SetVal(63,1); i.SetVal(41,1);
  concentration(i)->Intensity()[0][0] = -1.01;
  concentration(i)->Intensity()[0][1] = 1.0;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.01;
      concentration(i)->Intensity()[1][1] = -2.01;
      concentration(i)->Intensity()[1][2] = 2.0;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.01;
          concentration(i)->Intensity()[2][2] = -0.02;

  i.SetVal(63,0); i.SetVal(41,2);
  concentration(i)->Intensity()[0][0] = -0.02;
  concentration(i)->Intensity()[0][1] = 0.01;
  concentration(i)->Intensity()[0][2] = 0.01;
      concentration(i)->Intensity()[1][0] = 0.25;
      concentration(i)->Intensity()[1][1] = -0.26;
      concentration(i)->Intensity()[1][2] = 0.01;
          concentration(i)->Intensity()[2][0] = 0.01;
          concentration(i)->Intensity()[2][1] = 0.5;
          concentration(i)->Intensity()[2][2] = -0.51;

  i.SetVal(63,1); i.SetVal(41,2);
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
  pconcentration[0] = 1;
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
  pbarometer[1] = 1;
  barometerrv[0].SetDist(pbarometer);

  //------pain------
  Instantiation j(Context(barometerc, concentrationc));
  j.SetVal(94,0); j.SetVal(85,0);
  pain(j)->Intensity()[0][0] = -6.0;
  pain(j)->Intensity()[0][1] = 6.0;
      pain(j)->Intensity()[1][0] = 0.1;
      pain(j)->Intensity()[1][1] = -0.1;

  j.SetVal(94,1); j.SetVal(85,0);
  pain(j)->Intensity()[0][0] = -1.0;
  pain(j)->Intensity()[0][1] = 1.0;
      pain(j)->Intensity()[1][0] = 0.1;
      pain(j)->Intensity()[1][1] = -0.1;
  j.SetVal(94,2); j.SetVal(85,0);
  pain(j)->Intensity()[0][0] = -6.0;
  pain(j)->Intensity()[0][1] = 6.0;
      pain(j)->Intensity()[1][0] = 0.1;
      pain(j)->Intensity()[1][1] = -0.1;
  j.SetVal(94,0); j.SetVal(85,1);
  pain(j)->Intensity()[0][0] = -1.0;
  pain(j)->Intensity()[0][1] = 1.0;
      pain(j)->Intensity()[1][0] = 0.3;
      pain(j)->Intensity()[1][1] = -0.3;
  j.SetVal(94,1); j.SetVal(85,1);
  pain(j)->Intensity()[0][0] = -0.3;
  pain(j)->Intensity()[0][1] = 0.3;
      pain(j)->Intensity()[1][0] = 0.3;
      pain(j)->Intensity()[1][1] = -0.3;
  j.SetVal(94,2); j.SetVal(85,1);
  pain(j)->Intensity()[0][0] = -1.0;
  pain(j)->Intensity()[0][1] = 1.0;
      pain(j)->Intensity()[1][0] = 0.3;
      pain(j)->Intensity()[1][1] = -0.3;
  j.SetVal(94,0); j.SetVal(85,2);
  pain(j)->Intensity()[0][0] = -0.01;
  pain(j)->Intensity()[0][1] = 0.01;
      pain(j)->Intensity()[1][0] = 2.0;
      pain(j)->Intensity()[1][1] = -2.0;
  j.SetVal(94,1); j.SetVal(85,2);
  pain(j)->Intensity()[0][0] = -0.01;
  pain(j)->Intensity()[0][1] = 0.01;
      pain(j)->Intensity()[1][0] = 2.0;
      pain(j)->Intensity()[1][1] = -2.0;
  j.SetVal(94,2); j.SetVal(85,2);
  pain(j)->Intensity()[0][0] = -0.01;
  pain(j)->Intensity()[0][1] = 0.01;
      pain(j)->Intensity()[1][0] = 2.0;
      pain(j)->Intensity()[1][1] = -2.0;

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
  pdrowsy[0] = 1;
  drowsyrv[0].SetDist(pdrowsy);


  CTBNDyn ctbndyn;
  ctbndyn.AddNode(eat.Clone());
  ctbndyn.AddNode(full.Clone());
  ctbndyn.AddNode(hungry.Clone());
  ctbndyn.AddNode(uptake.Clone());
  ctbndyn.AddNode(concentration.Clone());
  ctbndyn.AddNode(barometer.Clone());
  ctbndyn.AddNode(pain.Clone());
  ctbndyn.AddNode(drowsy.Clone());

  BN bn;
  bn.AddNode(eatrv.Clone());
  bn.AddNode(fullrv.Clone());
  bn.AddNode(hungryrv.Clone());
  bn.AddNode(uptakerv.Clone());
  bn.AddNode(concentrationrv.Clone());
  bn.AddNode(barometerrv.Clone());
  bn.AddNode(painrv.Clone());
  bn.AddNode(drowsyrv.Clone());

  Markov m(bn.Clone(), ctbndyn.Clone());

 

	// ==========================
	// Building simpler structure

	MarkovDyn eatsimp(eatc, nullc);
	MarkovDyn fullsimp(fullc, eatc);
	MarkovDyn hungrysimp(hungryc, fullc);
	MarkovDyn uptakesimp(uptakec, nullc);
	MarkovDyn concentrationsimp(concentrationc, uptakec);
	MarkovDyn barometersimp(barometerc, nullc);
	MarkovDyn painsimp(painc, Context(barometerc, concentrationc) );
	MarkovDyn drowsysimp(drowsyc, nullc);
	
	CTBNDyn ctbndynsimp;
	ctbndynsimp.AddNode(eatsimp.Clone());
	ctbndynsimp.AddNode(fullsimp.Clone());
	ctbndynsimp.AddNode(hungrysimp.Clone());
	ctbndynsimp.AddNode(uptakesimp.Clone());
	ctbndynsimp.AddNode(concentrationsimp.Clone());
	ctbndynsimp.AddNode(barometersimp.Clone());
	ctbndynsimp.AddNode(painsimp.Clone());
	ctbndynsimp.AddNode(drowsysimp.Clone());

	Markov msimple(bn.Clone(), ctbndynsimp.Clone());

	//==============================
	//Building independent structure
	
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

	Markov mindep(bn.Clone(), ctbndynindep.Clone());

	//==============================
	//Building complex structure
	
	MarkovDyn eatcomplex(eatc, hungryc);
	MarkovDyn fullcomplex(fullc, eatc);
	MarkovDyn hungrycomplex(hungryc, fullc);
	MarkovDyn uptakecomplex(uptakec, eatc);
	Context fubhc(Context(fullc,uptakec),Context(barometerc,hungryc));
	MarkovDyn concentrationcomplex(concentrationc, fubhc);
	MarkovDyn barometercomplex(barometerc, nullc);
	Context bcdc(Context(barometerc,concentrationc),drowsyc);
	MarkovDyn paincomplex(painc, bcdc);
	MarkovDyn drowsycomplex(drowsyc, concentrationc);
	
	CTBNDyn ctbndyncomplex;
	ctbndyncomplex.AddNode(eatcomplex.Clone());
	ctbndyncomplex.AddNode(fullcomplex.Clone());
	ctbndyncomplex.AddNode(hungrycomplex.Clone());
	ctbndyncomplex.AddNode(uptakecomplex.Clone());
	ctbndyncomplex.AddNode(concentrationcomplex.Clone());
	ctbndyncomplex.AddNode(barometercomplex.Clone());
	ctbndyncomplex.AddNode(paincomplex.Clone());
	ctbndyncomplex.AddNode(drowsycomplex.Clone());

	Markov mcomplex(bn.Clone(), ctbndyncomplex.Clone());


	//==============================
	//Building reverse structure
	
	MarkovDyn eatrev(eatc, fullc);
	MarkovDyn fullrev(fullc, hungryc);
	MarkovDyn hungryrev(hungryc, eatc);
	MarkovDyn uptakerev(uptakec, concentrationc);
	MarkovDyn concentrationrev(concentrationc, Context(painc,drowsyc));
	MarkovDyn barometerrev(barometerc, painc);
	MarkovDyn painrev(painc, nullc);
	MarkovDyn drowsyrev(drowsyc, nullc);
	
	CTBNDyn ctbndynrev;
	ctbndynrev.AddNode(eatrev.Clone());
	ctbndynrev.AddNode(fullrev.Clone());
	ctbndynrev.AddNode(hungryrev.Clone());
	ctbndynrev.AddNode(uptakerev.Clone());
	ctbndynrev.AddNode(concentrationrev.Clone());
	ctbndynrev.AddNode(barometerrev.Clone());
	ctbndynrev.AddNode(painrev.Clone());
	ctbndynrev.AddNode(drowsyrev.Clone());

	Markov mrev(bn.Clone(), ctbndynrev.Clone());


	//=============
	//SCORE TESTING
	//=============
	
	vector <double> zeros;
	for(int i=0; i<=NUM_TRAJ; i++) zeros.push_back(0.0);

	double scores[5];

	for(int k=0; k<TRIALS; k++) {
		//Generate trajectories
		vector <Trajectory> trVec;
		for(int i=0; i<NUM_TRAJ; i++) {
			Trajectory tr; tr.SetEndTime(TRAJ_LEN);
			m.Sample(tr);
			trVec.push_back(tr);
		}

		score(m,trVec,scores[0]);
		score(msimple,trVec,scores[1]);
		score(mindep,trVec,scores[2]);
		score(mcomplex,trVec,scores[3]);
		score(mrev,trVec,scores[4]);
	}

	EXPECT(scores[0]>scores[1]);
	EXPECT(scores[0]>scores[2]);
	EXPECT(scores[0]>scores[3]);
	EXPECT(scores[0]>scores[4]);
	if (testing_result==0)
		cout << "PASS structurescore unittest." << endl;
	return testing_result!=0;
}

void normalizeVec(vector<double> &results) {
	for(unsigned int i=0; i<results.size();i++) results[i]=results[i]/TRIALS;
}


void score(const Markov& m, const vector <Trajectory>& trVec, 
		double &results) {
	SS *ss;
	vector<double> weights(trVec.size(),1.0);

	ss = m.SuffStats(trVec,weights);
	results += m.GetScore(PRIOR_TRANS,PRIOR_TIME,ss)/trVec.size();
	delete ss;
}
