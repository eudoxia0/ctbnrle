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
#include "famscore.h"
#include "markovdyn.h"
#include "ss.h"
#include "multirv.h"
#include <iostream>


namespace ctbn {

using namespace std;

FamScore::FamScore(double nTrans, double aTime, SuffStatsQuery* ssQ) : 
	numTrans(nTrans), amtTime(aTime), ssQuery(ssQ) {
}

FamScore::~FamScore() {
}

// Subclasses below
// ================

BNFamScore::BNFamScore(double nTrans, SuffStatsQuery* ssQ) :
	FamScore(nTrans,0,ssQ) {
}

CTBNFamScore::CTBNFamScore(double nTrans, double aTime, SuffStatsQuery* ssQ) :
	FamScore(nTrans,aTime,ssQ) {
}


//---

BNFamScore::~BNFamScore() {
}

CTBNFamScore::~CTBNFamScore() {
}

//---

double BNFamScore::GetScore(const Context& v, const Context& cv) {
	MultiRV* test = new MultiRV(v, cv);
	ssQuery->SetInitRV(test);
	SS* ss = ssQuery->GetP0SS();
	double ret = test->GetScore(numTrans,ss);
	delete ss;
	ssQuery->SetInitRV(NULL);
	delete test;
	return ret;
}

double CTBNFamScore::GetScore(const Context& v, const Context& cv) {
	MarkovDyn* test = new MarkovDyn(v, cv);
	ssQuery->SetDynamics(test);
	SS* ss = ssQuery->GetDSS();
	double ret = test->GetScore(numTrans,amtTime,ss);
	delete ss;
	ssQuery->SetDynamics(NULL);
	delete test;
	return ret;
}

} // end of ctbn namespace
