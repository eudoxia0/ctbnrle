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
#include "expsuffstatsquery.h"
#include "nullptr03.h"


namespace ctbn {

using namespace std;

ExpSuffStatsQuery::ExpSuffStatsQuery(Inference* i, bool cloneInf)
	: SuffStatsQuery(), inf(i), useInfCache(cloneInf) {
}

void ExpSuffStatsQuery::InitInfCache() {
	if(infCache.size() == tr->size()) return;
	for(unsigned int i = 0; i < tr->size(); i++) {
		infCache.push_back(inf->Clone());
		infCache[i]->SetTrajectory(&(*tr)[i]);
	}
}

void ExpSuffStatsQuery::SetProcess(const Process* p) {
	inf->SetProcess(p);
	for(unsigned int i=0;i<infCache.size();i++)
		delete infCache[i];
	infCache.resize(0);
}

SS* ExpSuffStatsQuery::GetP0SS() {
	if(p0 && tr) {
		SS* ss = p0->BlankSS();
		if(useInfCache) {
			InitInfCache();
			for(unsigned int i = 0; i < infCache.size(); i++) {
				infCache[i]->AddExpSuffStats(p0,ss);
			}
		}
		else {
			for(unsigned int i = 0; i < tr->size(); i++) {
				inf->SetTrajectory(&(*tr)[i]);
				inf->AddExpSuffStats(p0,ss);
			}
		}
		return ss;
	}
	else {
		cerr << "p0 and/or data are not initialized" << endl;
		exit(0);
	}	
}

SS* ExpSuffStatsQuery::GetDSS() {
	if(d && tr) {
		SS* ss = d->BlankSS();
		if(useInfCache) {
			InitInfCache();
			for(unsigned int i = 0; i < infCache.size(); i++) {
				infCache[i]->AddExpSuffStats(d,ss);
			}
		}
		else {
			const vector<Trajectory> &trj = *tr;
			for(unsigned int i = 0; i < trj.size(); i++) {
				inf->SetTrajectory(&trj[i]);
				inf->AddExpSuffStats(d,ss);
			}
		}
		return ss;
	}
	else {
		cerr << "d and/or data are not initialized" << endl;
		exit(0);
	}	
}

ExpSuffStatsQuery::~ExpSuffStatsQuery() {
	for(unsigned int i = 0; i < infCache.size(); i++)
		if (infCache[i] != nullptr03) delete infCache[i];
	infCache.clear();
	delete inf;
}

	
} // end of ctbn namespace
