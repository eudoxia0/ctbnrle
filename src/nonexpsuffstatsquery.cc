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
#include "nonexpsuffstatsquery.h"



namespace ctbn {

using namespace std;

NonExpSuffStatsQuery::NonExpSuffStatsQuery()
	: SuffStatsQuery() {
}

SS* NonExpSuffStatsQuery::GetP0SS() {
	if(p0 && tr) {
		return p0->SuffStats(*tr,vector<double>(tr->size(),1));
	}
	else {
		cerr << "p0 and/or data are not initialized" << endl;
		return NULL;
	}
}

SS* NonExpSuffStatsQuery::GetDSS() {
	if(d && tr) {
		return d->SuffStats(*tr,vector<double>(tr->size(),1));
	}
	else {
		cerr << "d and/or data are not initialized" << endl;
		return NULL;
	}
}

NonExpSuffStatsQuery::~NonExpSuffStatsQuery() {
}	

} // end of ctbn namespace
