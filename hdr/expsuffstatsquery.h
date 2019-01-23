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
#ifndef CTBNRLE_EXPSUFFSTATSQUERY_H
#define CTBNRLE_EXPSUFFSTATSQUERY_H

// A class for performing queries for sufficient statistics where we have 
// incomplete trajectories.

#include "suffstatsquery.h"
#include "trajectory.h"
#include "process.h"
#include "inference.h"
#include <vector>

namespace ctbn {

class ExpSuffStatsQuery : public SuffStatsQuery {
public:
	// Expected sufficient statistics need a process for its expectation,
	// along with an inference method. The pointer of the inference
	// object is owned by the caller.
	ExpSuffStatsQuery(Inference* i, bool cloneInf = false);

	virtual ~ExpSuffStatsQuery();

	// In the event of a high number of queries with the same process 
	// and data, caching the inference objects eliminates the need to 
	// set the same trajectories multiple times.
	void InitInfCache();

	// changes the process (associated with the Inference object
	// passed into the constructor) from which the expectaction is to
	// be taken
	// (SetInitRV and SetDynamics from SuffStatsQuery change the
	//  process who sufficient statistics are taken)
	void SetProcess(const Process* p);	

	// Calculate and return sufficient statistics.
	virtual SS* GetP0SS();
	virtual SS* GetDSS();
protected:
	Inference* inf;
	std::vector <Inference *> infCache;
	bool useInfCache;
};
} // end of ctbn namespace

#endif
