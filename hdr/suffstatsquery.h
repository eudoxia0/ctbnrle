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
#ifndef CTBNRLE_SUFFSTATSQUERY_H
#define CTBNRLE_SUFFSTATSQUERY_H

// A class to perform queries for sufficient statistics
// The object is to be passed into a function as an empty object. 
// The member variables are set by the function using this class.

#include "ss.h"
#include "rv.h"
#include "dynamics.h"
#include "trajectory.h"
#include <vector>

namespace ctbn {

class SuffStatsQuery {

public:
	SuffStatsQuery();
	virtual ~SuffStatsQuery();
	virtual SS* GetP0SS() = 0;
	virtual SS* GetDSS() = 0;

	// The caller of these functions own the pointers.
	void SetInitRV(const RV* p0) {this->p0 = p0;}
	void SetDynamics(const Dynamics* d) {this->d = d;}
	virtual void SetData(const std::vector<Trajectory>* data) {tr = data;}
protected:
	// This class does not own any of these pointers.
	const RV* p0;
	const Dynamics* d;
	const std::vector<Trajectory>* tr;
};

} // end of ctbn namespace

#endif
