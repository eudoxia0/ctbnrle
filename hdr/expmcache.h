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
#ifndef CTBNRLE_EXPMCACHE_H
#define CTBNRLE_EXPMCACHE_H

#include "contfunction.h"
#include "matrix.h"
#include<vector>

namespace ctbn {

class ExpMCache {
public:
	ExpMCache(const matrix &Q, double maxt, double eps=1e-6);

	void vexpmt(vectr &a, double t) const;
	// not yet implemented:
	//double expmtv(vectr &b, double t) const;

private:
/* VERSION 1:
	// one for each dimension of vectr
	// element i is the result of e_i exp(Qt)
	// where e_i is a vector of all 0s except a 1 in position i 
	std::vector<ContFunction<vectr> > expm;
*/
/* VERSION 2:
 */
	ContFunction<matrix> expm;
};

}
#endif

