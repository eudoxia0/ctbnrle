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
#ifndef CTBNRLE_BRUTESTRUCTURESEARCH_H
#define CTBNRLE_BRUTESTRUCTURESEARCH_H

#include "structuresearch.h"



namespace ctbn {

// Class to perform exhaustive "brute-force" structure search.

class BruteStructureSearch : public StructureSearch {
  public:
  	// The FamScore pointer is not owned by this class. It is owned by
	// the caller of the constructor.
	BruteStructureSearch(const Context &c, FamScore* fs, int mp = 1);
	inline virtual ~BruteStructureSearch() throw() {}
	virtual Structure LearnStructure();
	void SetMaxParents(int mp) {maxParents = mp;}
  protected:
	int maxParents;
};

} // end of ctbn namespace

#endif
