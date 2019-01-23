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
#ifndef CTBNRLE_GRAPHEDITSEARCH_H
#define CTBNRLE_GRAPHEDITSEARCH_H

#include "structuresearch.h"
#include "markovdyn.h"
#include "searchqueue.h"
#include <vector>
#include <map>



namespace ctbn {

// Performs the graph edit search algorithm (a greedy search)
// The user has the option to disallow cycles (generally for BNs)

class GraphEditSearch : public StructureSearch {
  public:
  	// The FamScore pointer is owned by the caller of the constructor
	  // but must remain valid for the lifetime of the GraphEditSearch object
	GraphEditSearch(const Context &c, FamScore* fs, bool noCyc=false);
	inline virtual ~GraphEditSearch() throw() { delete sQ; }
	virtual Structure LearnStructure();
	virtual void SetStructure(const Structure& s);
  protected:
	// Called by constructor to setup the nodes used to score
	// and determine initial scores.
	virtual void Initialize();


	// Enqueues the operators performed by the variable of nodeid i
	virtual void EnqueueOps(int i);

	// Checks if adding nodeid j to the parent set of i produces a cycle
	virtual bool ProducesCycle(int i, int j);
	// Checks if adding nodeid j to the parent set of i
	// and removing nodeidi from the parents set of j  produces a cycle
	virtual bool ProducesCycleRev(int i, int j);

	// Keeps the current scores
	std::vector<double> scores;
	
	// Priority queue designed for this algorithm
	// This class owns this pointer and creates the object for it upon
	// construction.
	SearchQueue* sQ;

	bool noCycles;
	
	bool learnCalled;
};

} // end of ctbn namespace

#endif
