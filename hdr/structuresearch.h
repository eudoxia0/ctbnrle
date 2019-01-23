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
#ifndef CTBNRLE_STRUCTURESEARCH_H
#define CTBNRLE_STRUCTURESEARCH_H

#include "context.h"
#include "structure.h"
#include "famscore.h"
#include <vector>



namespace ctbn {

// An abstract class for performing structure searches.
// Constructed using a Context containing the relevant variables
// and a FamScore object that takes care of the scoring.
// (see famscore.h for more details.)

class StructureSearch {
  public:
	// A generic node to keep track of the different 
	// Domains and CondDomains in the network
	class Node {
	  public:
		Node(const Context &var, const Context &cvar) :
			v(var), cv(cvar) {}
		Node* Clone() const {return new Node(*this);}
		const Context &Domain() const {return v;}
		const Context &CondDomain() const {return cv;}
	  private:
		Context v;
		Context cv;
	};

	// This assumes each variable in the context is a node in the network
	// Also, takes in a scoring object to choose between different
	// models (BN and CTBNDyn). This scoring object's pointer is
	// owned by the caller, but must remain valid for the lifetime of
	// the StructureSearch object.
	StructureSearch(const Context &c, FamScore* fs);

	inline virtual ~StructureSearch() throw() {}

	virtual Structure LearnStructure() = 0;

	// Use the node vector below to fill a structure object
	void GetStructure(Structure &s) const;
  protected:
	// Stores each variable of the Context passed in
	std::vector<Context> vars;

	// After running LearnStructure(), this contains the nodes that 
	// describe the learned network
	std::vector<Node *> nodes;

	// Mappings from nodeids to varids and vice versa
	mutable std::vector<int> node2var;
	mutable std::map<int,int> var2node;

	// This class does not own this pointer
	FamScore* fScore;
};

} // end of ctbn namespace

#endif
