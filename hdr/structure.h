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
#ifndef CTBNRLE_STRUCTURE_H
#define CTBNRLE_STRUCTURE_H

#include <vector>
#include "ctbndyn.h"
#include "context.h"


namespace ctbn {

//The Structure class copies node2var, var2node and children
//from class CTBNDyn. It can be used when the structure 
//cannot be accessed directly.
//It represents the graph structure of a factored model, like a BN
// or CTBNDyn
//It also understands that the variable IDs involved might not be
//consecutive and so it keeps a mapping from varid <-> array index
class Structure {
friend class CTBNDyn;
friend class BN;
friend class StructureSearch;
friend class GraphEditSearch;
public: 
	Structure();
	virtual ~Structure();
	// these return the parents or children for a variable
	// if isVarID is false (default), then child/parent is the node index
	// if isVarID is true, then child/parent is the variable number
	// the returned vector is of *node* indexes (not variable IDs)
	const std::vector<int>& GetParents(int child, bool isVarID=false) const;
	const std::vector<int>& GetChildren(int parent, bool isVarID=false) const;

	// these convert between the node index and the variable ID
	int Node2Var(int index) const;
	int Var2Node(int index) const;
	const std::vector<std::vector<int> >& GetAdjLists() const;
	int HammingDist(const Structure &other);
	std::vector<int> GetVarList() const;
	int GetNumParams(const Context& c) const;
  	bool IsCyclic() const;
	void Print(std::ostream& os) const;
protected:

	//mapping from index into "nodes" to the variable id represented
	std::vector<int> node2var; 	

	//parents of nodes
	//cache the parents variable id of each node index
	std::vector<std::vector<int> > parentlist;

	//mapping from variable id to the index into "nodes"
	std::map<int, int> var2node;

	// cache of
	// the indexes into "nodes" 
	// of the children for a particular variable id
	// (Note this one is slightly odd as it maps a variable id
	// into a list of node-indexes (not back into varids)
	std::map<int,std::vector<int> > children;

	SERIAL_START(Structure)
		SERIAL_VAR(std::vector<int>,node2var)
		SERIAL_VAR(std::vector<std::vector<int> >,parentlist)
	SERIAL_END
public:
	void serial_postload();
};
} // end of ctbn namespace

#endif
