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
#include "structuresearch.h"
#include "matrix.h"


namespace ctbn {

using namespace std;

StructureSearch::StructureSearch(const Context &c, FamScore* fs) 
 : fScore(fs) {
	vector<int> varL = c.VarList();
	for(unsigned int i=0; i<varL.size(); i++) {
		Context cur;
		cur.AddVar(varL[i], c.Cardinality(varL[i]));
		vars.push_back(cur);
	}

	unsigned int numNodes = vars.size();	
	for(unsigned int i=0; i<numNodes; i++)
		node2var.push_back(vars[i].MinVar());
	for(unsigned int i=0; i<numNodes; i++)
		var2node.insert(make_pair(node2var[i],i));
}

void StructureSearch::GetStructure(Structure &s) const {
	s.var2node = var2node;
	s.node2var = node2var;

	map<int, vector<int> > children;

	int n = nodes.size();
	for(int j=0;j<n;j++) {
		vector<int> par = nodes[j]->CondDomain().VarList();
		for(unsigned int k=0;k<par.size();k++)
			children[par[k]].push_back(j);
	}
	
	for(int j=0;j<n;j++) {
		if(children.find(node2var[j]) == children.end())
			children.insert(make_pair(node2var[j],
				vector<int>()));
	}

	s.children = children;

	s.parentlist.clear();
	for (int i=0; i<n; i++) {
		const Context &c = nodes[i]->CondDomain();
		s.parentlist.push_back(c.VarList());
	}	
}

} // end of ctbn namespace
