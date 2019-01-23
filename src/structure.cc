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
#include "structure.h"
#include "extramath.h"
#include <set>



namespace ctbn {

using namespace std;

Structure::Structure() {
}

Structure::~Structure() {
}

const vector<int> &Structure::GetChildren(int parent, bool isVarID) const {
	if(!isVarID) parent = node2var[parent];
	map<int, vector<int> >::const_iterator ci = children.find(parent);
	if (ci != children.end()) {
		return ci->second;
	} else {
		static vector<int> empty(0);
		return empty;
	}
}

const vector<int> &Structure::GetParents(int child, bool isVarID) const {
	if(isVarID) child = var2node.find(child)->second;
	if(static_cast<unsigned int>(child) >= parentlist.size()) {
		cerr << "incorrect node index" << endl;
		cerr << "child: " << child << endl;
		exit(1);
	}
	return parentlist[child];
}

const vector<vector<int> > &Structure::GetAdjLists() const {
	return parentlist;
}

int Structure::HammingDist(const Structure &other) {
	const vector<vector <int> > &listA = this->GetAdjLists();
	const vector<vector <int> > &listB = other.GetAdjLists();

	vector<int> vlA = this->GetVarList();
	vector<int> vlB = other.GetVarList();

	set<int> allvarid;
	for(unsigned int i=0;i<vlA.size();i++)
		allvarid.insert(vlA[i]);
	for(unsigned int i=0;i<vlB.size();i++)
		allvarid.insert(vlB[i]);
	

	int ret = 0;
	int count = 0;	

	set<int>::const_iterator it;
	for(it = allvarid.begin(); it != allvarid.end(); ++it) {
		int cur = *it;
		if(this->children.find(cur) == this->children.end()) {
			ret += listB[other.var2node.find(cur)->second].size();
			continue;
		}
		if(other.children.find(cur) == other.children.end()) {
			ret += listA[this->var2node.find(cur)->second].size();
			continue;
		}
		ret += listDiff(listA[this->var2node.find(cur)->second],
				listB[other.var2node.find(cur)->second]);
	}
	return ret;
}

vector<int> Structure::GetVarList() const {
	vector<int> ret;
	for(unsigned int i=0;i<node2var.size();i++)
		ret.push_back(node2var[i]);
	return ret;
}

int Structure::GetNumParams(const Context& c) const {
	int ret = 0;
	for(unsigned int i=0;i<parentlist.size();i++) {
		int numParInst = 1;
		for(unsigned int j=0;j<parentlist[i].size();j++)
			numParInst *= c.Cardinality(parentlist[i][j]);
		ret += numParInst * (c.Cardinality(node2var[i]) - 1);
	}
	return ret;
}

bool Structure::IsCyclic() const {
	vector<int> zeroInDeg;
	vector<vector<int> > pl(parentlist);
	map<int, vector<int> > cl(children);

	for(unsigned int i = 0; i < pl.size(); i++)
		if(pl[i].size() == 0) zeroInDeg.push_back(i);

	while(zeroInDeg.size() != 0) {
		int rem = zeroInDeg.back();
		rem = node2var[rem];
		zeroInDeg.pop_back();
		map<int, vector<int> >::iterator ci = cl.find(rem);
		if(ci == cl.end()) continue;
		
		for(unsigned int i = 0; i < ci->second.size(); i++) {
			int rem2 = ci->second[i];
			vector<int>::iterator it =
				find(pl[rem2].begin(),pl[rem2].end(), rem);
			pl[rem2].erase(it);
			if(pl[rem2].size() == 0)
				zeroInDeg.push_back(rem2);
		}
	}
	int totalsize = 0;
	for(unsigned int i = 0; i < pl.size(); i++)
		totalsize += pl[i].size();
	return totalsize != 0;
} 

int Structure::Node2Var(int index) const {
	return node2var[index];
}

int Structure::Var2Node(int varid) const {
	return var2node.find(varid)->second;
}

void Structure::Print(ostream &os) const {
	os << "Parent list: " << endl;
	for(unsigned int i=0; i<parentlist.size(); i++)
		os << "id " << Node2Var(i) << ": " << parentlist[i] << endl;
}

void Structure::serial_postload() {
	var2node.clear();
	for(unsigned int i=0;i<node2var.size();i++)
		var2node[node2var[i]] = i;
	children.clear();
	for(unsigned int i=0;i<parentlist.size();i++)
		for(unsigned int j=0;j<parentlist[i].size();j++)
			children[node2var[j]].push_back(i);
}

} // end of ctbn namespace
