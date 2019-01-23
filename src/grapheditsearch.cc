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
#include "grapheditsearch.h"


namespace ctbn {

using namespace std;

GraphEditSearch::GraphEditSearch(const Context& c, FamScore* fs, bool noCyc) 
 : StructureSearch(c, fs), noCycles(noCyc), learnCalled(false) {
	int numNodes = vars.size();
	sQ = new SearchQueue(numNodes*2*(numNodes-1), numNodes);
	Initialize();
}

void GraphEditSearch::SetStructure(const Structure& s) {
	for(unsigned int i=0; i<vars.size(); i++) {
		const vector<int>& curParList = s.parentlist[i];
		int curvarid = s.node2var[i];
		map<int,int>::iterator it = var2node.find(curvarid);
		int ssnodeid = it != var2node.end() ? it->second : -1;
		if(ssnodeid == -1) continue;
		delete nodes[ssnodeid];
		Context cv;
		for(unsigned int j=0; j<curParList.size(); j++) {
			it = var2node.find(curParList[j]);
			int pnodeid = it != var2node.end() ? it->second : -1;
			if(pnodeid == -1) continue;
			cv = cv + vars[pnodeid];
		}
		double s = fScore->GetScore(vars[ssnodeid], cv);
		nodes[ssnodeid] = new Node(vars[ssnodeid], cv);
		scores[ssnodeid] = s;
	}
}

void GraphEditSearch::Initialize() {
	scores.clear();
	for(unsigned int i=0; i<nodes.size(); i++)
		delete nodes[i];
	nodes.clear();
	Context nullc;
	for(unsigned int i=0; i<vars.size(); i++) {
		double s = fScore->GetScore(vars[i],nullc);
		nodes.push_back(new Node(vars[i],nullc));
		scores.push_back(s);
	}
}

void GraphEditSearch::EnqueueOps(int i) {
	vector<int> varLvid = nodes[i]->CondDomain().VarList();
	vector<int> varL;

	// Map the VarList into a new vector using node ids instead.
	for(unsigned int j=0; j<varLvid.size(); j++)
		varL.push_back(var2node[varLvid[j]]);

	// Enqueue all ADD operations
	for(unsigned int j=0; j<vars.size(); j++) {
		if (static_cast<unsigned int>(i) == j ||
         find(varL.begin(),varL.end(),j) != varL.end())
			continue;

		if(noCycles && ProducesCycle(i,j)) continue;

		double s = fScore->GetScore(nodes[i]->Domain(), 
					nodes[i]->CondDomain() + vars[j]);
		double deltas = s - scores[i];

		sQ->Add(SearchQueue::Action(i,j,
					SearchQueue::ADD,deltas,0));
	}

	// Enqueue all REMOVE and REVERSE operations
	for(unsigned int j=0; j<varL.size(); j++) {
		int k = varL[j];
		double srem = fScore->GetScore(nodes[i]->Domain(), 
					nodes[i]->CondDomain() - vars[k]);
		double deltasrem = srem - scores[i];

		sQ->Add(SearchQueue::Action(i,k,
				SearchQueue::REMOVE,deltasrem,0));		
		if(noCycles && ProducesCycleRev(k,i)) continue;
		
		double saddr = fScore->GetScore(nodes[k]->Domain(), 
					nodes[k]->CondDomain() + vars[i]);
		double deltasaddr = saddr - scores[k];

		sQ->Add(SearchQueue::Action(i,k,
				SearchQueue::REVERSE,
				deltasrem,deltasaddr));
	}
}

bool GraphEditSearch::ProducesCycle(int i, int j) {
	Node* temp = nodes[i];
	nodes[i] = new Node(nodes[i]->Domain(), 
				nodes[i]->CondDomain() + vars[j]);
	Structure sCycleCheck;
	GetStructure(sCycleCheck);
	delete nodes[i];
	nodes[i] = temp;
	return sCycleCheck.IsCyclic();
}

bool GraphEditSearch::ProducesCycleRev(int i, int j) {
	Node* tempi = nodes[i];
	Node* tempj = nodes[j];
	nodes[i] = new Node(nodes[i]->Domain(), 
				nodes[i]->CondDomain() + vars[j]);
	nodes[j] = new Node(nodes[j]->Domain(),
				nodes[j]->CondDomain() - vars[i]);
	Structure sCycleCheck;
	GetStructure(sCycleCheck);
	delete nodes[i];
	delete nodes[j];
	nodes[i] = tempi;
	nodes[j] = tempj;
	return sCycleCheck.IsCyclic();
}

Structure GraphEditSearch::LearnStructure() {
	Context nullc;
	int numNodes = vars.size();

	if(learnCalled) Initialize();

	for(int i=0; i<numNodes; i++)
		EnqueueOps(i);
	
	SearchQueue::Action nextAct;
	sQ->Head(nextAct);
	
	while(nextAct.gain + nextAct.gainrev > 0) {
		int ch = nextAct.child;
		int pa = nextAct.parent;
		SearchQueue::searchop op = nextAct.oper;
		if(op == SearchQueue::ADD) {
			Node* newNode = 
				new Node(nodes[ch]->Domain(),
						nodes[ch]->CondDomain() +
						vars[pa]);
			delete nodes[ch];
			nodes[ch] = newNode;
			scores[ch] += nextAct.gain;
		}

		else if(op == SearchQueue::REMOVE) {
			Node* newNode = 
				new Node(nodes[ch]->Domain(),
						nodes[ch]->CondDomain() -
						vars[pa]);
			delete nodes[ch];
			nodes[ch] = newNode;
			scores[ch] += nextAct.gain;
		}

		else if(op == SearchQueue::REVERSE) {
			Node* newNode = 
				new Node(nodes[ch]->Domain(),
						nodes[ch]->CondDomain() -
						vars[pa]);
			delete nodes[ch];
			nodes[ch] = newNode;
			scores[ch] += nextAct.gain;

			newNode = new Node(nodes[pa]->Domain(),
						nodes[pa]->CondDomain() +
						vars[ch]);
			delete nodes[pa];
			nodes[pa] = newNode;
			scores[pa] += nextAct.gainrev;
		}
		
		if(noCycles) {
			sQ->Clear(); // complicated to determine what would
				// now be allowed (or now disallowed) due to 
				// cycles... for now, just recalculate everything
			for(unsigned int i = 0; i < vars.size(); i++)
				EnqueueOps(i);
		}
		else {
			sQ->Remove(ch);
			if (op == SearchQueue::REVERSE)
				sQ->Remove(pa);
			EnqueueOps(ch);
			if(op == SearchQueue::REVERSE) 
				EnqueueOps(pa);
		}
			sQ->Head(nextAct);
	}
	Structure s;
	GetStructure(s);
	learnCalled = true;
	return s;
}

} // end of ctbn namespace
