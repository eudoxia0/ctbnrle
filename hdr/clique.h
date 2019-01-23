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
#ifndef CTBNRLE_CLIQUE_H_
#define CTBNRLE_CLIQUE_H_

#include "bn.h"
#include "dynamics.h"
#include "context.h"
#include "ctbndyn.h"
#include "markov.h"
#include "markovdyn.h"
#include "random.h"

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>




namespace ctbn {

// Class to build a clique tree from a set of factors,
// either defined from a BN, a CTBNDyn, or both (Markov process)
class CliqueTree: public StreamObj{
public:

	CliqueTree();

	//build the clique tree for a markov (both CTBNDyn and BN factors)
	void BuildClique(const Markov &process, Random &rand = randomizer);
	//build the clique tree for a CTBNDyn
	// caller retains ownership of dyn
	void BuildClique(const CTBNDyn *dyn, Random &rand = randomizer);
	//build the clique tree for a BN
	// caller retains ownership of bn
	void BuildClique(const BN *bn, Random &rand = randomizer);
    // build the clique truee for Contexts
    void BuildClique(const std::vector<Context> &, Random & rand = randomizer);

	//print the clique tree
	void CliquePrint(std::ostream &os);

	//save the current state
	virtual void SaveOld(std::ostream &os) const;

	//loading the value into the class so it will have all the data
	virtual void LoadOld(std::istream &is);

	//to use for testing:
	//if true is passed in random number generator is ignored
	//and factors are processed in a deterministic order
	void IgnoreRand(bool ignore) {
		unit_test = ignore;
	}

	//change the number of random restarts
	void SetNumRestarts(int num) {
		numrandrestarts = num;
	}

	// stores a node of the clique tree
	class Node {
	public:
		Context vars; // variables that are part of this factor
		std::vector<int> adj; // index of adjacent nodes in tree
	};

	//return the clique tree
	std::vector<Node> ReturnCliques() const;

	virtual ~CliqueTree();

private:
	//all the private functions with pointers own their pointers
	class CliqueNode {
		friend class CliqueTree;
	public:
		CliqueNode();
		void AddVar(const int &);
		void print (std::ostream &);
		~CliqueNode();
	private:
		std::set <int> connected;
		std::vector <CliqueNode *> next;
		CliqueNode * prev;

		SERIAL_START(CliqueNode)
			SERIAL_VAR(std::set<int>,connected)
			SERIAL_VAR(std::vector<CliqueNode*>,next)
		SERIAL_END
	public:
		void serial_preload();
		void serial_postload();
	};
	void CliqueTrees(Random &rand);
	void CliqueTrees(int);
	void CliqueNodePrint(CliqueNode *, int, std::ostream &);
	static void Destroy(CliqueNode *);
	void ReturnClique(CliqueNode *, int, std::vector<Node> &) const;
	void InsertCliqueNode(CliqueNode *, int, std::vector<Node> &) const;
	Context Clique2Context(CliqueNode *) const;
	void CliqueSave(std::ostream &os, CliqueNode * print, int size) const;
	void Eliminate();
	bool EliminateDups(CliqueNode *);
	void EliminatePrev(CliqueNode *);
	static void RemoveElem(std::vector <int> &, int);
	void prints(std::set <int>, std::ostream &) const;
	double MaxSize() const;
	double MaxSize(CliqueNode *) const;
	void Dequeue(std::vector <int> remove, CliqueNode * newitem);

	void BuildIt(Random &);
	void AddFromBN(const BN *);
	void AddFromCTBNDyn(const CTBNDyn *);
    void AddFromContext(const std::vector <Context> &);

	std::vector <CliqueNode> data;
	std::vector <CliqueNode> save_data;
	std::vector <CliqueNode *> clique_built;

	// these pointers at the bottom are owned by this class
	CliqueNode * Head_compare;
	int numrandrestarts;
	CliqueNode * Head;
	bool unit_test;
    Context vars;
/*
	int VarSize(int varid) const;
	void AddVarSize(int varid, int size);
	// maps varid -> size
	std::map<int,int> context_size;
*/
	SERIAL_START(CliqueTree)
		SERIAL_VAR(Context,vars)
		SERIAL_VAR(int,numrandrestarts)
		SERIAL_VAR(bool,unit_test)
		SERIAL_VAR(std::vector<CliqueNode>,data)
	SERIAL_END
public:
	void serial_preload();
};

} // end of ctbn namespace

#endif /* CLIQUE_H_ */
