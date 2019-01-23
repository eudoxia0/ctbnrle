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
#ifndef CTBNRLE_SEARCHQUEUE_H
#define CTBNRLE_SEARCHQUEUE_H

#include <iostream>
#include <deque>



namespace ctbn {

// A priority queue for the GraphEditSearch algorithm.
// All the actions possible in the search at each node in the 
// search tree are enqueued into this.

// There is no dequeue operation. Instead, a remove operation used will
// remove all actions associated with a specific child. We can do this 
// since taking a certain action with a node invalidates the scores
// associated with that node. The possible operations also change.

class SearchQueue {
  public:
	enum searchop {ADD,REMOVE,REVERSE};
	// a class to store an action
	// pretty much just a structure
	class Action {
	  public:
		Action(int ch, int pa, searchop op, double g, double gr) :
			child(ch), parent(pa), oper(op), 
			gain(g), gainrev(gr) {}
		Action() { }

		friend std::ostream& operator << (std::ostream& out, 
						const Action &rhs) {
			return out << "(C=" << rhs.child 
					<< ",P=" << rhs.parent
					<< ",O=" << rhs.oper 
					<< ",G=" << rhs.gain 
					<< ",GR=" << rhs.gainrev << ")";
		}
			
		int child;
		int parent;
		searchop oper;
		double gain;
		double gainrev;		
	  };

	// create a queue with a maximum of maxn actions
	SearchQueue(int maxn, int nvars);

	~SearchQueue();

	// place the earliest action in a
	// return false if there are no actions
	// Note:  this action is NOT removed from the queue
	bool Head(Action &a);

	// add an action to the queue
	void Add(const Action &a);

	// remove the actions of a specific node
	void Remove(int child);

	// empty the entire queue
	void Clear();

	void Print();

	void PrintPlaceList();

  private:
	int Heapify(int i, double tempk);

	Action *heap;
	std::deque <int>* places;
	int n;
	int nv;
};

} // end of ctbn namespace

#endif
