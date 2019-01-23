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
#ifndef CTBNRLE_SAMPLEQUEUE_H
#define CTBNRLE_SAMPLEQUEUE_H



namespace ctbn {

// a priority queue just for sampling quickly
// It keeps track of a set of event associated with a set of variables
// will return the next event (earliest event) and can remove an event
// associated with a particular variable even if that event is not next

// Different from STL's pqueue in that we need to be able
// to remove an event from the middle of the queue without having
// to search for it.  Also, we know the total number of events to be stored
// (this could be done by storing a structure with a copy operator that
//  kept track of its position.  However, this is simpler and faster.)

// Note:  the "var" stored in an event is assume to be one of a
//  set of consecutive integers starting at 0.  Therefore, you may need
//  a map to map your varid onto this value.
class SampleQueue {
  public:
	  // a class to store an event
	  // pretty much just a structure
	  class Event {
	    public:
		  Event(int v,int val, double t) { 
			  var = v; value = val; time = t; 
		  }
		  Event() { }
		  int var;
		  int value;
		  double time;
	  };

	  // create a queue with a maximum of maxn variables
	  SampleQueue(int maxn);

	  ~SampleQueue();

	  // place the earliest event in e
	  // return false if there are no events
	  // Note:  this event is NOT removed from the queue
	  bool Head(Event &e);

	  // add an event to the queue
	  void Add(const Event &e);

	  // remove the event for variable var from the queue
	  void Remove(int var);

  private:
	int Heapify(int i, double tempk);

	Event *heap;
	int *places;
	int n;
};

} // end of ctbn namespace

#endif
