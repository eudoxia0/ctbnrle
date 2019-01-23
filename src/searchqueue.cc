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
#include "searchqueue.h"
#include "defines.h"

#include <algorithm>
#include <cmath>



namespace ctbn {

using namespace std;

SearchQueue::SearchQueue(int maxn, int nvars) {
	heap = new Action[maxn];
	places = new deque<int>[nvars];
	n = 0;
	nv = nvars;
}

SearchQueue::~SearchQueue() {
	delete []heap;
	delete []places;
}

bool SearchQueue::Head(Action &a) {
	if (n<=0) return false;
	a = heap[0];
	return true;
}

void SearchQueue::Add(const Action &a) {
	heap[n] = a;
	++n;
	double key = a.gain + a.gainrev;
	int i;
	for(i=n-1;i>0 && key > heap[(i-1)/2].gain + heap[(i-1)/2].gainrev;
	i = (i-1)/2) {
		heap[i] = heap[(i-1)/2];
		deque<int> &list = places[heap[i].child];
		list.erase(find(list.begin(),list.end(),(i-1)/2));
		list.push_front(i);
	}
	heap[i] = a;
	places[a.child].push_front(i);
}

void SearchQueue::Remove(int child) {
	while(places[child].size() != 0) {
		int i = places[child].front();
		places[child].pop_front();
		if (i==0) {
			--n;
			if (n>0) {
				heap[0] = heap[n];
				deque<int> &list = places[heap[0].child];
				list.erase(find(list.begin(),list.end(),n));
				list.push_front(0);
				Heapify(0,heap[0].gain + heap[0].gainrev);
			}
		} else {
			i = Heapify(i,INFINITY);
			--n;
			deque<int> &list = places[heap[i].child];
			list.erase(find(list.begin(),list.end(),i));

			if (i==n) continue;
			Action temp = heap[n];
			while(i>0 && temp.gain + temp.gainrev > 
			heap[(i-1)/2].gain + heap[(i-1)/2].gainrev) {
				heap[i] = heap[(i-1)/2];

				deque<int> &list2 = places[heap[i].child];
				list2.erase(find(list2.begin(),list2.end(),
							(i-1)/2));
				list2.push_front(i);
				
				i = (i-1)/2;
			}
			heap[i] = temp;
			deque<int> &list2 = places[temp.child];
			list2.erase(find(list2.begin(),list2.end(),n));
			list2.push_front(i);
		}
	}
}

int SearchQueue::Heapify(int i, double tempk) {
	int l,r,large,oldlarge;
	double larget;
	oldlarge = -1;

	while(i<n) {
		l = i*2+1;
		r = i*2+2;
		if (l < n && (heap[l].gain + heap[l].gainrev > tempk || 
		isinf(tempk))) {
			large = l;
			larget = heap[l].gain + heap[l].gainrev;
		} else {
			large = i;
			larget = tempk;
		}
		if (r < n && heap[r].gain + heap[r].gainrev > larget) {
			large = r;
			larget = heap[r].gain + heap[r].gainrev;
		}
		if (large==i) {
			deque<int> &list = places[heap[i].child];
			deque<int>::iterator it = 
				find(list.begin(),list.end(),oldlarge);
			if(it != list.end()) list.erase(it);
			list.push_front(i);
			return i;
		}
		Action temp = heap[i];
		heap[i] = heap[large];
		heap[large] = temp;
	
		deque<int> &list = places[temp.child];
		deque<int>::iterator it =
				find(list.begin(),list.end(),i);
		if(it != list.end()) list.erase(it);
		list.push_front(large);

		deque<int> &list2 = places[heap[i].child];
		it = find(list2.begin(),list2.end(),large);
		if(it != list2.end()) list2.erase(it);
		list2.push_front(i);

		i = large;
		oldlarge = large;
	}
	return i; // shouldn't get here
}

void SearchQueue::Clear() {
	n = 0;
	for(int i = 0; i < nv; i++) {
		places[i].clear();
	}
}

void SearchQueue::Print() {
	for(int i = 0; i < n; i++) {
		cout << i << ": " << heap[i] << endl;
	}
}

void SearchQueue::PrintPlaceList() {
	for(int i = 0; i < nv; i++) {
		cout << i << ": ";
		for(unsigned int j = 0; j < places[i].size(); j++)
			cout << places[i][j] << " ";
		cout << endl;
	}
}

} // end of ctbn namespace
