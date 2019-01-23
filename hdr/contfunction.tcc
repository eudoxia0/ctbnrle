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
#include "contfunction.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <utility>

namespace ctbn {

using namespace std;

template <class VALUE>
ContFunction<VALUE>::ContFunction() : start(0), end(0), empty(true) {
}

template <class VALUE>
ContFunction<VALUE>::~ContFunction() {
	Clear();
}

/*
template <class VALUE>
ContFunction<VALUE>::ContFunction(const ContFunction &rhs)
: start(rhs.start), end(rhs.end), empty(rhs.empty) {
	values.clear();
	typename map<double,VALUE>::const_iterator it = rhs.values.begin();
	for(; it != rhs.values.end(); ++it) {
		AddVal(it->first,(it->second));
	}
}
*/

template <class VALUE>
VALUE ContFunction<VALUE>::GetVal(double t) const {
	if(empty) {
		cerr << "no values initialized" << endl;
		exit(0);
	}

	if(start == end)
		if(fabs(t-start) < 1e-9) return values.begin()->second;

	if(fabs(t-start) < 1e-9) return values.begin()->second;
	if(fabs(t-end) < 1e-9) {
		typename map<double,VALUE>::const_iterator it = values.end();
		it--;
		return it->second;
	}

	if(t < start || t > end) {

		cerr << "t=" << t << endl;
		cerr << "start=" << start << endl;
		cerr << "end=" << end << endl;
		cerr << "time out of bounds" << endl;
		cin.get();
		exit(0);
	}


/*
	typename map<double,VALUE>::const_iterator it = values.find(t);
	if(it != values.end())
		return it->second;

	it = values.begin();
	while(it != values.end() && t > it->first) ++it;
*/
	typename map<double,VALUE>::const_iterator it = values.lower_bound(t);
	typename map<double,VALUE>::const_iterator ite = it;
	--it;

	//VALUE *vs = new VALUE(*(it->second));
	//VALUE *ve = new VALUE(*(ite->second));

	double ts = it->first;
	double te = ite->first;

	double c = (t-ts) / (te-ts);
	return ite->second*c + it->second*(1-c);//*ve - *vs;
	//ret *= c;
	//ret += *vs;

	//delete vs,
	//delete ve;

	//return ret;
}

template<class VALUE>
VALUE ContFunction<VALUE>::GetDeriv(double t) {
	if(empty) {
		cerr << "no values initialized" << endl;
		exit(0);
	}

	else if(start == end ||
		fabs(t-start) < 1e-9 || fabs(t-end) < 1e-9) {
		cout << "Not differentiable." << endl;
		exit(0);
	}
	else if(t < start || t > end) {

		cerr << "t=" << t << endl;
		cerr << "start=" << start << endl;
		cerr << "end=" << end << endl;
		cerr << "time out of bounds" << endl;
		cin.get();
		exit(0);
	}
/*
	typename map<double,VALUE>::const_iterator it = values.find(t);
	if(it != values.end()) {
		cout << "Not differentiable." << endl;
		exit(0);
	}

	it = values.begin();
	while(it != values.end() && t > it->first) ++it;
*/

	typename map<double,VALUE>::const_iterator it = values.lower_bound(t);
	typename map<double,VALUE>::const_iterator ite = it;
	--it;

	return (ite->second - it->second)/(ite->first-it->first);

/*
	VALUE *vs = new VALUE(*(it->second));
	VALUE *ve = new VALUE(*(ite->second));

	double ts = it->first;
	double te = ite->first;

	VALUE ret = *ve - *vs;
	ret /= te - ts;

	delete vs,
	delete ve;

	return ret;
*/
}

template<class VALUE>
VALUE ContFunction<VALUE>::GetIntercept(double t) {
	if(empty) {
		cerr << "no values initialized" << endl;
		exit(0);
	}

	else if(start == end ||
		fabs(t-start) < 1e-9 || fabs(t-end) < 1e-9) {
		cout << "No intercept." << endl;
		exit(0);
	}
	else if(t < start || t > end) {

		cerr << "t=" << t << endl;
		cerr << "start=" << start << endl;
		cerr << "end=" << end << endl;
		cerr << "time out of bounds" << endl;
		cin.get();
		exit(0);
	}

	typename map<double,VALUE>::const_iterator it = values.lower_bound(t);
	typename map<double,VALUE>::const_iterator ite = it;
	--it;

	double den = ite->first - it->first;
	return it->second * (ite->first/den) + ite->second * (-it->first/den);

/*
	VALUE m = GetDeriv(t);
	typename map<double,VALUE>::const_iterator it = values.find(t);
	if(it != values.end()) {
		cout << "No intercept." << endl;
		exit(0);
	}

	it = values.begin();
	while(it != values.end() && t > it->first) ++it;
	it--;

	VALUE *vs = new VALUE(*(it->second));
	double ts = it->first;

	m *= ts;
	VALUE ret = *vs - m;

	delete vs;

	return ret;
*/
}

template <class VALUE>
vector<double> ContFunction<VALUE>::GetCutPoints() const {
	vector<double> ret;
	typename map<double,VALUE>::const_iterator it = values.begin();
	for(; it != values.end(); ++it) ret.push_back(it->first);
	return ret;
}

template <class VALUE>
void ContFunction<VALUE>::Print() const {
	cout << "start: " << start << endl;
	cout << "end: " << end << endl;
	typename map<double,VALUE>::const_iterator it = values.begin();
	for(; it != values.end(); ++it)
		cout << it->first << "\t" << it->second << endl;
	cout << endl;
}

template <class VALUE>
void ContFunction<VALUE>::AddVal(double t, const VALUE& val) {
	if(values.find(t) == values.end()) {
		values.insert(std::make_pair(t,val));
		if(empty) {
			start = t;
			end = t;
			empty = false;
		}
		if(t < start) start = t;
		if(t > end) end = t;
	}
}

template<class VALUE>
void ContFunction<VALUE>::Combine(const ContFunction<VALUE> &fun) {
	typename map<double,VALUE>::const_iterator it = fun.values.begin();
	for(; it!=fun.values.end(); it++)
		AddVal(it->first,it->second);
}

template<class VALUE>
void ContFunction<VALUE>::EraseRange(double t0, double t1) {
	//typename map<double,VALUE>::iterator it = values.begin();
	//while(it->first < t0) it++;
	typename map<double,VALUE>::iterator it = values.lower_bound(t0);
	//typename map<double,VALUE>::iterator itend = it;
	//while(itend->first < t1 && itend != values.end()) itend++;
	typename map<double,VALUE>::iterator itend = values.lower_bound(t1+1e-12);
	//while(itend != values.end() && abs(itend->first - t1) < 1e-12) itend++;
	values.erase(it,itend);
/*
	vector<double> times;
	for(;it != itend; it++) times.push_back(it->first);

	for(unsigned int i=0; i<times.size(); i++) EraseVal(times[i]);
*/
}

template<class VALUE>
void ContFunction<VALUE>::EraseVal(double t) {
	values.erase(t);
/*
	typename map<double,VALUE *>::iterator it = values.find(t);
	if(it != values.end()) {
		//delete it->second;
		values.erase(it);
	}
*/
}

template<class VALUE>
void ContFunction<VALUE>::Replace(const ContFunction<VALUE> &fun) {
	EraseRange(fun.GetStart(),fun.GetEnd());
	Combine(fun);
}

template <class VALUE>
void ContFunction<VALUE>::Clear() {
/*
	typename map<double,VALUE *>::iterator it = values.begin();
	for(; it != values.end(); ++it) {
		if(it->second) {
			delete it->second;
			it->second = NULL;
		}
	}
*/
	values.clear();
	start = 0;
	end = 0;
	empty = true;
}

} // end of ctbn namespace
