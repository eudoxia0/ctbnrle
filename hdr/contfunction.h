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
#ifndef CTBNRLE_CONTFUNCTION_H_
#define CTBNRLE_CONTFUNCTION_H_
#include <map>
#include <vector>

namespace ctbn {

// This class allows a function of the form double -> any numerical type
// to be represented by linear interpolation

template <class VALUE>
class ContFunction {
public:
	ContFunction();
	virtual ~ContFunction();
	// default okay now that map stores VALUE (not VALUE*)
	//ContFunction(const ContFunction &rhs);

	// Returns the value at time t
	VALUE GetVal(double t) const;

	// Returns the derivative at time t (not defined at defined t points)
	VALUE GetDeriv(double t);

	// Returns the y-intercept at time t (not defined at defined t points)
	VALUE GetIntercept(double t);

	// Adds a value to the function at time t
	void AddVal(double t, const VALUE& val);

	// Combines the values stored in fun to this function
	void Combine(const ContFunction<VALUE> &fun);

	// Erases the values from t0 to t1, inclusively
	void EraseRange(double t0, double t1);

	// Erases the value at time t
	void EraseVal(double t);

	// Replaces the values in the range of fun with the values of fun
	void Replace(const ContFunction<VALUE> &fun);

	double GetStart() const {return start;}
	double GetEnd() const {return end;}

	// Returns a vector with all the times that a value is stored
	std::vector<double> GetCutPoints() const;

	void Print() const;
	void Clear();

private:
	std::map<double,VALUE> values;
	double start, end;
	bool empty;
};


} // end of name ctbn namespace

#include "contfunction.tcc"
#endif
