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
#ifndef CTBNRLE_FACTORED_VECTOR_H
#define CTBNRLE_FACTORED_VECTOR_H

#include "ctbndyn.h"
#include "bn.h"
#include "rvsimple.h"

#include <iostream>
#include <vector>
#include <map>


namespace ctbn {

using namespace std;

class FactoredVector : public RVSimple {

public:
    FactoredVector();
    FactoredVector(const FactoredVector & v);
    FactoredVector(RV *rv);
	~FactoredVector();

	//copy constructor
	FactoredVector *Clone() const;
	
	void Restrict(const std::vector<int, std::allocator<int> >&);
	void RestrictById(int id, int val);
	void RestrictById(int id, int val, double x);
	void Condition(const map<int, int> &ev);
	double Normalize(); //(vector<double> &normf);
	double absnormalize(vectr &v);
	double maxnormalize(vectr &v);
	double NormalizeNeg();
	void GetDist(vectr &v, int varid, double &logfactor);
	// cshelton: made first arg const below
	void SetDist(const vectr &v, int varid, double logfactor=0.0);
	void SetLogf(double logfactor);
	double GetLogf() const;
	double GetDistInstNorm(int varid, int index);
	double GetDistInst(int varid, int index);
	void SetDistInst(int varid, int index, double value);
	int GetDistSize(int varid);
	int Size();
	std::ostream & Print(std::ostream & out, RV * const) const;
    std::ostream & PrintSimple(std::ostream & out) const;

	double GetMargMin();
	double GetJointMin();
	double GetMin(int varid);
	double GetJointMax();
	double GetMargMax();
	double GetMax(int varid);

	 void Load(std::istream &is);
	 void Save(std::ostream &os) const;
	 RVSimple *Condition(int ind) const;
	void GetDist(vectr &d, double &logfactor) const;
	void SetDist(const vectr &d, double logfactor=0.0);

	// returns the probability of an individual event (index)
	 double Prob(int ind, bool log=false) const;

	// Returns the sum (or log thereof)
	 double Sum(bool log=false) const;


	// Marginalizes, reindexes, or expands depending on the argument
	// ind is a list, for each new index, of the list of old indexes
	// which should be summed
	 void Reindex(const std::vector<std::vector<int> > &ind);

	// Multiplies by another measure (point-wise)
	 void MultBy(const RVSimple *x);

	void AddMult(const FactoredVector & v, double x);
	// Adds in another measure
	 void Add(const RVSimple *x);
	 void Add(const RVSimple *x, double w);

	// Multiplies measure by constant
	 void Mult(double x);

	// return suitable sufficient statistics object
	 SS *BlankSS() const;

	// add independent draw of x to ss
	 void AddSS(int x, SS *ss, double w=1.0) const;
	// add average of independent draws of x to ss
	 void AddExpSS(SS *ss, double w=1.0) const;

	 void AddSS(const SS *toadd, const RVSimple* rvs,
			const std::vector<std::vector<int> > &mapping, 
			SS *ss, double w=1.0) const;

	// returns a sample
	 int Sample(Random &rand = randomizer) const;

	// sets parameters to ML estimate
	 void Maximize(const SS *ss);

	// set measure to be 1 everywhere it has support
	 void MakeUniform();

	// Draws parameters from Dirichlet distribution
	// with all alphas the same (equal to alpha)
	// If degree equals 1, all parameters are replaced
	// If degree is between 0 and 1, the old parameters are
	// mixed with the new one (linear combination, degree controls
	// the mixing).  
	 void Scramble(double alpha=1.0, double degree=1.0, 
			Random &rand=randomizer); 
	// Returns the log-likelihood of a set of sufficient statistics
	 double LLH(const SS *ss) const;

	 double GetScore(double numTrans, const SS* ss) const;

private:
	std::vector<vectr> dists;
	double logf;

};

/*
inline FactoredVector operator*(const FactoredVector &v, const vector<double> &s) {
	FactoredVector newv = v;
	newv *= s;
	return newv;
}

inline  FactoredVector operator*(const vector<double> &s, const FactoredVector &v) {
	FactoredVector newv = v;
	newv *= s;
	return newv;
}

inline FactoredVector operator*(const FactoredVector &v, const double &s) {
	FactoredVector newv = v;
	newv *= s;
	return newv;
}

inline FactoredVector operator*(const double &s, const FactoredVector &v) {
	FactoredVector newv = v;
	newv *= s;
	return newv;
}
*/
}
#endif
