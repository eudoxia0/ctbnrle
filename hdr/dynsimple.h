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
#ifndef CTBNRLE_DYNSIMPLE_H
#define CTBNRLE_DYNSIMPLE_H

#include "random.h"
#include "ss.h"
#include "rvsimple.h"

#include <vector>

#include <vector>


namespace ctbn {

// DynSimple is the ABC for continuous-time semi-Markov dynamics
// Just like RVSimple, its domain is implied
// So this is an abstract representation of a semi-Markov process
// dynamics over the state space 0..n-1 with no conditioning
class DynSimple : public StreamObj {
public:
	DynSimple();
	DynSimple(std::istream &is);
	virtual ~DynSimple();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	// virtual copy constructor
	virtual DynSimple *Clone() const = 0;
	// virtual blank constructor
	virtual DynSimple *MakeNew(int nstates) const = 0;

	// Conditions the system being closed
	virtual void Normalize() = 0;

	// Sets to zero all probability of leaving the set of ind
	virtual void Restrict(const std::vector<int> &ind) = 0;

	// Changes the indexing.
	// newn is the new number of states
	// (let n be the old number of states)
	// ind is a list of size n
	// each element of ind is a list of size s (s is a constant)
	// each of these sublists is a set of states that were
	//   indistinguishable from the original dynamics point of view
	virtual void Expand(const std::vector<std::vector<int> > &ind,
			int newn) = 0;

	// Returns the conditional distribution of going from t0 to t1
	virtual RVCondSimple *Cond(double t0, double t1) const = 0;
	// same, but with the dynamics constrainted to the subsystem
	// described by ind
	virtual RVCondSimple *CondRestrict(double t0, double t1,
			const std::vector<int> &ind) const = 0;

			
	// Returns the density of instantaneous transitioning
	virtual RVCondSimple *Cond(double t) const = 0;
	// Same, but constrained (see same method in Dynamics)
	virtual RVCondSimple *CondRestrict(double t,
			const std::vector<int> &fromind,
			const std::vector<int> &toind, 
			bool transition=true) const = 0;

	// returns a suitable sufficient statistics object
	virtual SS *BlankSS() const = 0;

	// the amalgamation operation... assumes that indices
	// already agree...
	virtual void Mult(const DynSimple *x) = 0;

	// add interval to ss
	virtual void AddSS(int x, double t0, double t1,
			SS *ss, double w=1.0) const = 0;
	// add transition to ss
	virtual void AddTransSS(int x1, int x2, double t,
			SS *ss, double w=1.0) const = 0;
	// adds expected ss (assumes that alpha and beta have the
	//   same domain and that domain is the one to which the data
	//   restricts transitions)
	virtual void AddExpSS(const RVSimple *alpha, const RVSimple *beta,
			double t0, double deltat, SS *ss, 
			double w=1.0) const = 0;
	// adds a transition from x1 to x2
	virtual void AddExpTransSS(const RVSimple *x1, const RVSimple *x2,
			double t, SS *ss, double w=1.0) const = 0;

	// Add the sufficient statistics in "toadd" to those in "ss"
	// However, using "mapping" to map state-indexes for toadd to
	// those to ss.  "ss" should be for "this" while "toadd" should be
	// for "dyn."  Currently dyn needs to be of type MarkovSimple.
	// mapping[i] is a list of the states from toadd/dyn that correspond
	// to state i for ss/this.
	// [The code might also clarify -- see .cc file.]
	virtual void AddSS(const SS *toadd, const DynSimple *dyn,
			const std::vector<std::vector<int> > &mapping,
			SS *ss, double w=1.0) const = 0;

	// sample the next event - time: newt, state index: newind 
	virtual void SampleNextEvent(int ind, double t,
			double &newt, int &newind, 
			Random &rand=randomizer) const = 0;

	// sets parameters from ML estimate
	virtual void Maximize(const SS *ss) = 0;

	// Randomize the parameters
	virtual void Scramble(double a=1.0, double b=1.0, 
				double alpha=1.0, double degree=1.0, 
				Random &rand=randomizer) = 0; 
	// calculate the log-likelihood of the sufficient SS
	virtual double LLH(const SS *ss) const = 0;

	// Returns the structure search score where parentCard is
	// the number of parents this node has (this object represents
	// the process for just one parent instantiation, presumably)
	// numTrans and amtTime are the BDe prior parameters
	virtual double GetScore(int parentCard, 
				double numTrans, 
				double amtTime, 
				const SS* ss) const = 0;

	SERIAL_START_V(DynSimple)
	SERIAL_END
};

} // end of ctbn namespace

#endif
