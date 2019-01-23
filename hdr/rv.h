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

#ifndef CTBNRLE_RV_H
#define CTBNRLE_RV_H

#include "context.h"
#include "trajectory.h"
#include "random.h"
#include "streamserial.h"


namespace ctbn {

class RVSimple;
class SS;

// RV is a random variable.  It can be conditional and joint
// It keeps track of its context
class RV : public StreamObj {
public:
	// var = varaibles, cvar = conditioning variables
	RV(const Context &var, const Context &cvar);
	RV(std::istream &is);
	virtual ~RV();

	virtual void LoadOld(std::istream &is);
	virtual void SaveOld(std::ostream &os) const;

	// virtual copy constructor
	virtual RV *Clone() const = 0;

	const Context &Domain() const { return v; }
	const Context &CondDomain() const { return cv; }

	// makes sum to 1 for each conditional instantiation
	virtual void Normalize() = 0;

	// returns the probability (or log thereof)
	virtual double Prob(const Instantiation &x, bool log=false) const = 0;

	// Sets to zero all probabilities that are not consistent with v
	// (followed by a normalize, this is the same as conditioning
	//  except that it doesn't change the domain of the variable)
	// (if again followed by a Project, then you have conditioning)
	virtual void Restrict(const Instantiation &x) = 0;

	// Marginalizes to c
	// ids not mentioned in c are marginalized out
	// can also be used to transfer ids to/from conditioning
	virtual void Project(const Context &c, const Context &cc) {
		v = c; cv = cc;
	}

	// Multiplies by another measure (point-wise)
	// Should expand context to include union
	virtual void MultBy(const RV *x) = 0;

	// returns a suitable sufficient statistics object
	virtual SS *BlankSS() const =0;

	// adds the instantiation to ss
	virtual void AddSS(const Instantiation &x, SS *ss, double w=1.0) const = 0;

	// adds the expected sufficient statistic to ss
	// (ie x can be incomplete -- but not in the conditioning variables
	virtual void AddExpSS(const Instantiation &x, SS *ss, double w=1.0) const = 0;

	// adds the distribution to the sufficient statistics
	// (x only gives the conditional context.  rvs gives a distribution
	// over the domain of this variable)
	virtual void AddExpSS(const Instantiation &x, const RVSimple *rvs, 
			SS *ss, double w=1.0) const = 0;

	//This adds the sufficient statistics of a different random variable
	//to this random variable's sufficient statistics
	virtual void AddSS(const SS* toadd, const RV *rv,
			SS *ss, double w=1.0) const = 0;
	// Completes the instantiation v with a conditional sample
	virtual void Sample(Instantiation &x, Random &rand = randomizer) const = 0;

	// sets parameters to ML estimate
	virtual void Maximize(const SS *ss) = 0;

	// returns an RVSimple object (owned by caller) that represents
	// the distribution conditioned on x (note that x must include
	// a complete assignment to all variables in the CondDomain()
	// context) without explicit reference to the domain
	// (normalize indicates whether the RV should be renormalized
	//  -- in the case where x has assignments in Domain())
	virtual RVSimple * MakeSimple(Instantiation const & x, bool normalize = true) const = 0;

	virtual void Scramble(double alpha=1.0, double degree=1.0, 
			Random &rand=randomizer) = 0; 
	// computes the log-likelihood from the sufficient statistics passed in
	virtual double LLH(const SS *ss) const = 0;
	// scores this variable (and the implied parent assignment given
	//  by the conditional context) by the Bayesian structure score
	//  with a uniform BDe prior with numTrans strength.
	virtual double GetScore(double numTrans, const SS *aa) const = 0;

	SS* SuffStats(const std::vector<Trajectory> &tr, 
			const std::vector<double> &w) const;
protected:
	Context v,cv;

	SERIAL_START_V(RV)
		SERIAL_VAR(Context,cv)
		SERIAL_VAR(Context,v)
	SERIAL_END
};

} // end of ctbn namespace

#endif
