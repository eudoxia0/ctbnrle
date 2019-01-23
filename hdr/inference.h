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
#ifndef CTBNRLE_INFERENCE_H
#define CTBNRLE_INFERENCE_H

#include "context.h"
#include "bn.h"
#include "ss.h"
#include "process.h"



namespace ctbn {

// QueryCalculator is the base class to answer any query of 
// a process given a trajectory or sufficient statistics.
class QueryCalculator {
public:
	// virtual function that calculates queries, 
	// to be implemented by the inheriting classes 
	// QueryTime and QueryTransition
	// calculate query given a trajectory
	virtual double Calculate(const Trajectory &tr) = 0;
	// calculate query given sufficient statistics
	virtual double Calculate(const SS *ss) = 0;
	virtual ~QueryCalculator() throw() {};
};

// QueryTime is used to answer query that the total time
// a process stays in some state. The state is passed through
// constructor using an Instantiation c.  
class QueryTime : public QueryCalculator {
public:
	QueryTime(const Instantiation &c);
	double Calculate(const Trajectory &tr);
	double Calculate(const SS *ss);
	virtual ~QueryTime() throw() {};
private:
	Instantiation condi;
	int queryindex;
};

// QueryTransition is used to answer query that the total number
// of times a process transitions from one state to another. 
class QueryTransition : public QueryCalculator {
public:
	QueryTransition(const Context &c, int from, int to);
	double Calculate(const Trajectory &tr);
	double Calculate(const SS *ss);
	virtual ~QueryTransition() throw() {};
private:
	Context querycontext;
	int queryfrom;
	int queryto;
};

// this one only works from trajectories (cannot be calculated 
// from suffstats) -- useful for passing into sampling inferences
// (does the same thing as calling "Smooth")
class QueryProb : public QueryCalculator {
public:
	QueryProb(const Instantiation &val, double time);
	double Calculate(const Trajectory &tr);
	double Calculate(const SS *ss);
	virtual ~QueryProb() throw() {};
private:
	Context x;
	int xind;
	double t;
};

// Inference is not streamable for many reasons
// (although it could potentially be useful to be able
// to serialize it)
class Inference { // : public StreamObj {
	//SOBJCLASSDECL(Inference)
public:
	virtual ~Inference();
	virtual Inference *Clone() const = 0;

	// The Inference object owns neither of these
	// (they are still owned by the caller, but the
	// pointers passed in must remain valid (and unchanging)
	// for the life of the Inference object, or at least
	// until the next "SetProcess" or "SetTrajectory" (respectively)
	// is called
	virtual void SetProcess(const Process *p) = 0;
	virtual void SetTrajectory(const Trajectory *tr) = 0;

	// returns the probability of x at t given the trajectory up to time t
	virtual double Filter(const Instantiation &x, double t, bool log=false) = 0;

	// returns the probability of x at t given the whole trajectory
	virtual double Smooth(const Instantiation &x, double t, bool log=false) = 0;

	// returns the probability of the trajectory up to time t
	virtual double Prob(double t, bool log=false) = 0;

	// adds the Expected Sufficient Statistics for the process
	// given the trajectory
	virtual void AddExpSuffStats(SS *ss, double w=1.0) = 0;
	// adds the ESS for the dynamics given by dyn to ss
	// (note that the expectation is wrt to the process set 
	// in "SetProcess" but the sufficient statistics are those 
	// for the dynamics dyn)
	// used primarily for structure search
	virtual void AddExpSuffStats(const Dynamics *dyn, SS *ss,
					 double w=1.0) = 0;
	// same as above, but for an initial starting distribution p0
	virtual void AddExpSuffStats(const RV *p0, SS *ss,
					 double w=1.0) = 0;

	virtual double CalcQuery(QueryCalculator &calc) = 0;
};

} // end of ctbn namespace

#endif
