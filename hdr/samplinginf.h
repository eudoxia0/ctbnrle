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
#ifndef CTBNRLE_SAMPLINGINF_H
#define CTBNRLE_SAMPLINGINF_H

#include "inference.h"
#include "sampler.h"
#include "trajectory.h"
#include "context.h"
#include "markov.h"
#include "trajectory.h"
#include "ss.h"
#include <vector> 

namespace ctbn {

class SamplingInfAbs : public Inference {
public:
	SamplingInfAbs();
	virtual ~SamplingInfAbs();
	virtual SamplingInfAbs *Clone() const=0;
    //The caller of the function owns the pointers
	//  (but they must remain valid for the lifetime of the object)
	virtual void SetProcess(const Process *pr);
	virtual void SetTrajectory(const Trajectory *tr);
	virtual void SetSampler(const Sampler *s);
	void SetNumSample(int n);
	void SetTimeLimit(double t, int batchsize=10);

	virtual double Filter(const Instantiation &x, double t, 
				bool log=false) = 0;
	virtual double Smooth(const Instantiation &x, double t, 
				bool log=false) = 0;
	virtual double Prob(double t, bool log=false) = 0;
    //the caller of the function own the pointer
	virtual void AddExpSuffStats(SS *ss, double w=1.0) = 0;
    //The caller of the function own the pointers
	virtual void AddExpSuffStats(const Dynamics *dyn, SS *ss,
					double w=1.0) = 0;
    //The caller of the function own the pointers
	virtual void AddExpSuffStats(const RV *p0, SS *ss, 
					double w=1.0) = 0;

	//double LLH(const std::vector<Trajectory> &tr, 
	//		const std::vector<double> &w) const;
	double CalcQuery(QueryCalculator &calc) = 0;
protected:
	virtual void Sample();
	virtual void AddSamples(int nsamp, bool all) = 0;
	virtual void Reset() = 0;

	static double NormalizeSet(std::vector<double> &w);

	const Markov *p; //the class does not own this pointer
	Sampler *sampler; //SamplingInfAbs own this pointer

	int numofsamples;
	double timelimit;
	int bsize; // # of traj to grab at a time when time limited
};

// See inference.h for a description of the methods
class SamplingInf : public SamplingInfAbs {
public:
	SamplingInf();
	virtual ~SamplingInf();
	virtual SamplingInf *Clone() const;

	virtual double Filter(const Instantiation &x, double t, 
				bool log=false);
	virtual double Smooth(const Instantiation &x, double t, 
				bool log=false);
	virtual double Prob(double t, bool log=false);
    //the caller of the function own the pointer
	virtual void AddExpSuffStats(SS *ss, double w=1.0);
    //The caller of the function own the pointers
	virtual void AddExpSuffStats(const Dynamics *dyn, SS *ss,
					double w=1.0);
    //The caller of the function own the pointers
	virtual void AddExpSuffStats(const RV *p0, SS *ss, 
					double w=1.0);

	//double LLH(const std::vector<Trajectory> &tr, 
	//		const std::vector<double> &w) const;
	double CalcQuery(QueryCalculator &calc);
protected:
	virtual void Sample();
	virtual void AddSamples(int nsamp, bool all);
	virtual void Reset();

	std::vector<Trajectory> traj;
	std::vector<double> weights;
};

// Same as above, but the queries to be asked are known in
// advance, so the samples need not be stored.
class SamplingPreInf : public SamplingInfAbs {
public:
	// these pointers need to stay valid for the life of
	// this object.  Further, references to the exact same queries must be
	// used in CalcQuery.
	SamplingPreInf(const std::vector<QueryCalculator *> &queries,
			int workingsetsize=1000);
	virtual ~SamplingPreInf();
	virtual SamplingPreInf *Clone() const;

	virtual double Filter(const Instantiation &x, double t, 
				bool log=false);
	virtual double Smooth(const Instantiation &x, double t, 
				bool log=false);
	virtual double Prob(double t, bool log=false);
    //the caller of the function own the pointer
	virtual void AddExpSuffStats(SS *ss, double w=1.0);
    //The caller of the function own the pointers
	virtual void AddExpSuffStats(const Dynamics *dyn, SS *ss,
					double w=1.0);
    //The caller of the function own the pointers
	virtual void AddExpSuffStats(const RV *p0, SS *ss, 
					double w=1.0);

	//double LLH(const std::vector<Trajectory> &tr, 
	//		const std::vector<double> &w) const;
	double CalcQuery(QueryCalculator &calc);
protected:
	virtual void AddSamples(int nsamp, bool all);
	virtual void Reset();

	bool sampled;
	double maxlnw,ttlw;
	int setsize;
	std::vector<QueryCalculator *> qs;
	std::vector<double> ans;
};

} // end of ctbn namespace

#endif

