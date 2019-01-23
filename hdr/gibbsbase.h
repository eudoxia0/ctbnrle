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
#ifndef CTBNRLE_GIBBSBASE_H
#define CTBNRLE_GIBBSBASE_H

#include <map>
#include <vector>
#include "sampler.h"



namespace ctbn {

class matrix;
class MultiSimple;
class RV;
class MarkovDyn;

// abstract base class for Gibbs sampling.  Takes care of pulling each variable
// out in turn, splitting trajectory into regions where variable is free, etc.
// handles burnin, etc.
// Sample usage:

class GibbsBase : public Sampler
{
 public:
	// the caller of the function owns the process and evidence pointers
	 // but they must remain valid for the life of the GibbsBase object
	GibbsBase(const Process *process, const Trajectory *evidence, 
				 int burnin);
	virtual ~GibbsBase() { if (init_traj) delete init_traj; }

	// the caller of the function owns the process pointer 
	 // but it must remain valid for the life of the GibbsSampler object
	virtual void SetProcess(const Process *process) { 
		p = process; Initialize(); 
	}
	// the caller of the function owns the evidence pointer
	 // but it must remain valid for the life of the GibbsSampler object
	virtual void SetTrajectory(const Trajectory *evidence) { 
		Sampler::SetTrajectory(evidence); burntin=false;
	}
	virtual void SampleTrajectories(std::vector<Trajectory> &traj, std::vector<double> &w,
					int numsamples, Random &rand = randomizer);


 protected:
	// specific public interfaces for GibbsSampler
	// burn in, do numBurninIter steps of markov transitions
	void BurnIn(Random &rand = randomizer) const;

	// one step in the markov chain state transition; 
	void Next(int num_iter = 1, Random &rand = randomizer) const;

	// get the current trajectory
	const Trajectory &Get() const { return tr; }
	// get the start point - the initial trajectory
	const Trajectory &GetInitTraj() const { return *init_traj; }
	// set the initial trajectory
	void SetInitTraj(const Trajectory &t);
	// delete the initial trajetcory
	void ClearInitTraj();


	int numBurninIter; // # of burn in iterations
	Context context;
	std::vector<int> own_var_list;
	const RV *p0;
	const MarkovDyn *md;
	// mutable so that we can access parents[v] 
	// etc. easily in const member fuctions.
	// each variable's parents and children
	mutable std::map<int, Context> parents, children;
	// each variable's markovdyn
	mutable std::map<int, const MarkovDyn *> mds; 
	// each variable's Markov Blanket
	mutable std::map<int, Context> markovblanket;
	// initial trajectory
	mutable Trajectory *init_traj; 
	// the current sample (and previous one for MCMC samplers)
	mutable Trajectory tr,oldtr; 

	mutable bool burntin;

	// inherited functions from base class, currently do nothing
	double SampleInitial(Instantiation &inst, Random &rand=randomizer){return 0.0;};
	double SampleDyn(Trajectory &traj, Random &rand=randomizer){return 0.0;};

	// precompute context, p0, dynamics, parents/children and markovdyn
	// for each variable;
	void Initialize();
	// Find the {tau} splitting points (timing, full instantiation) for 
	// transitioned variables who are in v's Markov blanket
	static void GetTau(const Trajectory &t, double begin_time, double end_time, 
				int v, const Context &varset, std::vector<double> &tau, 
				std::vector<Instantiation> &inst);
	// Resample the entire trajectory of v given all the other variables'
	// full trajectory.
	void SampleVariable(int v, Random &rand = randomizer) const;
	// Sample initial trajectory that agrees with evidence *evid.
	void SampleInitialTrajectory(Random &rand = randomizer) const;
	
 protected:
	// find evidence on v, and split into "free" intervals for v
	// there is no v's evidence inside each [t0, tT] interval
	void SplitFreeIntervals(int v, 
				std::vector<int> &x0, 
				std::vector<int> &xT, 
				std::vector<double> &t0, 
				std::vector<double> &tT, 
				const Instantiation &inst0, 
				Random &rand = randomizer) const;
	// sampler v's trajectory for interval [t0, tT];
	virtual void SampleVariableInterval(int v, int x0, int xT, double t0, double tT,
					Random &rand = randomizer) const = 0;

	// returns the diagonal transition matrix for v when Markov Blanket
	// changes from oldinst to newinst
	// returns false iff no variable in child set changed
	// NOTE: assumes that ret is the right size
	// ret is just the diagonal of the T Matrix (which is a diag matrix)
	bool GetTMatrix(int v, const Instantiation &oldinst,
				const Instantiation &newinst, vectr &ret) const;
};

} // end of ctbn namespace

#endif
