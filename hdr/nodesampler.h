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
#ifndef CTBNRLE_NODE_SAMPLER_
#define CTBNRLE_NODE_SAMPLER_

#include "impbase.h"
#include "sampler.h"
#include "markovdyn.h"
#include "ctbn.h"
#include "trajectory.h"
#include "structure.h"
#include <vector> 
#include "varsample.h"
#include "process.h"

namespace ctbn {

class NodeSampler : public ImpBase{
public:
	NodeSampler(const Process *pr, const Trajectory *traj, const VarSample *m);
	virtual ~NodeSampler(){}
	NodeSampler *Clone() const;
	void SampleTrajectories(std::vector<Trajectory> &tr, std::vector<double> &w,
			int numsamples, Random &rand=randomizer){} ;
	double SampleNodeTrajectory(const Trajectory &tr, 
			Trajectory &nodetraj,
			int nodeid, Random &rand=randomizer);
	double NodeWeight(const Trajectory &tr, int nodeid) const {return 0.0;};

protected:
	void GetChildrenContext(); //get the context of each node's children and their parents
	double ChildrenSampleWeight(const Trajectory &tr, int nodeid) const;
	//sample the initial values
	virtual double SampleInitial(Instantiation &inst, Random &rand=randomizer){return 0.0;};
	//sample the dynamics
	virtual double SampleDyn(Trajectory &tr, Random &rand=randomizer){return 0.0;};

	std::vector<Context> childrencontext;

};

} // end of ctbn namespace
#endif
