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
#ifndef CTBNRLE_SAMPLER_H
#define CTBNRLE_SAMPLER_H

#include "process.h"
#include "trajectory.h"
#include "random.h" 
#include "varsample.h"
#include "structure.h"

namespace ctbn {

//Sampler abstract class for sampling method for dynamic process

class Sampler{
public:
	virtual ~Sampler();
	virtual Sampler *Clone() const = 0;
	//set the dynamic process, the caller of the function
    //owns the pointer of the process variable,
	// but it must remain valid for the lifetime of this object
	virtual void SetProcess(const Process *pr)=0;
	//set the evidence, the caller of the function owns the pointer 
	// but it must remain valid for the lifetime of this object
	virtual void SetTrajectory(const Trajectory *traj); 
	//set the sampling method, the caller of the function owns the pointer
	// but it must remain valid for the lifetime of this object
	virtual void SetMethod(const VarSample *m);
	//sampling trajectories given the evidence
	virtual void SampleTrajectories(std::vector<Trajectory> &tr, std::vector<double> &w, 
			int numsamples, Random &rand=randomizer) = 0;
protected:
	//sample the initial values
	virtual double SampleInitial(Instantiation &inst, Random &rand=randomizer) = 0;
	//sample the dynamics (complete the trajectory assuming the initial
	// state has been completed)
	virtual double SampleDyn(Trajectory &tr, Random &rand=randomizer) = 0;
    //The class does not own these pointers
	const Process *p;
	const Trajectory *evid;
	const VarSample *method;
	double begintime;
	double endtime;
};

} // end of ctbn namespace

#endif
